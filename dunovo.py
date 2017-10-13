#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import argparse
import collections
import multiprocessing
import consensus
import swalign
import shims
# There can be problems with the submodules, but none are essential.
# Try to load these modules, but if there's a problem, load a harmless dummy and continue.
simplewrap = shims.get_module_or_shim('utillib.simplewrap')
version = shims.get_module_or_shim('utillib.version')
phone = shims.get_module_or_shim('ET.phone')

SANGER_START = 33
SOLEXA_START = 64
USAGE = '$ %(prog)s [options] families.msa.tsv > duplex-consensuses.fa'
DESCRIPTION = """Build consensus sequences from read aligned families. Prints duplex consensus \
sequences in FASTA to stdout. The sequence ids are BARCODE.MATE, e.g. "CTCAGATAACATACCTTATATGCA.1", \
where "BARCODE" is the input barcode, and "MATE" is "1" or "2" as an arbitrary designation of the \
two reads in the pair. The id is followed by the count of the number of reads in the two families \
(one from each strand) that make up the duplex, in the format READS1/READS2. If the duplex is \
actually a single-strand consensus because the matching strand is missing, only one number is \
listed.
Rules for consensus building: Single-strand consensus sequences are made by counting how many of \
each base are at a given position. Bases with a PHRED quality score below the --qual threshold are \
not counted. If a majority of the reads (that pass the --qual threshold at that position) have one \
base at that position, then that base is used as the consensus base. If no base has a majority, then \
an N is used. Duplex consensus sequences are made by aligning pairs of single-strand consensuses, \
and comparing bases at each position. If they agree, that base is used in the consensus. Otherwise, \
the IUPAC ambiguity code for both bases is used (N + anything and gap + non-gap result in an N)."""


def make_argparser():

  wrapper = simplewrap.Wrapper()
  wrap = wrapper.wrap
  parser = argparse.ArgumentParser(usage=USAGE, description=wrap(DESCRIPTION), add_help=False,
                                   formatter_class=argparse.RawTextHelpFormatter)

  wrapper.width = wrapper.width - 24
  io = parser.add_argument_group('Inputs and outputs')
  io.add_argument('infile', metavar='families.msa.tsv', nargs='?', default=sys.stdin,
                  type=argparse.FileType('r'),
    help=wrap('The output of align_families.py. 6 columns:\n'
              '1. (canonical) barcode\n'
              '2. order ("ab" or "ba")\n'
              '3. mate ("1" or "2")\n'
              '4. read name\n'
              '5. aligned sequence\n'
              '6. aligned quality scores.'))
  io.add_argument('-1', '--dcs1', metavar='duplex_1.fa', type=argparse.FileType('w'),
    help=wrap('The file to output the first mates of the duplex consensus sequences into. '
              'Warning: This will be overwritten if it exists!'))
  io.add_argument('-2', '--dcs2', metavar='duplex_2.fa', type=argparse.FileType('w'),
    help=wrap('Same, but for mate 2.'))
  io.add_argument('--sscs1', metavar='sscs_1.fa', type=argparse.FileType('w'),
    help=wrap('Save the single-strand consensus sequences (mate 1) in this file (FASTA format). '
              'Warning: This will be overwritten if it exists!'))
  io.add_argument('--sscs2', metavar='sscs_2.fa', type=argparse.FileType('w'),
    help=wrap('Save the single-strand consensus sequences (mate 2) in this file (FASTA format). '
              'Warning: This will be overwritten if it exists!'))
  io.add_argument('--incl-sscs', action='store_true',
    help=wrap('When outputting duplex consensus sequences, include reads without a full duplex '
              '(missing one strand). The result will just be the single-strand consensus of the '
              'remaining read.'))
  params = parser.add_argument_group('Algorithm parameters')
  params.add_argument('-r', '--min-reads', type=int, default=3,
    help=wrap('The minimum number of reads (from each strand) required to form a single-strand '
              'consensus. Strands with fewer reads will be skipped. Default: %(default)s.'))
  params.add_argument('-q', '--qual', type=int, default=20,
    help=wrap('Base quality threshold. Bases below this quality will not be counted. '
              'Default: %(default)s.'))
  params.add_argument('-F', '--qual-format', choices=('sanger', 'solexa'), default='sanger',
    help=wrap('FASTQ quality score format. Sanger scores are assumed to begin at \'{}\' ({}). '
              'Default: %(default)s.'.format(SANGER_START, chr(SANGER_START))))
  params.add_argument('-c', '--cons-thres', type=float, default=0.5,
    help=wrap('The threshold to use when making consensus sequences. The consensus base must be '
              'present in more than this fraction of the reads, or N will be used. '
              'Default: %(default)s'))
  params.add_argument('-C', '--min-cons-reads', type=int, default=0,
    help=wrap('The minimum number of reads a base must appear in to be used as the consensus base. '
              'If no base at the position appears in at least this many reads, N will be used as '
              'the consensus base. Default: %(default)s'))
  phoning = parser.add_argument_group('Feedback')
  phoning.add_argument('--phone-home', action='store_true',
    help=wrap('Report helpful usage data to the developer, to better understand the use cases and '
              'performance of the tool. The only data which will be recorded is the name and '
              'version of the tool, the size of the input data, the time taken to process it, and '
              'the IP address of the machine running it. No filenames are sent, and the only '
              'parameter reported is the number of processes. All the reporting and recording code '
              'is available at https://github.com/NickSto/ET.'))
  phoning.add_argument('--galaxy', dest='platform', action='store_const', const='galaxy',
    help=wrap('Tell the script it\'s running on Galaxy. Currently this only affects data reported '
              'when phoning home.'))
  phoning.add_argument('--test', action='store_true',
    help=wrap('If reporting usage data, mark this as a test run.'))
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help=wrap('Print log messages to this file instead of to stderr. Warning: Will overwrite the '
              'file.'))
  log.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
                      default=logging.WARNING)
  log.add_argument('-V', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  log.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  misc = parser.add_argument_group('Miscellaneous')
  misc.add_argument('-p', '--processes', type=int, default=1,
    help=wrap('Number of worker subprocesses to use. Default: %(default)s.'))
  misc.add_argument('--max-results', type=int,
    help=wrap('How long to go accumulating responses from worker subprocesses before dealing '
              'with all of them. Default: the same as the number of workers.'))
  misc.add_argument('-v', '--version', action='version', version=str(version.get_version()),
    help=wrap('Print the version number and exit.'))
  misc.add_argument('-h', '--help', action='store_true',
    help='Print this text on usage and arguments.')

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  if args.help:
    parser.print_help()
    return 0

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  start_time = time.time()
  # If the user requested, report back some data about the start of the run.
  if args.phone_home:
    call = phone.Call(__file__, version.get_version(), platform=args.platform, test=args.test,
                      fail='warn')
    call.send_data('start')
    if args.infile is sys.stdin:
      data = {'stdin':True, 'input_size':None}
    else:
      data = {'stdin':False, 'input_size':os.path.getsize(args.infile.name)}
    call.send_data('prelim', run_data=data)

  # Process and validate arguments.
  if args.max_results is None:
    max_results = args.processes
  else:
    max_results = args.max_results
  if args.processes <= 0:
    fail('Error: --processes must be greater than zero.')
  if max_results <= 0:
    fail('Error: --max-results must be greater than zero.')
  if args.qual_format == 'sanger':
    qual_thres = chr(args.qual + SANGER_START)
  elif args.qual_format == 'solexa':
    qual_thres = chr(args.qual + SOLEXA_START)
  else:
    fail('Error: unrecognized --qual-format.')
  if args.min_cons_reads > args.min_reads:
    fail('Error: --min-reads must be greater than --min-cons-reads (or you\'ll have a lot of '
         'consensus sequences with only N\'s!). If you want to exclude families with fewer than X '
         'reads, give --min-reads X instead of --min-cons-reads X.')
  if not any((args.dcs1, args.dcs2, args.sscs1, args.sscs2)):
    fail('Error: must specify an output file!')
  # A dict of output filehandles.
  # 1-indexed so we can do filehandles['dcs'][mate].
  filehandles = {
    'dcs': (None, args.dcs1, args.dcs2),
    'sscs': (None, args.sscs1, args.sscs2),
  }

  # Open a pool of worker processes, and a list to cache results in.
  static_args = [args.incl_sscs, args.min_reads, args.cons_thres, args.min_cons_reads, qual_thres]
  pool = multiprocessing.Pool(args.processes)
  results = []

  try:

    stats = {'time':0, 'reads':0, 'runs':0, 'duplexes':0}
    all_reads = 0
    duplex = collections.OrderedDict()
    family = []
    barcode = None
    order = None
    mate = None
    for line in args.infile:
      # Allow comments (e.g. for test input files).
      if line.startswith('#'):
        continue
      fields = line.rstrip('\r\n').split('\t')
      if len(fields) != 6:
        continue
      this_barcode, this_order, this_mate, name, seq, qual = fields
      this_mate = int(this_mate)
      # If the barcode or order has changed, we're in a new single-stranded family.
      # Process the reads we've previously gathered as one family and start a new family.
      new_barcode = this_barcode != barcode
      new_order = this_order != order
      new_mate = this_mate != mate
      if new_barcode or new_order or new_mate:
        if order is not None and mate is not None:
          duplex[(order, mate)] = family
        # We're at the end of the duplex pair if the barcode changes or if the order changes without
        # the mate changing, or vice versa (the second read in each duplex comes when the barcode
        # stays the same while both the order and mate switch). Process the duplex and start
        # a new one. If the barcode is the same, we're in the same duplex, but we've switched strands.
        if barcode is not None and (new_barcode or not (new_order and new_mate)):
          assert len(duplex) <= 2, duplex.keys()
          results.append(pool.apply_async(process_duplex, args=[duplex, barcode]+static_args))
          results = process_results(results, max_results, filehandles, stats)
          stats['duplexes'] += 1
          duplex = collections.OrderedDict()
        barcode = this_barcode
        order = this_order
        mate = this_mate
        family = []
      read = {'name': name, 'seq':seq, 'qual':qual}
      family.append(read)
      all_reads += 1
    # Process the last family.
    duplex[(order, mate)] = family
    assert len(duplex) <= 2, duplex.keys()
    results.append(pool.apply_async(process_duplex, args=[duplex, barcode]+static_args))
    results = process_results(results, max_results, filehandles, stats)
    stats['duplexes'] += 1

    # Retrieve the remaining results.
    logging.info('Flushing remaining results from worker processes..')
    results = process_results(results, 0, filehandles, stats)

  finally:
    # If the root process encounters an exception and doesn't tell the workers to stop, it will
    # hang forever.
    pool.close()
    pool.join()
    # Close all open filehandles.
    if args.infile is not sys.stdin:
      args.infile.close()
    for fh_group in filehandles.values():
      for fh in fh_group:
        if fh:
          fh.close()

  end_time = time.time()
  run_time = int(end_time - start_time)

  # Final stats on the run.
  logging.info('Processed {} reads and {} duplexes in {} seconds.'
               .format(all_reads, stats['runs'], run_time))
  if stats['reads'] > 0 and stats['runs'] > 0:
    per_read = stats['time'] / stats['reads']
    per_run = stats['time'] / stats['runs']
    logging.info('{:0.3f}s per read, {:0.3f}s per run.'.format(per_read, per_run))

  if args.phone_home:
    run_data = stats.copy()
    run_data['consensus_time'] = run_data['time']
    del run_data['time']
    run_data['processes'] = args.processes
    call.send_data('end', run_time=run_time, run_data=run_data)


def process_duplex(duplex, barcode, incl_sscs=False, min_reads=1, cons_thres=0.5, min_cons_reads=0,
                   qual_thres=' '):
  """Create one duplex consensus sequence from a pair of single-stranded families."""
  # The code in the main loop ensures that "duplex" contains only reads belonging to one final
  # duplex consensus read: ab.1 and ba.2 reads OR ab.2 and ba.1 reads. (Of course, one half might
  # be missing).
  logging.info('Starting duplex {}'.format(barcode))
  start = time.time()
  # Construct consensus sequences.
  sscss, dcs_seq, duplex_mate = make_consensuses(duplex, min_reads, cons_thres, min_cons_reads,
                                                 qual_thres)
  reads_per_strand = [sscs['reads'] for sscs in sscss]
  # Format output.
  dcs_str, sscs_strs = format_outputs(sscss, dcs_seq, barcode, incl_sscs, reads_per_strand)
  # Calculate run statistics.
  elapsed = time.time() - start
  logging.info('{} sec for {} reads.'.format(elapsed, sum(reads_per_strand)))
  if len(sscss) > 0:
    run_stats = {'time':elapsed, 'runs':1, 'reads':sum(reads_per_strand)}
  else:
    run_stats = {'time':0, 'runs':0, 'reads':0}
  return dcs_str, sscs_strs, duplex_mate, run_stats


def make_consensuses(duplex, min_reads, cons_thres, min_cons_reads, qual_thres):
  # Make SSCSs.
  sscss, duplex_mate = make_sscs(duplex, min_reads, cons_thres, min_cons_reads, qual_thres)
  # Make DCS, if possible.
  dcs_seq = None
  if len(sscss) == 2:
    align = swalign.smith_waterman(sscss[0]['consensus'], sscss[1]['consensus'])
    #TODO: log error & return if len(align.target) != len(align.query)
    dcs_seq = consensus.build_consensus_duplex_simple(align.target, align.query)
  return sscss, dcs_seq, duplex_mate


def make_sscs(duplex, min_reads, cons_thres, min_cons_reads, qual_thres):
  """Create single-strand consensus sequences from families of raw reads."""
  sscss = []
  duplex_mate = None
  for (order, mate), family in duplex.items():
    # logging.info('\t{0}.{1}:'.format(order, mate))
    # for read in family:
    #   logging.info('\t\t{name}\t{seq}'.format(**read))
    nreads = len(family)
    if nreads < min_reads:
      continue
    # The mate number for the duplex consensus. It's arbitrary, but all that matters is that the
    # two mates have different numbers. This system ensures that:
    # Mate 1 is from the consensus of ab/1 and ba/2 families, while mate 2 is from ba/1 and ab/2.
    if (order == 'ab' and mate == 1) or (order == 'ba' and mate == 2):
      duplex_mate = 1
    else:
      duplex_mate = 2
    seqs = [read['seq'] for read in family]
    quals = [read['qual'] for read in family]
    consensus_seq = consensus.get_consensus(seqs, quals, cons_thres=cons_thres,
                                            min_reads=min_cons_reads, qual_thres=qual_thres)
    sscss.append({'consensus':consensus_seq, 'order':order, 'mate':mate, 'reads':nreads})
  return sscss, duplex_mate


def format_outputs(sscss, dcs_seq, barcode, incl_sscs, reads_per_strand):
  # SSCS
  sscs_strs = [None, None, None]
  for sscs in sscss:
    sscs_strs[sscs['mate']] = '>{0}.{order} {reads}\n{consensus}\n'.format(barcode, **sscs)
  # DCS
  dcs_str = ''
  if dcs_seq:
    dcs_str = format_consensus(dcs_seq, barcode, reads_per_strand)
  elif incl_sscs and sscss:
    dcs_str = format_consensus(sscss[0]['consensus'], barcode, reads_per_strand)
  return dcs_str, sscs_strs


def format_consensus(seq, barcode, reads_per_strand):
  reads_str = '-'.join(map(str, reads_per_strand))
  return '>{bar} {reads}\n{seq}\n'.format(bar=barcode, reads=reads_str, seq=seq)


def process_results(results, max_results, filehandles, stats):
  if len(results) < max_results:
    return results
  for result in results:
    dcs_str, sscs_strs, duplex_mate, run_stats = result.get()
    for key, value in run_stats.items():
      stats[key] += value
    if dcs_str and filehandles['dcs'][duplex_mate]:
      filehandles['dcs'][duplex_mate].write(dcs_str)
    for mate in (1, 2):
      if sscs_strs[mate] and filehandles['sscs'][mate]:
        filehandles['sscs'][mate].write(sscs_strs[mate])
  return []


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
