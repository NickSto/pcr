#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import argparse
import collections
import parallel_tools
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
USAGE = """$ %(prog)s [options] families.msa.tsv -1 duplexes_1.fa -2 duplexes_2.fa
       $ cat families.msa.tsv | %(prog)s [options] -1 duplexes_1.fa -2 duplexes_2.fa"""
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
              'parameters reported are the number of --processes and the --queue-size. All the '
              'reporting and recording code is available at https://github.com/NickSto/ET.'))
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
  misc.add_argument('-p', '--processes', default=0,
    help=wrap('Number of worker subprocesses to use. If 0, no subprocesses will be started and '
              'everything will be done inside one process. Give "auto" to use as many processes '
              'as there are CPU cores. Default: %(default)s.'))
  misc.add_argument('--queue-size', type=int,
    help=wrap('How long to go accumulating responses from worker subprocesses before dealing '
              'with all of them. Default: {} * the number of worker --processes.'
              .format(parallel_tools.QUEUE_SIZE_MULTIPLIER)))
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
  else:
    call = None

  # Execute as much of the script as possible in a try/except to catch any exception that occurs
  # and report it via ET.phone.
  try:
    # Process and validate arguments.
    if args.queue_size is not None and args.queue_size <= 0:
      fail('Error: --queue-size must be greater than zero.')
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
    # Indexed so we can do filehandles['dcs'][mate].
    filehandles = {
      'dcs': (args.dcs1, args.dcs2),
      'sscs': (args.sscs1, args.sscs2),
    }

    # Open a pool of worker processes.
    stats = {'time':0, 'reads':0, 'runs':0, 'duplexes':0}
    static_kwargs = {
      'min_reads': args.min_reads,
      'cons_thres': args.cons_thres,
      'min_cons_reads': args.min_cons_reads,
      'qual_thres': qual_thres,
    }
    pool = parallel_tools.SyncAsyncPool(process_duplex,
                                        processes=args.processes,
                                        static_kwargs=static_kwargs,
                                        queue_size=args.queue_size,
                                        callback=process_result,
                                        callback_args=[filehandles, stats],
                                       )
    try:
      total_reads = 0
      duplex = collections.OrderedDict()
      family = []
      barcode = None
      order = None
      # Note: mate is a 0-indexed integer ("mate 1" from the input file is mate 0 here).
      mate = None
      for line in args.infile:
        # Allow comments (e.g. for test input files).
        if line.startswith('#'):
          continue
        fields = line.rstrip('\r\n').split('\t')
        if len(fields) != 6:
          continue
        this_barcode, this_order, this_mate, name, seq, qual = fields
        this_mate = int(this_mate)-1
        # If the barcode, order, and mate are the same, we're just continuing the add reads to the
        # current family. Otherwise, store the current family, start a new one, and process the
        # duplex if we're at the end of one.
        new_barcode = this_barcode != barcode
        new_order = this_order != order
        new_mate = this_mate != mate
        if new_barcode or new_order or new_mate:
          if order is not None and mate is not None:
            duplex[(order, mate)] = family
          # If the barcode changed, process the last duplex and start a new one.
          if new_barcode and barcode is not None:
            assert len(duplex) <= 4, duplex.keys()
            pool.compute(duplex, barcode)
            stats['duplexes'] += 1
            duplex = collections.OrderedDict()
          barcode = this_barcode
          order = this_order
          mate = this_mate
          family = []
        read = {'name': name, 'seq':seq, 'qual':qual}
        family.append(read)
        total_reads += 1
      # Process the last family.
      if order is not None and mate is not None:
        duplex[(order, mate)] = family
      assert len(duplex) <= 4, duplex.keys()
      pool.compute(duplex, barcode)
      stats['duplexes'] += 1

      # Retrieve the remaining results.
      logging.info('Flushing remaining results from worker processes..')
      pool.flush()

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

    # Final stats on the run.
    run_time = int(time.time() - start_time)
    logging.info('Processed {} reads and {} duplexes in {} seconds.'
                 .format(total_reads, stats['runs'], run_time))
    if stats['reads'] > 0 and stats['runs'] > 0:
      per_read = stats['time'] / stats['reads']
      per_run = stats['time'] / stats['runs']
      logging.info('{:0.3f}s per read, {:0.3f}s per run.'.format(per_read, per_run))

  except (Exception, KeyboardInterrupt) as exception:
    if args.phone_home and call:
      exception_data = getattr(exception, 'child_context', parallel_tools.get_exception_data())
      run_time = int(time.time() - start_time)
      try:
        run_data = get_run_data(stats, pool)
      except (Exception, UnboundLocalError):
        run_data = {}
      run_data['failed'] = True
      run_data['exception'] = exception_data
      call.send_data('end', run_time=run_time, run_data=run_data)
      logging.critical(parallel_tools.format_traceback(exception_data))
      raise exception
    else:
      raise

  if args.phone_home and call:
    run_data = get_run_data(stats, pool)
    call.send_data('end', run_time=run_time, run_data=run_data)


def get_run_data(stats, pool):
  run_data = stats.copy()
  run_data['consensus_time'] = run_data['time']
  del run_data['time']
  run_data['processes'] = pool.processes
  run_data['queue_size'] = pool.queue_size
  return run_data


def process_duplex(duplex, barcode, min_reads=3, cons_thres=0.5, min_cons_reads=0, qual_thres=' '):
  """Create duplex consensus sequences for the reads from one barcode."""
  # The code in the main loop used to ensure that "duplex" contains only reads belonging to one final
  # duplex consensus read: ab.1 and ba.2 reads OR ab.2 and ba.1 reads. (Of course, one half might
  # be missing).
  logging.info('Starting duplex {}'.format(barcode))
  start = time.time()
  # Construct consensus sequences.
  try:
    sscss = make_sscss(duplex, min_reads, cons_thres, min_cons_reads, qual_thres)
    dcss = make_dcss(sscss)
  except AssertionError:
    logging.exception('While processing duplex {}:'.format(barcode))
    raise
  # Format output.
  dcs_strs, sscs_strs = format_outputs(dcss, sscss, barcode)
  # Calculate run statistics.
  elapsed = time.time() - start
  total_reads = sum([len(family) for family in duplex.values()])
  logging.debug('{} sec for {} reads.'.format(elapsed, total_reads))
  if len(sscss) > 0:
    run_stats = {'time':elapsed, 'runs':1, 'reads':total_reads}
  else:
    run_stats = {'time':0, 'runs':0, 'reads':0}
  return dcs_strs, sscs_strs, run_stats


def make_sscss(duplex, min_reads, cons_thres, min_cons_reads, qual_thres):
  """Create single-strand consensus sequences from families of raw reads."""
  sscss = {}
  for (order, mate), family in duplex.items():
    # logging.info('\t{0}.{1}:'.format(order, mate))
    # for read in family:
    #   logging.info('\t\t{name}\t{seq}'.format(**read))
    if len(family) < min_reads:
      logging.debug('\tnot enough reads ({} < {})'.format(len(family), min_reads))
      continue
    sscs = make_sscs(family, order, mate, qual_thres, cons_thres, min_cons_reads)
    sscss[(order, mate)] = sscs
  return sscss


def make_sscs(family, order, mate, qual_thres, cons_thres, min_cons_reads):
  seqs = [read['seq'] for read in family]
  quals = [read['qual'] for read in family]
  consensus_seq = consensus.get_consensus(seqs,
                                          quals,
                                          cons_thres=cons_thres,
                                          min_reads=min_cons_reads,
                                          qual_thres=qual_thres
                                         )
  return {'seq':consensus_seq, 'order':order, 'mate':mate, 'nreads':len(family)}


def make_dcss(sscss):
  # ordermates is the mapping between the duplex consensus mate number and the order/mates of the
  # SSCSs it's composed of. It's arbitrary but consistent, to make sure the duplex consensuses have
  # different mate numbers, and they're the same from run to run.
  ordermates = {
    0: (('ab', 0), ('ba', 1)),
    1: (('ab', 1), ('ba', 0)),
  }
  # Get the consensus of each pair of SSCSs.
  # dcss is indexed by (0-based) mate.
  dcss = []
  for duplex_mate in 0, 1:
    # Gather the pair of reads for this duplex consensus.
    sscs_pair = []
    for order, mate in ordermates[duplex_mate]:
      sscs = sscss.get((order, mate))
      if sscs:
        sscs_pair.append(sscs)
    if len(sscs_pair) < 2:
      # If we didn't find two SSCSs for this duplex mate, we can't make a complete pair of duplex
      # consensus sequences.
      break
    align = swalign.smith_waterman(sscs_pair[0]['seq'], sscs_pair[1]['seq'])
    if len(align.target) != len(align.query):
      message = '{} != {}:\n'.format(len(align.target), len(align.query))
      message += '\n'.join([repr(sscs) for sscs in sscs_pair])
      raise AssertionError(message)
    seq = consensus.build_consensus_duplex_simple(align.target, align.query)
    reads_per_strand = [sscs['nreads'] for sscs in sscs_pair]
    dcss.append({'seq':seq, 'nreads':reads_per_strand})
  assert len(dcss) == 0 or len(dcss) == 2, len(dcss)
  return dcss


def format_outputs(dcss, sscss, barcode):
  """Format the consensus sequences into FASTA-formatted strings ready for printing.
  sscs_strs is structured so that sscs_strs[order][mate] is the FASTA-formatted output string for
  one SSCS (including ending newline)."""
  # SSCS
  sscs_strs = {}
  for order in 'ab', 'ba':
    sscs_str_pair = []
    for mate in 0, 1:
      sscs = sscss.get((order, mate))
      if sscs:
        sscs_str_pair.append('>{bar}.{order} {nreads}\n{seq}\n'.format(bar=barcode, **sscs))
    if len(sscs_str_pair) == 2:
      sscs_strs[order] = sscs_str_pair
  # DCS
  dcs_strs = []
  if dcss:
    assert len(dcss) == 2, (barcode, len(dcss))
    for duplex_mate in 0, 1:
      dcs = dcss[duplex_mate]
      nreads_str = '-'.join([str(nreads) for nreads in dcs['nreads']])
      dcs_str = '>{bar} {nreads}\n{seq}\n'.format(bar=barcode, nreads=nreads_str, seq=dcs['seq'])
      dcs_strs.append(dcs_str)
  return dcs_strs, sscs_strs


def process_result(result, filehandles, stats):
  dcs_strs, sscs_strs, run_stats = result
  # Stats
  for key, value in run_stats.items():
    stats[key] += value
  # SSCSs
  for order, sscs_pair in sscs_strs.items():
    if sscs_pair:
      for mate in 0, 1:
        sscs_fh = filehandles['sscs'][mate]
        if sscs_fh:
          sscs_fh.write(sscs_pair[mate])
  # DCSs
  if dcs_strs:
    for duplex_mate in 0, 1:
      dcs_fh = filehandles['dcs'][duplex_mate]
      if dcs_fh:
        dcs_fh.write(dcs_strs[duplex_mate])


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
