#!/usr/bin/env python
from __future__ import division
import os
import sys
import time
import logging
import tempfile
import argparse
import subprocess
import collections
import distutils.spawn
import parallel_tools
import seqtools
import shims
# There can be problems with the submodules, but none are essential.
# Try to load these modules, but if there's a problem, load a harmless dummy and continue.
simplewrap = shims.get_module_or_shim('utillib.simplewrap')
version = shims.get_module_or_shim('utillib.version')
phone = shims.get_module_or_shim('ET.phone')

#TODO: Warn if it looks like the two input FASTQ files are the same (i.e. the _1 file was given
#      twice). Can tell by whether the alpha and beta (first and last 12bp) portions of the barcodes
#      are always identical. This would be a good thing to warn about, since it's an easy mistake
#      to make, but it's not obvious that it happened. The pipeline won't fail, but will just
#      produce pretty weird results.

USAGE = """$ %(prog)s [options] families.tsv > families.msa.tsv
       $ cat families.tsv | %(prog)s [options] > families.msa.tsv"""
DESCRIPTION = """Read in sorted FASTQ data and do multiple sequence alignments of each family."""

def make_argparser():

  wrapper = simplewrap.Wrapper()
  wrap = wrapper.wrap
  parser = argparse.ArgumentParser(usage=USAGE, description=wrap(DESCRIPTION),
                                   formatter_class=argparse.RawTextHelpFormatter)

  wrapper.width = wrapper.width - 24
  parser.add_argument('infile', metavar='read-families.tsv', nargs='?', default=sys.stdin,
                      type=argparse.FileType('r'),
    help=wrap('The input reads, sorted into families. One line per read pair, 8 tab-delimited '
              'columns:\n'
              '1. canonical barcode\n'
              '2. barcode order ("ab" for alpha+beta, "ba" for beta-alpha)\n'
              '3. read 1 name\n'
              '4. read 1 sequence\n'
              '5. read 1 quality scores\n'
              '6. read 2 name\n'
              '7. read 2 sequence\n'
              '8. read 2 quality scores'))
  parser.add_argument('-a', '--aligner', choices=('mafft', 'kalign'), default='kalign',
    help=wrap('The multiple sequence aligner to use. Default: %(default)s'))
  parser.add_argument('-p', '--processes', default=0,
    help=wrap('Number of worker subprocesses to use. If 0, no subprocesses will be started and '
              'everything will be done inside one process. Give "auto" to use as many processes '
              'as there are CPU cores. Default: %(default)s.'))
  parser.add_argument('--queue-size', type=int,
    help=wrap('How long to go accumulating responses from worker subprocesses before dealing '
              'with all of them. Default: {} * the number of worker --processes.'
              .format(parallel_tools.QUEUE_SIZE_MULTIPLIER)))
  parser.add_argument('--phone-home', action='store_true',
    help=wrap('Report helpful usage data to the developer, to better understand the use cases and '
              'performance of the tool. The only data which will be recorded is the name and '
              'version of the tool, the size of the input data, the time taken to process it, and '
              'the IP address of the machine running it. No filenames are sent, and the only '
              'parameters reported are --aligner, --processes, and --queue-size, which are '
              'necessary to evaluate performance. All the reporting and recording code is '
              'available at https://github.com/NickSto/ET.'))
  parser.add_argument('--galaxy', dest='platform', action='store_const', const='galaxy',
    help=wrap('Tell the script it\'s running on Galaxy. Currently this only affects data reported '
              'when phoning home.'))
  parser.add_argument('--test', action='store_true',
    help=wrap('If reporting usage data, mark this as a test run.'))
  parser.add_argument('--version', action='version', version=str(version.get_version()),
    help=wrap('Print the version number and exit.'))
  parser.add_argument('-L', '--log-file', type=argparse.FileType('w'), default=sys.stderr,
    help=wrap('Print log messages to this file instead of to stderr. NOTE: Will overwrite the file.'))
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
                      default=logging.WARNING)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log_file, level=args.volume, format='%(message)s')
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

  # Execute as much of the script as possible in a try/except to catch any exception that occurs
  # and report it via ET.phone.
  try:
    if args.queue_size is not None and args.queue_size <= 0:
      fail('Error: --queue-size must be greater than zero.')

    # If we're using mafft, check that we can execute it.
    if args.aligner == 'mafft' and not distutils.spawn.find_executable('mafft'):
      fail('Error: Could not find "mafft" command on $PATH.')

    # Open a pool of worker processes.
    stats = {'duplexes':0, 'time':0, 'pairs':0, 'runs':0, 'failures':0, 'aligned_pairs':0}
    pool = parallel_tools.SyncAsyncPool(process_duplex,
                                        processes=args.processes,
                                        static_kwargs={'aligner':args.aligner},
                                        queue_size=args.queue_size,
                                        callback=process_result,
                                        callback_args=[stats],
                                       )
    """Now the main loop.
    This processes whole duplexes (pairs of strands) at a time for a future option to align the
    whole duplex at a time.
    duplex data structure:
    duplex = {
      'ab': [
        {'name1': 'read_name1a',
         'seq1':  'GATT-ACA',
         'qual1': 'sc!0 /J*',
         'name2': 'read_name1b',
         'seq2':  'ACTGACTA',
         'qual2': '34I&SDF)'
        },
        {'name1': 'read_name2a',
         ...
        },
        ...
      ],
      'ba': [
        ...
      ]
    }
    e.g.:
    seq = duplex[order][pair_num]['seq1']"""

    try:
      duplex = collections.OrderedDict()
      family = []
      barcode = None
      order = None
      for line in args.infile:
        fields = line.rstrip('\r\n').split('\t')
        if len(fields) != 8:
          continue
        (this_barcode, this_order, name1, seq1, qual1, name2, seq2, qual2) = fields
        # If the barcode or order has changed, we're in a new family.
        # Process the reads we've previously gathered as one family and start a new family.
        if this_barcode != barcode or this_order != order:
          duplex[order] = family
          # If the barcode is different, we're at the end of the whole duplex. Process the it and start
          # a new one. If the barcode is the same, we're in the same duplex, but we've switched strands.
          if this_barcode != barcode:
            # logging.debug('processing {}: {} orders ({})'.format(barcode, len(duplex),
            #               '/'.join([str(len(duplex[o])) for o in duplex])))
            pool.compute(duplex, barcode)
            stats['duplexes'] += 1
            duplex = collections.OrderedDict()
          barcode = this_barcode
          order = this_order
          family = []
        pair = {'name1': name1, 'seq1':seq1, 'qual1':qual1, 'name2':name2, 'seq2':seq2, 'qual2':qual2}
        family.append(pair)
        stats['pairs'] += 1
      # Process the last family.
      duplex[order] = family
      # logging.debug('processing {}: {} orders ({}) [last]'.format(barcode, len(duplex),
      #               '/'.join([str(len(duplex[o])) for o in duplex])))
      pool.compute(duplex, barcode)
      stats['duplexes'] += 1

      # Retrieve the remaining results.
      logging.info('Flushing remaining results from worker processes..')
      pool.flush()

    finally:
      # If an exception occurs in the parent without stopping the child processes, this will hang.
      # Make sure to kill the children in all cases.
      pool.close()
      pool.join()

    if args.infile is not sys.stdin:
      args.infile.close()

    # Final stats on the run.
    run_time = int(time.time() - start_time)
    logging.error('Processed {pairs} read pairs in {duplexes} duplexes, with {failures} alignment '
                  'failures.'.format(**stats))
    if stats['aligned_pairs'] > 0 and stats['runs'] > 0:
      per_pair = stats['time'] / stats['aligned_pairs']
      per_run = stats['time'] / stats['runs']
      logging.error('{:0.3f}s per pair, {:0.3f}s per run.'.format(per_pair, per_run))
    logging.error('in {}s total time.'.format(run_time))

  except (Exception, KeyboardInterrupt) as exception:
    if args.phone_home and call:
      exception_data = getattr(exception, 'child_context', parallel_tools.get_exception_data())
      run_time = int(time.time() - start_time)
      try:
        run_data = get_run_data(stats, pool, args.aligner)
      except (Exception, UnboundLocalError):
        run_data = {}
      run_data['failed'] = True
      run_data['exception'] = exception_data
      call.send_data('end', run_time=run_time, run_data=run_data)
      logging.critical(parallel_tools.format_traceback(exception_data))
      raise exception
    else:
      raise

  if args.phone_home:
    run_data = get_run_data(stats, pool, args.aligner)
    call.send_data('end', run_time=run_time, run_data=run_data)


def get_run_data(stats, pool, aligner):
  run_data = stats.copy()
  run_data['align_time'] = run_data['time']
  del run_data['time']
  run_data['processes'] = pool.processes
  run_data['queue_size'] = pool.queue_size
  run_data['aligner'] = aligner
  return run_data


def process_duplex(duplex, barcode, aligner='mafft'):
  output = ''
  run_stats = {'time':0, 'runs':0, 'aligned_pairs':0, 'failures':0}
  orders = duplex.keys()
  if len(duplex) == 0 or None in duplex:
    return '', {}
  elif len(duplex) == 1:
    # If there's only one strand in the duplex, just process the first mate, then the second.
    combos = ((1, orders[0]), (2, orders[0]))
  elif len(duplex) == 2:
    # If there's two strands, process in a criss-cross order:
    # strand1/mate1, strand2/mate2, strand1/mate2, strand2/mate1
    combos = ((1, orders[0]), (2, orders[1]), (2, orders[0]), (1, orders[1]))
  else:
    raise AssertionError('More than 2 orders in duplex {}: {}'.format(barcode, orders))
  for mate, order in combos:
    family = duplex[order]
    start = time.time()
    try:
      alignment = align_family(family, mate, aligner=aligner)
    except AssertionError as error:
      logging.exception('While processing duplex {}, order {}, mate {}:'.format(barcode, order, mate))
      raise
    except (OSError, subprocess.CalledProcessError) as error:
      logging.warning('{} on family {}, order {}, mate {}:\n{}'
                      .format(type(error).__name__, barcode, order, mate, error))
      alignment = None
    # Compile statistics.
    elapsed = time.time() - start
    pairs = len(family)
    logging.info('{} sec for {} read pairs.'.format(elapsed, pairs))
    if pairs > 1:
      run_stats['time'] += elapsed
      run_stats['runs'] += 1
      run_stats['aligned_pairs'] += pairs
    if alignment is None:
      logging.warning('Error aligning family {}/{} (read {}).'.format(barcode, order, mate))
      run_stats['failures'] += 1
    else:
      output += format_msa(alignment, barcode, order, mate)
  return output, run_stats


def align_family(family, mate, aligner='mafft'):
  """Do a multiple sequence alignment of the reads in a family and their quality scores."""
  mate = str(mate)
  assert mate == '1' or mate == '2'
  if len(family) == 0:
    return None
  elif len(family) == 1:
    # If there's only one read pair, there's no alignment to be done (and MAFFT won't accept it).
    aligned_seqs = [family[0]['seq'+mate]]
  else:
    # Do the multiple sequence alignment.
    aligned_seqs = make_msa(family, mate, aligner=aligner)
  # Transfer the alignment to the quality scores.
  ## Get a list of all quality scores in the family for this mate.
  quals_raw = [pair['qual'+mate] for pair in family]
  qual_alignment = seqtools.transfer_gaps_multi(quals_raw, aligned_seqs, gap_char_out=' ')
  # Package them up in the output data structure.
  alignment = []
  for pair, aligned_seq, aligned_qual in zip(family, aligned_seqs, qual_alignment):
    alignment.append({'name':pair['name'+mate], 'seq':aligned_seq, 'qual':aligned_qual})
  return alignment


def make_msa(family, mate, aligner='mafft'):
  if aligner == 'mafft':
    return make_msa_mafft(family, mate)
  elif aligner == 'kalign':
    return make_msa_kalign(family, mate)


def make_msa_kalign(family, mate):
  logging.info('Aligning with kalign.')
  from kalign import kalign
  seqs = [pair['seq'+mate] for pair in family]
  aln_struct = kalign.align(seqs)
  return [aln_struct.seqs[i] for i in range(aln_struct.nseqs)]


def make_msa_mafft(family, mate):
  """Perform a multiple sequence alignment on a set of sequences and parse the result.
  Uses MAFFT."""
  logging.info('Aligning with mafft.')
  #TODO: Replace with tempfile.mkstemp()?
  with tempfile.NamedTemporaryFile('w', delete=False, prefix='align.msa.') as family_file:
    for pair in family:
      name = pair['name'+mate]
      seq = pair['seq'+mate]
      family_file.write('>'+name+'\n')
      family_file.write(seq+'\n')
  with open(os.devnull, 'w') as devnull:
    try:
      command = ['mafft', '--nuc', '--quiet', family_file.name]
      output = subprocess.check_output(command, stderr=devnull)
    except (OSError, subprocess.CalledProcessError):
      raise
    finally:
      # Make sure we delete the temporary file.
      os.remove(family_file.name)
  return read_fasta(output)


def read_fasta(fasta):
  """Quick and dirty FASTA parser. Return the sequences and their names.
  Returns a list of sequences.
  Warning: Reads the entire contents of the file into memory at once."""
  sequences = []
  sequence = ''
  for line in fasta.splitlines():
    if line.startswith('>'):
      if sequence:
        sequences.append(sequence.upper())
      sequence = ''
      continue
    sequence += line.strip()
  if sequence:
    sequences.append(sequence.upper())
  return sequences


def format_msa(align, barcode, order, mate, outfile=sys.stdout):
  output = ''
  for sequence in align:
    output += '{bar}\t{order}\t{mate}\t{name}\t{seq}\t{qual}\n'.format(bar=barcode, order=order,
                                                                       mate=mate, **sequence)
  return output


def process_result(result, stats):
  """Process the outcome of a duplex run.
  Print the aligned output and sum the stats from the run with the running totals."""
  output, run_stats = result
  for key, value in run_stats.items():
    stats[key] += value
  if output:
    sys.stdout.write(output)


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
