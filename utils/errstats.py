#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import os
import sys
import errno
import random
import logging
import argparse
import pyBamParser.bam
# sys.path hack to access lib package in root directory.
sys.path.insert(1, os.path.dirname(sys.path[0]))
sys.path.insert(2, os.path.join(sys.path[1], 'lib'))
from lib import simplewrap
import consensus

ARG_DEFAULTS = {'input':sys.stdin, 'qual_thres':0, 'seed':0, 'log':sys.stderr,
                'volume':logging.ERROR}
DESCRIPTION = """Tally statistics on errors in reads, compared to the rest of their (single-\
stranded) families.
Output columns without --all-repeats:
1. barcode
2. order
3. mate
4. number of reads
5. number of errors that occurred more than once
6-end. number of errors in each read.
With --all-repeats:
1. barcode
2. order
3. mate
4. number of reads
5-end. count of how many times each error occurred in the family."""


def make_argparser():

  # Need to use argparse.RawDescriptionHelpFormatter to preserve formatting in the
  # description of columns in the tsv output. But to still accommodate different
  # terminal widths, dynamic wrapping with simplewrap will be necessary.
  wrap = simplewrap.Wrapper().wrap

  parser = argparse.ArgumentParser(description=wrap(DESCRIPTION),
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('input', metavar='families.msa.tsv', nargs='?', type=argparse.FileType('r'),
    help='')
  parser.add_argument('-a', '--alignment', action='store_true',
    help='Print the full alignment, with consensus bases masked (to highlight errors).')
  parser.add_argument('-r', '--all-repeats', action='store_true',
    help='Output the full count of how many times each error recurred in each single-strand '
         'alignment.')
  parser.add_argument('-q', '--qual-thres', type=int,
    help='PHRED quality threshold for consensus making. NOTE: This should be the same as was used '
         'for producing the reads in the bam file, if provided!')
  parser.add_argument('-o', '--overlap', action='store_true',
    help='Figure out whether there is overlap between mates in read pairs and deduplicate errors '
         'that appear twice because of it. Requires --bam.')
  parser.add_argument('-b', '--bam',
    help='The final duplex consensus reads, aligned to a reference. Used to find overlaps.')
  parser.add_argument('-s', '--seed', type=int,
    help='The random seed. Used to choose which error to keep when deduplicating errors in '
         'overlaps. Default: %(default)s')
  parser.add_argument('-d', '--dedup-log', type=argparse.FileType('w'),
    help='Log overlap error deduplication to this file. Warning: Will overwrite any existing file.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  logging.info('Calculating consensus sequences and counting errors..')
  family_stats = {}
  for family in parse_families(args.input):
    barcode = family['bar']
    consensi = get_family_stat(family, get_consensus, args.qual_thres)
    errors = get_family_stat(family, get_family_errors, consensi)
    num_seqs = get_family_stat(family, get_num_seqs)
    if args.overlap:
      family_stats[barcode] = collate_stats(consensi, errors, num_seqs)
    else:
      print_errors(errors, barcode, num_seqs, args.all_repeats, args.alignment, consensi, family)

  if args.overlap:
    logging.info('Deduplicating errors in overlaps..')
    dedup_all_errors(args.bam, family_stats, args.dedup_log)
    for barcode in family_stats:
      consensi, errors, num_seqs, duplicates = uncollate_stats(family_stats[barcode])
      print_errors(errors, barcode, num_seqs, args.all_repeats)


def parse_families(infile):
  """Parse a families.msa.tsv file.
  Yields a data structure for each family:
  family = {
    'bar': barcode,                                     # family 'AAACCGACACAGGACTAGGGATCA'
    'ab': (                                               # order ab
            ([seq1, seq2, seq3], [quals1, quals2, quals3]), # mate 1
            ([seq1, seq2], [quals1, quals2]),               # mate 2
          ),
    'ba': (                                               # order ba
            ([seq1, seq2], [quals1, quals2]),               # mate 1
            ([], []),                                       # mate 2
          )
  }
  That is, each family is a dict with the 'bar' key giving the barcode sequence, and a key for both
  orders ('ab' and 'ba'). The value for each order is a tuple of 2 values, one for each mate. Each
  value in the tuple is itself a 2-tuple containing the aligned bases and quality scores.
  Examples:
  Getting the sequences for mate 1 of order "ab":
  seq_align = family['ab'][0][0]
  Getting the quality scores:
  seq_align = family['ab'][0][1]
  Getting the sequences for mate 2 of order "ba":
  seq_align = family['ba'][1][0]
  """
  last_barcode = None
  family = {'bar':None, 'ab':(([],[]), ([],[])), 'ba':(([],[]), ([],[]))}
  for line in infile:
    fields = line.rstrip('\r\n').split('\t')
    barcode = fields[0]
    order = fields[1]
    mate = int(fields[2])-1
    seq = fields[4]
    quals = fields[5]
    if barcode != last_barcode:
      if last_barcode is not None:
        yield family
      family = {'bar':barcode, 'ab':(([],[]), ([],[])), 'ba':(([],[]), ([],[]))}
      last_barcode = barcode
    family[order][mate][0].append(seq)
    family[order][mate][1].append(quals)
  yield family


def get_family_stat(family, stat_fxn, *args, **kwargs):
  """Common function for going through the four alignments in each family (+/- strand, 1st/2nd mate)
  and calculating a statistic on each.
  Pass in a function which takes the fixed arguments seq_align, qual_align, order, and mate,
  plus whatever custom ones you want after that. The arguments to this function that come after the
  stat function will be passed directly to the stat function.
  Returns a data structure like that from parse_families(), but family[order][mate] == the stat you
  requested (the output of the given function)."""
  stats = {'ab':[None, None], 'ba':[None, None]}
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      seq_align, qual_align = family[order][mate]
      stats[order][mate] = stat_fxn(seq_align, qual_align, order, mate, *args, **kwargs)
  return stats


def get_consensus(seq_align, qual_align, order, mate, qual_thres):
  """Wrapper around consensus.get_consensus().
  When running under Python 3, this encodes strings passed to it as bytes and decodes its return
  value into str."""
  if not (seq_align and qual_align):
    return None
  if sys.version_info.major == 3:
    seqs_bytes = [bytes(seq, 'utf8') for seq in seq_align]
    quals_bytes = [bytes(qual, 'utf8') for qual in qual_align]
    qual_thres_byte = qual_thres
  else:
    seqs_bytes = seq_align
    quals_bytes = qual_align
    qual_thres_byte = chr(qual_thres)
  cons_bytes = consensus.get_consensus(seqs_bytes,
                                       quals_bytes,
                                       qual_thres=qual_thres_byte,
                                       gapped=True)
  if sys.version_info.major == 3:
    cons_seq = str(cons_bytes, 'utf8')
  else:
    cons_seq = cons_bytes
  return cons_seq


def get_num_seqs(seq_align, qual_align, order, mate):
  return len(seq_align)


def get_family_errors(seq_align, qual_align, order, mate, consensi):
  if not (seq_align and qual_align):
    return None
  consensus_seq = consensi[order][mate]
  errors = get_alignment_errors(consensus_seq, seq_align)
  error_types = group_errors(errors)
  return error_types


def collate_stats(consensi, errors, num_seqs):
  stats = {'ab':[None, None], 'ba':[None, None]}
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      stats[order][mate] = {
        'consensus': consensi[order][mate],
        'errors': errors[order][mate],
        'num_seqs': num_seqs[order][mate],
        'duplicates': 0,
      }
  return stats


def uncollate_stats(stats):
  consensi = {'ab':[None, None], 'ba':[None, None]}
  errors = {'ab':[None, None], 'ba':[None, None]}
  num_seqs = {'ab':[None, None], 'ba':[None, None]}
  duplicates = {'ab':[None, None], 'ba':[None, None]}
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      consensi[order][mate] = stats[order][mate]['consensus']
      errors[order][mate] = stats[order][mate]['errors']
      num_seqs[order][mate] = stats[order][mate]['num_seqs']
      duplicates[order][mate] = stats[order][mate]['duplicates']
  return consensi, errors, num_seqs, duplicates


def print_errors(family_errors, barcode, num_seqs, all_repeats, print_alignment=False,
                 consensi=None, family=None):
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      num_seq = num_seqs[order][mate]
      if num_seq == 0:
        continue
      error_types = family_errors[order][mate]
      errors_per_seq, repeated_errors, error_repeat_counts = tally_errors(error_types, num_seq)
      if print_alignment:
        consensus_seq = consensi[order][mate]
        seq_align, qual_align = family[order][mate]
        masked_alignment = mask_alignment(seq_align, error_types)
        for seq, seq_errors in zip(masked_alignment, errors_per_seq):
          print('{} errors: {}'.format(seq, seq_errors))
        if all_repeats:
          print('{} errors: {}, repeat errors: {}\n'.format(consensus_seq,
                                                            sum(errors_per_seq),
                                                            ', '.join(map(str, error_repeat_counts))))
        else:
          print('{} errors: {}, repeat errors: {}\n'.format(consensus_seq,
                                                            sum(errors_per_seq),
                                                            repeated_errors))
      elif all_repeats:
        print(barcode, order, mate, num_seq, *error_repeat_counts, sep='\t')
      else:
        print(barcode, order, mate, num_seq, repeated_errors, *errors_per_seq, sep='\t')


def get_alignment_errors(consensus_seq, alignment):
  errors = []
  for coord, (cons_base, bases) in enumerate(zip(consensus_seq, zip(*alignment))):
    for seq_num, base in enumerate(bases):
      if base != cons_base:
        errors.append((seq_num, coord+1, base))
  return errors


def group_errors(errors):
  """Group errors by coordinate and base."""
  error_types = []
  last_error = None
  current_types = []
  for error in sorted(errors, key=lambda error: error[1]):
    if last_error is not None and last_error[1] == error[1] and last_error[2] == error[2]:
      current_types.append(error)
    else:
      if current_types:
        error_types.append(tuple(current_types))
      current_types = [error]
    last_error = error
  if current_types:
    error_types.append(tuple(current_types))
  return error_types


def tally_errors(error_types, num_seqs):
  errors_per_seq = [0] * num_seqs
  repeated_errors = 0
  error_repeat_counts = []
  for error_type in error_types:
    error_repeat_counts.append(len(error_type))
    if len(error_type) > 1:
      repeated_errors += 1
    for error in error_type:
      errors_per_seq[error[0]] += 1
  return errors_per_seq, repeated_errors, error_repeat_counts


def mask_alignment(seq_alignment, error_types):
  masked_alignment = [['.'] * len(seq) for seq in seq_alignment]
  for error_type in error_types:
    for error in error_type:
      seq_num = error[0]
      coord = error[1]
      base = error[2]
      masked_alignment[seq_num][coord-1] = base
  return [''.join(seq) for seq in masked_alignment]


def dedup_all_errors(bam_path, family_stats, dedup_log):
  pair = [None, None]
  for read in pyBamParser.bam.Reader(bam_path):
    barcode, order, mate = get_read_identifiers(read)
    # Skip if it's a secondary alignment or a supplementary alignment.
    if read.get_flag() & (256+2048):
      continue
    if pair[mate]:
      # We already have this mate for this pair.
      # We must be on a new pair now, and the matching mate for the last one is missing.
      logging.debug('Failed to complete the pair for {}'.format(pair[mate].get_read_name()))
      pair = [None, None]
      pair[mate] = read
    else:
      other_mate = mate ^ 1
      if pair[other_mate]:
        barcode2, order2, mate2 = get_read_identifiers(pair[other_mate])
        if barcode2 == barcode and order2 == order:
          # It's a matching pair.
          pair[mate] = read
          try:
            dedup_pair(pair, family_stats[barcode][order], dedup_log)
          except KeyError:
            fail('Read pair found in BAM but not in alignment:\nbar: {}, order: {}'
                 .format(barcode, order))
          pair = [None, None]
        else:
          # The reads do not match; they're from different pairs.
          # We must be on a new pair now, and the matching mate for the last one is missing.
          logging.debug('Failed to complete the pair for {}.{}'.format(barcode2, order2))
          pair = [None, None]
          pair[mate] = read
      else:
        # The pair is empty ([None, None]).
        # This happens on the first loop, and after a pair has been completed on the previous loop.
        logging.debug('Pair for {} empty ([None, None]).'.format(barcode))
        pair[mate] = read
  if pair[0] and pair[1]:
    try:
      dedup_pair(pair, family_stats[barcode][order], dedup_log)
    except KeyError:
      fail('Read pair found in BAM but not in alignment:\nbar: {}, order: {}'
           .format(barcode, order))
  elif pair[0] or pair[1]:
    read = pair[0] or pair[1]
    logging.debug('Failed to complete the pair for {}'.format(read.get_read_name()))


def get_read_identifiers(read):
  name = read.get_read_name()
  barcode, order = name.split('.')
  flags = read.get_flag()
  if flags & 64:
    mate = 0
  elif flags & 128:
    mate = 1
  else:
    raise ValueError('Neither flag 64 nor 128 are set: {}'.format(read.get_flag()))
  return barcode, order, mate


def dedup_pair(pair, pair_stats, dedup_log=None):
  dedup_log and dedup_log.write('{} ({} read pairs)\n'.format(pair[0].get_read_name(),
                                                              pair_stats[0]['num_seqs']))
  errors_by_ref_coord, nonref_errors = convert_pair_errors(pair, pair_stats)
  new_errors_lists = null_duplicate_errors(errors_by_ref_coord, nonref_errors, pair_stats, dedup_log)
  if dedup_log and (nonref_errors[0] or nonref_errors[1]):
    log_nonref_errors(nonref_errors, dedup_log)
  dedup_log and dedup_log.write('\n')
  for mate in (0, 1):
    new_errors_lists[mate].extend(nonref_errors[mate])
    pair_stats[mate]['errors'] = filter(lambda e: e is not None, new_errors_lists[mate])


def convert_pair_errors(pair, pair_stats):
  errors_by_ref_coord = [{}, {}]
  nonref_errors = [[], []]
  for mate in (0, 1):
    read = pair[mate]
    error_types = pair_stats[mate]['errors']
    for i, error_type in enumerate(error_types):
      error = error_type[0]
      read_coord = error[1]
      base = error[2]
      ref_coord = read.to_ref_coord(read_coord)
      if ref_coord is None:
        nonref_errors[mate].append(error_type)
      else:
        errors_by_ref_coord[mate][(ref_coord, base)] = error_type
  return errors_by_ref_coord, nonref_errors


def null_duplicate_errors(errors_by_ref_coord, nonref_errors, pair_stats, dedup_log=None):
  new_errors_lists = ([], [])
  errors1 = set(errors_by_ref_coord[0].keys())
  errors2 = set(errors_by_ref_coord[1].keys())
  all_errors = list(errors1.union(errors2))
  for ref_coord, base in sorted(all_errors):
    these_error_types = [None, None]
    for mate in (0, 1):
      these_error_types[mate] = errors_by_ref_coord[mate].get((ref_coord, base))
    if these_error_types[0] and these_error_types[1]:
      # The same error occurred at the same position in both reads. Keep only one of them.
      mate = random.choice((0, 1))
      these_error_types[mate] = None
      pair_stats[mate]['duplicates'] += 1
      dedup_log and dedup_log.write('omitting error {} {} from mate {}\n'
                                    .format(ref_coord, base, mate+1))
    dedup_log and dedup_log.write('{:5d} {:1s}:  '.format(ref_coord, base))
    for mate in (0, 1):
      if dedup_log:
        error_type = these_error_types[mate]
        if error_type is None:
          dedup_log.write('             ')
        else:
          dedup_log.write('  {:2d} errors  '.format(len(error_type)))
      new_errors_lists[mate].append(these_error_types[mate])
    dedup_log and dedup_log.write('\n')
  return new_errors_lists


def log_nonref_errors(nonref_errors, dedup_log):
  dedup_log.write('NO REF COORD:\n')
  for mate in (0, 1):
    for error_type in nonref_errors[mate]:
      error = error_type[0]
      dedup_log.write('{:5d} {:1s}:  '.format(error[1], error[2]))
      if mate == 1:
        dedup_log.write('             ')
      dedup_log.write('  {:2d} errors\n'.format(len(error_type)))


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
