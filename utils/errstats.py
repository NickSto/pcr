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
import collections
# sys.path hack to access lib package in root directory.
sys.path.insert(1, os.path.dirname(sys.path[0]))
sys.path.insert(2, os.path.join(sys.path[1], 'lib'))
# sys.path hack to allow overriding installed pyBamParser.
if os.environ.get('PYTHONPATH'):
  sys.path.insert(1, os.environ.get('PYTHONPATH'))
import pyBamParser.bam
from lib import simplewrap
import consensus

DESCRIPTION = """Tally statistics on errors in reads, compared to their (single-stranded) \
consensus sequences. Output is one tab-delimited line per single-read alignment (one mate within \
one strand (order) within one family (barcode)).
A "unique error" is a class of error defined by its reference coordinate and the erroneous base.
A single unique error may occur several times in the same family, if it happened on multiple reads.
The default columns are:
1. barcode
2. order
3. mate
4. number of reads
5. number of unique errors that were observed in more than one read
6-end. number of errors in each read
With --all-repeats, these columns are changed:
5-end. count of how many times each unique error was observed in the reads
Format of --overlap-stats is tab-delimited statistics on each mate:
1. barcode
2. order
3. mate
4. "True"/"False": did we find this read's opposite mate in the alignment?
5. length of the overlap region
6. length of the non-overlap region (in this mate)
7. number of unique errors in the overlap region
8. number of unique errors outside the overlap, but aligned to the reference
9. number of unique errors with no reference coordinate
10. number of unique errors that appeared on both mates in the pair (duplicates)"""


def make_argparser():

  # Need to use argparse.RawDescriptionHelpFormatter to preserve formatting in the
  # description of columns in the tsv output. But to still accommodate different
  # terminal widths, dynamic wrapping with simplewrap will be necessary.
  wrap = simplewrap.Wrapper().wrap
  parser = argparse.ArgumentParser(description=wrap(DESCRIPTION),
                                   formatter_class=argparse.RawDescriptionHelpFormatter)

  parser.add_argument('input', metavar='families.msa.tsv', nargs='?', type=argparse.FileType('r'),
    default=sys.stdin,
    help='Aligned families (output of align_families.py). Omit to read from stdin.')
  parser.add_argument('-a', '--alignment', action='store_true',
    help='Print the full alignment, with consensus bases masked (to highlight errors).')
  parser.add_argument('-R', '--all-repeats', action='store_true',
    help='Output the full count of how many times each error recurred in each single-strand '
         'alignment.')
  parser.add_argument('-r', '--min-reads', type=int, default=1,
    help='Minimum number of reads to form a consensus (and thus get any statistics). '
         'Default: %(default)s')
  parser.add_argument('-q', '--qual-thres', type=int, default=0,
    help='PHRED quality score threshold for consensus making. NOTE: This should be the same as was '
         'used for producing the reads in the bam file, if provided! Default: %(default)s')
  parser.add_argument('-Q', '--qual-errors', action='store_true',
    help='Don\'t count errors with quality scores below the --qual-thres in the error counts.')
  parser.add_argument('-d', '--dedup', action='store_true',
    help='Figure out whether there is overlap between mates in read pairs and deduplicate errors '
         'that appear twice because of it. Requires --bam.')
  parser.add_argument('-b', '--bam',
    help='The final duplex consensus reads, aligned to a reference. Used to find overlaps.')
  parser.add_argument('-s', '--seed', type=int, default=0,
    help='The random seed. Used to choose which error to keep when deduplicating errors in '
         'overlaps. Default: %(default)s')
  parser.add_argument('-o', '--overlap-stats', type=argparse.FileType('w'),
    help='Write statistics on overlaps and errors in overlaps to this file. Warning: will '
         'overwrite any existing file.')
  parser.add_argument('-L', '--dedup-log', type=argparse.FileType('w'),
    help='Log overlap error deduplication to this file. Warning: Will overwrite any existing file.')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-S', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.ERROR)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  if args.qual_errors:
    error_qual_thres = args.qual_thres
  else:
    error_qual_thres = 0

  logging.info('Calculating consensus sequences and counting errors..')
  family_stats = {}
  for family in parse_families(args.input):
    barcode = family['bar']
    num_seqs = get_family_stat(family, get_num_seqs)
    consensi = get_family_stat(family, get_consensus, args.qual_thres)
    errors = get_family_stat(family, get_family_errors, consensi, error_qual_thres)
    overlap = get_family_stat(family, lambda a, b, c, d: collections.defaultdict(int))
    if args.dedup:
      family_stats[barcode] = collate_stats(consensi, errors, num_seqs, overlap)
    else:
      print_errors(errors, barcode, num_seqs, args.all_repeats, args.min_reads, args.alignment,
                   consensi, family)

  if args.dedup:
    logging.info('Deduplicating errors in overlaps..')
    dedup_all_errors(args.bam, family_stats, args.dedup_log)
    for barcode in family_stats:
      consensi, errors, num_seqs, overlap = uncollate_stats(family_stats[barcode])
      print_errors(errors, barcode, num_seqs, args.all_repeats, args.min_reads)
      print_overlap_stats(args.overlap_stats, overlap, barcode, num_seqs, args.min_reads)


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
    qual_thres_byte = qual_thres+32
  else:
    seqs_bytes = seq_align
    quals_bytes = qual_align
    qual_thres_byte = chr(qual_thres+32)
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


def get_family_errors(seq_align, qual_align, order, mate, consensi, qual_thres):
  if not (seq_align and qual_align):
    return None
  consensus_seq = consensi[order][mate]
  errors = get_alignment_errors(consensus_seq, seq_align, qual_align, qual_thres)
  error_types = group_errors(errors)
  return list(error_types)


def collate_stats(consensi, errors, num_seqs, overlap):
  stats = {'ab':[None, None], 'ba':[None, None]}
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      stats[order][mate] = {
        'consensus': consensi[order][mate],
        'errors': errors[order][mate],
        'num_seqs': num_seqs[order][mate],
        'overlap': overlap[order][mate],
      }
  return stats


def uncollate_stats(stats):
  consensi = {'ab':[None, None], 'ba':[None, None]}
  errors = {'ab':[None, None], 'ba':[None, None]}
  num_seqs = {'ab':[None, None], 'ba':[None, None]}
  overlap = {'ab':[None, None], 'ba':[None, None]}
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      consensi[order][mate] = stats[order][mate]['consensus']
      errors[order][mate] = stats[order][mate]['errors']
      num_seqs[order][mate] = stats[order][mate]['num_seqs']
      overlap[order][mate] = stats[order][mate]['overlap']
  return consensi, errors, num_seqs, overlap


def print_errors(family_errors, barcode, num_seqs, all_repeats, min_reads=1, print_alignment=False,
                 consensi=None, family=None):
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      num_seq = num_seqs[order][mate]
      if num_seq < min_reads:
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


def print_overlap_stats(stats_fh, family_stats, barcode, num_seqs, min_reads):
  if not stats_fh:
    return
  for order in ('ab', 'ba'):
    for mate in (0, 1):
      num_seq = num_seqs[order][mate]
      if num_seq < min_reads:
        continue
      stats = family_stats[order][mate]
      columns = [barcode, order, mate]
      columns.append(stats['paired'])
      columns.append(stats['overlap_len'])
      columns.append(stats['non_overlap_len'])
      columns.append(stats['overlap_errors'])
      non_overlap_errors = stats['total_errors'] - stats['overlap_errors'] - stats['nonref_errors']
      columns.append(non_overlap_errors)
      columns.append(stats['nonref_errors'])
      columns.append(stats['duplicates'])
      stats_fh.write('\t'.join(map(str, columns))+'\n')


def get_alignment_errors(consensus_seq, seq_align, qual_align, qual_thres):
  qual_thres_char = chr(qual_thres+32)
  errors = []
  for coord, (cons_base, bases, quals) in enumerate(zip(consensus_seq, zip(*seq_align), zip(*qual_align))):
    for seq_num, (base, qual) in enumerate(zip(bases, quals)):
      if base != cons_base and qual >= qual_thres_char:
        errors.append((seq_num, coord+1, base))
  return errors


def group_errors(errors):
  """Group errors by coordinate and base."""
  last_error = None
  current_types = []
  for error in sorted(errors, key=lambda error: error[1]):
    if last_error is not None and last_error[1] == error[1] and last_error[2] == error[2]:
      current_types.append(error)
    else:
      if current_types:
        yield tuple(current_types)
      current_types = [error]
    last_error = error
  if current_types:
    yield tuple(current_types)


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
    try:
      pair_stats = family_stats[barcode][order]
    except KeyError:
      fail('Read pair found in BAM but not in alignment:\nbar: {}, order: {}'
           .format(barcode, order))
    pair_stats[mate]['overlap']['found'] = True
    pair_stats[mate]['overlap']['paired'] = False
    # Skip if it's a secondary alignment or a supplementary alignment, or if it's not mapped in
    # the proper pair.
    flags = read.get_flag()
    if flags & (256+2048) or not flags & 2:
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
          pair_stats[0]['overlap']['paired'] = True
          pair_stats[1]['overlap']['paired'] = True
          pair[mate] = read
          dedup_pair(pair, pair_stats, dedup_log)
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
  if (pair[0] and not pair[1]) or (pair[1] and not pair[0]):
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
  """We've gathered a pair of reads and the errors in them. Now correlate the data between them."""
  dedup_log and dedup_log.write('{} ({} read pairs)\n'.format(pair[0].get_read_name(),
                                                              pair_stats[0]['num_seqs']))
  edges = get_edges(pair)
  overlap_len, non_overlap_lens = get_overlap_len(edges)
  errors_by_ref_coord, nonref_errors = convert_pair_errors(pair, pair_stats)
  count_errors_by_location(errors_by_ref_coord, nonref_errors, edges, pair_stats)
  new_errors_lists = null_duplicate_errors(errors_by_ref_coord, pair_stats, dedup_log)
  if dedup_log and (nonref_errors[0] or nonref_errors[1]):
    log_nonref_errors(nonref_errors, dedup_log)
  dedup_log and dedup_log.write('\n')
  for mate in (0, 1):
    pair_stats[mate]['overlap']['total_errors'] = len(pair_stats[mate]['errors'])
    new_errors_lists[mate].extend(nonref_errors[mate])
    pair_stats[mate]['errors'] = filter(lambda e: e is not None, new_errors_lists[mate])
    pair_stats[mate]['overlap']['non_overlap_len'] = non_overlap_lens[mate]
    pair_stats[mate]['overlap']['overlap_len'] = overlap_len


def convert_pair_errors(pair, pair_stats):
  """Convert the coordinate of each error type into reference coordinates.
  Returns two data structures:
  1. errors which have a reference coordinate
  - A list of two dicts, one per mate.
    - Each dict is indexed by a 2-tuple: (the reference coordinate, the erroneous base)
      - Each value of the dicts is an error_type, as created by group_errors().
  2. errors with no reference coordinate (failed conversion: read.to_ref_coord() returns None)
  - A list of two lists, one per mate.
    - Each value of each list is an error_type."""
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


def count_errors_by_location(errors_by_ref_coord, nonref_errors, edges, pair_stats):
  for mate in (0, 1):
    for error_type in nonref_errors[mate]:
      pair_stats[mate]['overlap']['nonref_errors'] += 1
    for ref_coord, base in errors_by_ref_coord[mate].keys():
      if (edges[mate]['overlap_start'] and edges[mate]['overlap_end'] and
          edges[mate]['overlap_start'] <= ref_coord <= edges[mate]['overlap_end']):
        pair_stats[mate]['overlap']['overlap_errors'] += 1
      else:
        pair_stats[mate]['overlap']['non_overlap_errors'] += 1


def null_duplicate_errors(errors_by_ref_coord, pair_stats, dedup_log=None):
  """Correlate errors between the two mates of the pair & replace one of each duplicate pair with None.
  Whenever the same error (identified by reference coordinate and base) appears on both reads in the
  pair, randomly choose one and replace it with None.
  Returns a new data structure, a tuple containing two lists, one for each mate.
  Each list is simply a list of error_types (from group_errors()), except some elements which are
  None."""
  new_errors_lists = ([], [])
  errors1 = set(errors_by_ref_coord[0].keys())
  errors2 = set(errors_by_ref_coord[1].keys())
  all_errors = list(errors1.union(errors2))
  for ref_coord, base in sorted(all_errors):
    these_error_types = [None, None]
    for mate in (0, 1):
      error_type = errors_by_ref_coord[mate].get((ref_coord, base))
      these_error_types[mate] = error_type
    if these_error_types[0] and these_error_types[1]:
      # The same error occurred at the same position in both reads. Keep only one of them.
      mate = random.choice((0, 1))
      these_error_types[mate] = None
      pair_stats[mate]['overlap']['duplicates'] += 1
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


def get_edges(pair):
  edges = [{}, {}]
  for mate in (0, 1):
    edges[mate]['start'] = pair[mate].get_position()
    edges[mate]['end'] = pair[mate].get_end_position()
  overlap = {}
  overlap['start'] = max(edges[0]['start'], edges[1]['start'])
  overlap['end'] = min(edges[0]['end'], edges[1]['end'])
  for mate in (0, 1):
    for edge in ('start', 'end'):
      if edges[mate]['start'] < overlap[edge] < edges[mate]['end']:
        edges[mate]['overlap_'+edge] = overlap[edge]
      else:
        edges[mate]['overlap_'+edge] = None
  return edges


def get_overlap_len(edges):
  """Get the total number of reference base pairs that the two alignments span, as well as the
  number of base pairs in the overlap between the two alignments."""
  overlap_start = edges[0]['overlap_start'] or edges[1]['overlap_start']
  overlap_end = edges[0]['overlap_end'] or edges[1]['overlap_end']
  if overlap_start is None or overlap_end is None:
    overlap_len = 0
  else:
    overlap_len = overlap_end - overlap_start + 1
  non_overlap_lens = [None, None]
  for mate in (0, 1):
    non_overlap_lens[mate] = edges[mate]['end'] - edges[mate]['start'] + 1 - overlap_len
  return overlap_len, non_overlap_lens


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
