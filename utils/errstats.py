#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import os
import sys
import errno
import logging
import argparse
# sys.path hack to access lib package in root directory.
sys.path.insert(1, os.path.dirname(sys.path[0]))
sys.path.insert(2, os.path.join(sys.path[1], 'lib'))
from lib import simplewrap
import consensus

ARG_DEFAULTS = {'input':sys.stdin, 'qual_thres':0, 'log':sys.stderr,
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
  parser.add_argument('-q', '--qual-thres', type=int)
  parser.add_argument('-o', '--overlap', action='store_true',
    help='Figure out whether there is overlap between mates in read pairs and deduplicate errors '
         'that appear twice because of it. Requires --bam.')
  parser.add_argument('-b', '--bam',
    help='The final duplex consensus reads, aligned to a reference. Used to find overlaps.')
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

  if args.overlap:
    raise NotImplementedError

  for family in parse_families(args.input):
    barcode = family['bar']
    consensi = get_family_stat(family, get_consensus, args.qual_thres)
    errors = get_family_stat(family, get_family_errors, consensi)
    num_seqs = get_family_stat(family, get_num_seqs)
    if not args.overlap:
      print_errors(errors, barcode, num_seqs, args.alignment, args.all_repeats, consensi, family)


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


def print_errors(family_errors, barcode, num_seqs, print_alignment, all_repeats, consensi=None, family=None):
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
        print(family['bar'], order, mate, num_seq, *error_repeat_counts, sep='\t')
      else:
        print(family['bar'], order, mate, num_seq, repeated_errors, *errors_per_seq, sep='\t')


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
        error_types.append(current_types)
      current_types = [error]
    last_error = error
  if current_types:
    error_types.append(current_types)
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
