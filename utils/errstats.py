#!/usr/bin/env python3
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

ARG_DEFAULTS = {'input':sys.stdin, 'qual_thres':0, 'qual_format':'sanger', 'log':sys.stderr,
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
  parser.add_argument('-H', '--human', action='store_true')
  parser.add_argument('-r', '--all-repeats', action='store_true')
  parser.add_argument('-q', '--qual-thres', type=int)
  parser.add_argument('-F', '--qual-format', choices=('sanger', 'solexa', 'illumina1.3',
                                                      'illumina1.5', 'illumina1.8'))
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

  print_family_errors(args.input, args.qual_thres, args.qual_format, args.human, args.all_repeats)


def print_family_errors(infile, qual_thres, qual_format, human, all_repeats):
  for family in parse_families(infile):
    for order in ('ab', 'ba'):
      for mate in (0, 1):
        seq_align, qual_align = family[order][mate]
        if not seq_align:
          continue
        consensus_seq = get_consensus_wrapper(seq_align, qual_align, qual_thres)
        errors = get_family_errors(consensus_seq, seq_align)
        error_types = group_errors(errors)
        errors_per_seq, repeated_errors, error_repeat_counts = tally_errors(error_types, len(seq_align))
        masked_alignment = mask_alignment(seq_align, errors)
        if human:
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
          print(family['bar'], order, mate, len(seq_align), *error_repeat_counts, sep='\t')
        else:
          print(family['bar'], order, mate, len(seq_align), repeated_errors, *errors_per_seq, sep='\t')


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


def get_consensus_wrapper(seqs, quals, qual_thres):
  """Encode outgoing strings as bytes and decode return value into str."""
  seqs_bytes = [bytes(seq, 'utf8') for seq in seqs]
  quals_bytes = [bytes(qual, 'utf8') for qual in quals]
  cons_bytes = consensus.get_consensus(seqs_bytes, quals_bytes, qual_thres=qual_thres)
  return str(cons_bytes, 'utf8')


def compare(consensus_seq, seq_align, all_repeats=False):
  num_seqs = len(seq_align)
  new_align = [''] * num_seqs
  errors = [0] * num_seqs
  if all_repeats:
    repeat_errors = []
  else:
    repeat_errors = 0
  # Walk along the alignment, one position at a time.
  for bases, cons_base in zip(zip(*seq_align), consensus_seq):
    # Tally how many of each base there are at this position.
    votes = {'A':0, 'C':0, 'G':0, 'T':0, '-':0, 'N':0}
    # Count votes and make new alignment without consensus bases.
    for i, base in enumerate(bases):
      votes[base] += 1
      if base == cons_base:
        new_align[i] += '.'
      else:
        new_align[i] += base
        errors[i] += 1
    # How often did each error occur at this position?
    # (counting each unique non-consensus base as an error)
    for base, vote in votes.items():
      if base != cons_base:
        if all_repeats:
          if vote > 0:
            repeat_errors.append(vote)
        elif vote > 1:
          repeat_errors += 1
    i += 1
  return errors, repeat_errors, new_align


def get_family_errors(consensus_seq, alignment):
  for coord, (cons_base, bases) in enumerate(zip(consensus_seq, zip(*alignment))):
    for seq_num, base in enumerate(bases):
      if base != cons_base:
        yield (seq_num, coord+1, base)


def group_errors(errors):
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


def mask_alignment(seq_alignment, errors):
  masked_alignment = [['.'] * len(seq) for seq in seq_alignment]
  for error in sorted(errors, key=lambda error: error[1]):
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
