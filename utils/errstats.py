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
        consensus, errors, repeat_errors, new_align = compare(seq_align,
                                                              qual_align,
                                                              thres=qual_thres,
                                                              qual_format=qual_format,
                                                              all_repeats=all_repeats,
                                                              make_new_align=human)
        if human:
          for seq, seq_errors in zip(new_align, errors):
            print('{} errors: {}'.format(seq, seq_errors))
          if all_repeats:
            print('{} errors: {}, repeat errors: {}\n'.format(consensus, sum(errors),
                                                              ', '.join(map(str, repeat_errors))))
          else:
            print('{} errors: {}, repeat errors: {}\n'.format(consensus, sum(errors), repeat_errors))
        elif all_repeats:
          print(family['bar'], order, mate, len(seq_align), *repeat_errors, sep='\t')
        else:
          print(family['bar'], order, mate, len(seq_align), repeat_errors, *errors, sep='\t')


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


def compare(seq_align, qual_align=None, thres=None, qual_format='sanger', all_repeats=False,
            make_new_align=False):
  consensus = ''
  thres_char = phred_to_qual_char(thres, qual_format)
  num_seqs = len(seq_align)
  if qual_align is None:
    qual_align = [''] * num_seqs
  new_align = [''] * num_seqs
  errors = [0] * num_seqs
  last_qual = chr(126)
  if all_repeats:
    repeat_errors = []
  else:
    repeat_errors = 0
  for bases, quals in zip(zip(*seq_align), zip(*qual_align)):
    # Tally how many of each base there are at this position, and find the "winner" (the consensus).
    votes = {'A':0, 'C':0, 'G':0, 'T':0, '-':0, 'N':0}
    cons = 'N'
    max_vote = 0
    for base, qual in zip(bases, quals):
      #TODO: Get gap quality score using algorithm in consensus.c.
      if qual == ' ':
        qual = last_qual
      last_qual = qual
      if qual < thres_char:
        continue
      vote = votes[base]
      vote += 1
      if vote == max_vote:
        cons = 'N'
      elif vote > max_vote:
        max_vote = vote
        cons = base
      votes[base] = vote
    # Make new alignment without consensus bases.
    last_qual = chr(126)
    for i, (base, qual) in enumerate(zip(bases, quals)):
      if qual == ' ':
        qual = last_qual
      last_qual = qual
      if qual < thres_char:
        new_align[i] += ' '
      elif base == cons:
        new_align[i] += '.'
      else:
        new_align[i] += base
        errors[i] += 1
    # How often did each error occur at this position?
    # (counting each unique non-consensus base as an error)
    for base, vote in votes.items():
      if base != cons:
        if all_repeats:
          if vote > 0:
            repeat_errors.append(vote)
        elif vote > 1:
          repeat_errors += 1
    consensus += cons
    i += 1
  return consensus, errors, repeat_errors, new_align


def phred_to_qual_char(phred, qual_format):
  # Based on the trusty "standard": https://en.wikipedia.org/wiki/Fastq#Encoding
  if qual_format == 'solexa' or qual_format == 'illumina1.3' or qual_format == 'illumina1.5':
    return chr(phred+64)
  if qual_format == 'sanger' or qual_format == 'illumina1.8':
    return chr(phred+33)
  else:
    raise ValueError('Unrecognized qual_format "{}"'.format(qual_format))


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
