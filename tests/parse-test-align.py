#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import random
import logging
import argparse
assert sys.version_info.major >= 3, 'Python 3 required'

REVCOMP_TABLE = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')

DESCRIPTION = """Generate test files from a human-readable alignment."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.add_argument('alignment', type=argparse.FileType('r'), nargs='?', default=sys.stdin,
    help='The input alignment file. Omit to read from stdin.')
  parser.add_argument('-1', '--fq1', type=argparse.FileType('w'),
    help='Write the first-mate reads to this fastq file. Warning: will overwrite any existing file.')
  parser.add_argument('-2', '--fq2', type=argparse.FileType('w'),
    help='Write the first-mate reads to this fastq file. Warning: will overwrite any existing file.')
  parser.add_argument('-r', '--ref', type=argparse.FileType('w'),
    help='Write the reference sequence to this file. Warning: will overwrite any existing file.')
  parser.add_argument('-b', '--bar-len', type=int, default=12)
  parser.add_argument('-c', '--const-len', type=int, default=5)
  parser.add_argument('-q', '--default-qual', type=int, default=40)
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  const = rand_seq(args.const_len)

  qual_char = chr(args.default_qual+32)

  if args.fq1 and args.fq2:
    fq_files = [args.fq1, args.fq2]
  else:
    fq_files = []

  ref_seq = None
  family_num = 0
  pair_num = 0
  first_mate = None
  last_family_char = None
  for line_raw in args.alignment:
    prefix = line_raw[:3]
    line = line_raw[3:].rstrip()
    if not line:
      continue
    if prefix[0] == 'f':
      ref_seq = line
      if args.ref:
        args.ref.write('>ref\n')
        args.ref.write(ref_seq+'\n')
    elif prefix[0] == 'r':
      assert ref_seq is not None, line_raw
      mate = int(prefix[1])
      if first_mate is None:
        first_mate = mate
      family_char = prefix[2]
      if family_char != last_family_char:
        family_num += 1
        pair_num = 0
      if mate == first_mate:
        pair_num += 1
      raw_seq, pos, direction = get_raw_seq(line)
      barcodes = get_barcodes(family_num, args.bar_len)
      final_seq = substitute_ref_bases(raw_seq, pos, ref_seq)
      if fq_files:
        for line in format_read(final_seq, direction, mate, family_num, pair_num, barcodes, const,
                                qual_char):
          fq_files[mate-1].write(line+'\n')
      last_family_char = family_char


def get_raw_seq(line):
  direction = None
  seq = line.lstrip(' ')
  if seq.endswith('+'):
    direction = 'forward'
    seq = seq.rstrip('+')
    pos = len(line) - len(seq)
  elif seq.startswith('-'):
    direction = 'reverse'
    seq = seq.lstrip('-')
    pos = len(line) - len(seq) + 1
  else:
    fail('A +/- direction is required at the 3\' end of the read.')
  return seq, pos, direction


def substitute_ref_bases(raw_seq, pos, ref_seq):
  """Replace dots in the raw sequence with reference bases."""
  final_seq = ''
  for i, raw_char in enumerate(raw_seq):
    if raw_char == '.':
      final_seq += ref_seq[pos+i-1]
    else:
      final_seq += raw_char
  return final_seq


def get_barcodes(family_num, bar_len):
  random.seed(family_num)
  barcodes = []
  for i in range(2):
    barcodes.append(rand_seq(bar_len))
  return barcodes


def format_read(seq, direction, mate, family_num, pair_num, barcodes, const, qual_char):
  # Read name (line 1):
  if (direction == 'forward' and mate == 1) or (direction == 'reverse' and mate == 2):
    order = 'ab'
  else:
    order = 'ba'
  yield '@fam{}.{}.pair{} mate{}'.format(family_num, order, pair_num, mate)
  # Read sequence (line 2):
  if order == 'ab':
    barcodes_ordered = barcodes
  if order == 'ba':
    barcodes_ordered = list(reversed(barcodes))
  barcode = barcodes_ordered[mate-1]
  if direction == 'reverse':
    yield barcode + const + revcomp(seq)
  else:
    yield barcode + const + seq
  # Plus (line 3):
  yield '+'
  # Quality scores (line 4):
  read_len = len(barcode + const + seq)
  yield qual_char * read_len


def revcomp(seq):
  return seq.translate(REVCOMP_TABLE)[::-1]


def rand_seq(length):
  barcode = ''
  for i in range(length):
    barcode += random.choice('ACGT')
  return barcode


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
