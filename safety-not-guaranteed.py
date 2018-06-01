#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import argparse
import errno
import logging
import sys
from bfx import getreads

# The ascii values that represent a 0 PHRED score.
QUAL_OFFSETS = {'sanger':33, 'solexa':64}
USAGE = '$ %(prog)s -q [phred score] reads.fa > reads.fq'
DESCRIPTION = """Convert FASTA to FASTQ."""
EPILOG = """WARNING: There is no quality information in a FASTA, so this must use a dummy value
for the quality scores. This is inherently false/artifactual information, so keep that in mind when
using the output in downstream tools."""


def make_argparser():
  parser = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)
  parser.add_argument('input', metavar='reads.fa', nargs='?', type=argparse.FileType('r'),
    default=sys.stdin,
    help='Input FASTA reads. Omit to read from stdin.')
  parser.add_argument('-o', '--output', metavar='reads.fq', type=argparse.FileType('w'),
    default=sys.stdout,
    help='Write output FASTQ reads to this file. If not given, will print to stdout.')
  parser.add_argument('-q', '--phred-score', required=True, type=int,
    help='The quality score to give to all bases. There is no meaningful quality score we can '
         'automatically give, so you will have to specify an artificial one. A good choice is 40, '
         'the maximum score normally output by sequencers.')
  parser.add_argument('-F', '--qual-format', choices=('sanger', 'solexa'), default='sanger',
    help='FASTQ quality score format. Sanger scores are assumed to begin at \'{}\' ({}). '
         'Default: %(default)s.'.format(QUAL_OFFSETS['sanger'], chr(QUAL_OFFSETS['sanger'])))
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  qual_start = QUAL_OFFSETS[args.qual_format]
  if qual_start+args.phred_score > 126:
    fail('Error: PHRED score ({}) is too large.'.format(args.fastq_out))
  qual_char = chr(qual_start+args.phred_score)

  fasta_to_fastq(args.input, args.output, qual_char)

  if args.input is not sys.stdin and not args.input.closed:
    args.input.close()
  if args.output is not sys.stdout and not args.output.closed:
    args.output.close()


def fasta_to_fastq(fasta_file, fastq_file, qual_char):
  for read in getreads.getparser(fasta_file, filetype='fasta'):
    quals = qual_char * len(read.seq)
    fastq_file.write('@{0}\n{1}\n+\n{2}\n'.format(read.name, read.seq, quals))


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
