#!/usr/bin/env python3
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import math
import errno
import logging
import argparse
import collections
assert sys.version_info.major >= 3, 'Python 3 required'

ARG_DEFAULTS = {'families':sys.stdin, 'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """Check bias in the split of how many reads per family are on strand ab vs ba."""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('families', metavar='families.tsv', type=argparse.FileType('r'), nargs='?',
    help='')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()

  counts = collections.defaultdict(dict)
  strand1 = 0
  strand2 = 0
  last_bar = None
  last_order = None
  for line in args.families:
    fields = line.split('\t')
    bar = fields[0]
    order = fields[1]
    if bar != last_bar:
      if last_bar is not None:
        add_counts(counts, strand1, strand2)
      strand1 = 0
      strand2 = 0
    elif order != last_order:
      strand2 = strand1
      strand1 = 0
    strand1 += 1
    last_bar = bar
    last_order = order
  add_counts(counts, strand1, strand2)

  for famsize in sorted(counts.keys()):
    family_counts = counts[famsize]
    i = math.ceil(famsize/2)
    strand_counts = []
    while i <= famsize:
      j = famsize - i
      try:
        strand_counts.append(family_counts[(i,j)])
      except KeyError:
        strand_counts.append(0)
      i += 1
    print(famsize, *strand_counts, sep='\t')


def add_counts(counts, strand1, strand2):
  if strand1 < strand2:
    strand1, strand2 = strand2, strand1
  try:
    counts[strand1+strand2][(strand1, strand2)] += 1
  except KeyError:
    counts[strand1+strand2][(strand1, strand2)] = 1


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
