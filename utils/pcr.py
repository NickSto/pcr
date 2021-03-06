#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import errno
import logging
import argparse

ARG_DEFAULTS = {'single':False, 'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """Calculate the probability that an error that occurs somewhere during the a PCR
process with k cycles ends up in x reads out of n in a duplex family. Assumes a simple PCR model
where every fragment is doubled every cycle. It prints four tab-delimited columns: k, n, x, and the
probability of that x."""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('-x', type=int,
    help='Number of reads with the error. If omitted, it will output the probability of every '
         'x from 1 to n/2 (or n-1 if --one-sided).')
  parser.add_argument('-n', type=int, required=True,
    help='Total number of reads.')
  parser.add_argument('-k', type=int, required=True,
    help='Number of PCR cycles.')
  parser.add_argument('-1', '--1-sided', dest='single', action='store_true',
    help='Only calculate the literal probability of x errors in n reads. Contrast with '
         '--double-sided.')
  parser.add_argument('-2', '--2-sided', dest='single', action='store_false',
    help='Output P(x/n) + P((n-x)/n) (default). This is how real errors will appear in families, '
         'since errors over 50%% will be considered the "correct", consensus base.')
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

  if args.x is not None and (args.x < 1 or args.x >= args.n):
    fail('-x must be between 0 and -n (you gave {})'.format(args.x))

  if args.x:
    x_values = [args.x]
  else:
    if args.single:
      x_values = range(1, args.n)
    else:
      x_values = range(1, args.n//2+1)

  for x in x_values:
    if args.single or x/args.n == 0.5:
      p = get_maf_prob(args.k, args.n, x)
    else:
      p1 = get_maf_prob(args.k, args.n, x)
      p2 = get_maf_prob(args.k, args.n, args.n-x)
      p = p1 + p2
    print(args.k, args.n, x, p, sep='\t')


def get_maf_prob(k, n, x):
  """Calculate the equation:
  $\frac{\sum_{i=1}^k 2^i {n \choose x} \frac{1}{2^i}^x (1 - \frac{1}{2^i})^{n - x} }
  {\sum_{y=1}^{n-1} \sum_{i=1}^k 2^i {n \choose y} \frac{1}{2^i}^y (1 - \frac{1}{2^i})^{n-y}}$
  Where n is the total number of reads in the family, x is the number of reads containing a given
  error, and k is the number of PCR cycles used. x/n should then be the frequency of the error in
  the family.
  """
  numerator = summation(equation1, 1, k, n, x)
  denominator = 0
  for y in range(1, n):
    denominator += summation(equation1, 1, k, n, y)
  return numerator/denominator


def summation(function, start, end, *args):
  sum = 0
  for i in range(start, end+1):
    sum += function(i, *args)
  return sum


def equation1(i, n, x):
  two_i = 2**i
  mult1 = two_i
  mult2 = n_choose_k(n, x)
  mult3 = (1/two_i)**x
  mult4 = (1-(1/two_i))**(n-x)
  return mult1 * mult2 * mult3 * mult4


def n_choose_k(n, k):
  return factorial(n)/(factorial(k)*factorial(n-k))


def factorial(n):
  """A non-recursive factorial function. Because why not."""
  product = 1
  for i in range(n, 1, -1):
    product *= i
  return product


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
