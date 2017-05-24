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
import scipy.stats

ARG_DEFAULTS = {'min_obs':20, 'min_value':5, 'output':'stats',
                'log_file':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """"""


def make_argparser():

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('expected', metavar='expected.tsv', type=argparse.FileType('r'))
  parser.add_argument('observed', metavar='observed.tsv', type=argparse.FileType('r'))
  parser.add_argument('-m', '--min-value', type=float,
    help='The minimum acceptable value for an individual observation. Default: %(default)s')
  parser.add_argument('-M', '--min-obs', type=int,
    help='The minimum number of errors required to perform the test on. Default: %(default)s')
  parser.add_argument('-d', '--dists', dest='output', action='store_const', const='dists',
    help='Instead of printing the calculated test statistic for each family size, print the raw '
         'expected and observed distributions. The output is '
         'tab-delimited, with three columns: # reads in family, # of reads the errors occurred in, '
         'expected number of errors occurring in this number of reads, observed number.')
  parser.add_argument('-l', '--log', action='store_true',
    help='Log10 transform the chi statistics and p-values.')
  parser.add_argument('-L', '--log-file', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
  parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log_file, level=args.volume, format='%(message)s')
  tone_down_logger()

  dists = read_expected(args.expected)
  counts = read_observed(args.observed)
  totals = sum_counts(counts)

  for famsize, exp_freqs in dists.items():
    # Calculate the observed and expected distributions.
    total_errors = totals[famsize]
    observations = counts[famsize]
    exp_counts = []
    obs_counts = []
    for i, exp_freq in enumerate(exp_freqs):
      repeat_num = i+1
      obs_count = float(observations[repeat_num])
      if obs_count < args.min_value:
        logging.info('Observed value {} in family {} is below --min-value {}.'
                     .format(int(obs_count), famsize, args.min_value))
        continue
      exp_count = exp_freq * total_errors
      if exp_count < args.min_value:
        logging.info('Expected value {:5.2f} in family {} is below --min-value {}.'
                     .format(exp_count, famsize, args.min_value))
        continue
      obs_counts.append(obs_count)
      exp_counts.append(exp_count)
    # Check filters.
    if len(exp_counts) <= 1:
      logging.info('{} in family {} is too few frequencies'.format(len(exp_freqs), famsize))
      continue
    if total_errors < args.min_obs:
      logging.info('{} in family {} is too few errors to test against.'
                   .format(total_errors, famsize))
      continue
    # Do the Chi Square test.
    chi = scipy.stats.chisquare(obs_counts, exp_counts)
    if args.log:
      chi_stat = math.log(chi.statistic, 10)
      chi_p = math.log(chi.pvalue, 10)
    else:
      chi_stat = chi.statistic
      chi_p = chi.pvalue
    if args.output == 'dists':
      i = 0
      for exp, obs in zip(exp_counts, obs_counts):
        i += 1
        print('{}\t{}\t{}\t{}'.format(famsize, i, smart_round(exp), smart_round(obs)))
    else:
      print('{}\t{}\t{:8.5f}\t{}'.format(famsize, total_errors, chi_stat, chi_p))


def read_expected(expected_file):
  dists = collections.defaultdict(list)
  for line in expected_file:
    fields = line.rstrip('\r\n').split('\t')
    famsize = int(fields[1])
    prob = float(fields[3])
    dists[famsize].append(prob)
  return dists


def read_observed(observed_file):
  # Each line represents one unique error in a family (e.g. a C->T at position 52).
  # There are two columns:
  # 1. famsize: The size of the family
  # 2. errors:  The number of reads the error appeared in
  counts = collections.defaultdict(lambda: collections.defaultdict(int))
  for line in observed_file:
    fields = line.rstrip('\r\n').split('\t')
    famsize = int(fields[0])
    errors = int(fields[1])
    counts[famsize][errors] += 1
  return counts


def sum_counts(counts):
  totals = {}
  # "observations" is a dict mapping x to c, where c is the number of times you see an error
  # repeated x times.
  for famsize, observations in counts.items():
    totals[famsize] = sum(observations.values())
  # "totals" is a dict mapping family sizes to error totals (the total number of unique errors
  # observed in all families of that size).
  return totals


def calc_freqs(counts, totals):
  freqs = collections.defaultdict(list)
  for famsize, observations in counts.items():
    total = totals[famsize]
    for num_errors in range(1, famsize//2+1):
      freqs[famsize].append(observations[num_errors] / total)
  return freqs


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


def smart_round(x):
  if x == 0:
    return x
  power = int(math.log(x, 10))
  if x < 1:
    power -= 1
  decimals = max(2-power, 0)
  rounded = round(x, decimals)
  if decimals == 0:
    return int(rounded)
  else:
    return rounded


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
