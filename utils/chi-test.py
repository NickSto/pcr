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
from matplotlib import pyplot

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
  parser.add_argument('-p', '--plots', dest='output', action='store_const', const='plots',
    help='Use matplotlib to plot the distributions.')
  parser.add_argument('-o', '--plots-outpath',
    help='Instead of displaying the plots live, save them to image files. Use this path base to '
         'create the filenames, with the format "[plots-outpath][famsize].png"')
  parser.add_argument('-T', '--plots-title')
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
    exp_counts, obs_counts, exp_counts_filt, obs_counts_filt = calc_obs_exp(famsize,
                                                                            total_errors,
                                                                            counts,
                                                                            exp_freqs,
                                                                            args.min_value)
    if args.output == 'dists':
      i = 0
      for exp, obs in zip(exp_counts, obs_counts):
        i += 1
        print('{}\t{}\t{}\t{}'.format(famsize, i, sig_round(exp), sig_round(obs)))
      continue
    # Is there enough data for a valid chi square test?
    if passes_thresholds(exp_counts, total_errors, args.min_obs, famsize):
      chi = scipy.stats.chisquare(exp_counts_filt, obs_counts_filt)
      if args.log:
        chi_stat, chi_p = log_transform((chi.statistic, chi.pvalue), default=float('nan'))
      else:
        chi_stat = chi.statistic
        chi_p = chi.pvalue
    else:
      chi = None
    if args.output == 'stats':
      if chi is not None:
        print('{}\t{}\t{:8.5f}\t{}'.format(famsize, total_errors, chi_stat, chi_p))
    elif args.output == 'plots':
      plot_dists(exp_counts, obs_counts, famsize, args.log, args.plots_outpath, args.plots_title)


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


def calc_obs_exp(famsize, total_errors, counts, exp_freqs, min_value=5):
  observations = counts[famsize]
  exp_counts = []
  obs_counts = []
  exp_counts_filt = []
  obs_counts_filt = []
  for i, exp_freq in enumerate(exp_freqs):
    repeat_num = i+1
    obs_count = float(observations[repeat_num])
    exp_count = exp_freq * total_errors
    obs_counts.append(obs_count)
    exp_counts.append(exp_count)
    if obs_count < min_value:
      logging.info('Observed value {} in family {} is below --min-value {}.'
                   .format(int(obs_count), famsize, min_value))
      continue
    if exp_count < min_value:
      logging.info('Expected value {:5.2f} in family {} is below --min-value {}.'
                   .format(exp_count, famsize, min_value))
      continue
    exp_counts_filt.append(obs_count)
    obs_counts_filt.append(exp_count)
  return exp_counts, obs_counts, exp_counts_filt, obs_counts_filt


def passes_thresholds(exp_counts, total_errors, min_obs, famsize):
  if len(exp_counts) <= 1:
    logging.info('{} in family {} is too few frequencies'.format(len(exp_counts), famsize))
    return False
  if total_errors < min_obs:
    logging.info('{} in family {} is too few errors to test against.'
                 .format(total_errors, famsize))
    return False
  return True


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def chisquare(obs_dist, exp_dist):
  chi = 0
  for obs, exp in zip(obs_dist, exp_dist):
    chi += ((obs-exp)**2)/exp
  return chi


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


def sig_round(n, figs=3):
  """Round a number "n" to "figs" significant figures."""
  if n == 0 or math.isnan(n) or math.isinf(n):
    return n
  elif n < 0:
    n = -n
    sign = -1
  else:
    sign = 1
  magnitude = int(math.floor(math.log10(n)))
  return sign * round(n, figs - 1 - magnitude)


def log_transform(raw, default=-2):
  log = []
  for num in raw:
    if num == 0:
      log.append(default)
    else:
      log.append(math.log(num, 10))
  return log


def plot_dists(exp_counts, obs_counts, famsize, log, outpath, title):
  figure = pyplot.figure(dpi=120, figsize=(6,4.5))
  axis = figure.add_subplot(1, 1, 1)
  if log:
    exp_counts = log_transform(exp_counts)
    obs_counts = log_transform(obs_counts)
  x_values = list(range(1, len(exp_counts)+1))
  axis.scatter(x_values, exp_counts, s=10, color='blue', label='expected')
  axis.scatter(x_values, obs_counts, s=10, color='red', label='observed')
  set_plot_display(axis, x_values, exp_counts, obs_counts, famsize, title, log)
  show_plot(outpath, famsize)


def set_plot_display(axis, x_values, exp_counts, obs_counts, famsize, title, log):
  # Set the limits of the graph.
  MARGIN_SIZE = 0.05
  X_MAX = 5
  Y_MAX = 4
  X_MIN = 0
  if log:
    Y_MIN = -2
  else:
    Y_MIN = 0
  x_max = max(x_values[-1], X_MAX)
  x_margin = x_max * MARGIN_SIZE
  pyplot.xlim(X_MIN-x_margin, x_max+x_margin)
  y_max = max(max(exp_counts + obs_counts), Y_MAX)
  y_margin = y_max * MARGIN_SIZE
  pyplot.ylim(Y_MIN-y_margin, y_max+y_margin)
  # Draw a grid line across the 0 of both axes.
  axis.set_xticks([0], minor=True)
  if log:
    axis.set_yticks([-2], minor=True)
    # Replace the tick label at -2 with "Zero".
    yticks = axis.get_yticks()
    labels = axis.get_yticklabels()
    i = 0
    while i < len(yticks):
      if yticks[i] == -2.0:
        labels[i] = 'Zero'
      else:
        labels[i] = int(yticks[i])
      i += 1
    axis.set_yticklabels(labels)
  else:
    axis.set_yticks([0], minor=True)
  pyplot.grid(axis='both', which='minor')
  # Set the figure text.
  pyplot.legend(loc='upper right')
  pyplot.xlabel('x (# of reads errors occur in)')
  if log:
    pyplot.ylabel('Log10(Number of errors)')
  else:
    pyplot.ylabel('Number of errors')
  title += '\nfamilies with {} reads'.format(famsize)
  pyplot.title(title)


def show_plot(plots_outpath, famsize):
  try:
    if plots_outpath:
      plot_path = '{}{}.png'.format(plots_outpath, famsize)
      pyplot.savefig(plot_path)
    else:
      pyplot.show()
  except KeyboardInterrupt:
    if __name__ == '__main__':
      sys.exit()
    else:
      raise
  pyplot.close()


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except IOError as ioe:
    if ioe.errno != errno.EPIPE:
      raise
