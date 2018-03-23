#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import argparse
import collections
import getreads

STATS = collections.OrderedDict(
  (
    ('pairs',       'Total pairs of reads in input'),
    ('barcodes',    'Unique barcodes (double-standed families)'),
    ('bar_orders',  'Unique barcode/order combinations (single-stranded families)'),
    ('avg_pairs',   'Average # of read pairs per unique barcode'),
    ('stranded_singletons', 'Single-strand singletons (only 1 read with that barcode/order combination)'),
    ('duplex_singletons',   'Double-strand singletons (only 1 read with that barcode)'),
    ('duplexes',    'Families with both strands present'),
    ('passed_sscs', 'Single-strand families with > {min_reads} reads'),
    ('passed_duplexes',     'Families with both strands with > {min_reads} reads'),
  )
)
OPT_DEFAULTS = {'tag_len':12, 'const_len':5, 'min_reads':3, 'human':True}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Print statistics on the raw duplex sequencing reads."""
EPILOG = """Warning: This tracks all barcodes in a dict, so it can take a lot of memory. A guideline
is about 200 bytes per (12bp) tag. For example, it took about 800MB for a 10GB, 32 million read
dataset with an average of 4 pairs per barcode."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
  parser.set_defaults(**OPT_DEFAULTS)

  parser.add_argument('infile1', metavar='reads_1.fq', nargs='?',
    help='The first mates in the read pairs.')
  parser.add_argument('infile2', metavar='reads_2.fq', nargs='?',
    help='The second mates in the read pairs.')
  parser.add_argument('-f', '--families', metavar='families.tsv')
  parser.add_argument('-t', '--tag-length', dest='tag_len', type=int)
  parser.add_argument('-c', '--constant-length', dest='const_len', type=int)
  parser.add_argument('-C', '--computer', dest='human', action='store_false',
    help='Print results in computer-readable format. This will be a tab-delimited version of the '
         'output, in the same order, but with two columns: value and stat name.')
  parser.add_argument('-m', '--min-reads', type=int,
    help='The minimum number of reads required in each single-stranded family. Default: '
         '%(default)s')
  parser.add_argument('-v', '--validate', action='store_true',
    help='Check the id\'s of the reads to make sure the correct reads are mated into pairs (the '
         'id\'s of mates must be identical).')

  args = parser.parse_args(argv[1:])

  if args.families:
    barcodes = {}
    with open(args.families) as infile:
      barcodes = read_families(infile)
  else:
    if not (args.infile1 and args.infile2):
      fail('Error: You must provide either a --families file or two fastq files.')
    with open(args.infile1) as infileh1:
      with open(args.infile2) as infileh2:
        barcodes = read_fastqs(infileh1, infileh2, tag_len=args.tag_len, validate=args.validate)

  stats = get_stats(barcodes, stats_key=STATS, min_reads=args.min_reads)
  print_stats(stats, stats_key=STATS, human=args.human)


def read_families(infile):
  barcodes = collections.Counter()
  for line in infile:
    fields = line.split('\t')
    barcode = fields[0]
    order = fields[1]
    barcodes[(barcode, order)] += 1
  return barcodes


def read_fastqs(infileh1, infileh2, tag_len=12, validate=False):
  reader1 = getreads.getparser(infileh1, filetype='fastq').parser()
  reader2 = getreads.getparser(infileh2, filetype='fastq').parser()
  barcodes = collections.Counter()
  while True:
    try:
      read1 = next(reader1)
      read2 = next(reader2)
    except StopIteration:
      break
    if validate and read1.id != read2.id:
      raise getreads.FormatError('Read pair mismatch: "{}" and "{}"'.format(read1.id, read2.id))
    alpha = read1.seq[:tag_len]
    beta  = read2.seq[:tag_len]
    if alpha < beta:
      order = 'ab'
      barcode = alpha + beta
    else:
      order = 'ba'
      barcode = beta + alpha
    barcodes[(barcode, order)] += 1
  return barcodes


def get_stats(barcodes, stats_key=STATS, min_reads=3):
  stats = {'min_reads':min_reads}
  for stat_name in stats_key:
    stats[stat_name] = 0
  barcode_counts = collections.Counter()
  for (barcode, order), count in barcodes.items():
    stats['pairs'] += count
    barcode_counts[barcode] += count
    if count == 1:
      stats['stranded_singletons'] += 1
    if count >= min_reads:
      stats['passed_sscs'] += 1
    if order == 'ab':
      other_order = 'ba'
    else:
      other_order = 'ab'
    if (barcode, other_order) in barcodes:
      other_count = barcodes[(barcode, other_order)]
      # The loop encounters both strands of each duplex, so to only count each duplex once,
      # only increment on one of the strands (ab).
      if order == 'ab':
        stats['duplexes'] += 1
        if count >= min_reads and other_count >= min_reads:
          stats['passed_duplexes'] += 1
  stats['barcodes'] = len(barcode_counts)
  stats['bar_orders'] = len(barcodes)
  stats['avg_pairs'] = stats['pairs']/stats['barcodes']
  stats['duplex_singletons'] = 0
  for count in barcode_counts.values():
    if count == 1:
      stats['duplex_singletons'] += 1
  return stats


def print_stats(stats, stats_key=STATS, human=True):
  if human:
    for stat_name, stat_description in stats_key.items():
      print(stats[stat_name], stat_description.format(**stats), sep='\t')
  else:
    for stat_name in stats_key:
      print(stats[stat_name], stat_name, sep='\t')


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
