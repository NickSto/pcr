#!/usr/bin/env python
from __future__ import division
import sys
import argparse
import getreads

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
         'output, in the same order, but with two columns: stat name and value.')
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

  stats = get_stats(barcodes, tag_len=args.tag_len, min_reads=args.min_reads)
  print_stats(stats, min_reads=args.min_reads, human=args.human)


def read_families(infile):
  barcodes = {}
  for line in infile:
    fields = line.split('\t')
    barcode = fields[0]
    order = fields[1]
    try:
      barcodes[(barcode, order)] += 1
    except KeyError:
      barcodes[(barcode, order)] = 1
  return barcodes


def read_fastqs(infileh1, infileh2, tag_len=12, validate=False):
  reader1 = getreads.getparser(infileh1, filetype='fastq').parser()
  reader2 = getreads.getparser(infileh2, filetype='fastq').parser()
  barcodes = {}
  while True:
    try:
      read1 = reader1.next()
      read2 = reader2.next()
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
    try:
      barcodes[(barcode, order)] += 1
    except KeyError:
      barcodes[(barcode, order)] = 1
  return barcodes


def get_stats(barcodes, tag_len=12, min_reads=3):
  passed_sscs = 0
  duplexes = 0
  passed_duplexes = 0
  singletons = 0
  total_pairs = 0
  for (barcode, order), count in barcodes.items():
    total_pairs += count
    if count == 1:
      singletons += 1
    if count >= min_reads:
      passed_sscs += 1
    if order == 'ab':
      other_order = 'ba'
    else:
      other_order = 'ab'
    if (barcode, other_order) in barcodes:
      duplexes += 1
      other_count = barcodes[(barcode, other_order)]
      if count >= min_reads and other_count >= min_reads:
        passed_duplexes += 1
  # Each full duplex ends up being counted twice. Halve it to get the real total.
  stats = {
    'pairs':total_pairs,
    'barcodes':len(barcodes),
    'avg_pairs':total_pairs/len(barcodes),
    'singletons':singletons,
    'duplexes':duplexes//2,
    'passed_sscs':passed_sscs*2,
    'passed_duplexes':passed_duplexes,
  }
  return stats


def print_stats(stats, min_reads=3, human=True):
  all_stats = stats.copy()
  all_stats.update({'min_reads':min_reads})
  if human:
    print """Total read pairs:\t{pairs}
Unique barcodes:\t{barcodes}
Avg # of read pairs per barcode:\t{avg_pairs}
Singletons:\t{singletons}
Barcodes with reverse (other strand) present:\t{duplexes}
Passing threshold of {min_reads} reads per single-strand consensus:
\tSingle-strand consensus sequences:\t{passed_sscs}
\tDuplex consensus sequences:\t{passed_duplexes}""".format(**all_stats)
  else:
    for stat in ('pairs', 'barcodes', 'avg_pairs', 'singletons', 'duplexes', 'min_reads',
                 'passed_sscs', 'passed_duplexes'):
      print '{}\t{}'.format(stat, all_stats[stat])


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
