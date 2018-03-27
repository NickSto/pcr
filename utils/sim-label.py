#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import argparse
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(root_dir)
from makrutenko import getreads

ARG_DEFAULTS = {}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Convert the results of the duplex pipeline on simulated data to a single-line tsv
format, and label them with their original fragments."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('reads', type=getreads.FastaReader,
    help='Output fasta file from duplex pipeline. This script depends on their names being the '
         'exact format produced by sim.py.')
  parser.add_argument('families', type=argparse.FileType('r'),
    help='Families tsv file (output of make-barcodes.awk).')

  args = parser.parse_args(argv[1:])

  bars_to_frags = {}
  for line in args.families:
    barcode, order, read_name = line.rstrip('\r\n').split('\t')[:3]
    if read_name.startswith('@'):
      read_name = read_name[1:]
    chrom, frag_id, read_num = read_name.split('-')
    bars_to_frags[barcode] = (chrom, frag_id)

  reads = iter(args.reads)
  while True:
    try:
      read = next(reads)
    except StopIteration:
      break
    barcode = read.id
    try:
      chrom, frag_id = bars_to_frags[barcode]
    except KeyError:
      sys.stderr.write('Missing barcode: {}\n'.format(barcode))
      continue
    try:
      frag_num = int(frag_id, 16)
    except ValueError:
      sys.stderr.write('Invalid fragment id: {}\n'.format(frag_id))
      continue
    fam1size, fam2size = get_famsizes(read.name)
    print(chrom, frag_num, frag_id, barcode, fam1size, fam2size, read.seq, sep='\t')


def get_famsizes(read_name):
  try:
    faminfo = read_name.split()[1]
  except IndexError:
    faminfo = ''
  famsizes = faminfo.split('-')
  try:
    fam1size = famsizes[0]
  except IndexError:
    fam1size = ''
  try:
    fam2size = famsizes[1]
  except IndexError:
    fam2size = ''
  return fam1size, fam2size


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
