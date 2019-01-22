#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import re
import os
import sys
import argparse
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(root_dir)
import swalign
from bfx import getreads

CANON = 'ACGT-'

WGSIM_ID_REGEX = r'^(.+)_\d+_\d+_\d+:\d+:\d+_\d+:\d+:\d+_([0-9a-f]+)/[12]$'
ARG_DEFAULTS = {'print_stats':True}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Correlate (labeled) reads from duplex pipeline with truth from simulator input,
and print the number of errors."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('reads', type=argparse.FileType('r'),
    help='Output from duplex pipeline. Should be the tsv produced by sim-label.py.')
  parser.add_argument('frags', type=getreads.FastqReader,
    help='--frag-file from sim.py.')
  parser.add_argument('-i', '--ignore-ambiguous', action='store_true',
    help='Don\'t consider ambiguous bases ("N", "R", etc.) in SNV errors. Specifically, it will '
         'ignore any mismatch between a non-gap base in the fragment and read base that isn\'t '
         'one of "'+CANON+'".')
  parser.add_argument('-a', '--print-alignments', action='store_true',
    help='Print the alignments of each read with each fragment. Mostly for debug purposes.')
  parser.add_argument('-S', '--no-stats', dest='print_stats', action='store_false',
    help='Don\'t print the normal output of statistics on differences.')

  args = parser.parse_args(argv[1:])

  reads = iter(args.reads)
  frags = iter(args.frags)
  while True:
    # Read in the next output read.
    try:
      read_line = next(reads)
    except StopIteration:
      break
    fields = read_line.rstrip('\r\n').split('\t')
    assert len(fields) == 7, fields
    read = dict(zip(('chrom', 'frag_num', 'frag_id', 'bar', 'reads+', 'reads-', 'seq'), fields))
    # Read in fragments until we find the one corresponding to the output read.
    frag_chrom = None
    frag_frag_id = None
    while not (read['chrom'] == frag_chrom and read['frag_id'] == frag_frag_id):
      try:
        frag = next(frags)
      except StopIteration:
        break
      match = re.search(WGSIM_ID_REGEX, frag.id)
      if match:
        frag_chrom = match.group(1)
        frag_frag_id = match.group(2)
      else:
        sys.stderr.write('Invalid wgsim read name: {}\n'.format(frag.id))
    if frag_chrom is None and frag_frag_id is None:
      break
    # Align the output read to the fragment.
    align = swalign.smith_waterman_duplex(frag.seq, read['seq'])
    assert len(align.target) == len(align.query)
    if args.print_alignments:
      print(align.target)
    diffs = get_diffs(align.target, align.query, print_mid=args.print_alignments,
                      ignore_ambig=args.ignore_ambiguous)
    if args.print_alignments:
      print(align.query)
    read_len = len(read['seq'])
    snvs = ins = dels = 0
    for diff in diffs:
      if diff['type'] == 'snv':
        snvs += 1
      elif diff['type'] == 'ins':
        ins += 1
      elif diff['type'] == 'del':
        dels += 1
    match_rate = round(align.matches/read_len, 2)
    if args.print_stats:
      print(read['bar'], read['frag_id'], read['reads+'], read['reads-'], read_len,
            read_len-align.matches, match_rate, len(diffs), snvs, ins, dels, sep='\t')


def get_diffs(target, query, print_mid=False, ignore_ambig=False):
  diffs = []
  diff = None
  coord = 0
  for base1, base2 in zip(target, query):
    if base1 != '-':
      coord += 1
    if base1 == base2:
      # Finish ongoing indels and add them to the list.
      if diff is not None:
        # But omit the "indel" that's just the unaligned portion at the start.
        if diff['coord'] > 1:
          diffs.append(diff)
        diff = None
      if print_mid:
        sys.stdout.write('|')
    elif ignore_ambig and base1 != '-' and base2 not in CANON:
      if print_mid:
        sys.stdout.write(' ')
    elif base1 == '-':
      if diff is None:
        diff = {'coord':coord, 'type':'ins', 'alt':base2}
      else:
        diff['alt'] += base2
      if print_mid:
        sys.stdout.write(' ')
    elif base2 == '-':
      if diff is None:
        diff = {'coord':coord-1, 'type':'del', 'alt':1}
      else:
        diff['alt'] += 1
      if print_mid:
        sys.stdout.write(' ')
    else:
      diffs.append({'coord':coord, 'type':'snv', 'alt':base2})
      if print_mid:
        sys.stdout.write(' ')
  if diff is not None:
    diffs.append(diff)
  if print_mid:
    print()
  return diffs


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
