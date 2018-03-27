#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import os
import sys
import random
import argparse
import numpy
import sim
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(root_dir)
from makrutenko import getreads

ARG_DEFAULTS = {'spacing':600, 'indel_rate':0.15, 'ext_rate':0.3}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Generate a version of the input reference with randomly added variants."""


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('ref', type=getreads.FastaReader,
    help='The original reference genome (FASTA).')
  parser.add_argument('-m', '--mutations', type=argparse.FileType('w'),
    help='Write inserted mutations here.')
  parser.add_argument('-s', '--spacing', type=int,
    help='Average (or exact) spacing between variants, in bp. Note: must be less than the line '
         'width in the reference file. Default: %(default)s')
  parser.add_argument('-r', '--random', action='store_true',
    help='Randomly distribute the variants instead of evenly spacing them.')
  parser.add_argument('-N', '--ref-name',
    help='Name the output sequence this.')
  parser.add_argument('-i', '--indel-rate', type=float,
    help='Default: %(default)s')
  parser.add_argument('-E', '--extension-rate', dest='ext_rate', type=float,
    help='Default: %(default)s')
  parser.add_argument('-S', '--seed', type=int,
    help='Default: random.')

  args = parser.parse_args(argv[1:])

  if args.seed is None:
    seed = random.randint(0, 2**31-1)
    sys.stderr.write('seed: {}\n'.format(seed))
  else:
    seed = args.seed
  random.seed(seed)

  var_prob = 1/args.spacing

  next_var = args.spacing
  coord = 0
  chr_name = None
  for read in args.ref:
    if read.name != chr_name:
      chr_name = read.name
      if args.ref_name:
        chr_id = args.ref_name.split()[0]
        print('>'+args.ref_name)
      else:
        chr_id = read.id
        print('>'+chr_name)
    end_coord = coord + len(read.seq)
    new_line = read.seq
    if args.random:
      n_vars = numpy.random.binomial(len(read.seq), var_prob)
      for i in range(n_vars):
        vcoord = random.randint(coord+1, end_coord)
        vtype, alt = sim.make_mutation(args.indel_rate, args.ext_rate)
        var = {'coord':vcoord-coord-1, 'type':vtype, 'alt':alt}
        # sys.stderr.write('Adding var at {}: {} ({})\n'.format(vcoord, vtype, alt))
        new_line = sim.apply_mutation(var, new_line)
        if args.mutations:
          var['coord'] = vcoord
          sim.log_mutations(args.mutations, [var], '.', chr_id, '.', '.')
    else:
      if next_var <= end_coord:
        vtype, alt = sim.make_mutation(args.indel_rate, args.ext_rate)
        var = {'coord':next_var-coord-1, 'type':vtype, 'alt':alt}
        # sys.stderr.write('Adding var at {}: {} ({})\n'.format(next_var, vtype, alt))
        new_line = sim.apply_mutation(var, new_line)
        if args.mutations:
          var['coord'] = next_var
          sim.log_mutations(args.mutations, [var], '.', chr_id, '.', '.')
        next_var += args.spacing
    print(new_line)
    coord += len(read.seq)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
