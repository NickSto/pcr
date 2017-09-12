#!/usr/bin/env python
import os
import sys
import errno
import ctypes
import argparse
import tempfile

# Locate the library file.
LIBFILE = 'kalign/kalign.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

kalign = ctypes.cdll.LoadLibrary(library_path)


def make_argparser():
  parser = argparse.ArgumentParser(description='Align a set of sequences.')
  parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
    help='Input sequences.')
  # parser.add_argument('infile', metavar='seqs.fa')
  return parser


def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  seqs = []
  for line_raw in args.input:
    line = line_raw.rstrip('\r\n')
    if line.startswith('>'):
      continue
    else:
      seqs.append(line)
  alignment = align(seqs)
  print
  for i, seq in enumerate(alignment):
    print '>seq{}\n{}'.format(i, seq)


def align(seqs):
  """Perform a multiple sequence alignment on a set of sequences and parse the result.
  Uses MAFFT."""
  i = 0
  try:
    with tempfile.NamedTemporaryFile('w', delete=False, prefix='align.msa.') as input_file:
      for seq in seqs:
        i += 1
        input_file.write('>seq{}\n'.format(i))
        input_file.write(seq+'\n')
    output_file = tempfile.NamedTemporaryFile('r', delete=False, prefix='align.msa.')
    output_file.close()
    argc, argv = make_args(input_file.name, output_file.name)
    # print '$ '+' '.join([arg for arg in argv])
    kalign.main(argc, argv)
    return read_fasta(output_file.name)
  finally:
    # Make sure we delete the temporary files.
    try:
      os.remove(input_file.name)
    except OSError:
      pass
    try:
      os.remove(output_file.name)
    except OSError:
      pass


def make_args(infile, outfile):
  argv_list = ('kalign', infile, '-o', outfile)
  argc = len(argv_list)
  argv_c = strlist_to_c(argv_list)
  return argc, argv_c


def strlist_to_c(strlist):
  c_strs = (ctypes.c_char_p * len(strlist))()
  for i, s in enumerate(strlist):
    c_strs[i] = ctypes.c_char_p(s)
  return c_strs


def read_fasta(fasta):
  """Quick and dirty FASTA parser. Return the sequences and their names.
  Returns a list of sequences.
  Warning: Reads the entire contents of the file into memory at once."""
  sequences = []
  sequence = ''
  with open(fasta) as fasta_file:
    for line in fasta_file:
      if line.startswith('>'):
        if sequence:
          sequences.append(sequence.upper())
        sequence = ''
        continue
      sequence += line.strip()
  if sequence:
    sequences.append(sequence.upper())
  return sequences


if __name__ == '__main__':
  sys.exit(main(sys.argv))
