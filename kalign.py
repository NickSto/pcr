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


class AlnStrs(ctypes.Structure):
  _fields_ = [
    ('nseqs', ctypes.c_int),
    ('seqlen', ctypes.c_int),
    ('names', ctypes.POINTER(ctypes.c_char_p)),
    ('seqs', ctypes.POINTER(ctypes.c_char_p)),
  ]

kalign.main.restype = ctypes.POINTER(AlnStrs)


def make_argparser():
  parser = argparse.ArgumentParser(description='Align a set of sequences.')
  parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
    help='Input sequences.')
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
  for i in range(alignment.nseqs):
    print alignment.seqs[i]


def align(seqs):
  """Perform a multiple sequence alignment on a set of sequences and parse the result."""
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
    aln = kalign.main(argc, argv)
    return aln.contents
  finally:
    # Make sure we delete the temporary files.
    try:
      os.remove(input_file.name)
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


if __name__ == '__main__':
  sys.exit(main(sys.argv))
