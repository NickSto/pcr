#!/usr/bin/env python
import os
import sys
import errno
import ctypes
import argparse

# Locate the library file.
LIBFILE = 'libconsensus.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

consensus = ctypes.cdll.LoadLibrary(library_path)
consensus.get_consensus.restype = ctypes.c_char_p
consensus.get_consensus_duplex.restype = ctypes.c_char_p
consensus.build_consensus_duplex_simple.restype = ctypes.c_char_p

ARG_DEFAULTS = {'alignment':sys.stdin}
DESCRIPTION = "Get the consensus of a set of aligned sequences."


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)
  parser.add_argument('alignment', type=argparse.FileType('r'),
    help='The aligned sequences, in FASTA format (but no multi-line sequences).')
  return parser


def main(argv):
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  sequences = []
  line_num = 0
  for line in args.alignment:
    line_num += 1
    if line_num % 2 == 0:
      sequences.append(line.rstrip('\r\n'))
  cons = get_consensus(sequences)
  print(cons)


# N.B.: The quality scores must be aligned with their accompanying sequences.
def get_consensus(align, quals=[], cons_thres=-1.0, min_reads=0, qual_thres=' ', gapped=False):
  cons_thres_c = ctypes.c_double(cons_thres)
  qual_thres_c = ctypes.c_char(qual_thres)
  n_seqs = len(align)
  if gapped:
    gapped_c = 1
  else:
    gapped_c = 0
  assert not quals or len(quals) == n_seqs, 'Different number of sequences and quals.'
  seq_len = None
  for seq in (align + quals):
    if seq_len is None:
      seq_len = len(seq)
    else:
      if seq_len != len(seq):
        raise AssertionError('All sequences in the alignment must be the same length: {}bp != {}bp.'
                             '\nAlignment:\n{}'.format(seq_len, len(seq), '\n'.join(align)))
  align_c = (ctypes.c_char_p * n_seqs)()
  for i, seq in enumerate(align):
    align_c[i] = ctypes.c_char_p(seq)
  quals_c = (ctypes.c_char_p * n_seqs)()
  for i, qual in enumerate(quals):
    quals_c[i] = ctypes.c_char_p(qual)
  if not quals:
    quals_c = 0
  return consensus.get_consensus(align_c, quals_c, n_seqs, seq_len, cons_thres_c, min_reads,
                                 qual_thres_c, gapped_c)


# N.B.: The quality scores must be aligned with their accompanying sequences.
def get_consensus_duplex(align1, align2, quals1=[], quals2=[], cons_thres=-1.0, min_reads=0,
                         qual_thres=' ', method='iupac'):
  assert method in ('iupac', 'freq')
  cons_thres_c = ctypes.c_double(cons_thres)
  qual_thres_c = ctypes.c_char(qual_thres)
  n_seqs1 = len(align1)
  n_seqs2 = len(align2)
  assert (not quals1 and not quals2) or (quals1 and quals2)
  assert not quals1 or len(quals1) == n_seqs1
  assert not quals2 or len(quals2) == n_seqs2
  seq_len = None
  for seq in (align1 + align2 + quals1 + quals2):
    if seq_len is None:
      seq_len = len(seq)
    else:
      assert seq_len == len(seq), 'All sequences in the alignment must be the same length.'
  align1_c = (ctypes.c_char_p * n_seqs1)()
  for i, seq in enumerate(align1):
    align1_c[i] = ctypes.c_char_p(seq)
  align2_c = (ctypes.c_char_p * n_seqs1)()
  for i, seq in enumerate(align2):
    align2_c[i] = ctypes.c_char_p(seq)
  quals1_c = (ctypes.c_char_p * n_seqs1)()
  for i, seq in enumerate(quals1):
    quals1_c[i] = ctypes.c_char_p(seq)
  quals2_c = (ctypes.c_char_p * n_seqs1)()
  for i, seq in enumerate(quals2):
    quals2_c[i] = ctypes.c_char_p(seq)
  if not quals1:
    quals1_c = 0
  if not quals2:
    quals2_c = 0
  return consensus.get_consensus_duplex(align1_c, align2_c, quals1_c, quals2_c, n_seqs1, n_seqs2,
                                        seq_len, cons_thres_c, min_reads, qual_thres_c, method)


def build_consensus_duplex_simple(cons1, cons2, gapped=False):
  assert len(cons1) == len(cons2)
  cons1_c = ctypes.c_char_p(cons1)
  cons2_c = ctypes.c_char_p(cons2)
  if gapped:
    gapped_c = 1
  else:
    gapped_c = 0
  return consensus.build_consensus_duplex_simple(cons1_c, cons2_c, gapped_c)


if __name__ == '__main__':
  sys.exit(main(sys.argv))
