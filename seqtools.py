import os
import errno
import ctypes

# Locate the library file.
LIBFILE = 'libseqtools.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

seqtools = ctypes.cdll.LoadLibrary(library_path)
seqtools.get_revcomp.restype = ctypes.c_char_p
seqtools.transfer_gaps.restype = ctypes.c_char_p


def get_revcomp(seq):
  return seqtools.get_revcomp(seq)


def get_diffs_frac_simple(consensus, family):
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    c_family[i] = ctypes.c_char_p(seq)
  seqtools.get_diffs_frac_simple.restype = ctypes.POINTER(ctypes.c_double * len(c_family))
  diffs = seqtools.get_diffs_frac_simple(c_consensus, c_family, len(c_family))
  return diffs.contents


def get_diffs_frac_binned(consensus, family, bins):
  seq_len = None
  c_consensus = ctypes.c_char_p(consensus)
  c_family = (ctypes.c_char_p * len(family))()
  for i, seq in enumerate(family):
    if seq_len:
      if seq_len != len(seq):
        return None
    else:
      seq_len = len(seq)
    c_family[i] = ctypes.c_char_p(seq)
  double_array_pointer = ctypes.POINTER(ctypes.c_double * bins)
  seqtools.get_diffs_frac_binned.restype = ctypes.POINTER(double_array_pointer * len(c_family))
  diffs_binned_c = seqtools.get_diffs_frac_binned(c_consensus, c_family, len(c_family), seq_len, bins)
  diffs_binned = []
  for diffs_c in diffs_binned_c.contents:
    diffs_binned.append(diffs_c.contents)
  return diffs_binned


def transfer_gaps(aligned, seq, gap_char_in='-', gap_char_out='-'):
  gap_char_in_c = ctypes.c_char(gap_char_in)
  gap_char_out_c = ctypes.c_char(gap_char_out)
  return seqtools.transfer_gaps(aligned, seq, gap_char_in_c, gap_char_out_c)


def transfer_gaps_multi(seqs, aligned, gap_char_in='-', gap_char_out='-'):
  gap_char_in_c = ctypes.c_char(gap_char_in)
  gap_char_out_c = ctypes.c_char(gap_char_out)
  n_seqs = len(seqs)
  assert n_seqs == len(aligned), ('Error: Unequal number of gapped and ungapped sequences ({} vs '
                                  '{} sequences, respectively).'.format(len(aligned), n_seqs))
  seqs_c = (ctypes.c_char_p * n_seqs)()
  for i, seq in enumerate(seqs):
    seqs_c[i] = ctypes.c_char_p(seq)
  aligned_c = (ctypes.c_char_p * n_seqs)()
  for i, seq in enumerate(aligned):
    aligned_c[i] = ctypes.c_char_p(seq)
  seqtools.transfer_gaps_multi.restype = ctypes.POINTER(ctypes.c_char_p * n_seqs)
  output_c = seqtools.transfer_gaps_multi(n_seqs, aligned_c, seqs_c, gap_char_in_c, gap_char_out_c)
  output = []
  for seq in output_c.contents:
    output.append(seq)
  return output
