import os
import errno
import ctypes
import string

# Locate the library file.
LIBFILE = 'libswalign.so'
script_dir = os.path.dirname(os.path.realpath(__file__))
library_path = os.path.join(script_dir, LIBFILE)
if not os.path.isfile(library_path):
  ioe = IOError('Library file "'+LIBFILE+'" not found.')
  ioe.errno = errno.ENOENT
  raise ioe

swalign = ctypes.cdll.LoadLibrary(library_path)

REVCOMP_TABLE = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')


# C struct for ctypes
class SeqPairC(ctypes.Structure):
  _fields_ = [
    ('a', ctypes.c_char_p),
    ('alen', ctypes.c_uint),
    ('b', ctypes.c_char_p),
    ('blen', ctypes.c_uint),
  ]


# C struct for ctypes
class AlignC(ctypes.Structure):
  _fields_ = [
    ('seqs', ctypes.POINTER(SeqPairC)),
    ('start_a', ctypes.c_int),
    ('start_b', ctypes.c_int),
    ('end_a', ctypes.c_int),
    ('end_b', ctypes.c_int),
    ('matches', ctypes.c_int),
    ('score', ctypes.c_double),
  ]


# The Python version
class Align(object):
  def __init__(self, align_c):
    self.target = align_c.seqs.contents.a
    self.query = align_c.seqs.contents.b
    # Where the first base of the target aligns on the query, in query coordinates (or 1, if <= 0).
    self.start_target = align_c.start_a
    # Where the first base of the query aligns on the target, in target coordinates (or 1, if <= 0).
    self.start_query = align_c.start_b
    # Where the last base of the target aligns on the query, in query coordinates.
    self.end_target = align_c.end_a
    # Where the last base of the query aligns on the target, in target coordinates.
    self.end_query = align_c.end_b
    self.matches = align_c.matches
    self.score = align_c.score

  # Provide this common function.
  def __str__(self):
    """Print a human-readable representation of the alignment."""
    start_query = str(self.start_query)
    start_target = str(self.start_target)
    start_width = str(max(len(start_query), len(start_target)))
    line_format = '{:'+start_width+'} {} {}'
    output = line_format.format(start_target, self.target, self.end_target) + '\n'
    output += line_format.format(start_query, self.query, self.end_query)
    return output


# Initialize functions (define types).
swalign.smith_waterman.restype = ctypes.POINTER(AlignC)
swalign.revcomp.restype = ctypes.c_char_p


def smith_waterman(target, query):
  seq_pair = SeqPairC(target, len(target), query, len(query))
  align_c = swalign.smith_waterman(ctypes.pointer(seq_pair), 1).contents
  return Align(align_c)


def smith_waterman_duplex(target, query):
  """Smith-Waterman align query to target in both orientations and return the best.
  Convenience function that calls smith_waterman() twice, and returns the
  alignment with the highest score."""
  align = smith_waterman(target, query)
  query_rc = revcomp(query)
  align_rc = smith_waterman(target, query_rc)
  if align_rc.score > align.score:
    return align_rc
  else:
    return align


def revcomp(seq):
  """Return the reverse complement of the input sequence.
  Leaves the input string unaltered."""
  return seq.translate(REVCOMP_TABLE)[::-1]


def revcomp_inplace(seq):
  """Convert the input sequence to its reverse complement.
  WARNING: This will alter the input string in-place!"""
  swalign.revcomp(seq)
