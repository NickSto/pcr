#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
import sys
import copy
import random
import string
import logging
from bfx import getreads

PY3 = sys.version_info[0] >= 3

if PY3:
  REVCOMP_TABLE = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
else:
  REVCOMP_TABLE = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
USAGE = """%(prog)s [options]"""
DESCRIPTION = """Simulate a PCR experiment."""


def main(argv):
  return 0


#TODO: Clean up "mutation" vs "error" terminology.
def generate_mutations(seq_len, error_rate, indel_rate, extension_rate):
  """Generate all the mutations that occur over the length of a sequence."""
  i = 0
  while i <= seq_len:
    if random.random() < error_rate:
      mtype, alt = make_mutation(indel_rate, extension_rate)
      # Allow mutation after the last base only if it's an insertion.
      if i < seq_len or mtype == 'ins':
        yield {'coord':i, 'type':mtype, 'alt':alt}
      # Compensate for length variations to keep i tracking the original read's base coordinates.
      if mtype == 'ins':
        i += len(alt)
      elif mtype == 'del':
        i -= alt
    i += 1


def make_mutation(indel_rate, extension_rate):
  """Simulate a random mutation."""
  # Is it an indel?
  rand = random.random()
  if rand < indel_rate:
    # Is it an insertion or deletion? Decide, then initialize it.
    # Re-use the random number from above. Just check if it's in the lower or upper half of the
    # range from 0 to indel_rate.
    if rand < indel_rate/2:
      mtype = 'del'
      alt = 1
    else:
      mtype = 'ins'
      alt = get_rand_base()
    # Extend the indel as long as the extension rate allows.
    while random.random() < extension_rate:
      if mtype == 'ins':
        alt += get_rand_base()
      else:
        alt += 1
  else:
    # What is the new base for the SNV?
    mtype = 'snv'
    alt = get_rand_base()
  return mtype, alt


def get_rand_base(bases='ACGT'):
  return random.choice(bases)


def get_rand_seq(seq_len):
  return ''.join([get_rand_base() for i in range(seq_len)])


def get_revcomp(seq):
  return seq.translate(REVCOMP_TABLE)[::-1]


def apply_mutation(mut, seq):
  i = mut['coord']
  if mut['type'] == 'snv':
    # Replace the base at "coord".
    new_seq = seq[:i] + mut['alt'] + seq[i+1:]
  else:
    # Indels are handled by inserting or deleting bases starting *before* the base at "coord".
    # This goes agains the VCF convention, but it allows deleting the first and last base, as well
    # as inserting before and after the sequence without as much special-casing.
    if mut['type'] == 'ins':
      # Example: 'ACGTACGT' + ins 'GC' at 4 = 'ACGTGCACGT'
      new_seq = seq[:i] + mut['alt'] + seq[i:]
    else:
      # Example: 'ACGTACGT' + del 2 at 4 = 'ACGTGT'
      new_seq = seq[:i] + seq[i+mut['alt']:]
  return new_seq


def get_mutations_subset(mutations_old, start, length, revcomp=False, seqlen=None):
  """Get a list of the input mutations which are within a certain region.
  The output list maintains the order in the input list, only filtering out
  mutations outside the specified region.
  "start" is the start of the region (0-based). If revcomp, this start should be
  in the coordinate system of the reverse-complemented sequence.
  "length" is the length of the region.
  "revcomp" causes the mutations to be converted to their reverse complements, and
  the "start" to refer to the reverse complement sequence's coordinates. The order
  of the mutations is unchanged, though.
  "seqlen" is the length of the sequence the mutations occurred in. This is only
  needed when revcomp is True, to convert coordinates to the reverse complement
  coordinate system."""
  stop = start + length
  mutations_new = []
  for mutation in mutations_old:
    if revcomp:
      mutation = get_mutation_revcomp(mutation, seqlen)
    if start <= mutation['coord'] < stop:
      mutations_new.append(mutation)
    elif mutation['coord'] == stop and mutation['type'] == 'ins':
      # Allow insertions at the last coordinate.
      mutations_new.append(mutation)
  return mutations_new


def get_mutation_revcomp(mut, seqlen):
  """Convert a mutation to its reverse complement.
  "seqlen" is the length of the sequence the mutation is being applied to. Needed
  to convert the coordinate to a coordinate system starting at the end of the
  sequence."""
  mut_rc = {'type':mut['type']}
  if mut['type'] == 'snv':
    mut_rc['coord'] = seqlen - mut['coord'] - 1
    mut_rc['alt'] = get_revcomp(mut['alt'])
  elif mut['type'] == 'ins':
    mut_rc['coord'] = seqlen - mut['coord']
    mut_rc['alt'] = get_revcomp(mut['alt'])
  elif mut['type'] == 'del':
    mut_rc['coord'] = seqlen - mut['coord'] - mut['alt']
    mut_rc['alt'] = mut['alt']
  return mut_rc


def log_mutations(mutfile, mutations, read_id, chrom, start, stop):
  for mutation in mutations:
    mutfile.write('{read_id}\t{chrom}\t{start}\t{stop}\t{coord}\t{type}\t{alt}\n'
                  .format(read_id=read_id, chrom=chrom, start=start, stop=stop, **mutation))


def add_pcr_errors(subtree, strand, read_len, error_rate, indel_rate, extension_rate):
  """Add simulated PCR errors to a node in a tree and all its descendants."""
  # Note: The errors are intended as "errors made in creating this molecule", so don't apply this to
  # the root node, since that is supposed to be the original, unaltered molecule.
  # Go down the subtree and simulate errors in creating each fragment.
  # Process all the first-child descendants of the original node in a loop, and recursively call
  # this function to process all second children.
  node = subtree
  while node:
    node.strand = strand
    node.errors = list(generate_mutations(read_len, error_rate, indel_rate, extension_rate))
    add_pcr_errors(node.child2, strand, read_len, error_rate, indel_rate, extension_rate)
    node = node.child1


def apply_pcr_errors(subtree, seq):
  node = subtree
  while node:
    for error in node.errors:
      seq = apply_mutation(error, seq)
    if not node.child1:
      node.seq = seq
    apply_pcr_errors(node.child2, seq)
    node = node.child1


def get_final_fragments(tree):
  """Walk to the leaf nodes of the tree and get the post-PCR sequences of all the fragments.
  Returns a dict mapping fragment id number to a dict representing the fragment. Its only two keys
  are 'seq' (the final sequence) and 'strand' ('+' or '-')."""
  fragments = {}
  nodes = [tree]
  while nodes:
    node = nodes.pop()
    if node.child1:
      nodes.append(node.child1)
    else:
      fragments[node.leaf_id] = {'seq':node.seq, 'strand':node.strand}
    if node.child2:
      nodes.append(node.child2)
  return fragments


def add_mutation_lists(subtree, fragments, mut_list1):
  """Compile the list of mutations that each fragment has undergone in PCR.
  To call from the root, give [] as "mut_list1" and a dict mapping all existing node id's to a dict
  as "fragments". Instead of returning the data, this will add a 'mutations' key to the dict for
  each fragment, mapping it to a list of PCR mutations that occurred in the lineage of the fragment,
  in chronological order."""
  node = subtree
  while node:
    mut_list1.extend(node.errors)
    if not node.child1:
      fragments[node.leaf_id]['mutations'] = mut_list1
    if node.child2:
      mut_list2 = copy.deepcopy(mut_list1)
      add_mutation_lists(node.child2, fragments, mut_list2)
    node = node.child1


def build_good_pcr_tree(n_cycles, final_reads, efficiency_decline, max_tries=100):
  """Try to get a PCR tree with final_reads leaf nodes.
  Retries up to max_tries times until success. If unsuccessful, returns the tree with the most leaf
  nodes."""
  best_tree = None
  success = False
  tries = 0
  while not success and tries < max_tries:
    tries += 1
    tree = build_pcr_tree(n_cycles, final_reads, efficiency_decline)
    if best_tree is None or tree.leaves > best_tree.leaves:
      best_tree = tree
    if best_tree.leaves == final_reads:
      success = True
  logging.debug('Tries to get a good PCR tree: {}. Success: {}'
                .format(tries, best_tree.leaves == final_reads))
  return best_tree


def build_pcr_tree(n_cycles, final_reads, efficiency_decline):
  """Create a simulated descent lineage of how all the final PCR fragments are related.
  Each node represents a fragment molecule at one stage of PCR.
  Returns the root node.
  n_cycles is the number of PCR cycles to simulate.
  final_reads is the target number of observed PCR fragments (the duplex family size).
    This is not guaranteed to return a tree with this many leaf nodes. Lower n_cycles, higher
    final_reads, and higher efficiency_decline will make it more likely to fail to reach the
    intended final_reads target (always resulting in fewer than intended).
    For example, with 25 cycles and an efficiency_decline of 1.05, a final_reads target of 2 will
    fail to be met about 0.015% of the time. A target of 3 fails about 0.04%, and for 20, it's 2%.
  efficiency_decline is the rate at which replication efficiency declines with each cycle.
    That is, the efficiency is divided by this number every cycle. So it should usually be a little
    greater than 1.
  """
  efficiency = 2
  root = Node(skipped_branches=0, taken_branches=0, reads_left=final_reads)
  skipped_buffer = 0
  leaves = [root]
  for i in range(n_cycles):
    new_leaves = []
    for leaf in leaves:
      if random.random() * 2 > efficiency:
        # Molecule didn't replicate this cycle.
        new_leaves.append(leaf)
        continue
      # Decide which branch each final read follows.
      child1s = 0
      child2s = 0
      for read in range(leaf.reads_left):
        if random.random() < 0.5:
          child1s += 1
        else:
          child2s += 1
      if child1s == 0:
        child1s, child2s = child2s, child1s
      if child1s:
        leaf.child1 = Node(parent=leaf, reads_left=child1s)
        new_leaves.append(leaf.child1)
      if child2s:
        leaf.child2 = Node(parent=leaf, reads_left=child2s)
        new_leaves.append(leaf.child2)
        root.taken_branches += 1
        root.skipped_branches += skipped_buffer
        skipped_buffer = 0
      else:
        skipped_buffer += 1
    leaves = new_leaves
    efficiency = efficiency / efficiency_decline
  for i, leaf in enumerate(leaves):
    leaf.leaf_id = i
  return root


class Node(object):
  __slots__ = ('parent', 'child1', 'child2', 'seq', 'errors', 'strand', 'branch', 'leaf_id',
               'skipped_branches', 'taken_branches', '_level', '_leaves', 'reads_left')

  def __init__(self, parent=None, child1=None, child2=None, seq=None, leaves=None,
               skipped_branches=None, taken_branches=None, reads_left=None):
    self.parent = parent
    self.child1 = child1
    self.child2 = child2
    self.seq = seq
    self.skipped_branches = skipped_branches
    self.taken_branches = taken_branches
    self.branch = None
    self._leaves = leaves
    self._level = None
    self.reads_left = reads_left
    self.errors = []
    self.strand = None
    self.leaf_id = None

  @property
  def leaves(self):
    """How many leaves are (at or) below this node?
    If this is a leaf, return 1."""
    if self._leaves is None:
      self._leaves = 0
      if self.child1:
        self._leaves += self.child1.leaves
      if self.child2:
        self._leaves += self.child2.leaves
      if not (self.child1 or self.child2):
        self._leaves = 1
    return self._leaves

  @leaves.setter
  def leaves(self, value):
    self._leaves = value

  @property
  def level(self):
    """How long is our branch, from the root to this node?
    The root's level is 0."""
    if self._level is None:
      if self.parent is None:
        self._level = 0
      else:
        self._level = self.parent.level + 1
    return self._level

  @property
  def compactness(self):
    total = self.taken_branches + self.skipped_branches
    if total > 0:
      return self.taken_branches / (self.taken_branches + self.skipped_branches)
    elif self.leaves == 1:
      return 1

  def print_tree(self):
    # We "write" strings to an output buffer instead of directly printing, so we can post-process
    # the output. The buffer is a matrix of cells, each holding a string representing one element.
    lines = [[]]
    # Add some bookkeeping data.
    self.label_branches()
    branches = [self]
    while branches:
      line = lines[-1]
      branch = branches.pop()
      level = branch.level
      while level > 0:
        line.append('  ')
        level -= 1
      node = branch
      while node:
        # Is it the root node? (Have we written anything yet?)
        if lines[0]:
          # Are we at the start of the line? (Is it only spaces so far?)
          if line[-1] == '  ':
            line.append('\-')
          elif line[-1].endswith('-'):
            line.append('=-')
        else:
          line.append('*-')
        if node.child2:
          branches.append(node.child2)
        parent = node
        node = node.child1
        if not node:
          line.append(' {}'.format(parent.branch))
          lines.append([])
    # Post-process output: Add lines connecting branches to parents.
    x = 0
    done = False
    while not done:
      # Draw vertical lines upward from branch points.
      drawing = False
      for line in reversed(lines):
        done = True
        if x < len(line):
          done = False
          cell = line[x]
          if cell == '\-':
            drawing = True
          elif cell == '  ' and drawing:
            line[x] = '| '
          elif cell == '=-' and drawing:
            drawing = False
      x += 1
    # Print the final output.
    for line in lines:
      print(''.join(line))

  def label_branches(self):
    """Label each vertical branch (line of 'child1's) with an id number."""
    counter = 1
    self.branch = counter
    nodes = [self]
    while nodes:
      node = nodes.pop(0)
      if node.child1:
        node.child1.branch = node.branch
        nodes.append(node.child1)
      if node.child2:
        counter += 1
        node.child2.branch = counter
        nodes.append(node.child2)


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')

if __name__ == '__main__':
  sys.exit(main(sys.argv))
