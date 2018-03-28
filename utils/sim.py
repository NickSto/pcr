#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import re
import os
import sys
import copy
import numpy
import bisect
import random
import string
import numbers
import tempfile
import argparse
import subprocess
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(root_dir)
from bfx import getreads

REVCOMP_TABLE = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
WGSIM_ID_REGEX = r'^(.+)_(\d+)_(\d+)_\d+:\d+:\d+_\d+:\d+:\d+_([0-9a-f]+)/[12]$'
ARG_DEFAULTS = {'read_len':100, 'frag_len':400, 'n_frags':1000, 'out_format':'fasta',
                'seq_error':0.001, 'pcr_error':0.001, 'cycles':25, 'indel_rate':0.15,
                'ext_rate':0.3, 'seed':None, 'invariant':'TGACT', 'bar_len':12, 'fastq_qual':'I'}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Simulate a duplex sequencing experiment."""

RAW_DISTRIBUTION = (
  #  0     1     2     3     4     5     6     7     8     9
  # Low singletons, but then constant drop-off. From pML113 (see 2015-09-28 report).
  #  0,  100,   36,   31,   27,   22,   17,   12,    7,  4.3,
  #2.4,  1.2,  0.6,  0.3,  0.2, 0.15,  0.1, 0.07, 0.05, 0.03,
  # High singletons, but then a second peak around 10. From Christine plasmid (2015-10-06 report).
  #    0,  100, 5.24, 3.67, 3.50, 3.67, 3.85, 4.02, 4.11, 4.20,
  # 4.17, 4.10, 4.00, 3.85, 3.69, 3.55, 3.38, 3.15, 2.92, 2.62,
  # 2.27, 2.01, 1.74, 1.56, 1.38, 1.20, 1.02, 0.85,
  # Same as above, but low singletons, 2's, and 3's (rely on errors to fill out those).
     0,    1,    2,    3, 3.50, 3.67, 3.85, 4.02, 4.11, 4.20,
  4.17, 4.10, 4.00, 3.85, 3.69, 3.55, 3.38, 3.15, 2.92, 2.62,
  2.27, 2.01, 1.74, 1.56, 1.38, 1.20, 1.02, 0.85,
)


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('ref', metavar='ref.fa', nargs='?',
    help='Reference sequence. Omit if giving --frag-file.')
  parser.add_argument('out1', type=argparse.FileType('w'),
    help='Write final mate 1 reads to this file.')
  parser.add_argument('out2', type=argparse.FileType('w'),
    help='Write final mate 2 reads to this file.')
  parser.add_argument('-o', '--out-format', choices=('fastq', 'fasta'))
  parser.add_argument('--stdout', action='store_true',
    help='Print interleaved output reads to stdout.')
  parser.add_argument('-m', '--mutations', type=argparse.FileType('w'),
    help='Write a log of the PCR and sequencing errors introduced to this file. Will overwrite any '
         'existing file at this path.')
  parser.add_argument('-b', '--barcodes', type=argparse.FileType('w'),
    help='Write a log of which barcodes were ligated to which fragments. Will overwrite any '
         'existing file at this path.')
  parser.add_argument('--frag-file',
    help='The path of the FASTQ file of fragments. If --ref is given, these will be generated with '
         'wgsim and kept (normally a temporary file is used, then deleted). Note: the file will be '
         'overwritten! If --ref is not given, then this should be a file of already generated '
         'fragments, and they will be used instead of generating new ones.')
  parser.add_argument('-Q', '--fastq-qual',
    help='The quality score to assign to all bases in FASTQ output. Give a character or PHRED '
         'score (integer). A PHRED score will be converted using the Sanger offset (33). Default: '
         '"%(default)s"')
  parser.add_argument('-S', '--seed', type=int,
    help='Random number generator seed. By default, a random, 32-bit seed will be generated and '
         'logged to stdout.')
  params = parser.add_argument_group('simulation parameters')
  params.add_argument('-n', '--n-frags', type=int,
    help='The number of original fragment molecules to simulate. The final number of reads will be '
         'this multiplied by the average number of reads per family. If you provide fragments with '
         '--frag-file, the script will still only read in the number specified here. Default: '
         '%(default)s')
  params.add_argument('-r', '--read-len', type=int,
    help='Default: %(default)s')
  params.add_argument('-f', '--frag-len', type=int,
    help='Default: %(default)s')
  params.add_argument('-s', '--seq-error', type=float,
    help='Sequencing error rate per base (0-1 proportion, not percent). Default: %(default)s')
  params.add_argument('-p', '--pcr-error', type=float,
    help='PCR error rate per base (0-1 proportion, not percent). Default: %(default)s')
  params.add_argument('-c', '--cycles', type=int,
    help='Number of PCR cycles to simulate. Default: %(default)s')
  params.add_argument('-i', '--indel-rate', type=float,
    help='Fraction of errors which are indels. Default: %(default)s')
  params.add_argument('-E', '--extension-rate', dest='ext_rate', type=float,
    help='Probability an indel is extended. Default: %(default)s')
  params.add_argument('-B', '--bar-len', type=int,
    help='Length of the barcodes to generate. Default: %(default)s')
  params.add_argument('-I', '--invariant',
    help='The invariant linker sequence between the barcode and sample sequence in each read. '
         'Default: %(default)s')

  # Parse and interpret arguments.
  args = parser.parse_args(argv[1:])
  assert args.ref or args.frag_file, 'You must provide either a reference or fragments file.'
  if args.seed is None:
    seed = random.randint(0, 2**31-1)
    sys.stderr.write('seed: {}\n'.format(seed))
  else:
    seed = args.seed
  random.seed(seed)
  if args.stdout:
    out1 = sys.stdout
    out2 = sys.stdout
  else:
    out1 = args.out1
    out2 = args.out2
  if isinstance(args.fastq_qual, numbers.Integral):
    assert args.fastq_qual >= 0, '--fastq-qual cannot be negative.'
    fastq_qual = chr(args.fastq_qual + 33)
  elif isinstance(args.fastq_qual, basestring):
    assert len(args.fastq_qual) == 1, '--fastq-qual cannot be more than a single character.'
    fastq_qual = args.fastq_qual
  else:
    raise AssertionError('--fastq-qual must be a positive integer or single character.')
  qual_line = fastq_qual * args.read_len

  invariant_rc = get_revcomp(args.invariant)

  # Create a temporary directory to do our work in. Then work inside a try so we can finally remove
  # the directory no matter what exceptions are encountered.
  tmpfile = tempfile.NamedTemporaryFile(prefix='wgdsim.frags.')
  tmpfile.close()
  try:
    # Step 1: Use wgsim to create fragments from the reference.
    if args.frag_file:
      frag_path = args.frag_file
    else:
      frag_path = tmpfile.name
    if args.ref and os.path.isfile(args.ref) and os.path.getsize(args.ref):
      #TODO: Check exit status
      #TODO: Check for wgsim on the PATH.
      # Set error and mutation rates to 0 to just slice sequences out of the reference without
      # modification.
      run_command('wgsim', '-e', '0', '-r', '0', '-d', '0', '-R', args.indel_rate, '-S', seed,
                  '-N', args.n_frags, '-X', args.ext_rate, '-1', args.frag_len,
                  args.ref, frag_path, os.devnull)

    # NOTE: Coordinates here are 0-based (0 is the first base in the sequence).
    extended_dist = extend_dist(RAW_DISTRIBUTION)
    proportional_dist = compile_dist(extended_dist)
    n_frags = 0
    for raw_fragment in getreads.getparser(frag_path, filetype='fastq'):
      n_frags += 1
      if n_frags > args.n_frags:
        break
      chrom, id_num, start, stop = parse_read_id(raw_fragment.id)
      barcode1 = get_rand_seq(args.bar_len)
      barcode2 = get_rand_seq(args.bar_len)
      barcode2_rc = get_revcomp(barcode2)
      raw_frag_full = barcode1 + args.invariant + raw_fragment.seq + invariant_rc + barcode2

      # Step 2: Determine how many reads to produce from each fragment.
      # - Use random.random() and divide the range 0-1 into segments of sizes proportional to
      #   the likelihood of each family size.
      # bisect.bisect() finds where an element belongs in a sorted list, returning the index.
      # proportional_dist is just such a sorted list, with values from 0 to 1.
      n_reads = bisect.bisect(proportional_dist, random.random())

      # Step 3: Introduce PCR errors.
      # - Determine the mutations and their frequencies.
      #   - Could get frequency from the cycle of PCR it occurs in.
      #     - Important to have PCR errors shared between reads.
      # - For each read, determine which mutations it contains.
      #   - Use random.random() < mut_freq.
      tree = get_good_pcr_tree(n_reads, args.cycles, 1000, max_diff=1)
      # Add errors to all children of original fragment.
      subtree1 = tree.get('child1')
      subtree2 = tree.get('child2')
      #TODO: Only simulate errors on portions of fragment that will become reads.
      add_pcr_errors(subtree1, '+', len(raw_frag_full), args.pcr_error, args.indel_rate, args.ext_rate)
      add_pcr_errors(subtree2, '-', len(raw_frag_full), args.pcr_error, args.indel_rate, args.ext_rate)
      apply_pcr_errors(tree, raw_frag_full)
      fragments = get_final_fragments(tree)
      add_mutation_lists(tree, fragments, [])

      # Step 4: Introduce sequencing errors.
      for fragment in fragments.values():
        for mutation in generate_mutations(args.read_len, args.seq_error, args.indel_rate,
                                           args.ext_rate):
          fragment['mutations'].append(mutation)
          fragment['seq'] = apply_mutation(mutation, fragment['seq'])

      # Print barcodes to log file.
      if args.barcodes:
        args.barcodes.write('{}-{}\t{}\t{}\n'.format(chrom, id_num, barcode1, barcode2_rc))
      # Print family.
      for frag_id in sorted(fragments.keys()):
        fragment = fragments[frag_id]
        read_id = '{}-{}-{}'.format(chrom, id_num, frag_id)
        # Print mutations to log file.
        if args.mutations:
          read1_muts = get_mutations_subset(fragment['mutations'], 0, args.read_len)
          read2_muts = get_mutations_subset(fragment['mutations'], 0, args.read_len, revcomp=True,
                                            seqlen=len(fragment['seq']))
          if fragment['strand'] == '-':
            read1_muts, read2_muts = read2_muts, read1_muts
          log_mutations(args.mutations, read1_muts, read_id+'/1', chrom, start, stop)
          log_mutations(args.mutations, read2_muts, read_id+'/2', chrom, start, stop)
        frag_seq = fragment['seq']
        read1_seq = frag_seq[:args.read_len]
        read2_seq = get_revcomp(frag_seq[len(frag_seq)-args.read_len:])
        if fragment['strand'] == '-':
          read1_seq, read2_seq = read2_seq, read1_seq
        if args.out_format == 'fasta':
          out1.write('>{}\n{}\n'.format(read_id, read1_seq))
          out2.write('>{}\n{}\n'.format(read_id, read2_seq))
        elif args.out_format == 'fastq':
          out1.write('@{}\n{}\n+\n{}\n'.format(read_id, read1_seq, qual_line))
          out2.write('@{}\n{}\n+\n{}\n'.format(read_id, read2_seq, qual_line))

  finally:
    try:
      os.remove(tmpfile.name)
    except OSError:
      pass


def run_command(*command, **kwargs):
  """Run a command and return the exit code.
  run_command('echo', 'hello')
  Will print the command to stderr before running, unless "silent" is set to True."""
  command_strs = map(str, command)
  if not kwargs.get('silent'):
    sys.stderr.write('$ '+' '.join(command_strs)+'\n')
  devnull = open(os.devnull, 'w')
  try:
    exit_status = subprocess.call(map(str, command), stderr=devnull)
  except OSError:
    exit_status = None
  finally:
    devnull.close()
  return exit_status


def extend_dist(raw_dist, exponent=1.25, min_prob=0.00001, max_len_mult=2):
  """Add an exponentially decreasing tail to the distribution.
  It takes the final value in the distribution and keeps dividing it by
  "exponent", adding each new value to the end. It will not add probabilities
  smaller than "min_prob" or extend the length of the list by more than
  "max_len_mult" times."""
  extended_dist = list(raw_dist)
  final_sum = sum(raw_dist)
  value = raw_dist[-1]
  value /= exponent
  while value/final_sum >= min_prob and len(extended_dist) < len(raw_dist)*max_len_mult:
    extended_dist.append(value)
    final_sum += value
    value /= exponent
  return extended_dist


def compile_dist(raw_dist):
  """Turn the human-readable list of probabilities defined at the top into
  proportional probabilities.
  E.g. [10, 5, 5] -> [0.5, 0.75, 1.0]"""
  proportional_dist = []
  final_sum = sum(raw_dist)
  current_sum = 0
  for magnitude in raw_dist:
    current_sum += magnitude
    proportional_dist.append(current_sum/final_sum)
  return proportional_dist


def parse_read_id(read_id):
  match = re.search(WGSIM_ID_REGEX, read_id)
  if match:
    chrom = match.group(1)
    start = match.group(2)
    stop = match.group(3)
    id_num = match.group(4)
  else:
    chrom, id_num, start, stop = read_id, None, None, None
  return chrom, id_num, start, stop


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
    node['strand'] = strand
    node['errors'] = list(generate_mutations(read_len, error_rate, indel_rate, extension_rate))
    add_pcr_errors(node.get('child2'), strand, read_len, error_rate, indel_rate, extension_rate)
    node = node.get('child1')


def apply_pcr_errors(subtree, seq):
  node = subtree
  while node:
    for error in node.get('errors', ()):
      seq = apply_mutation(error, seq)
    if 'child1' not in node:
      node['seq'] = seq
    apply_pcr_errors(node.get('child2'), seq)
    node = node.get('child1')


def get_final_fragments(tree):
  """Walk to the leaf nodes of the tree and get the post-PCR sequences of all the fragments.
  Returns a dict mapping fragment id number to a dict representing the fragment. Its only two keys
  are 'seq' (the final sequence) and 'strand' ('+' or '-')."""
  fragments = {}
  nodes = [tree]
  while nodes:
    node = nodes.pop()
    child1 = node.get('child1')
    if child1:
      nodes.append(child1)
    else:
      fragments[node['id']] = {'seq':node['seq'], 'strand':node['strand']}
    child2 = node.get('child2')
    if child2:
      nodes.append(child2)
  return fragments


def add_mutation_lists(subtree, fragments, mut_list1):
  """Compile the list of mutations that each fragment has undergone in PCR.
  To call from the root, give [] as "mut_list1" and a dict mapping all existing node id's to a dict
  as "fragments". Instead of returning the data, this will add a 'mutations' key to the dict for
  each fragment, mapping it to a list of PCR mutations that occurred in the lineage of the fragment,
  in chronological order."""
  node = subtree
  while node:
    mut_list1.extend(node.get('errors', ()))
    if 'child1' not in node:
      fragments[node['id']]['mutations'] = mut_list1
    if 'child2' in node:
      mut_list2 = copy.deepcopy(mut_list1)
      add_mutation_lists(node.get('child2'), fragments, mut_list2)
    node = node.get('child1')


def check_tree_balance(subtree):
  """Find all points in the tree where the cycles of sibling nodes is unequal, and
  return the maximum difference."""
  node = subtree
  if node:
    child1 = node.get('child1')
    child2 = node.get('child2')
    if child1 and child2:
      diff = abs(child1['cycle'] - child2['cycle'])
    else:
      diff = 0
    diff_child1 = check_tree_balance(child1)
    diff_child2 = check_tree_balance(child2)
    return max(diff, diff_child1, diff_child2)
  else:
    return 0


def get_good_pcr_tree(n_reads, n_cycles, max_tries, max_diff=1):
  """Return a single, balanced PCR tree from build_pcr_tree(), or fail if one cannot
  be found in max_tries.
  Compensate for bugs in build_pcr_tree() that sometimes result in multiple trees,
  or trees with siblings from different cycles."""
  tries = 0
  while tries <= max_tries:
    trees = build_pcr_tree(n_reads, n_cycles)
    if len(trees) == 1 and check_tree_balance(trees[0]) <= max_diff:
      return trees[0]
    tries += 1
  raise AssertionError('Could not generate a single, balanced tree! (tried {} times)'
                       .format(max_tries))


def build_pcr_tree(n_reads, n_cycles):
  """Create a simulated descent lineage of how all the final PCR fragments are related.
  Each node represents a fragment molecule at one stage of PCR. Each node is a dict containing the
  fragment's children (other nodes) ('child1' and 'child2'), the PCR cycle number ('cycle'), and,
  at the leaves, a unique id number for each final fragment.
  Returns a list of root nodes. Usually there will only be one, but about 1-3% of the time it fails
  to unify the subtrees and results in a broken tree.
  """
  #TODO: Make it always return a single tree.
  # Begin a branch for each of the fragments. These are the leaf nodes. We'll work backward from
  # these, simulating the points at which they share ancestors, eventually coalescing into the
  # single (root) ancestor.
  branches = []
  for frag_id in range(n_reads):
    branches.append({'cycle':n_cycles-1, 'id':frag_id})
  # Build up all the branches in parallel. Start from the second-to-last PCR cycle.
  for cycle in reversed(range(n_cycles-1)):
    # Probability of 2 fragments sharing an ancestor at cycle c is 1/2^c.
    prob = 1/2**cycle
    frag_i = 0
    while frag_i < len(branches):
      current_root = branches[frag_i]
      # Does the current fragment share this ancestor with any of the other fragments?
      # numpy.random.binomial() is a fast way to simulate going through every other fragment and
      # asking if random.random() < prob.
      shared = numpy.random.binomial(len(branches)-1, prob)
      if shared == 0:
        # No branch point here. Just add another level to the lineage.
        branches[frag_i] = {'cycle':cycle, 'child1':current_root}
      else:
        # Pick a random other fragment to share this ancestor with.
        # Make a list of candidates to pick from.
        candidates = []
        for candidate_i, candidate in enumerate(branches):
          # Don't include ourselves.
          if candidate is current_root:
            continue
          # If it's at a cycle above us and it already has a child, skip it.
          if candidate['cycle'] == cycle and candidate.get('child2'):
            continue
          candidates.append(candidate_i)
        if candidates:
          relative_i = random.choice(candidates)
          relative = branches[relative_i]
          # Have we already passed this fragment on this cycle?
          if relative['cycle'] == cycle:
            # If we've already passed it, we're looking at the fragment's parent. We want the child.
            relative = relative['child1']
          # Join the lineages of our current fragment and the relative to a new parent.
          #TODO: Sometimes, we end up matching up subtrees of different depths. But the discrepancy
          #      is rarely greater than 1. Figure out why.
          # assert abs(current_root['cycle'] - relative['cycle']) < 3, ('cycle: {}, current_root: {},'
          #   ' relative: {}, frag_i: {}, relative_i: {}, branches: {}, candidates: {}, shared: {}'
          #   .format(cycle, current_root['cycle'], relative['cycle'], frag_i, relative_i,
          #           len(branches), len(candidates), shared))
          branches[frag_i] = {'cycle':cycle, 'child1':current_root, 'child2':relative}
          # Remove the relative from the list of lineages.
          del(branches[relative_i])
          if relative_i < frag_i:
            frag_i -= 1
      frag_i += 1
  return branches


def get_depth(tree):
  depth = 0
  node = tree
  while node:
    depth += 1
    node = node.get('child1')
  return depth


def convert_tree(tree_orig):
  # Let's operate on a copy only.
  tree = copy.deepcopy(tree_orig)
  # Turn the tree vertical.
  tree['line'] = 1
  tree['children'] = 0
  levels = [[tree]]
  done = False
  while not done:
    last_level = levels[-1]
    this_level = []
    done = True
    for node in last_level:
      for child_name in ('child1', 'child2'):
        child = node.get(child_name)
        if child:
          done = False
          child['parent'] = node
          child['branch'] = child['parent']['branch']
          if child_name == 'child2':
            child['branch'] += 1
          this_level.append(child)
    this_level.sort(key=lambda node: node['branch'])
    levels.append(this_level)
  return levels


def print_levels(levels):
  last_level = []
  for level in levels:
    for node in level:
      child = 1
      for parent in last_level:
        if parent.get('child2') is node:
          child = 2
      if child == 1:
        sys.stdout.write('| ')
      else:
        sys.stdout.write('\ ')
    last_level = level
    print()


def label_branches(tree):
  """Label each vertical branch (line of 'child1's) with an id number."""
  counter = 1
  tree['branch'] = counter
  nodes = [tree]
  while nodes:
    node = nodes.pop(0)
    child1 = node.get('child1')
    if child1:
      child1['branch'] = node['branch']
      nodes.append(child1)
    child2 = node.get('child2')
    if child2:
      counter += 1
      child2['branch'] = counter
      nodes.append(child2)


def print_tree(tree_orig):
  # We "write" strings to an output buffer instead of directly printing, so we can post-process the
  # output. The buffer is a matrix of cells, each holding a string representing one element.
  lines = [[]]
  # Let's operate on a copy only.
  tree = copy.deepcopy(tree_orig)
  # Add some bookkeeping data.
  label_branches(tree)
  tree['level'] = 0
  branches = [tree]
  while branches:
    line = lines[-1]
    branch = branches.pop()
    level = branch['level']
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
      child2 = node.get('child2')
      if child2:
        child2['level'] = node['level'] + 1
        branches.append(child2)
      parent = node
      node = node.get('child1')
      if node:
        node['level'] = parent['level'] + 1
      else:
        line.append(' {}'.format(parent['branch']))
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


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == '__main__':
  sys.exit(main(sys.argv))
