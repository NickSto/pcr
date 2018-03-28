#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import re
import os
import sys
import copy
import bisect
import random
import string
import logging
import numbers
import tempfile
import argparse
import subprocess
script_path = os.path.realpath(__file__)
root_dir = os.path.dirname(os.path.dirname(script_path))
sys.path.append(root_dir)
from bfx import getreads

PY3 = sys.version_info[0] == 3

if PY3:
  REVCOMP_TABLE = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
else:
  REVCOMP_TABLE = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
WGSIM_ID_REGEX = r'^(.+)_(\d+)_(\d+)_\d+:\d+:\d+_\d+:\d+:\d+_([0-9a-f]+)/[12]$'
USAGE = """%(prog)s [options] ref.fa [--frag-file frags.fq] -1 reads_1.fa -2 reads_2.fa
or     %(prog)s [options] ref.fa --stdout > reads.fa
or     %(prog)s [options] --frag-file frags.fq -1 reads_1.fa -2 reads_2.fa"""
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


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION, usage=USAGE)
  io = parser.add_argument_group('I/O')
  io.add_argument('ref', metavar='ref.fa', nargs='?',
    help='Reference sequence. Omit if giving --frag-file.')
  io.add_argument('-1', '--reads1', type=argparse.FileType('w'),
    help='Write final mate 1 reads to this file.')
  io.add_argument('-2', '--reads2', type=argparse.FileType('w'),
    help='Write final mate 2 reads to this file.')
  io.add_argument('-o', '--out-format', choices=('fastq', 'fasta'), default='fasta',
    help='Default: %(default)s')
  io.add_argument('--stdout', action='store_true',
    help='Print interleaved output reads to stdout.')
  io.add_argument('-m', '--mutations', type=argparse.FileType('w'),
    help='Write a log of the PCR and sequencing errors introduced to this file. Will overwrite any '
         'existing file at this path.')
  io.add_argument('-b', '--barcodes', type=argparse.FileType('w'),
    help='Write a log of which barcodes were ligated to which fragments. Will overwrite any '
         'existing file at this path.')
  io.add_argument('--frag-file',
    help='The path of the FASTQ file of fragments. If --ref is given, these will be generated with '
         'wgsim and kept (normally a temporary file is used, then deleted). Note: the file will be '
         'overwritten! If --ref is not given, then this should be a file of already generated '
         'fragments, and they will be used instead of generating new ones.')
  io.add_argument('-Q', '--fastq-qual', default='I',
    help='The quality score to assign to all bases in FASTQ output. Give a character or PHRED '
         'score (integer). A PHRED score will be converted using the Sanger offset (33). Default: '
         '"%(default)s"')
  params = parser.add_argument_group('Simulation Parameters')
  params.add_argument('-n', '--n-frags', type=int, default=1000,
    help='The number of original fragment molecules to simulate. The final number of reads will be '
         'this multiplied by the average number of reads per family. If you provide fragments with '
         '--frag-file, the script will still only read in the number specified here. Default: '
         '%(default)s')
  params.add_argument('-r', '--read-len', type=int, default=100,
    help='Default: %(default)s')
  params.add_argument('-f', '--frag-len', type=int, default=400,
    help='Default: %(default)s')
  params.add_argument('-s', '--seq-error', type=float, default=0.001,
    help='Sequencing error rate per base (0-1 proportion, not percent). Default: %(default)s')
  params.add_argument('-p', '--pcr-error', type=float, default=0.001,
    help='PCR error rate per base (0-1 proportion, not percent). Default: %(default)s')
  params.add_argument('-c', '--cycles', type=int, default=25,
    help='Number of PCR cycles to simulate. Default: %(default)s')
  params.add_argument('-e', '--efficiency-decline', type=float, default=1.05,
    help='Rate at which the PCR replication efficiency declines.')
  params.add_argument('-i', '--indel-rate', type=float, default=0.15,
    help='Fraction of errors which are indels. Default: %(default)s')
  params.add_argument('-E', '--extension-rate', dest='ext_rate', type=float, default=0.3,
    help='Probability an indel is extended. Default: %(default)s')
  params.add_argument('-B', '--bar-len', type=int, default=12,
    help='Length of the barcodes to generate. Default: %(default)s')
  params.add_argument('-I', '--invariant', default='TGACT',
    help='The invariant linker sequence between the barcode and sample sequence in each read. '
         'Default: %(default)s')
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  log.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  log.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  log.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  misc = parser.add_argument_group('Misc')
  misc.add_argument('-S', '--seed', type=int,
    help='Random number generator seed. By default, a random, 32-bit seed will be generated and '
         'logged to stdout.')
  return parser


def main(argv):
  # Parse and interpret arguments.
  parser = make_argparser()
  args = parser.parse_args(argv[1:])
  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
  tone_down_logger()
  if not (args.ref or args.frag_file):
    parser.print_usage()
    fail('You must provide either a reference or fragments file.')
  if args.ref:
    if not os.path.isfile(args.ref):
      fail('Error: reference file {!r} not found.'.format(args.ref))
    if not os.path.getsize(args.ref):
      fail('Error: reference file {!r} empty (0 bytes).'.format(args.ref))
  else:
    if not (args.reads1 and args.reads2):
      fail('Error: must provide output --reads1 and --reads2 files.')
  if args.seed is None:
    seed = random.randint(0, 2**31-1)
    logging.info('seed: {}\n'.format(seed))
  else:
    seed = args.seed
  random.seed(seed)
  if args.stdout:
    reads1 = sys.stdout
    reads2 = sys.stdout
  else:
    reads1 = args.reads1
    reads2 = args.reads2
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
  tmpfile = tempfile.NamedTemporaryFile(prefix='wgdsim.frags.', delete=False)
  tmpfile.close()
  try:
    # Step 1: Use wgsim to create fragments from the reference.
    if args.frag_file:
      frag_path = args.frag_file
    else:
      frag_path = tmpfile.name
    if args.ref:
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
      #TODO: Vary the size of the fragment.
      #      Could add ~100bp to frag_len arg to wgsim, then randomly select a subsequence here.
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
      tree = build_good_pcr_tree(args.cycles, n_reads, args.efficiency_decline, 1000)
      # Add errors to all children of original fragment.
      subtree1 = tree.child1
      subtree2 = tree.child2
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
          reads1.write('>{}\n{}\n'.format(read_id, read1_seq))
          reads2.write('>{}\n{}\n'.format(read_id, read2_seq))
        elif args.out_format == 'fastq':
          reads1.write('@{}\n{}\n+\n{}\n'.format(read_id, read1_seq, qual_line))
          reads2.write('@{}\n{}\n+\n{}\n'.format(read_id, read2_seq, qual_line))

  finally:
    try:
      os.remove(tmpfile.name)
    except OSError:
      pass


def run_command(*command, **kwargs):
  """Run a command and return the exit code.
  run_command('echo', 'hello')
  If "echo" keyword argument is set to True, this will print the command to stdout first."""
  command_strs = map(str, command)
  command_line = '$ '+' '.join(command_strs)+'\n'
  logging.info(command_line)
  if kwargs.get('echo'):
    print(command_line)
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


def build_pcr_tree_old(n_cycles, efficiency_decline, branch_rate):
  """Create a simulated descent lineage of how all the final PCR fragments are related.
  Each node represents a fragment molecule at one stage of PCR.
  Returns the root node.
  efficiency_decline is the rate at which replication efficiency declines with each cycle.
    That is, the efficiency is divided by this number every cycle. So it should usually be a little
    greater than 1.
  branch_rate is the probability of the tree branching in the first cycle (it declines thereafter).
    More specifically, it's the probability that two of the final, observed reads are derived from
    different children of the original molecule. This probability is halved every cycle.
    Note to self: This might be easy to derive from the desired number of final reads. E.g. if we
    want 12 reads, what is 1 - the probability that all 12 of the reads are derived from only one of
    the first two children?
    Other note to self: it looks like this formulation, with the probability halving every time,
    might not work out. Even with an initial probability of 1, the median number of final fragments
    is 4. Increasing it beyond 4 can get a higher median, but that breaks the theory of this value.
  """
  # Note: Remember that each child includes one of the strands of the parent. Only mutate an
  #       appropriate proportion (half, but check that) of the children.
  efficiency = 2
  root = Node(skipped_branches=0, taken_branches=0)
  skipped_buffer = 0
  leaves = [root]
  for i in range(n_cycles):
    # print('Cycle {}:'.format(i+1))
    # print('  branch_rate: {}'.format(branch_rate))
    # print('  efficiency:  {}'.format(efficiency))
    new_leaves = []
    for leaf in leaves:
      if random.random() * 2 > efficiency:
        # Molecule didn't replicate this cycle.
        new_leaves.append(leaf)
        continue
      leaf.child1 = Node(parent=leaf)
      new_leaves.append(leaf.child1)
      if random.random() < branch_rate:
        root.taken_branches += 1
        root.skipped_branches += skipped_buffer
        skipped_buffer = 0
        leaf.child2 = Node(parent=leaf)
        new_leaves.append(leaf.child2)
      else:
        skipped_buffer += 1
    leaves = new_leaves
    #TODO: Determine how exactly the branch rate should decline.
    branch_rate = branch_rate / 2
    efficiency = efficiency / efficiency_decline
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
