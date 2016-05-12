#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import logging
import argparse

ARG_DEFAULTS = {'sam':sys.stdin, 'qual':20, 'pos':2, 'dist':1, 'log':sys.stderr,
                'log_level':logging.ERROR}
USAGE = "%(prog)s [options]"
DESCRIPTION = """Correct barcodes using an alignment of all barcodes to themselves. Reads the
alignment in SAM format and """


def main(argv):

  parser = argparse.ArgumentParser(description=DESCRIPTION)
  parser.set_defaults(**ARG_DEFAULTS)

  parser.add_argument('families', type=argparse.FileType('r'),
    help='The sorted output of make-barcodes.awk. The important part is that it\'s a tab-delimited '
         'file with at least one column, which is the barcode sequence, and it must be sorted in '
         'the same order as the "reads" in the SAM file.')
  parser.add_argument('sam', type=argparse.FileType('r'), nargs='?',
    help='Barcode alignment, in SAM format. Omit to read from stdin. The read names must be '
         'integers.')
  parser.add_argument('-c', '--add-column', action='store_true',
    help='Include corrected barcodes as an additional column in ')
  parser.add_argument('-d', '--dist', type=int,
    help='NM edit distance threshold. Default: %(default)s')
  parser.add_argument('-q', '--qual', type=int,
    help='MAPQ threshold. Default: %(default)s')
  parser.add_argument('-p', '--pos', type=int,
    help='POS tolerance. Alignments will be ignored if abs(POS - 1) is greater than this value. '
         'Set to greater than the barcode length for no threshold. Default: %(default)s')
  parser.add_argument('-l', '--log', type=argparse.FileType('w'),
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  parser.add_argument('-D', '--debug', dest='log_level', action='store_const', const=logging.DEBUG,
    help='Print debug messages (very verbose).')

  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.log_level, format='%(message)s')
  tone_down_logger()

  sys.stderr.write('Starting to build groups from SAM alignment..\n')
  groups = read_alignment(args.sam, args.pos, args.qual, args.dist)

  #TODO: Gather orders (ab, ba) too. Erroneous barcodes can also end up with the wrong order.
  sys.stderr.write('Starting to read barcode sequences from families file..\n')
  votes, barcodes = count_barcodes(args.families)

  #TODO: Instead of going through all the data again, instead at this point we could stream through
  #      the lines in the families file, using the "member_to_group map" to look up each barcode in
  #      "groups" and determine the correction. At that point we could compute all the corrections
  #      for the group for later (essentially compute the "corrections" table on the fly).
  sys.stderr.write('Starting to build corrections table from collected data..\n')
  corrections = make_correction_table(groups, votes, barcodes)

  sys.stderr.write('Starting to write corrected output..\n')
  line_num = 0
  barcode_num = 0
  barcode_last = None
  with open(args.families.name) as families:
    for line in families:
      line_num += 1
      fields = line.rstrip('\r\n').split('\t')
      raw_barcode = fields[0]
      if raw_barcode != barcode_last:
        barcode_num += 1
        barcode_last = raw_barcode
      try:
        corrected_barcode = corrections[raw_barcode]
      except KeyError:
        logging.debug('Barcode number {} not in corrections table. Families file line {}, barcode: '
                      '{}'.format(barcode_num, line_num, raw_barcode))
        corrected_barcode = raw_barcode
      if args.add_column:
        fields[1:1] = [corrected_barcode]
      else:
        fields[0] = corrected_barcode
      print(*fields, sep='\t')


def read_alignment(sam, pos_thres, qual_thres, dist_thres):
  # Group reads (barcodes) into sets of reads linked by alignments to each other.
  # Each group of reads is a dict mapping read names to read sequences. "groups" is a list of these
  # dicts.
  groups = []
  # "member_to_group" is a dict mapping read names to their group's index in "groups".
  member_to_group = {}
  line_num = 0
  for line in sam:
    line_num += 1
    if line.startswith('@'):
      logging.debug('Header line ({})'.format(line_num))
      continue
    fields = line.split('\t')
    # if fields[0] == '18190' or fields[2] == '18190':
    #   logging.getLogger().setLevel(logging.DEBUG)
    logging.debug('Ref {}, read {} ({}):'.format(fields[0], fields[2], fields[9]))
    try:
      read_name = int(fields[0])
      ref_name = int(fields[2])
    except ValueError:
      if fields[2] == '*':
        logging.debug('\tRead unmapped (reference == "*")')
        continue
      else:
        logging.error('Non-integer read name(s) on line {}: "{}", "{}".'
                      .format(line_num, read_name, ref_name))
        raise
    # Apply alignment quality filters.
    try:
      flags = int(fields[1])
      pos = int(fields[3])
      mapq = int(fields[4])
    except ValueError:
      logging.warn('\tNon-integer flag ({}), pos ({}), or mapq ({})'
                   .format(fields[1], fields[3], fields[4]))
      continue
    if flags & 4:
      # Read unmapped.
      logging.debug('\tRead unmapped (flag & 4 == True)')
      continue
    if abs(pos - 1) > pos_thres:
      logging.debug('\tAlignment failed pos filter: abs({} - 1) > {}'.format(pos, pos_thres))
      continue
    if mapq < qual_thres:
      logging.debug('\tAlignment failed mapq filter: {} > {}'.format(mapq, qual_thres))
      continue
    nm = None
    for tag in fields[11:]:
      if tag.startswith('NM:i:'):
        try:
          nm = int(tag[5:])
        except ValueError:
          logging.error('Invalid NM tag "{}" on line {}.'.format(tag, line_num))
          raise
        break
    assert nm is not None, line_num
    if nm > dist_thres:
      logging.debug('\tAlignment failed NM distance filter: {} > {}'.format(nm, dist_thres))
      continue
    # It's a good alignment. Add it to a group.
    group_i_ref = member_to_group.get(ref_name)
    group_i_read = member_to_group.get(read_name)
    if group_i_ref is None:
      group_i = group_i_read
      logging.debug('\tgroup_i_ref  is None, using {}.'.format(group_i_read))
    elif group_i_read is None:
      group_i = group_i_ref
      logging.debug('\tgroup_i_read is None, using {}.'.format(group_i_ref))
    elif group_i_ref == group_i_read:
      group_i = group_i_ref
      logging.debug('\tgroup_i_ref == group_i_read ({}).'.format(group_i_ref))
    else:
      # Both the read and ref are already in a group, but they're different groups. We need to join
      # them.
      logging.debug('\tgroup_i_ref  is {}, group_i_read is {}.'.format(group_i_ref, group_i_read))
      logging.debug('\tJoining groups {} and {}.'.format(group_i_ref, group_i_read))
      join_groups(groups, member_to_group, ref_name, read_name)
      group_i = member_to_group[ref_name]
    if group_i is None:
      # No group exists yet. Create one and add it.
      new_group = set((ref_name, read_name))
      groups.append(new_group)
      group_i = len(groups) - 1
      member_to_group[ref_name] = group_i
      member_to_group[read_name] = group_i
      logging.debug('\tAdding new group ({}): {} and {}.'.format(group_i, ref_name, read_name))
    else:
      # Add these reads to the group.
      logging.debug('\tAdding {} and {} to group {}.'.format(ref_name, read_name, group_i))
      group = groups[group_i]
      group.add(ref_name)
      group.add(read_name)
      member_to_group[ref_name] = group_i
      member_to_group[read_name] = group_i
    # logging.getLogger().setLevel(logging.ERROR)
  return groups


def join_groups(groups, member_to_group, ref_name, read_name):
  group_i_ref = member_to_group[ref_name]
  group_i_read = member_to_group[read_name]
  logging.debug('\t\tmember_to_group[{}] == {}'.format(ref_name, group_i_ref))
  logging.debug('\t\tmember_to_group[{}] == {}'.format(read_name, group_i_read))
  group_ref = groups[group_i_ref]
  group_read = groups[group_i_read]
  # Get a union of all elements in both groups.
  group_union = group_ref.union(group_read)
  # Put the new, union group back in place of the ref group, and make that the canonical group.
  groups[group_i_ref] = group_union
  groups[group_i_read] = None
  for name in group_union:
    member_to_group[name] = group_i_ref


def count_barcodes(families):
  """Tally how many times each barcode appears in the raw data.
  Returns two lists: "votes" and "barcodes". In both, each element corresponds to a different
  barcode, in the order it appears in families. But the first element is None, to correspond with
  the 1-based read names in the alignment. So, you should be able to look up info on each barcode
  by taking its read name in the alignment and using it as an index into these lists.
  Elements in "votes" are integers representing how many times each barcode appears.
  Elements in "barcodes" are the barcode strings themselves."""
  votes = []
  barcodes = []
  vote = -999999
  last_barcode = None
  for line in families:
    fields = line.rstrip('\r\n').split('\t')
    barcode = fields[0]
    if barcode != last_barcode:
      votes.append(vote)
      barcodes.append(last_barcode)
      vote = 0
      last_barcode = barcode
    vote += 1
  votes.append(vote)
  barcodes.append(last_barcode)
  return votes, barcodes


def make_correction_table(groups, votes, barcodes):
  """Take all the gathered data and create a mapping of raw barcode sequences to corrected barcodes.
  """
  corrections = {}
  for i, group in enumerate(groups):
    # if i == 14835:
    #   sys.stderr.write('Found group 14835: {}.\n'.format(group))
    #   logging.getLogger().setLevel(logging.DEBUG)
    if group is None:
      logging.debug('{}:\tNone'.format(i))
      continue
    logging.debug('{}:\t{}'.format(i, ', '.join(map(str, group))))
    best_read = choose_read(group, votes)
    corrected_barcode = barcodes[best_read]
    logging.debug('\tChose read number {} ({})'.format(best_read, corrected_barcode))
    for read_name in group:
      raw_barcode = barcodes[read_name]
      logging.debug('\tMapping {} -> {} ({} -> {})'
                    .format(read_name, best_read, raw_barcode, corrected_barcode))
      corrections[raw_barcode] = corrected_barcode
    # logging.getLogger().setLevel(logging.ERROR)
  return corrections


def choose_read(group, votes):
  max_votes = -1
  best_read = None
  for read_name in group:
    if votes[read_name] > max_votes:
      max_votes = votes[read_name]
      best_read = read_name
  return best_read


def tone_down_logger():
  """Change the logging level names from all-caps to capitalized lowercase.
  E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
  for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
    level_name = logging.getLevelName(level)
    logging.addLevelName(level, level_name.capitalize())


if __name__ == '__main__':
  sys.exit(main(sys.argv))
