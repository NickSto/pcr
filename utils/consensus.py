#!/usr/bin/env python3
import argparse
import logging
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

WIN_LEN = 4
GAP_CHAR = ' '
DESCRIPTION = """The algorithms from consensus.c, in Python. For testing purposes."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  volume = parser.add_mutually_exclusive_group()
  volume.add_argument('-Q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.INFO)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
  parser.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  subparsers = parser.add_subparsers(dest='command', title='subcommands',
    help='Available subcommands')
  window = subparsers.add_parser('window', help='Test quality score gap window algorithm.')
  window.add_argument('-s', '--seq', type=lstrip,
    help='The literal sequence, as a string. If it starts with gaps, quote and put a space before '
         'the sequence to keep argparse from thinking it\'s an argument. Whitespace will be '
         'stripped before processing.')
  window.add_argument('-q', '--qual',
    help='The literal quality score string.')
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  if args.command == 'window':
    if not (args.seq and args.qual):
      fail('Error: --seq and --qual are required for "window" command.')
    test_window(args.seq, args.qual)


def test_window(seq, quals):
  seq_len = len(seq)
  assert seq_len == len(quals), (seq_len, len(quals))
  window = [None]*(WIN_LEN*2)
  win_edge = init_gap_qual_window(window, quals, seq_len)
  logging.debug('Init: '+format_window(window, win_edge))
  for i, (base, qual) in enumerate(zip(seq, quals)):
    if qual == GAP_CHAR:
      gap_qual = get_gap_qual(window)
      logging.info('gap {:2d}: "{}"/{:<2d} ({:2d})'
                   .format(i+1, gap_qual, ord(gap_qual)-33, ord(gap_qual)))
    else:
      win_edge = push_qual(window, win_edge, quals, seq_len)
      logging.debug('{:4d}: {}'.format(i+1, format_window(window, win_edge)))


def init_gap_qual_window(window, quals, seq_len):
  """This does the initial fill of the window array, adding the first WIN_LEN bases to the right
  side and filling the left side with -1's."""
  # Fill left side with -1's (no quality information).
  i = 0
  while i < WIN_LEN:
    window[i] = -1
    i += 1
  # Fill right side with first WIN_LEN quality scores. Skip gaps, and if you run out of quality
  # scores (if seq_len < WIN_LEN), fill the rest with -1. Leave win_edge at the last base we added.
  i = WIN_LEN
  win_edge = -1
  quals_added = 0
  while quals_added < WIN_LEN:
    win_edge += 1
    if win_edge >= seq_len:
      win_edge = seq_len
      window[i] = -1
      i += 1
      quals_added += 1
    elif quals[win_edge] != GAP_CHAR:
      window[i] = ord(quals[win_edge])
      i += 1
      quals_added += 1
  return win_edge


def push_qual(window, win_edge, quals, seq_len):
  """Push the next non-gap quality score onto the right side of the window."""
  # Find the next quality score that's not a gap.
  next_qual = GAP_CHAR
  while next_qual == GAP_CHAR:
    win_edge += 1
    if win_edge < seq_len:
      next_qual = ord(quals[win_edge])
    else:
      win_edge = seq_len
      next_qual = -1
  # Shift all the quality scores left add the new one.
  i = WIN_LEN*2 - 1
  while i >= 0:
    last_qual = window[i]
    window[i] = next_qual
    next_qual = last_qual
    i -= 1
  return win_edge


def get_gap_qual(window):
  """Compute the quality of the gap based on a weighted average of the quality scores in the window.
  The scores near the center of the window are weighted higher than the ones further away."""
  score_sum = 0
  weight_sum = 0
  weight = 1
  i = 0
  while i < WIN_LEN*2:
    if window[i] != -1:
      score_sum += window[i] * weight
      weight_sum += weight
    # Increase the weight until we get to the middle of the window (at WIN_LEN), then decrease it.
    if i < WIN_LEN - 1:
      weight += 1
    elif i > WIN_LEN - 1:
      weight -= 1
    i += 1
  if weight_sum > 0:
    return chr(int(score_sum//weight_sum))
  else:
    return chr(0)


def format_window(window, win_edge):
  output = '['
  output += '|'.join([format_score(score) for score in window[:4]])
  output += '||'
  output += '|'.join([format_score(score) for score in window[4:]])
  output += '] {:>2d}'.format(win_edge)
  return output


def format_score(score):
  if score == -1:
    return '-1'
  else:
    return str(score-33)


def lstrip(s):
  return s.lstrip()


def fail(message):
  logging.critical(message)
  if __name__ == '__main__':
    sys.exit(1)
  else:
    raise Exception('Unrecoverable error')


if __name__ == '__main__':
  try:
    sys.exit(main(sys.argv))
  except BrokenPipeError:
    pass
