#!/usr/bin/env python3
import argparse
import gzip
import logging
import os
import subprocess
import sys
assert sys.version_info.major >= 3, 'Python 3 required'

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DESCRIPTION = """Run the entire Du Novo pipeline."""


def make_argparser():
  parser = argparse.ArgumentParser(description=DESCRIPTION)
  io = parser.add_argument_group('Inputs and outputs')
  io.add_argument('fastq1', metavar='reads_1.fq', type=open_as_text_or_gzip,
    help='Input reads (mate 1). Can be gzipped.')
  io.add_argument('fastq2', metavar='reads_2.fq', type=open_as_text_or_gzip,
    help='Input reads (mate 2). Can be gzipped.')
  io.add_argument('-o', '--outdir')
  io.add_argument('-s', '--suffix',
    help='A string to use in naming the output files. If given, will be put just before the file '
         'extension (like "families.suffix.tsv").')
  io.add_argument('-d', '--dunovo-dir', default=SCRIPT_DIR,
    help='The directory containing the version of Du Novo you want to run. Default: The directory '
         'containing this script (%(default)s).')
  io.add_argument('-T', '--tempdir')
  params = parser.add_argument_group('Algorithm parameters')
  params.add_argument('-D', '--dist', type=int,
    help='correct.py --dist. Default: the correct.py default.')
  params.add_argument('-m', '--mapq', type=int,
    help='correct.py --mapq. Default: the correct.py default.')
  params.add_argument('-P', '--pos', type=int,
    help='correct.py --pos. Default: the correct.py default.')
  params.add_argument('-a', '--aligner', choices=('mafft', 'kalign'), default='kalign',
    help='align-families.py --aligner. Default: %(default)s')
  params.add_argument('-r', '--min-reads', type=int,
    help='make-consensi.py --min-reads. Default: the make-consensi.py default.')
  params.add_argument('-q', '--qual', type=int,
    help='make-consensi.py --qual. Default: the make-consensi.py default.')
  params.add_argument('-c', '--cons-thres', type=float, default=0.7,
    help='make-consensi.py --cons-thres. Default: %(default)s.')
  params.add_argument('-C', '--min-cons-reads', type=int,
    help='make-consensi.py --min-cons-reads. Default: the make-consensi.py default.')
  params.add_argument('-Q', '--fake-phred', type=int, default=40,
    help='The PHRED score to assign to all output consensus bases.')
  params.add_argument('-I', '--no-check-ids', action='store_true',
    help='Pass --no-check-ids to correct.py and align-families.py.')
  params.add_argument('-p', '--processes', type=int,
    help='align-families.py --processes. Default: the align-families.py default.')
  params.add_argument('-t', '--threads', type=int,
    help='baralign.sh -t. Default: the baralign.sh default.')
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log', type=argparse.FileType('w'), default=sys.stderr,
    help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
  volume = log.add_mutually_exclusive_group()
  volume.add_argument('--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')

  # Create and check output paths.
  if args.outdir:
    outdir = args.outdir
  else:
    outdir = os.path.dirname(args.fastq1.name) or '.'
  if not os.path.isdir(outdir):
    fail('Error: Output directory {!r} does not exist.'.format(outdir))
  if args.suffix:
    suffix = '.'+args.suffix
  else:
    suffix = ''
  paths = {
    'families': os.path.join(outdir, 'families{}.tsv'.format(suffix)),
    'refdir': os.path.join(outdir, 'refdir{}'.format(suffix)),
    'correct_sam': os.path.join(outdir, 'correct{}.sam'.format(suffix)),
    'families_corrected': os.path.join(outdir, 'families.corrected{}.tsv'.format(suffix)),
    'msa': os.path.join(outdir, 'families.msa{}.tsv'.format(suffix)),
    'sscs1': os.path.join(outdir, 'sscs{}_1.fq'.format(suffix)),
    'sscs2': os.path.join(outdir, 'sscs{}_2.fq'.format(suffix)),
    'duplex1': os.path.join(outdir, 'duplex{}_1.fq'.format(suffix)),
    'duplex2': os.path.join(outdir, 'duplex{}_2.fq'.format(suffix)),
  }
  existing_paths = False
  for path in paths.values():
    if os.path.exists(path):
      existing_paths = True
      logging.critical('Error: {!r} already exists.'.format(path))
  if existing_paths:
    return 1

  input_size = estimate_filesize(args.fastq1)
  mem_req = get_mem_requirement(input_size)

  # The 1st pipeline.
  # $ make-barcodes.awk
  command = ['awk', '-f', os.path.join(args.dunovo_dir, 'make-barcodes.awk')]
  logging.warning('$ {} \\'.format(' '.join(command)))
  make_barcodes = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
  with open(paths['families'], 'w') as families_file:
    # $ sort
    command = ['sort'] + get_sort_args(mem_req, args.tempdir)
    logging.warning('  | {} > {}'.format(' '.join(command), paths['families']))
    sort = subprocess.Popen(command, stdin=make_barcodes.stdout, stdout=families_file)
    # Execute it by starting to feed data into the front of the pipeline.
    make_barcodes.stdout.close()
    for columns in paste_magic(args.fastq1, args.fastq2):
      line = '\t'.join(columns)+'\n'
      make_barcodes.stdin.write(bytes(line, 'utf8'))
    make_barcodes.stdin.close()
  for process in (make_barcodes, sort):
    result = process.wait()
    if result != 0:
      fail('Error: Process exited with code {}: $ {}'.format(result, ' '.join(process.args)))

  # The 2nd pipeline.
  # $ baralign.sh
  command = [os.path.join(args.dunovo_dir, 'baralign.sh'),
             paths['families'], paths['refdir'], paths['correct_sam']]
  if args.threads:
    command.insert(1, '-t')
    command.insert(2, str(args.threads))
  logging.warning('$ '+' '.join(command))
  baralign = subprocess.Popen(command)
  result = baralign.wait()
  if result != 0:
    fail('Error: Process exited with code {}: $ {}'.format(result, ' '.join(baralign.args)))

  # The 3rd pipeline.
  # $ correct.py
  command = ([os.path.join(args.dunovo_dir, 'correct.py')] + get_correct_args(**vars(args))
             + [paths['families'], paths['refdir']+'/barcodes.fa', paths['correct_sam']])
  logging.warning('$ {} \\'.format(' '.join(command)))
  correct = subprocess.Popen(command, stdout=subprocess.PIPE)
  # $ sort
  command = ['sort'] + get_sort_args(mem_req, args.tempdir)
  logging.warning('  | {} \\'.format(' '.join(command)))
  sort = subprocess.Popen(command, stdin=correct.stdout, stdout=subprocess.PIPE)
  # $ tee
  command = ['tee', '-a', paths['families_corrected']]
  logging.warning('  | {} \\'.format(' '.join(command)))
  tee1 = subprocess.Popen(command, stdin=sort.stdout, stdout=subprocess.PIPE)
  # $ align-families.py
  command = ([os.path.join(args.dunovo_dir, 'align-families.py')]
             + get_align_families_args(**vars(args)))
  logging.warning('  | {} \\'.format(' '.join(command)))
  align_families = subprocess.Popen(command, stdin=tee1.stdout, stdout=subprocess.PIPE)
  # $ tee
  command = ['tee', '-a', paths['msa']]
  logging.warning('  | {} \\'.format(' '.join(command)))
  tee2 = subprocess.Popen(command, stdin=align_families.stdout, stdout=subprocess.PIPE)
  # $ make-consensi.py
  command = ([os.path.join(args.dunovo_dir, 'make-consensi.py')]
             + get_make_consensi_args(**vars(args)) + ['--sscs1', paths['sscs1'], '--sscs2',
             paths['sscs2'], '-1', paths['duplex1'], '-2', paths['duplex2']])
  logging.warning('  | '+' '.join(command))
  make_consensi = subprocess.Popen(command, stdin=tee2.stdout)
  # Kick it off by closing the first stdout (don't ask me).
  correct.stdout.close()
  # Check return codes.
  for process in (correct, sort, tee1, align_families, tee2, make_consensi):
    result = process.wait()
    if result != 0:
      fail('Error: Process exited with code {}: $ {}'.format(result, ' '.join(process.args)))


def open_as_text_or_gzip(path):
  """Return an open file-like object reading the path as a text file or a gzip file, depending on
  which it looks like."""
  if detect_gzip(path):
    return gzip.open(path, 'rt')
  else:
    return open(path, 'rU')


def detect_gzip(path):
  """Return True if the file looks like a gzip file: ends with .gz or contains non-ASCII bytes."""
  ext = os.path.splitext(path)[1]
  if ext == '.gz':
    return True
  elif ext in ('.txt', '.tsv', '.csv'):
    return False
  with open(path) as fh:
    is_not_ascii = detect_non_ascii(fh.read(100))
  if is_not_ascii:
    return True


def detect_non_ascii(bytes, max_test=100):
  """Return True if any of the first "max_test" bytes are non-ASCII (the high bit set to 1).
  Return False otherwise."""
  for i, char in enumerate(bytes):
    # Is the high bit a 1?
    if ord(char) & 128:
      return True
    if i >= max_test:
      return False
  return False


def paste_magic(reads1, reads2):
  # Emulate $ paste reads_1.fq reads_2.fq | paste - - - -
  columns = [None] * 8
  for i, (line1, line2) in enumerate(zip(reads1, reads2)):
    line_type = i % 4
    columns[line_type*2] = line1.rstrip('\r\n')
    columns[(line_type*2)+1] = line2.rstrip('\r\n')
    if line_type == 3:
      yield columns


def estimate_filesize(file_obj):
  size = os.path.getsize(file_obj.name)
  if isinstance(file_obj, gzip.GzipFile):
    size = size * 3
  return size


def get_mem_requirement(filesize):
  gigs = filesize/1024/1024/1024
  return round(gigs*3.5)


def get_sort_args(mem_req, tempdir):
  args = []
  if mem_req >= 2:
    logging.info('Telling sort to use a {}GB memory space.'.format(mem_req))
    args.extend(['-S', '{}G'.format(mem_req)])
  if tempdir:
    args.extend(['-T', tempdir])
  return args


def get_correct_args(**kwargs):
  arg_list = ('dist', 'mapq', 'pos')
  flag_list = ('no_check_ids',)
  return get_generic_args(arg_list, flag_list, kwargs)


def get_align_families_args(**kwargs):
  arg_list = ('aligner', 'processes')
  flag_list = ('no_check_ids',)
  return get_generic_args(arg_list, flag_list, kwargs)


def get_make_consensi_args(fake_phred=40, **kwargs):
  arg_list = ('min_reads', 'qual', 'cons_thres', 'min_cons_thres')
  args = ['--fastq-out', str(fake_phred)]
  return args + get_generic_args(arg_list, (), kwargs)


def get_generic_args(arg_list, flag_list, kwargs):
  args = []
  for arg_name in arg_list:
    value = kwargs.get(arg_name)
    if value is not None:
      arg_str = '--'+arg_name.replace('_', '-')
      args.extend([arg_str, str(value)])
  for flag_name in flag_list:
    value = kwargs.get(flag_name)
    if value is True:
      flag_str = '--'+flag_name.replace('_', '-')
      args.append(flag_str)
  return args


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
