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
  params.add_argument('-D', '--dist', type=int, default=3,
    help='correct.py --dist. Default: %(default)s')
  params.add_argument('-m', '--mapq', type=int, default=25,
    help='correct.py --mapq. Default: %(default)s')
  params.add_argument('-P', '--pos', type=int,
    help='correct.py --pos. Default: the correct.py default.')
  params.add_argument('-a', '--aligner', choices=('mafft', 'kalign'), default='kalign',
    help='align-families.py --aligner. Default: %(default)s')
  params.add_argument('-r', '--min-reads', type=int,
    help='make-consensi.py --min-reads. Default: the make-consensi.py default.')
  params.add_argument('-q', '--qual', type=int, default=25,
    help='make-consensi.py --qual. Default: %(default)s')
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
  params.add_argument('-b', '--filt-bases', default='N',
    help='trimmer.py --filt-bases. Default: %(default)s')
  params.add_argument('-e', '--thres', type=float, default=0.3,
    help='trimmer.py --thres. Default: %(default)s')
  params.add_argument('-w', '--window', type=int, default=10,
    help='trimmer.py --window. Default: %(default)s')
  params.add_argument('-M', '--min-length', type=int, default=75,
    help='trimmer.py --min-length. Default: %(default)s')
  log = parser.add_argument_group('Logging')
  log.add_argument('-l', '--log-dir',
    help='Write log output to files in this directory instead of to stderr.')
  volume = log.add_mutually_exclusive_group()
  volume.add_argument('--quiet', dest='volume', action='store_const', const=logging.CRITICAL,
    default=logging.WARNING)
  volume.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
  volume.add_argument('--debug', dest='volume', action='store_const', const=logging.DEBUG)
  return parser


def main(argv):

  parser = make_argparser()
  args = parser.parse_args(argv[1:])

  # Configure logging.
  make_log_dir(args.log_dir)
  stream = make_main_log_stream(args.log_dir, args.suffix)
  logging.basicConfig(stream=stream, level=args.volume, format='%(message)s')

  # Create and check output paths.
  paths, log_paths = make_paths(args.outdir, args.log_dir, args.fastq1.name, args.suffix)
  if invalid_paths(paths, log_paths, args.log_dir):
    return 1
  logs = open_logs(log_paths)

  input_size = estimate_filesize(args.fastq1)
  mem_req = get_mem_requirement(input_size)

  # The 1st pipeline.
  # $ paste
  stdin = {'function':paste_magic, 'fxn_args':(args.fastq1, args.fastq2)}
  steps = [
    {  # $ make-barcodes.awk
      'command': ('awk', '-f', os.path.join(args.dunovo_dir, 'make-barcodes.awk')),
      'stderr': logs['make-barcodes']
    },
    {  # $ sort
      'command': ['sort'] + get_sort_args(mem_req, args.tempdir),
      'stderr': logs['sort1']
    }
  ]
  run_pipeline(steps, stdin=stdin, stdout=paths['families'])

  # The 2nd pipeline.
  steps = [
    {  # $ baralign.sh
      'command': (os.path.join(args.dunovo_dir, 'baralign.sh'), paths['families'], paths['refdir'],
                  paths['correct_sam']),
      'stderr': logs['baralign']
    }
  ]
  run_pipeline(steps)


  # The 3rd pipeline.
  steps = [
    {  # $ correct.py
      'command': ([os.path.join(args.dunovo_dir, 'correct.py')] + get_correct_args(**vars(args))
                  + [paths['families'], paths['refdir']+'/barcodes.fa', paths['correct_sam']]),
      'stderr': logs['correct']
    },
    {  # $ sort
      'command': ['sort'] + get_sort_args(mem_req, args.tempdir),
      'stderr': logs['sort2']
    },
    {  # $ tee
      'command': ('tee', '-a', paths['families_corrected']),
      'stderr': None
    },
    {  # $ align-families.py
      'command': ([os.path.join(args.dunovo_dir, 'align-families.py')]
                  + get_align_families_args(**vars(args))),
      'stderr': logs['align-families']
    },
    {  # $ tee
      'command': ('tee', '-a', paths['msa']),
      'stderr': None
    },
    {  # $ make-consensi.py
      'command': ([os.path.join(args.dunovo_dir, 'make-consensi.py')]
                  + get_make_consensi_args(**vars(args)) + ['--sscs1', paths['sscs1'], '--sscs2',
                  paths['sscs2'], '-1', paths['duplex1'], '-2', paths['duplex2']]),
      'stderr': logs['make-consensi']
    }
  ]
  run_pipeline(steps)

  # The 4th pipeline.
  steps = [
    {  # $ trimmer.py
      'command': ([os.path.join(args.dunovo_dir, 'bfx/trimmer.py'), '--format', 'fastq']
                  + get_trimmer_args(**vars(args)) + [paths['duplex1'], paths['duplex2'],
                   paths['dupfilt1'], paths['dupfilt2']]),
      'stderr': logs['trimmer']
    }
  ]
  run_pipeline(steps)


def open_as_text_or_gzip(path):
  """Return an open file-like object reading the path as a text file or a gzip file, depending on
  which it looks like."""
  if detect_gzip(path):
    file_obj = gzip.open(path, 'rt')
    file_obj.type = 'gzip'
  else:
    file_obj = open(path, 'rU')
    file_obj.type = 'raw'
  return file_obj


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


def make_log_dir(log_dir):
  if log_dir is None or os.path.isdir(log_dir):
    return False
  elif os.path.exists(log_dir):
    fail('Error: {!r} exists but is not a directory.'.format(log_dir))
  else:
    parent_dir = os.path.dirname(log_dir)
    if not os.path.isdir(parent_dir):
      fail('Error: {!r} does not exist or is not a directory.'.format(parent_dir))
    os.mkdir(log_dir, mode=0o775)
    return True


def make_main_log_stream(log_dir, suffix_arg):
  if log_dir is None:
    return sys.stderr
  else:
    if suffix_arg:
      suffix = '.'+suffix_arg
    else:
      suffix = ''
    main_log_path = os.path.join(log_dir, 'dunovo{}.log'.format(suffix))
    return open(main_log_path, 'w')


def make_paths(outdir_arg, log_dir, fastq1_path, suffix_arg):
  if outdir_arg:
    outdir = outdir_arg
  else:
    outdir = os.path.dirname(fastq1_path) or '.'
  if suffix_arg:
    suffix = '.'+suffix_arg
  else:
    suffix = ''
  output_templates = {
    'families': 'families{}.tsv',
    'refdir': 'refdir{}',
    'correct_sam': 'correct{}.sam',
    'families_corrected': 'families.corrected{}.tsv',
    'msa': 'families.msa{}.tsv',
    'sscs1': 'sscs{}_1.fq',
    'sscs2': 'sscs{}_2.fq',
    'duplex1': 'duplex{}_1.fq',
    'duplex2': 'duplex{}_2.fq',
    'dupfilt1': 'duplex.filt{}_1.fq',
    'dupfilt2': 'duplex.filt{}_2.fq',
  }
  log_templates = {
    'make-barcodes': 'make-barcodes{}.log',
    'sort1': 'sort1{}.log',
    'baralign': 'baralign{}.log',
    'correct': 'correct{}.log',
    'sort2': 'sort2{}.log',
    'align-families': 'align-families{}.log',
    'make-consensi': 'make-consensi{}.log',
    'trimmer': 'trimmer{}.log',
  }
  output_paths = {}
  for name, template in output_templates.items():
    output_paths[name] = os.path.join(outdir, template.format(suffix))
  log_paths = {}
  for name, template in log_templates.items():
    if log_dir:
      log_paths[name] = os.path.join(log_dir, template.format(suffix))
    else:
      log_paths[name] = None
  return output_paths, log_paths


def invalid_paths(paths, log_paths, log_dir):
  invalid_paths = False
  for path in list(paths.values()) + list(log_paths.values()):
    if path is None:
      continue
    if not os.path.isdir(os.path.dirname(path)):
      logging.critical('Error: Directory {!r} doesn\'t exist.'.format(os.path.dirname(path)))
      invalid_paths = True
    if os.path.exists(path):
      logging.critical('Error: {!r} already exists.'.format(path))
      invalid_paths = True
  return invalid_paths


def open_logs(log_paths):
  # Open the log files.
  logs = {}
  for name, path in log_paths.items():
    if path is None:
      logs[name] = sys.stderr
    else:
      logs[name] = open(path, 'w')
  return logs


def paste_magic(reads1, reads2):
  # Emulate $ paste reads_1.fq reads_2.fq | paste - - - -
  columns = [None] * 8
  for i, (line1, line2) in enumerate(zip(reads1, reads2)):
    line_type = i % 4
    columns[line_type*2] = line1.rstrip('\r\n')
    columns[(line_type*2)+1] = line2.rstrip('\r\n')
    if line_type == 3:
      line = '\t'.join(columns)+'\n'
      yield bytes(line, 'utf8')


def estimate_filesize(file_obj):
  size = os.path.getsize(file_obj.name)
  if file_obj.type == 'gzip':
    logging.info('Input is gzip. Tripling estimated filesize from {}.'.format(size, size*3))
    size = size * 3
  logging.info('Estimated filesize: {} bytes'.format(size))
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


def get_trimmer_args(**kwargs):
  arg_list = ('filt_bases', 'thres', 'window', 'min_length')
  return get_generic_args(arg_list, (), kwargs)


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


def run_pipeline(steps, stdin=None, stdout=None):
  processes = []
  fxn_stdin = False
  for i, step in enumerate(steps):
    cmd_str = ' '.join(step['command'])
    kwargs = {}
    # Determine the inputs and outputs.
    #   stderr
    if step['stderr']:
      kwargs['stderr'] = step['stderr']
      if step['stderr'] is not sys.stderr:
        cmd_str += ' 2> '+step['stderr'].name
    #  stdin
    if i == 0:
      if isinstance(stdin, str):
        kwargs['stdin'] = open(stdin)
        cmd_str = '$ < '+stdin+' '+cmd_str
      elif isinstance(stdin, dict):
        assert 'function' in stdin, stdin
        kwargs['stdin'] = subprocess.PIPE
        fxn_stdin = True
        cmd_str = '$ <{}> | {}'.format(stdin['function'].__name__, cmd_str)
      else:
        kwargs['stdin'] = stdin
        if hasattr(stdin, 'name'):
          cmd_str = '< '+stdin.name+' '+cmd_str
        cmd_str = '$ '+cmd_str
    else:
      kwargs['stdin'] = processes[i-1].stdout
      cmd_str = '  | '+cmd_str
    #   stdout
    if i == len(steps)-1:
      if isinstance(stdout, str):
        kwargs['stdout'] = open(stdout, 'w')
        cmd_str += ' > '+stdout
      else:
        kwargs['stdout'] = stdout
        if hasattr(stdout, 'name'):
          cmd_str += ' > '+stdout.name
    else:
      kwargs['stdout'] = subprocess.PIPE
      cmd_str += ' \\'
    # Create the actual process and add it to the list.
    processes.append(subprocess.Popen(step['command'], **kwargs))
    logging.warning(cmd_str)
  # Kick it off by closing the first stdout (don't ask me).
  if processes[0].stdout:
    processes[0].stdout.close()
  # If the input is from a function, start feeding it in.
  if fxn_stdin:
    function = stdin['function']
    args = stdin['fxn_args']
    for line in function(*args):
      processes[0].stdin.write(line)
    processes[0].stdin.close()
  # Wait for all processes to finish and check their exit codes.
  for process in processes:
    result = process.wait()
    if result != 0:
      fail('Error: Process exited with code {}: $ {}'.format(result, ' '.join(process.args)))


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
