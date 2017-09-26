import time
import logging
import traceback
import collections
import multiprocessing

LOG_LEVELS = {
  logging.DEBUG: logging.DEBUG,
  'debug': logging.DEBUG,
  logging.INFO: logging.INFO,
  'info': logging.INFO,
  logging.WARNING: logging.WARNING,
  'warn': logging.WARN,
  'warning': logging.WARNING,
  logging.ERROR: logging.ERROR,
  'error': logging.ERROR,
  logging.CRITICAL: logging.CRITICAL,
  'critcal': logging.CRITICAL,
}


ExceptionEvent = collections.namedtuple('ExceptionEvent', ('type', 'exception', 'traceback'))
LogEvent = collections.namedtuple('LogEvent', ('type', 'level', 'message', 'args', 'kwargs'))


class PipeLogger(object):
  """A drop-in replacement for a logging.Logger which will send logs back to the root process.
  This implements only the methods for actually logging a message: log, debug, info, warning,
  error, and critical.
  Instead of sending the message through the logging module, this will send the message through
  the log pipe back to the parent process which will then log it using the logging module.
  The message will be prefixed with the worker's name."""

  def __init__(self, child_log_pipe):
    self.child_log_pipe = child_log_pipe

  def log(self, level, message, *args, **kwargs):
    event = LogEvent(type='log', level=level, message=message, args=args, kwargs=kwargs)
    self.child_log_pipe.send(event)

  def debug(self, message, *args, **kwargs):
    self.log(logging.DEBUG, message, *args, **kwargs)

  def info(self, message, *args, **kwargs):
    self.log(logging.INFO, message, *args, **kwargs)

  def warning(self, message, *args, **kwargs):
    self.log(logging.WARNING, message, *args, **kwargs)

  def error(self, message, *args, **kwargs):
    self.log(logging.ERROR, message, *args, **kwargs)

  def critical(self, message, *args, **kwargs):
    self.log(logging.CRITICAL, message, *args, **kwargs)


class Sentinel(object):
  pass


class WorkerDiedError(Exception):
  pass


class RotatingPool(list):

  def __init__(self, num_workers, function, static_args=None, pause=0.1, give_logger=False):
    """Open a pool of worker processes.
    num_workers is the number of worker processes to create.
    function is the function to execute on each set of input data.
    static_args is a list of arguments to always give to the function. These will go at the end of
    the arguments (appended to the input data arguments).
    If give_logger is true, the function will be given an object as the last static argument which
    it can use to send messages which will be logged. The object implements methods like log() and
    warning(), so it can often be used as a drop-in replacement for a logging.Logger."""
    list.__init__(self)
    self.function = function
    self.pause = pause
    self.give_logger = give_logger
    self.jobs_submitted = 0
    if static_args is None:
      self.static_args = []
    else:
      self.static_args = list(static_args)
    for i in range(num_workers):
      self.append(self.make_worker())

  def make_worker(self):
    """Convenience method to spawn and start a new worker with the Pool's parameters."""
    worker = Worker(self.function, static_args=self.static_args, pause=self.pause,
                    give_logger=self.give_logger)
    worker.start()
    while not worker.is_alive():
      time.sleep(self.pause)
    worker.state = 'idle'
    return worker

  @property
  def states(self):
    return [worker.state for worker in self]

  def compute(self, *input_data):
    """Process the input data by handing it as arguments to the pool's function.
    Returns a generator where each element is the result of a previous computation,
    or None if none has finished yet.
    If the current worker is still computing, this will block until it's done."""
    current_worker_index = self.jobs_submitted % len(self)
    worker = self[current_worker_index]
    worker.log_events()
    try:
      result = worker.get_result()
    except WorkerDiedError:
      logging.warning('{} died. Replacing it with a new one..'.format(worker.name))
      worker = self.make_worker()
      self[current_worker_index] = worker
      result = Sentinel
    worker.parent_data_pipe.send(input_data)
    worker.state = 'working'
    self.jobs_submitted += 1
    return result

  def flush(self):
    """Get the results of all pending computations.
    Returns a generator where each element is the result of a computation.
    This will block until all pending computations finish."""
    current_worker_index = self.jobs_submitted % len(self)
    workers_processed = 0
    while workers_processed < len(self):
      worker = self[current_worker_index]
      worker.log_events()
      try:
        result = worker.get_result()
        if result is not Sentinel:
          yield result
      except WorkerDiedError:
        logging.warning('{} died.'.format(worker.name))
      workers_processed += 1
      current_worker_index += 1
      if current_worker_index >= len(self):
        current_worker_index = 0

  def stop(self):
    """End all worker processes. Returns no data."""
    for worker in self:
      worker.parent_data_pipe.send(Sentinel)


class Worker(multiprocessing.Process):

  ##### Parent methods #####
  def __init__(self, function, static_args=None, pause=0.1, give_logger=False, *args, **kwargs):
    multiprocessing.Process.__init__(self, target=self.worker_loop, *args, **kwargs)
    self.parent_data_pipe, self.child_data_pipe = multiprocessing.Pipe()
    self.parent_log_pipe, self.child_log_pipe = multiprocessing.Pipe()
    self.function = function
    self.pause = pause
    if static_args is None:
      self.static_args = []
    else:
      self.static_args = list(static_args)
    if give_logger:
      self.static_args.append(PipeLogger(self.child_log_pipe))
    self.state = 'initialized'

  def get_result(self):
    if not self.is_alive():
      raise WorkerDiedError
    if self.state == 'working':
      result = self.parent_data_pipe.recv()
      self.state = 'idle'
      return result
    else:
      return Sentinel

  def log_events(self):
    while self.parent_log_pipe.poll():
      event = self.parent_log_pipe.recv()
      if event.type == 'exception':
        error = event.exception
        logging.warning('{} threw a {}: {}'.format(self.name, type(error).__name__, error))
        self.state = 'error'
      elif event.type == 'log':
        message = '({}) {}'.format(self.name, event.message)
        logging.log(LOG_LEVELS[event.level], message, *event.args, **event.kwargs)

  ##### Child methods #####
  def run(self):
    try:
      multiprocessing.Process.run(self)
    except Exception as error:
      trace = traceback.format_exc()
      self.child_log_pipe.send(ExceptionEvent(type='exception', exception=error, traceback=trace))
      raise

  def worker_loop(self):
    input_data = self.child_data_pipe.recv()
    while input_data is not Sentinel:
      args = list(input_data) + self.static_args
      result = self.function(*args)
      self.child_data_pipe.send(result)
      input_data = self.child_data_pipe.recv()
