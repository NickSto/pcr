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


LogEvent = collections.namedtuple('LogEvent', ('type', 'data1', 'data2'))


class Sentinel(object):
  pass


class WorkerDiedError(Exception):
  pass


class RotatingPool(list):

  def __init__(self, num_workers, function, static_args=None, pause=0.1):
    """Open a pool of worker processes.
    num_workers is the number of worker processes to create.
    function is the function to execute on each set of input data.
    static_args is a list of arguments to always give to the function. These will go at the end of
    the arguments (appended to the input data arguments)."""
    list.__init__(self)
    self.function = function
    self.pause = pause
    self.jobs_submitted = 0
    self.started = False
    if static_args is None:
      self.static_args = []
    else:
      self.static_args = static_args
    for i in range(num_workers):
      self.append(self.make_worker())

  def make_worker(self):
    """Convenience method to spawn and start a new worker with the Pool's parameters."""
    worker = Worker(self.function, static_args=self.static_args, pause=self.pause)
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
  def __init__(self, function, static_args=None, pause=0.1, *args, **kwargs):
    multiprocessing.Process.__init__(self, target=self.worker_loop, *args, **kwargs)
    self.parent_data_pipe, self.child_data_pipe = multiprocessing.Pipe()
    self.parent_log_pipe, self.child_log_pipe = multiprocessing.Pipe()
    self.function = function
    self.pause = pause
    if static_args is None:
      self.static_args = []
    else:
      self.static_args = static_args
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
        error = event.data1
        logging.warning('{} threw a {}: {}'.format(self.name, type(error).__name__, error))
        self.state = 'error'
      elif event.type == 'log':
        level = event.data1
        message = event.data2
        logging.log(LOG_LEVELS[level], message)

  ##### Child methods #####
  def run(self):
    try:
      multiprocessing.Process.run(self)
    except Exception as error:
      trace = traceback.format_exc()
      self.child_log_pipe.send(LogEvent('exception', error, trace))
      raise

  def worker_loop(self):
    input_data = self.child_data_pipe.recv()
    while input_data is not Sentinel:
      args = list(input_data) + self.static_args
      result = self.function(*args)
      self.child_data_pipe.send(result)
      input_data = self.child_data_pipe.recv()
