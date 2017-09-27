import time
import logging
import traceback
import collections
import multiprocessing

LOG_LEVELS = {
  'debug': logging.DEBUG,
  'info': logging.INFO,
  'warn': logging.WARN,
  'warning': logging.WARNING,
  'error': logging.ERROR,
  'critical': logging.CRITICAL,
}


class RotatingPool(list):
  """A pool of workers which will execute one function on a series of input data and directly return
  the result.
  This is a conceptually simplified variation on multiprocessing.Pool. Instead of interacting with
  the workers asynchronously with callbacks, this pool will directly return the result of the
  function. It won't give it to you right after you hand it the input data, but a few calls later.
  This is useful for times when you want to process a stream of data in parallel, but want the
  results back in the same loop, and don't care to match inputs to outputs. That can be done, since
  this keeps outputs in the same order as the inputs, but it takes a little work to sync them back
  up.
  As the name suggests, this is implemented with a rotating pool of workers. It creates N workers,
  hands the first input to the first one, then the next, and so on. At the end, it loops back
  to the start and starts getting a result from each worker before handing it new input.
  So on the N+1th call to .compute(), it'll return the result from the 1st input.
  If you haven't called .compute() N+1 times yet, it'll return parallel.Sentinel to let you know
  it's not actually a real result.
  This is actually a subclass of list, where each member is a Worker."""

  def __init__(self, num_workers, function, static_args=None, give_logger=False):
    """Open a pool of worker processes.
    num_workers is the number of worker processes to create.
    function is the function to execute on each set of input data.
    static_args is a list of arguments to always give to the function. These will go at the end of
    the arguments (appended to the input data arguments).
    If give_logger is true, the function will be given a parallel.PipeLogger object as the last
    static argument which it can use to send messages which will be logged. The object implements
    methods like log() and warning(), so it can often be used as a drop-in replacement for a
    logging.Logger."""
    list.__init__(self)
    self.function = function
    self.give_logger = give_logger
    self.jobs_submitted = 0
    if static_args is None:
      self.static_args = []
    else:
      self.static_args = list(static_args)
    for i in range(num_workers):
      self.append(self.make_worker())

  def make_worker(self):
    """Convenience method to spawn and start a new worker with the pool's parameters."""
    worker = Worker(self.function, static_args=self.static_args, give_logger=self.give_logger)
    worker.start()
    while not worker.is_alive():
      time.sleep(0.01)
    worker.state = 'idle'
    return worker

  @property
  def states(self):
    return [worker.state for worker in self]

  def compute(self, *input_data):
    """Process the input data by handing it as arguments to the pool's function.
    Returns the result of a previous computation, or parallel.Sentinel (the class, not an instance
    of it) if none has finished yet.
    The computation it will return on the nth call will be from the n-num_workers input.
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
  def __init__(self, function, static_args=None, give_logger=False, *args, **kwargs):
    multiprocessing.Process.__init__(self, target=self.worker_loop, *args, **kwargs)
    self.parent_data_pipe, self.child_data_pipe = multiprocessing.Pipe()
    self.parent_log_pipe, self.child_log_pipe = multiprocessing.Pipe()
    self.function = function
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
        level = LOG_LEVELS.get(event.level, event.level)
        logging.log(level, message, *event.args, **event.kwargs)

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


class Sentinel(object):
  pass


class WorkerDiedError(Exception):
  pass


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
