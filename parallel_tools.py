import logging
import traceback
import multiprocessing.pool

QUEUE_SIZE_MULTIPLIER = 8

#TODO: Replace apply_async() implementation with map_async(). It returns an AsyncResult object whose
#      .get() returns a list of the result values. It returns after all jobs finish.
#      The external API can be the same, and the performance should be about the same, except
#      we'll have to wait until queue_size inputs are gathered before starting the first job,
#      and we can't start gathering results until all are done. But that might be negligble.

class SyncAsyncPool(multiprocessing.pool.Pool):
  """A wrapper around multiprocessing.Pool which allows parallel processing but ordered results.
  This offers a compromise between synchronous and asynchronous processing, trying to get the
  benefits of both.
  It issues jobs asynchronously, but syncs up the results so they're processed in the same order
  the inputs were given. It does this by chunking the jobs, periodically stopping to wait for all
  jobs in the chunk to finish.
  It allows giving a callback which will be executed in the parent process. It will also be given
  results in the same order they were submitted."""

  def __init__(self,
               function,
               processes=None,
               queue_size=None,
               static_args=(),
               static_kwargs=None,
               callback=None,
               callback_args=()
              ):
    """Create a new SyncAsyncPool.
    processes can be None, "auto", an integer 0 or greater, or something that produces an integer
      >= 0 when converted with int(). None or "auto" mean to use as many subprocesses as there are
      cpu cores (reported by multiprocessing.cpu_count()). 0 means don't use any subprocesses. It
      will execute the function directly in the main process (and won't actually create a
      multiprocessing.Pool).
    queue_size can be None or an integer greater than 0. If it's None, the queue_size will be set
      to QUEUE_SIZE_MULTIPLIER * the number of processes."""
    # Validate arguments.
    if processes is not None and processes != 'auto':
      try:
        processes = int(processes)
      except (ValueError, TypeError):
        raise ValueError('processes must be an integer, None, or "auto" (received {!r})'.format(processes))
      if processes < 0:
        raise ValueError('processes must be 0 or greater (received {!r})'.format(processes))
    if processes == 'auto':
      processes = None
    if queue_size is not None and queue_size <= 0:
      raise ValueError('queue_size must be > 0 (received {!r})'.format(queue_size))
    # Are we actually doing multiprocessing, or should we do everything directly in one process?
    if processes == 0:
      self.multiproc = False
    else:
      self.multiproc = True
    if self.multiproc:
      multiprocessing.pool.Pool.__init__(self, processes=processes)
    # Determine the number of processes.
    if processes is None or processes == 'auto':
      try:
        processes = multiprocessing.cpu_count()
      except NotImplementedError:
        processes = self._processes
    self.processes = processes
    # Determine the queue size.
    if queue_size is None:
      if self.processes == 0:
        queue_size = 1
      else:
        queue_size = self.processes * QUEUE_SIZE_MULTIPLIER
    self.queue_size = queue_size
    self.function = function
    self.static_args = list(static_args)
    if static_kwargs is None:
      self.static_kwargs = {}
    else:
      self.static_kwargs = static_kwargs
    self.callback = callback
    self.callback_args = callback_args
    self.results = []

  def compute(self, *args, **kwargs):
    # Combine the static arguments with the args for this invocation.
    all_args = list(args) + self.static_args
    all_kwargs = self.static_kwargs.copy()
    all_kwargs.update(kwargs)
    # Send args to multiprocessing pool worker, or execute directly in this process if we're not
    # multiprocessing.
    if self.multiproc:
      result = self.apply_async(with_context, [self.function]+all_args, all_kwargs)
    else:
      result = FakeResult(self.function(*all_args, **all_kwargs))
    self.results.append(result)
    if len(self.results) >= self.queue_size:
      self.flush()
      self.results = []

  def flush(self):
    if self.callback:
      for result in self.results:
        self.callback(result.get(), *self.callback_args)

  def close(self):
    if self.multiproc:
      multiprocessing.pool.Pool.close(self)

  def join(self):
    if self.multiproc:
      multiprocessing.pool.Pool.join(self)


class FakeResult(object):
  """A dummy version of multiprocessing.pool.AsyncResult to hold a result.
  It's convenient to have an object with the same API when we're not doing multiprocessing."""
  def __init__(self, result_data, timeout=None):
    self.result_data = result_data
  def get(self):
    return self.result_data


def with_context(fxn, *args, **kwargs):
  """Execute fxn, logging child process' stack trace for any Exceptions that are raised.
  When Exceptions are raised in a multiprocessing subprocess, the stack trace it gives ends where
  you call .get() on the .apply_async() return value.
  This logs the real stack trace and re-raises it.
  Usage:
  To execute real_fxn(arg1, arg2) through this, do:
    result = pool.apply_async(with_context, args=(real_fxn, arg1, arg2))
  Note: This would be a decorator, but that doesn't work with multiprocessing.
  Functions below the top-level of a module aren't picklable and so can't be passed through
  a Queue to subprocesses:
  https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error/8805244#8805244
  NOTE: This must execute in the child process.
  """
  try:
    return fxn(*args, **kwargs)
  except Exception as exception:
    tb = traceback.format_exc()
    logging.exception(tb)
    raise exception
