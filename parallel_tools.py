import logging
import traceback


class FakePool(object):
  """A dummy version of multiprocessing.Pool which isn't actually parallelized.
  Use this instead of multiprocessing.Pool to do all work inside one process."""
  def apply_async(self, fxn, args=[], kwargs={}):
    result_data = fxn(*args, **kwargs)
    return FakeResult(result_data)
  def close(self):
    pass
  def join(self):
    pass


class FakeResult(object):
  """A dummy version of multiprocessing.pool.AsyncResult to hold a result from FakePool."""
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
  """
  try:
    return fxn(*args, **kwargs)
  except Exception as exception:
    tb = traceback.format_exc()
    logging.exception(tb)
    raise exception
