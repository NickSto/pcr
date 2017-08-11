"""Stub versions of optional submodules which may fail to clone."""

class version(object):
  def __init__(self):
    self.is_shim = True
  def get_version(self):
    return None


class simplewrap(object):
  def __init__(self):
    self.is_shim = True
  class Wrapper(object):
    def __init__(self):
      self.width = 80
    def wrap(self, string):
      return string


class phone(object):
  def __init__(self):
    self.is_shim = True
  def send_start(script_path,
                 version,
                 domain=None,
                 secure=None,
                 platform=None,
                 test=False):
    pass
  def send_end(script_path,
               version,
               run_id,
               run_time,
               optional_run_data={},
               domain=None,
               secure=None,
               platform=None,
               test=False):
    pass
