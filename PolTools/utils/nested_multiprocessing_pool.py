import multiprocessing.pool

# Code taken from https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic

class _NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class NestedPool(multiprocessing.pool.Pool):
    Process = _NoDaemonProcess