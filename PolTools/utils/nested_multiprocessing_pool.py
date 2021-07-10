import multiprocessing.pool

class NestedPool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwds):
        proc = super(NestedPool, self).Process(*args, **kwds)

        class NonDaemonProcess(proc.__class__):
            """Monkey-patch process to ensure it is never daemonized"""
            @property
            def daemon(self):
                return False

            @daemon.setter
            def daemon(self, val):
                pass

        proc.__class__ = NonDaemonProcess
        return proc