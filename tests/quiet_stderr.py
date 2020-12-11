import sys

class NullWriter:
    def write(self, s):
        pass


class Quieter:
    def __enter__(self):
        self.old_stderr = sys.stderr
        sys.stderr = NullWriter()

    def __exit__(self, type, value, traceback):
        sys.stderr = self.old_stderr
