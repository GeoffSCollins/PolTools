import sys

class NullWriter:
    def write(self, s):
        pass

    def flush(self):
        pass


class Quieter:
    def __enter__(self):
        self.old_stderr = sys.stderr
        self.old_stdout = sys.stdout
        sys.stderr = NullWriter()
        sys.stdout = NullWriter()

    def __exit__(self, type, value, traceback):
        sys.stderr = self.old_stderr
        sys.stdout = self.old_stdout
