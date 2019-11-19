import os
from clockwork import utils


class Error(Exception):
    pass


class LockFile:
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)
        if os.path.exists(self.filename):
            raise Error("Lock file already exists: " + self.filename)

        utils.make_empty_file(self.filename)

    def stop(self):
        try:
            os.unlink(self.filename)
        except:
            raise Error("Error deleting lock file " + self.filename)
