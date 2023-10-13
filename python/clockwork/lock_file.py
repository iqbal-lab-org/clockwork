import os
from clockwork import utils


class LockFile:
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)
        if os.path.exists(self.filename):
            raise Exception("Lock file already exists: " + self.filename)

        utils.make_empty_file(self.filename)

    def stop(self):
        try:
            os.unlink(self.filename)
        except:
            raise Exception("Error deleting lock file " + self.filename)
