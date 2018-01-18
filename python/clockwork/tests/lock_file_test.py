import unittest
import os
from clockwork import lock_file, utils

class TestDbMaker(unittest.TestCase):
    def test_init_and_stop(self):
        '''test init and stop'''
        tmp_file = 'test.lock_file'
        utils.make_empty_file(tmp_file)
        with self.assertRaises(lock_file.Error):
            lock_file.LockFile(tmp_file)
        os.unlink(tmp_file)

        lock = lock_file.LockFile(tmp_file)
        self.assertTrue(os.path.exists(tmp_file))

        os.unlink(tmp_file)
        with self.assertRaises(lock_file.Error):
            lock.stop()

        lock = lock_file.LockFile(tmp_file)
        self.assertTrue(os.path.exists(tmp_file))
        lock.stop()
        self.assertFalse(os.path.exists(tmp_file))
