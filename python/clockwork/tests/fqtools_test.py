import os
import datetime
import filecmp
import unittest
from clockwork import fqtools

modules_dir = os.path.dirname(os.path.abspath(fqtools.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "fqtools")


class TestFqtools(unittest.TestCase):
    def test_validate(self):
        """test validate"""
        bad_file = os.path.join(data_dir, "validate.bad.fq")
        with self.assertRaises(Exception):
            fqtools.validate([bad_file])

        bad_pair_1 = os.path.join(data_dir, "validate.bad.pair.1.fq")
        bad_pair_2 = os.path.join(data_dir, "validate.bad.pair.2.fq")
        with self.assertRaises(Exception):
            fqtools.validate([bad_pair_1, bad_pair_2])

        ok_pair_1 = os.path.join(data_dir, "validate.ok.pair.1.fq")
        ok_pair_2 = os.path.join(data_dir, "validate.ok.pair.2.fq")
        fqtools.validate([ok_pair_1, ok_pair_2])

    def test_count(self):
        """test count"""
        file1 = os.path.join(data_dir, "count.1.fq")
        file2 = os.path.join(data_dir, "count.2.fq")
        got = fqtools.count([file1, file2])
        self.assertEqual(3, got)
