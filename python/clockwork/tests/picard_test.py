import unittest
import os
from clockwork import picard

modules_dir = os.path.dirname(os.path.abspath(picard.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'picard')

class TestPicard(unittest.TestCase):
    def test_mark_duplicates(self):
        '''test mark_duplicates'''
        bam_in = os.path.join(data_dir, 'mark_duplicates.in.bam')
        tmp_out = 'tmp.picard_mark_duplicates.bam'
        picard.mark_duplicates(bam_in, tmp_out)
        self.assertTrue(os.path.exists(tmp_out))
        os.unlink(tmp_out)
