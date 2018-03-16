import unittest
import os
import shutil
from clockwork import mykrobe

modules_dir = os.path.dirname(os.path.abspath(mykrobe.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'mykrobe')

class TestMykrobe(unittest.TestCase):
    def test_run_predict(self):
        '''test run_predict'''
        reads_file = os.path.join(data_dir, 'run_predict.reads.fq.gz')
        tmp_out = 'tmp.mykrobe_run_predict'
        mykrobe.run_predict(reads_file, tmp_out, 'sample_name', 'tb')
        shutil.rmtree(tmp_out)

