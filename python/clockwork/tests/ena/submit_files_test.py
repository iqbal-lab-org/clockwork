import unittest
import filecmp
import os
from clockwork.ena import submit_files, submission_receipt

modules_dir = os.path.join(os.path.dirname(os.path.abspath(submit_files.__file__)), os.pardir)
data_dir = os.path.normpath(os.path.join(modules_dir, 'tests', 'data', 'ena', 'submit_files'))


class TestSubmitFiles(unittest.TestCase):
    def test_make_dummy_success_receipt(self):
        '''test _make_dummy_success_receipt'''
        tmp_receipt = 'tmp.ena.submit_files.make_dummy_success_receipt.xml'
        submit_files._make_dummy_success_receipt(tmp_receipt, 'project')
        receipt = submission_receipt.SubmissionReceipt(tmp_receipt)
        self.assertTrue(receipt.successful)
        os.unlink(tmp_receipt)


    def test_make_dummy_fail_receipt(self):
        '''test _make_dummy_fail_receipt'''
        tmp_receipt = 'tmp.ena.submit_files.make_dummy_fail_receipt.xml'
        submit_files._make_dummy_fail_receipt(tmp_receipt)
        receipt = submission_receipt.SubmissionReceipt(tmp_receipt)
        self.assertFalse(receipt.successful)
        os.unlink(tmp_receipt)


    def test_parse_config_file(self):
        '''test test_parse_config_file'''
        good_file = os.path.join(data_dir, 'conf.good.ini')
        user, password = submit_files.parse_config_file(good_file)
        self.assertEqual('username', user)
        self.assertEqual('correcthorsebatterystaple', password)

        for x in ['no_ena_login', 'no_user', 'no_password']:
            filename = os.path.join(data_dir, 'conf.bad.' + x + '.ini')
            with self.assertRaises(Exception):
                submit_files.parse_config_file(filename)

