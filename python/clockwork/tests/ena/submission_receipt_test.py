import unittest
import filecmp
import os
from clockwork.ena import submission_receipt

modules_dir = os.path.join(os.path.dirname(os.path.abspath(submission_receipt.__file__)), os.pardir)
data_dir = os.path.normpath(os.path.join(modules_dir, 'tests', 'data', 'ena', 'submission_receipt'))


class TestSubmissionReceipt(unittest.TestCase):
    def test_init_success_receipt(self):
        '''test __init__ with successful receipt'''
        receipt_xml = os.path.join(data_dir, 'init.good_receipt.xml')
        receipt = submission_receipt.SubmissionReceipt(receipt_xml)
        self.assertTrue(receipt.successful)
        expected_accessions = {
            'SAMPLE': 'ERS1234567',
            'SUBMISSION': 'ERA1234567',
        }
        self.assertEqual(expected_accessions, receipt.accessions)


    def test_init_fail_receipt(self):
        '''test __init__ with failed receipt'''
        receipt_xml = os.path.join(data_dir, 'init.bad_receipt.xml')
        receipt = submission_receipt.SubmissionReceipt(receipt_xml)
        self.assertFalse(receipt.successful)
        self.assertEqual({}, receipt.accessions)

