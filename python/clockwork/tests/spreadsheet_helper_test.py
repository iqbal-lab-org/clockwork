import unittest
import copy
import datetime
import shutil
import filecmp
import os
import re
from clockwork import spreadsheet_helper

modules_dir = os.path.dirname(os.path.abspath(spreadsheet_helper.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'spreadsheet_helper')
db_ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')


class TestSpreadsheetHelper(unittest.TestCase):
    def test_looks_like_date(self):
        '''Test _looks_like_date'''
        self.assertFalse(spreadsheet_helper.looks_like_date('not even close!'))
        self.assertFalse(spreadsheet_helper.looks_like_date('20001-23-11'))
        self.assertFalse(spreadsheet_helper.looks_like_date('2000-13-01'))
        self.assertFalse(spreadsheet_helper.looks_like_date('2000-04-31'))
        self.assertFalse(spreadsheet_helper.looks_like_date('20000430'))
        self.assertTrue(spreadsheet_helper.looks_like_date('2000-04-30'))


    def test_load_data_from_spreadsheet_tsv(self):
        '''test load_data_from_spreadsheet tsv file'''
        expected = [
           {'subject_id': 'p1',
            'site_id': 's1',
            'lab_id': 'l1',
            'isolate_number': '42',
            'sequence_replicate_number': '43',
            'submission_date': datetime.date(2017, 12, 25),
            'reads_file_1': 'reads_1_1.fq',
            'reads_file_1_md5': 'abcdefghijklmnopqrstuvwyx123456',
            'reads_file_2': 'reads_1_2.fq',
            'reads_file_2_md5': 'abcdefghijklmnopqrstuvwyx123457',
            'dataset_name': 'g1',
            'instrument_model': 'Illumina HiSeq 2000',
            'ena_center_name': 'Center 1',
            'submit_to_ena': '0',
            'ena_on_hold': '0',
            'ena_run_accession': 'ERR123456',
            'ena_sample_accession': 'ERS123456',
           },
           {'subject_id': 'p2',
            'site_id': 's2',
            'lab_id': 'l2',
            'isolate_number': '44',
            'sequence_replicate_number': '45',
            'submission_date': datetime.date(2017, 12, 26),
            'reads_file_1': 'reads_2_1.fq',
            'reads_file_1_md5': None,
            'reads_file_2': 'reads_2_2.fq',
            'reads_file_2_md5': None,
            'dataset_name': 'g2',
            'instrument_model': 'Illumina HiSeq 2000',
            'ena_center_name': 'Center 1',
            'submit_to_ena': '1',
            'ena_on_hold': '1',
            'ena_run_accession': None,
            'ena_sample_accession': None,
           },
        ]

        filename = os.path.join(data_dir, 'load_data_from_spreadsheet.tsv')
        got = spreadsheet_helper.load_data_from_spreadsheet(filename)
        self.maxDiff = None
        self.assertEqual(expected, got)


    def test_load_data_from_spreadsheet_xlsx(self):
        '''test load_data_from_spreadsheet with good input xlsx file'''
        expected = [
           {'subject_id': 'p1',
            'site_id': 's1',
            'lab_id': 'l1',
            'isolate_number': '42',
            'sequence_replicate_number': '43',
            'submission_date': datetime.date(2017, 12, 25),
            'reads_file_1': 'reads_1_1.fq',
            'reads_file_1_md5': 'abcdefghijklmnopqrstuvwyx123456',
            'reads_file_2': 'reads_1_2.fq',
            'reads_file_2_md5': 'abcdefghijklmnopqrstuvwyx123457',
            'dataset_name': 'g1',
            'instrument_model': 'Illumina HiSeq 2000',
            'ena_center_name': 'Center 1',
            'submit_to_ena': '0',
            'ena_on_hold': '0',
            'ena_run_accession': 'ERR123456',
            'ena_sample_accession': 'ERS123456',
           },
           {'subject_id': 'p2',
            'site_id': 's2',
            'lab_id': 'l2',
            'isolate_number': '44',
            'sequence_replicate_number': '45',
            'submission_date': datetime.date(2017, 12, 26),
            'reads_file_1': 'reads_2_1.fq',
            'reads_file_1_md5': None,
            'reads_file_2': 'reads_2_2.fq',
            'reads_file_2_md5': None,
            'dataset_name': 'g2',
            'instrument_model': 'Illumina HiSeq 2000',
            'ena_center_name': 'Center 1',
            'submit_to_ena': '1',
            'ena_on_hold': '1',
            'ena_run_accession': None,
            'ena_sample_accession': None,
           },
        ]

        filename = os.path.join(data_dir, 'load_data_from_spreadsheet.xlsx')
        got = spreadsheet_helper.load_data_from_spreadsheet(filename)
        self.assertEqual(expected, got)


    def test_load_data_from_spreadsheet_bad_input(self):
        '''Test load_data_from_spreadsheet with bad input files'''
        filenames = [
            'load_data_from_spreadsheet_bad_column_names.xlsx',
            'load_data_from_spreadsheet_wrong_field_number.xlsx',
        ]

        for filename in filenames:
            with self.assertRaises(spreadsheet_helper.Error):
                xlsx = os.path.join(data_dir, filename)
                spreadsheet_helper.load_data_from_spreadsheet(xlsx)


