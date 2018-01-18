import unittest
import os
from clockwork import spreadsheet_validator, utils

modules_dir = os.path.dirname(os.path.abspath(spreadsheet_validator.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'spreadsheet_validator')
db_ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')


class TestSpreadsheetValidator(unittest.TestCase):
    def test_check_no_blank_values(self):
        '''test _check_no_blank_values'''
        dict_list = [{'a': 1, 'b':2}]
        expected = []
        self.assertEqual([], spreadsheet_validator.SpreadsheetValidator._check_no_blank_values(dict_list))
        dict_list.append({'a': None, 'b':2})
        expected.append('Empty_field\ta\t3')
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_no_blank_values(dict_list))
        dict_list.append({'a': '', 'b':None})
        expected.append('Empty_field\ta\t4')
        expected.append('Empty_field\tb\t4')
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_no_blank_values(dict_list))


    def test_get_non_unique_dict_list_values(self):
        '''test _get_non_unique_dict_list_values'''
        dict_list = [
            {'a': 1, 'b': 2, 'c': 3},
            {'a': 2, 'b': 2, 'c': 42},
            {'a': 3, 'b': 3, 'c': 43},
            {'a': 4, 'b': 3, 'c': 44},
            {'a': 5, 'b': 3, 'c': 44},
        ]

        self.assertEqual({}, spreadsheet_validator.SpreadsheetValidator._get_non_unique_dict_list_values(dict_list, key='a'))
        self.assertEqual({2: [0,1], 3: [2,3,4]}, spreadsheet_validator.SpreadsheetValidator._get_non_unique_dict_list_values(dict_list, key='b'))
        self.assertEqual({44: [3,4]}, spreadsheet_validator.SpreadsheetValidator._get_non_unique_dict_list_values(dict_list, key='c'))
        self.assertEqual({}, spreadsheet_validator.SpreadsheetValidator._get_non_unique_dict_list_values(dict_list, key_list=['a', 'b']))
        self.assertEqual({(3,44): [3,4]}, spreadsheet_validator.SpreadsheetValidator._get_non_unique_dict_list_values(dict_list, key_list=['b', 'c']))
        self.assertEqual({(44,3): [3,4]}, spreadsheet_validator.SpreadsheetValidator._get_non_unique_dict_list_values(dict_list, key_list=['c', 'b']))


    def test_check_uniqueness_of_values(self):
        '''test _check_uniqueness_of_values'''
        dict_list = [
            {
              'reads_file_1': 'f1.1',
              'reads_file_1_md5': 'md51.1',
              'reads_file_2': 'f1.2',
              'reads_file_2_md5': 'md51.2',
              'subject_id': 'sub1',
              'site_id': 'site1',
              'lab_id': 'l1',
              'isolate_number': 'i1',
              'sequence_replicate_number': '1',
            }
        ]

        expected = []
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_uniqueness_of_values(dict_list))

        dict_list.append({
              'reads_file_1': 'f2.1',
              'reads_file_1_md5': 'md52.1',
              'reads_file_2': 'f2.2',
              'reads_file_2_md5': 'md52.2',
              'subject_id': 'sub1',
              'site_id': 'site1',
              'lab_id': 'l1',
              'isolate_number': 'i1',
              'sequence_replicate_number': '2',
        })
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_uniqueness_of_values(dict_list))

        dict_list[1]['reads_file_1'] = 'f1.1'
        expected.append('Non_unique\treads_file_1\tf1.1\t2,3')
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_uniqueness_of_values(dict_list))

        dict_list[1]['sequence_replicate_number'] = '1'
        expected.append('Non_unique\tsubject_id,site_id,lab_id,isolate_number,sequence_replicate_number\tsub1,site1,l1,i1,1\t2,3')
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_uniqueness_of_values(dict_list))


    def test_check_global_file_and_md5_column_intersection(self):
        '''test _check_global_file_and_md5_column_intersection'''
        dict_list = [
            {
              'reads_file_1': 'f1.1',
              'reads_file_1_md5': 'md51.1',
              'reads_file_2': 'f1.2',
              'reads_file_2_md5': 'md51.2',
            },
            {
              'reads_file_1': 'f2.1',
              'reads_file_1_md5': 'md52.1',
              'reads_file_2': 'f2.2',
              'reads_file_2_md5': 'md52.2',
            },
        ]
        expected = []
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_global_file_and_md5_column_intersection(dict_list))
        dict_list[1]['reads_file_2'] = 'f1.1'
        expected.append('Filename_in_both_columns\tf1.1')
        dict_list[0]['reads_file_1_md5'] = 'md52.2'


    def test_check_files_exist_and_md5(self):
        '''test _check_files_exist_and_md5'''
        file1 = 'tmp.spreadsheet_validator.file1'
        file2 = 'tmp.spreadsheet_validator.file2'
        with open(file1, 'w') as f1, open(file2, 'w') as f2:
            print('foo', file=f1)
            print('bar', file=f2)
        reads_file_1_md5 = utils.md5(file1)
        dict_list = [{
            'reads_file_1': file1,
            'reads_file_1_md5': reads_file_1_md5,
            'reads_file_2': file2,
            'reads_file_2_md5': utils.md5(file2),
        }]
        self.assertEqual([], spreadsheet_validator.SpreadsheetValidator._check_files_exist_and_md5(dict_list, os.getcwd(), md5_threads=1))
        self.assertEqual([], spreadsheet_validator.SpreadsheetValidator._check_files_exist_and_md5(dict_list, os.getcwd(), md5_threads=2))
        dict_list[0]['reads_file_1_md5'] = '12345'
        expected = ['md5_error\t' + file1 + '\t12345\t' + reads_file_1_md5]
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_files_exist_and_md5(dict_list, os.getcwd(), md5_threads=1))
        self.assertEqual([], spreadsheet_validator.SpreadsheetValidator._check_files_exist_and_md5(dict_list, os.getcwd(), md5_threads=0))
        os.unlink(file1)
        expected = ['File_not_found\t' + file1 + '\t2']
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_files_exist_and_md5(dict_list, os.getcwd()))
        os.unlink(file2)
        expected.append('File_not_found\t' + file2 + '\t2')
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_files_exist_and_md5(dict_list, os.getcwd()))


    def test_check_integers(self):
        '''test _check_integers'''
        dict_list = [{'a': 'oops'}, {'a': '-1'}, {'a': '0'}, {'a': '1'}, {'a': '2'}]
        expected = [
            'Not_an_integer\ta\toops\t2',
            'Integer_out_of_range\ta\t-1\t3',
            'Integer_out_of_range\ta\t2\t6',
        ]
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_integers(dict_list, 'a', min_value=0, max_value=1))


    def test_check_instrument_model(self):
        '''test _check_instrument_model'''
        dict_list = [{'instrument_model': 'oops'}, {'instrument_model': 'NextSeq 500'}]
        expected = ['Unknown_instrument_model\toops\t2']
        self.assertEqual(expected, spreadsheet_validator.SpreadsheetValidator._check_instrument_model(dict_list))


    def test_run(self):
        '''test run'''
        # We won't exhaustively check we get all errors, because that
        # is tested for each individual function. Just check it reads a file
        # and write a plausible output file
        xlsx_file_good = os.path.join(data_dir, 'run.in.good.xlsx')
        xlsx_file_bad = os.path.join(data_dir, 'run.in.bad.xlsx')
        tmp_out = 'tmp.spreadsheet_validator.out.txt'
        if os.path.exists(tmp_out):
            os.unlink(tmp_out)
        validator = spreadsheet_validator.SpreadsheetValidator(xlsx_file_good, data_dir, tmp_out)
        validator.run()
        with open(tmp_out) as f:
            got_lines = [x.rstrip() for x in f]
        self.assertEqual([], got_lines)

        validator = spreadsheet_validator.SpreadsheetValidator(xlsx_file_bad, data_dir, tmp_out)
        validator.run()
        with open(tmp_out) as f:
            got_lines = [x.rstrip() for x in f]
        self.assertEqual(['File_not_found\treads_2_3.fq\t3'], got_lines)
        os.unlink(tmp_out)

