import os
import datetime
import filecmp
import unittest
from clockwork import utils

modules_dir = os.path.dirname(os.path.abspath(utils.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'utils')


class TestUtils(unittest.TestCase):
    def test_syscall_no_error(self):
        '''test syscall with no error'''
        got = utils.syscall("echo testing 123")
        self.assertEqual("testing 123\n", got.stdout)


    def test_syscall_with_error(self):
        '''test syscall when there is an error'''
        with self.assertRaises(utils.Error):
            utils.syscall('notacommandunlessyoumadeitone')


    def test_md5(self):
        '''test md5'''
        filename = os.path.join(data_dir, 'md5.txt')
        expected = 'f738cdb99c0c6a670aa1e5742066b0cf'
        got = utils.md5(filename)
        self.assertEqual(expected, got)


    def test_load_md5_from_file(self):
        '''test load_md5_from_file'''
        expected = '43247f482b82e38a190c4d3243f97ea8'
        prefix = os.path.join(data_dir, 'load_md5_from_file.')
        self.assertEqual(expected, utils.load_md5_from_file(prefix + 'good_mac'))
        self.assertEqual(expected, utils.load_md5_from_file(prefix + 'good_linux'))
        with self.assertRaises(utils.Error):
            utils.load_md5_from_file(prefix + 'bad_mac')
        with self.assertRaises(utils.Error):
            utils.load_md5_from_file(prefix + 'bad_linux')


    def test_rsync_and_md5(self):
        '''test rsync_and_md5'''
        old_name = os.path.join(data_dir, 'rsync_and_md5.txt')
        new_name = 'tmp.test_rsync_and_md5'
        got = utils.rsync_and_md5(old_name, new_name)
        expected = 'a00096ee7316167c6ceef09f43433667'
        self.assertEqual(expected, got)
        self.assertTrue(filecmp.cmp(old_name, new_name, shallow=False))
        os.unlink(new_name)


    def test_date_string_from_file_mtime(self):
        '''test date_string_from_file_mtime'''
        tmpfile = 'tmp.test_date_string_from_file_mtime'
        with open(tmpfile, 'w'):
            pass

        today = datetime.date.today()
        got = utils.date_string_from_file_mtime(tmpfile)
        def int_to_str(x):
            if x < 10:
                return '0' + str(x)
            else:
                return str(x)

        self.assertEqual(str(today.year), got[0:4])
        self.assertEqual(int_to_str(today.month), got[4:6])
        self.assertEqual(int_to_str(today.day), got[6:])
        os.unlink(tmpfile)


    def test_make_empty_file(self):
        '''test make_empty_file'''
        tmpfile = 'tmp.test_make_empty_file'''
        self.assertFalse(os.path.exists(tmpfile))
        utils.make_empty_file(tmpfile)
        self.assertTrue(os.path.exists(tmpfile))
        with open(tmpfile) as f:
            lines = f.readlines()
            self.assertEqual([], lines)

        with open(tmpfile, 'w') as f:
            print('spam spam wonderful spam', file=f)

        utils.make_empty_file(tmpfile)
        self.assertTrue(os.path.exists(tmpfile))
        with open(tmpfile) as f:
            lines = f.readlines()
            self.assertEqual([], lines)

        os.unlink(tmpfile)


    def test_sam_record_count(self):
        '''test sam_record_count'''
        tmpfile = 'tmp.utils.sam_record_count'
        with open(tmpfile, 'w') as f:
            print('@foo bar', file=f)
            print('@PG baz', file=f)
            for i in range(42):
                print('zaphod', file=f)

        self.assertEqual(42, utils.sam_record_count(tmpfile))
        os.unlink(tmpfile)

