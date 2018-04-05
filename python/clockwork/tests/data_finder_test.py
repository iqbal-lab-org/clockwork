import datetime
import filecmp
import os
import shutil
import unittest

from clockwork import db, db_connection, db_maker, data_finder

modules_dir = os.path.dirname(os.path.abspath(data_finder.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'data_finder')
ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')


# Pipeline root directory is always absolute path, which we only
# know at run time. Thie function replaces every PIPEROOT in the
# input file with the actual pipeline root dir, then can use
# the output file as the correct expected file
def fix_pipeline_root_in_file(infile, outfile, pipeline_root):
    with open(infile) as f_in, open(outfile, 'w') as f_out:
        for line in f_in:
            print(line.replace('PIPEROOT', pipeline_root), end='', file=f_out)


class TestDataFinder(unittest.TestCase):
    def setUp(self):
        self.pipeline_root = os.path.abspath('piperoot')
        os.mkdir(self.pipeline_root)

        try:
            db_connection.DbConnection(ini_file, destroy=True)
        except:
            pass

        dbm = db_maker.DbMaker(ini_file)
        dbm.run()
        self.db = db.Db(ini_file)

        sample_dicts = [
          {
            'subject_id': 'subject_1',
            'site_id': '01',
            'lab_id': 'lab_id_1',
            'isolate_number': '1',
            'sequence_replicate_number': 1,
            'submission_date': datetime.date(2018, 4, 4),
            'reads_file_1': 'reads_1_1.fq',
            'reads_file_1_md5': 'md5_1_1',
            'reads_file_2_md5': 'md5_1_2',
            'reads_file_2': 'reads_1_2.fq',
            'dataset_name': 'set1',
            'submit_to_ena': '0',
            'instrument_model': 'Illumina HiSeq 2500',
            'ena_center_name': 'Centre 1',
            'ena_on_hold': '0',
            'ena_run_accession': 'ERR123456',
            'ena_sample_accession': 'ERS123456',
          },
          {
            'subject_id': 'subject_2',
            'site_id': '01',
            'lab_id': 'lab_id_2',
            'isolate_number': '1',
            'sequence_replicate_number': 1,
            'submission_date': datetime.date(2018, 4, 4),
            'reads_file_1': 'reads_2_1.fq',
            'reads_file_1_md5': 'md5_2_1',
            'reads_file_2_md5': 'md5_2_2',
            'reads_file_2': 'reads_2_2.fq',
            'dataset_name': 'set1',
            'submit_to_ena': '0',
            'instrument_model': 'Illumina HiSeq 2500',
            'ena_center_name': 'Centre 1',
            'ena_on_hold': '0',
            'ena_run_accession': 'ERR123457',
            'ena_sample_accession': 'ERS123457',
          },
          {
            'subject_id': 'subject_3',
            'site_id': '02',
            'lab_id': 'lab_id_3',
            'isolate_number': '1',
            'sequence_replicate_number': 1,
            'submission_date': datetime.date(2018, 4, 4),
            'reads_file_1': 'reads_3_1.fq',
            'reads_file_1_md5': 'md5_3_1',
            'reads_file_2_md5': 'md5_3_2',
            'reads_file_2': 'reads_3_2.fq',
            'dataset_name': 'set2',
            'submit_to_ena': '0',
            'instrument_model': 'Illumina HiSeq 2500',
            'ena_center_name': 'Centre 2',
            'ena_on_hold': '0',
            'ena_run_accession': None,
            'ena_sample_accession': None,
          },
          {
            'subject_id': 'subject_3',
            'site_id': '02',
            'lab_id': 'lab_id_3',
            'isolate_number': '1',
            'sequence_replicate_number': 2,
            'submission_date': datetime.date(2018, 4, 4),
            'reads_file_1': 'reads_4_1.fq',
            'reads_file_1_md5': 'md5_4_1',
            'reads_file_2_md5': 'md5_4_2',
            'reads_file_2': 'reads_4_2.fq',
            'dataset_name': 'set2',
            'submit_to_ena': '0',
            'instrument_model': 'Illumina HiSeq 2500',
            'ena_center_name': 'Centre 2',
            'ena_on_hold': '0',
            'ena_run_accession': None,
            'ena_sample_accession': None,
          },
        ]

        for d in sample_dicts:
            self.db.add_one_seqrep(d)
            where_dict = {'original_reads_file_1_md5': d['reads_file_1_md5']}
            update_dict = {
                'remove_contam_reads_file_1_md5': d['reads_file_1_md5'] + '.remove_contam',
                'remove_contam_reads_file_2_md5': d['reads_file_2_md5'] + '.remove_contam',
            }
            self.db.update_row('Seqrep', where_dict, update_dict)


        seqrep_to_isolate = {1: 1, 2: 2, 3: 3, 4: 3}
        for seqrep, isolate in seqrep_to_isolate.items():
            ref_id = 1 if seqrep in {1,2} else 2
            version = '0.1.1' if seqrep in {1,2} else '0.1.3'
            d = {'isolate_id': isolate, 'seqrep_id': seqrep, 'seqrep_pool': None,
                 'version': version, 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': ref_id}
            self.db.add_row_to_table('Pipeline', d)
            d = {'isolate_id': isolate, 'seqrep_id': seqrep, 'seqrep_pool': None,
                 'version': version, 'pipeline_name': 'qc', 'status': 1, 'reference_id': ref_id + 2}
            self.db.add_row_to_table('Pipeline', d)

        var_call_rows = [
            {'isolate_id': 1, 'seqrep_id': None, 'seqrep_pool': '1',
             'version': '1.2.3', 'pipeline_name': 'variant_call', 'status': 1, 'reference_id': 10},
            {'isolate_id': 2, 'seqrep_id': None, 'seqrep_pool': '2',
             'version': '1.2.3', 'pipeline_name': 'variant_call', 'status': 1, 'reference_id': 10},
            {'isolate_id': 3, 'seqrep_id': None, 'seqrep_pool': '1_2',
             'version': '1.2.3', 'pipeline_name': 'variant_call', 'status': 1, 'reference_id': 10},
        ]
        for d in var_call_rows:
            self.db.add_row_to_table('Pipeline', d)


        self.db.commit()


    def tearDown(self):
        self.db.commit_and_close()
        shutil.rmtree(self.pipeline_root)


    def test_write_seqrep_data_to_file(self):
        '''test write_seqrep_data_to_file'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root)
        tmpfile = 'tmp.data_finder.write_seqrep_data_to_file.out'
        expected_file = os.path.join(data_dir, 'write_seqrep_data_to_file.tsv')
        tmp_expected = 'tmp.data_finder.write_seqrep_data_to_file.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_seqrep_data_to_file(tmpfile)
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_seqrep_data_to_file_dataset_option(self):
        '''test write_seqrep_data_to_file limit to one dataset'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root, dataset_name='set2')
        tmpfile = 'tmp.data_finder.write_seqrep_data_to_file.set2.out'
        expected_file = os.path.join(data_dir, 'write_seqrep_data_to_file.set2.tsv')
        tmp_expected = 'tmp.data_finder.write_seqrep_data_to_file.set2.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_seqrep_data_to_file(tmpfile)
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_remove_contam_pipeline(self):
        '''test write_pipeline_data_to_file for remove_contam pipeline'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root)
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.out'
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.remove_contam.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_pipeline_data_to_file(tmpfile, 'remove_contam')
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_remove_contam_pipeline_with_internal_ids(self):
        '''test write_pipeline_data_to_file for remove_contam pipeline with internal ids'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root, include_internal_ids=True)
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.with_internal_ids.out'
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.remove_contam.with_internal_ids.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_pipeline_data_to_file(tmpfile, 'remove_contam')
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_remove_contam_pipeline_dataset_option(self):
        '''test write_pipeline_data_to_file for remove_contam pipeline limit to one dataset'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root, dataset_name='set2')
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.set2.out'
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.remove_contam.set2.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.set2.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_pipeline_data_to_file(tmpfile, 'remove_contam')
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_remove_contam_pipeline_version_option(self):
        '''test write_pipeline_data_to_file for remove_contam pipeline use version option'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root)
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.v0.1.1.out'
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.remove_contam.v0.1.1.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.v0.1.1.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_pipeline_data_to_file(tmpfile, 'remove_contam', pipeline_version='0.1.1')
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_remove_contam_pipeline_reference_option(self):
        '''test write_pipeline_data_to_file for remove_contam pipeline use reference option'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root)
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.ref.1.out'
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.remove_contam.ref.1.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.remove_contam.ref.1.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        finder.write_pipeline_data_to_file(tmpfile, 'remove_contam', reference_id=1)
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_qc_pipeline(self):
        '''test write_pipeline_data_to_file for qc pipeline'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root)
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.qc.out'
        finder.write_pipeline_data_to_file(tmpfile, 'qc')
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.qc.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.qc.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)


    def test_write_pipeline_data_to_file_variant_call_pipeline(self):
        '''test write_pipeline_data_to_file for variant_call pipeline'''
        finder = data_finder.DataFinder(ini_file, self.pipeline_root)
        tmpfile = 'tmp.data_finder.write_pipeline_data_to_file.variant_call.out'
        finder.write_pipeline_data_to_file(tmpfile, 'variant_call')
        expected_file = os.path.join(data_dir, 'write_pipeline_data_to_file.variant_call.tsv')
        tmp_expected = 'tmp.data_finder.write_pipeline_data_to_file.variant_call.expect'
        fix_pipeline_root_in_file(expected_file, tmp_expected, self.pipeline_root)
        self.assertTrue(filecmp.cmp(tmp_expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
        os.unlink(tmp_expected)
