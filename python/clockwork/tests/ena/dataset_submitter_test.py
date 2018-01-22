import unittest
import copy
import datetime
import os
import shutil
from clockwork import db, db_connection, db_maker, utils
from clockwork.ena import dataset_submitter

modules_dir = os.path.join(os.path.dirname(os.path.abspath(dataset_submitter.__file__)), os.pardir)
data_dir = os.path.normpath(os.path.join(modules_dir, 'tests', 'data', 'ena', 'dataset_submitter'))
ini_file = os.path.join(data_dir, 'config.ini')


class TestDatasetSubmitter(unittest.TestCase):
    def test_get_centres_from_ini_file(self):
        '''test _get_centres_from_ini_file'''
        tmp_ini_file = 'tmp.ena.dataset_submitter.ini'
        with open(tmp_ini_file, 'w') as f:
            print('[foo]', file=f)
            print('bar=baz', file=f)

        self.assertEqual({}, dataset_submitter.DatasetSubmitter._get_centres_from_ini_file(tmp_ini_file))

        with open(tmp_ini_file, 'a') as f:
            print('[sequencing_centres]', file=f)
            print('01=centre1', file=f)
            print('02=centre2', file=f)

        expected = {'01': 'centre1', '02': 'centre2'}
        self.assertEqual(expected, dataset_submitter.DatasetSubmitter._get_centres_from_ini_file(tmp_ini_file))
        os.unlink(tmp_ini_file)


    def test_get_key_from_ini_file(self):
        '''test _get_key_from_ini_file'''
        tmp_ini_file = 'tmp.ena.dataset_submitter.ini'
        with open(tmp_ini_file, 'w') as f:
            print('[foo]', file=f)
            print('bar=baz', file=f)

        self.assertEqual(None, dataset_submitter.DatasetSubmitter._get_key_from_ini_file(tmp_ini_file, 'spam', 'bar'))
        self.assertEqual(None, dataset_submitter.DatasetSubmitter._get_key_from_ini_file(tmp_ini_file, 'foo', 'baz'))
        self.assertEqual('baz', dataset_submitter.DatasetSubmitter._get_key_from_ini_file(tmp_ini_file, 'foo', 'bar'))
        os.unlink(tmp_ini_file)


    def test_get_broker_name_from_ini_file(self):
        '''test _get_broker_name_from_ini_file'''
        tmp_ini_file = 'tmp.ena.dataset_submitter.ini'
        with open(tmp_ini_file, 'w') as f:
            print('[foo]', file=f)
            print('bar=baz', file=f)

        self.assertEqual(None, dataset_submitter.DatasetSubmitter._get_broker_name_from_ini_file(tmp_ini_file))

        with open(tmp_ini_file, 'a') as f:
            print('[ena_login]', file=f)
            print('spam=eggs', file=f)

        self.assertEqual(None, dataset_submitter.DatasetSubmitter._get_broker_name_from_ini_file(tmp_ini_file))

        with open(tmp_ini_file, 'a') as f:
            print('broker_name = Broker Name', file=f)
        self.assertEqual('Broker Name', dataset_submitter.DatasetSubmitter._get_broker_name_from_ini_file(tmp_ini_file))
        os.unlink(tmp_ini_file)


    def test_ena_center_name_from_db_data(self):
        '''test _ena_center_name_from_db_data'''
        test_data = [{'ena_center_name': '01', 'foo': 'bar'}]
        self.assertEqual('01', dataset_submitter.DatasetSubmitter._ena_center_name_from_db_data(test_data))
        test_data.append({'ena_center_name': '01', 'spam': 'eggs'})
        self.assertEqual('01', dataset_submitter.DatasetSubmitter._ena_center_name_from_db_data(test_data))
        number_to_name_dict = {'02': 'centre2'}
        self.assertEqual('01', dataset_submitter.DatasetSubmitter._ena_center_name_from_db_data(test_data, number_to_name_dict=number_to_name_dict))
        number_to_name_dict['01'] = 'centre1'
        self.assertEqual('centre1', dataset_submitter.DatasetSubmitter._ena_center_name_from_db_data(test_data, number_to_name_dict=number_to_name_dict))
        test_data.append({'ena_center_name': '02', 'spam': 'eggs'})
        with self.assertRaises(dataset_submitter.Error):
            dataset_submitter.DatasetSubmitter._ena_center_name_from_db_data(test_data)


    def test_get_data_from_db(self):
        '''test _get_data_from_db'''
        pipeline_test_dir = os.path.join(data_dir, 'get_data_from_db')
        mysql_dump = os.path.join(pipeline_test_dir, 'mysql.dump')
        mysql_config_file = os.path.join(data_dir, 'db.cnf')
        db_config_data = db_connection.DbConnection._parse_config_file(ini_file)
        utils.syscall('mysql --defaults-file=' + mysql_config_file + ' -e "DROP DATABASE IF EXISTS ' + db_config_data['db'] + '; CREATE DATABASE ' + db_config_data['db'] + '"')
        utils.syscall('mysql --defaults-file=' + mysql_config_file + ' ' + db_config_data['db'] + ' < ' + mysql_dump)

        gsub = dataset_submitter.DatasetSubmitter(ini_file, 'g1', 'pipeline_root', 42)
        got = gsub._get_data_from_db()

        expected = [
           {'ena_center_name': 'Center 42',
            'ena_experiment_accession': None,
            'ena_run_accession': None,
            'ena_sample_accession': None,
            'ena_study_accession': None,
            'instrument_model': 'Illumina HiSeq 2000',
            'isolate_id': 1,
            'isolate_number_from_lab': '42',
            'remove_contam_reads_file_1_md5': 'c40e0ce9e03810a2fd81fe7edf86e4c6',
            'remove_contam_reads_file_2_md5': 'b464216441da0d3dd539b583a93c635d',
            'sample_id': 1,
            'subject_id': 'p1',
            'sample_id_from_lab': 'l1',
            'seqrep_id': 1,
            'sequence_replicate_number': 43,
            'site_id': 's1'},
           {'ena_center_name': 'Center 42',
            'ena_experiment_accession': None,
            'ena_run_accession': None,
            'ena_sample_accession': None,
            'ena_study_accession': None,
            'instrument_model': 'Illumina HiSeq 2000',
            'isolate_id': 2,
            'isolate_number_from_lab': '44',
            'remove_contam_reads_file_1_md5': '42b6872fac4f221325466583de464919',
            'remove_contam_reads_file_2_md5': '247e5d7cb8feefc95045b3e598dae245',
            'sample_id': 2,
            'subject_id': 'p2',
            'sample_id_from_lab': 'l1',
            'seqrep_id': 2,
            'sequence_replicate_number': 1,
            'site_id': 's1'},
           {'ena_center_name': 'Center 42',
            'ena_experiment_accession': None,
            'ena_run_accession': None,
            'ena_sample_accession': None,
            'ena_study_accession': None,
            'instrument_model': 'Illumina HiSeq 2000',
            'isolate_id': 2,
            'isolate_number_from_lab': '44',
            'remove_contam_reads_file_1_md5': 'ec82f1f67913b8d817119c304d2ef638',
            'remove_contam_reads_file_2_md5': '4cbd5d2980198ad9e87bf2968cddfd5a',
            'sample_id': 2,
            'subject_id': 'p2',
            'sample_id_from_lab': 'l1',
            'seqrep_id': 3,
            'sequence_replicate_number': 2,
            'site_id': 's1'}
        ]
        self.maxDiff = None
        self.assertEqual(expected, got)


    def test_run(self):
        '''test run'''
        original_pipeline_root = os.path.join(data_dir, 'run', 'Pipeline_root')
        tmp_pipeline_root = 'tmp.dataset_submitter.pipeline_root'
        shutil.copytree(original_pipeline_root, tmp_pipeline_root)
        pipeline_test_dir = os.path.join(data_dir, 'run')
        mysql_dump = os.path.join(pipeline_test_dir, 'mysql.dump')
        mysql_config_file = os.path.join(data_dir, 'db.cnf')
        db_config_data = db_connection.DbConnection._parse_config_file(ini_file)
        utils.syscall('mysql --defaults-file=' + mysql_config_file + ' -e "DROP DATABASE IF EXISTS ' + db_config_data['db'] + '; CREATE DATABASE ' + db_config_data['db'] + '"')
        utils.syscall('mysql --defaults-file=' + mysql_config_file + ' ' + db_config_data['db'] + ' < ' + mysql_dump)

        gsub = dataset_submitter.DatasetSubmitter(ini_file, 'g1', tmp_pipeline_root, 42, unit_test='success')
        gsub.run()
        columns = 'Seqrep.seqrep_id, sequence_replicate_number, remove_contam_reads_file_1_md5, remove_contam_reads_file_2_md5, ena_center_name, ena_run_accession, Isolate.isolate_id, ena_experiment_accession, Sample.sample_id, site_id, instrument_model, ena_sample_accession, ena_study_accession'
        join = 'Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id JOIN Sample ON Isolate.sample_id = Sample.sample_id'
        where = 'submit_to_ena=1 AND import_status=1 AND dataset_name="g1"'
        query = ' '.join(['SELECT', columns, 'FROM', '(' + join + ')', 'WHERE', '(' + where + ')'])
        database = db.Db(ini_file)
        got_data = database.query_to_dict(query)

        for row in got_data:
            accessions = {row[x] for x in row if x.endswith('accession')}
            self.assertNotIn(None, accessions)

        run_accessions = {x['ena_run_accession'] for x in got_data}
        self.assertEqual(5, len(run_accessions))
        study_accessions = {x['ena_study_accession'] for x in got_data}
        self.assertEqual(1, len(study_accessions))

        # hash the rows by md5 of file 1, since we don't know the auto
        # generated IDs in the DB.
        data_by_md5_1 = {x['remove_contam_reads_file_1_md5']: x for x in got_data}
        md51 = '83d842db2d9ea84faa747cefa4b2f1b4'
        md52 = '67ff4c03bd637e027f372b4b5a833935'
        md53 = 'bfde82c3a5ec16ffefb32fdfcfd4cf53'
        md54 = 'be5c2e07716c119a2e86f6421df5f63b'
        md55 = '21544f51d9d620ca99bc445219b1018d'
        self.assertNotEqual(data_by_md5_1[md51]['ena_sample_accession'], data_by_md5_1[md52]['ena_sample_accession'])
        self.assertEqual(data_by_md5_1[md52]['ena_sample_accession'], data_by_md5_1[md53]['ena_sample_accession'])
        self.assertEqual(data_by_md5_1[md53]['ena_sample_accession'], data_by_md5_1[md54]['ena_sample_accession'])
        self.assertNotEqual(data_by_md5_1[md54]['ena_sample_accession'], data_by_md5_1[md55]['ena_sample_accession'])
        self.assertNotEqual(data_by_md5_1[md51]['ena_experiment_accession'], data_by_md5_1[md52]['ena_experiment_accession'])
        self.assertEqual(data_by_md5_1[md52]['ena_experiment_accession'], data_by_md5_1[md53]['ena_experiment_accession'])
        self.assertNotEqual(data_by_md5_1[md53]['ena_experiment_accession'], data_by_md5_1[md54]['ena_experiment_accession'])
        self.assertNotEqual(data_by_md5_1[md54]['ena_experiment_accession'], data_by_md5_1[md55]['ena_experiment_accession'])
        shutil.rmtree(tmp_pipeline_root)
