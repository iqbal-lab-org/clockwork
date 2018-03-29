import unittest
import glob
import sys
import shutil
import os
from operator import itemgetter

from clockwork import db, db_connection, isolate_dir, utils
from clockwork import __version__ as clockwork_version

sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
import nextflow_helper
data_dir = os.path.join(nextflow_helper.data_root_dir, 'nextflow_mykrobe_predict')
modules_dir = os.path.dirname(os.path.abspath(db.__file__))
db_ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')



class TestNextflowMykrobe(unittest.TestCase):
    def test_nextflow_mykrobe_predict(self):
        '''test nextflow_mykrobe using database'''
        tmp_data_dir = 'tmp.nextflow_mykrobe_db_input.data'
        if os.path.exists(tmp_data_dir):
            shutil.rmtree(tmp_data_dir)
        shutil.copytree(data_dir, tmp_data_dir)
        nextflow_helper.write_config_file()
        mysql_config_file = os.path.join(data_dir, 'db.cnf')
        mysql_dump = os.path.join(data_dir, 'mysql.dump')
        db_config_data = db_connection.DbConnection._parse_config_file(db_ini_file)
        utils.syscall('mysql --defaults-file=' + mysql_config_file + ' -e "DROP DATABASE IF EXISTS ' + db_config_data['db'] + '; CREATE DATABASE ' + db_config_data['db'] + '"')
        utils.syscall('mysql --defaults-file=' + mysql_config_file + ' ' + db_config_data['db'] + ' < ' + mysql_dump)
        pipeline_root = os.path.join(tmp_data_dir, 'Pipeline_root')
        references_root = os.path.join(tmp_data_dir, 'Pipeline_refs')
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'mykrobe_predict.nf')
        work_dir = 'tmp.nextflow_mykrobe_db_input.work'
        dag_file = 'nextflow.mykrobe.dag.db.pdf'
        try:
            os.unlink(dag_file)
        except:
            pass

        command = ' '.join([
            'nextflow run',
            '--dataset_name g1', # one read pair is from group 2 and should get ignored
            '--ref_id 2',
            '--references_root', os.path.abspath(references_root),
            '--pipeline_root', pipeline_root,
            '--db_config_file', db_ini_file,
            '--testing',
            '-with-dag', dag_file,
            '-c', nextflow_helper.config_file,
            '-w', work_dir,
            nextflow_file
        ])
        utils.syscall(command)
        os.unlink(nextflow_helper.config_file)
        shutil.rmtree(work_dir)

        # check database Pipeline table updated as expected.
        # The --testing option is set up so that the pooled
        # sample fails, hence it gets a status of -1.
        database = db.Db(db_ini_file)
        got_rows = database.get_rows_from_table('Pipeline')
        got_rows.sort(key=itemgetter('isolate_id', 'pipeline_name'))
        expected_rows = [
            {'isolate_id': 1, 'seqrep_id': None, 'seqrep_pool': '1_2', 'version': clockwork_version, 'pipeline_name': 'mykrobe_predict', 'status': -1, 'reference_id': 2},
            {'isolate_id': 1, 'seqrep_id': 1, 'seqrep_pool': None, 'version': '0.4.0', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 1, 'seqrep_id': 2, 'seqrep_pool': None, 'version': '0.4.0', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 3, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'mykrobe_predict', 'status': 1, 'reference_id': 2},
            {'isolate_id': 2, 'seqrep_id': 3, 'seqrep_pool': None, 'version': '0.4.0', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 4, 'seqrep_pool': None, 'version': '0.4.0', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 4, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'mykrobe_predict', 'status': 1, 'reference_id': 2},
            {'isolate_id': 3, 'seqrep_id': None, 'seqrep_pool': '1', 'version': clockwork_version, 'pipeline_name': 'mykrobe_predict', 'status': 1, 'reference_id': 2},
            {'isolate_id': 3, 'seqrep_id': 5, 'seqrep_pool': None, 'version': '0.4.0', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 4, 'seqrep_id': 6, 'seqrep_pool': None, 'version': '0.4.0', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
        ]
        expected_rows.sort(key=itemgetter('isolate_id', 'pipeline_name'))
        self.assertEqual(expected_rows, got_rows)

        # check mykrobe output files etc got written. No need to check contents, trust the tools
        # We're just checking nextflow runs OK here.
        ids = [
            {'sample': 1, 'seqrep_id': '1_2', 'isolate_id': 1, 'seq_repl': '1_2', 'sample_name': 'site.s1.iso.42.subject.p1.lab_id.l1.seq_reps.1_2'},
            {'sample': 2, 'seqrep_id': 3, 'isolate_id': 2, 'seq_repl': '1', 'sample_name': 'site.s2.iso.43.subject.p2.lab_id.l2.seq_reps.1'},
            {'sample': 2, 'seqrep_id': 4, 'isolate_id': 2, 'seq_repl': '2', 'sample_name': 'site.s2.iso.43.subject.p2.lab_id.l2.seq_reps.2'},
        ]
        for id_dict in ids:
            iso_dir = isolate_dir.IsolateDir(pipeline_root, id_dict['sample'], id_dict['isolate_id'])
            pipeline_dir = iso_dir.pipeline_dir(id_dict['seq_repl'], 'mykrobe_predict', clockwork_version, reference_id=2)
            self.assertTrue(os.path.exists(pipeline_dir))
            log = os.path.join(pipeline_dir, 'log.txt')
            json_file = os.path.join(pipeline_dir, 'out.json')

            if id_dict['sample_name'].endswith('1_2'):
                self.assertFalse(os.path.exists(log))
                self.assertFalse(os.path.exists(json_file))
            else:
                self.assertTrue(os.path.exists(log))
                self.assertTrue(os.path.exists(json_file))

        shutil.rmtree(tmp_data_dir)
        nextflow_helper.clean_files()

