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
data_dir = os.path.join(nextflow_helper.data_root_dir, 'nextflow_remove_contam')
modules_dir = os.path.dirname(os.path.abspath(db.__file__))
db_ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')


class TestNextflowRemoveContam(unittest.TestCase):
    def test_nextflow_remove_contam_using_database(self):
        '''test nextflow_remove_contam using database'''
        tmp_data_dir = 'tmp.nextflow_remove_contam'
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
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'remove_contam.nf')
        work_dir = 'tmp.nextflow_remove_contam.work'
        dag_file = 'nextflow.remove_contam.dag.db.pdf'
        try:
            os.unlink(dag_file)
        except:
            pass

        command = ' '.join([
            'nextflow run',
            '--dataset_name g1', # one read pair has group g2, so should get ignored
            '--ref_id 1',
            '--references_root', os.path.abspath(references_root),
            '--pipeline_root', os.path.abspath(pipeline_root),
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

        # check database Pipeline table updated as expected
        database = db.Db(db_ini_file)
        got_rows = database.get_rows_from_table('Pipeline')
        got_rows.sort(key=itemgetter('seqrep_id'))
        expected_rows = [
                {'isolate_id': 1, 'seqrep_id': 1, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
                {'isolate_id': 2, 'seqrep_id': 2, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
                {'isolate_id': 3, 'seqrep_id': 3, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'remove_contam', 'status': -1, 'reference_id': 1},
        ]
        self.assertEqual(expected_rows, got_rows)

        # check database Read_counts table updated
        got_rows = database.get_rows_from_table('Read_counts')
        got_rows.sort(key=itemgetter('seqrep_id'))
        expected_rows = [
          {
            'seqrep_id': 1,
            'original_total': 198,
            'contamination': 40,
            'not_contamination': 132,
            'unmapped': 26,
            'total_after_remove_contam': 158,
          },
          {
            'seqrep_id': 2,
            'original_total': 156,
            'contamination': 12,
            'not_contamination': 132,
            'unmapped': 12,
            'total_after_remove_contam': 144,
          },
        ]
        self.assertEqual(expected_rows, got_rows)

        # check FASTQ files got written. No need to check contents, as that is done
        # elsewhere. We're just checking nextflow runs OK here.
        ids = [
            {'sample': 1, 'isolate_id': 1, 'seq_repl': 43},
            {'sample': 2, 'isolate_id': 2, 'seq_repl': 45},
        ]
        for id_dict in ids:
            iso_dir = isolate_dir.IsolateDir(pipeline_root, id_dict['sample'], id_dict['isolate_id'])
            for read_type in ('original', 'remove_contam', 'contam'):
                for i in (1, 2):
                    self.assertTrue(os.path.exists(iso_dir.reads_filename(read_type, id_dict['seq_repl'], i)))

        shutil.rmtree(tmp_data_dir)
        nextflow_helper.clean_files()


    def test_nextflow_remove_contam_using_fastq_input(self):
        '''test nextflow_remove_contam using fastq input'''
        reads1 = os.path.join(data_dir, 'Reads', 'reads.1.1.fq.gz')
        reads2 = os.path.join(data_dir, 'Reads', 'reads.1.2.fq.gz')
        outprefix = 'tmp.test_nextflow_remove_contam_using_fastq_input'
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'remove_contam.nf')
        nextflow_helper.write_config_file()
        work_dir = 'tmp.nextflow_remove_contam.work'
        dag_file = 'nextflow.remove_contam.dag.no_db.pdf'
        try:
            os.unlink(dag_file)
        except:
            pass

        command = ' '.join([
            'nextflow run',
            '--reads_in1', reads1,
            '--reads_in2', reads2,
            '--outprefix', outprefix,
            '--ref_metadata_tsv', os.path.join(data_dir, 'Reference', 'remove_contam_metadata.tsv'),
            '--ref_fasta', os.path.join(data_dir, 'Reference', 'ref.fa'),
            '--testing',
            '-with-dag', dag_file,
            '-c', nextflow_helper.config_file,
            '-w', work_dir,
            nextflow_file
        ])
        utils.syscall(command)
        os.unlink(nextflow_helper.config_file)
        shutil.rmtree(work_dir)

        for reads_type in ('contam', 'remove_contam'):
            for i in ('1', '2'):
                filename = outprefix + '.' + reads_type + '.' + i + '.fq.gz'
                self.assertTrue(os.path.exists(filename))
                os.unlink(filename)


        expected_counts_lines = [
            'Name\tIs_contam\tReads\n',
            'contam\t1\t40\n',
            'ref\t0\t132\n',
            'Unmapped\t0\t26\n',
            'Reads_kept_after_remove_contam\t0\t158\n',
        ]

        counts_tsv = outprefix + '.counts.tsv'
        with open(counts_tsv) as f:
            got_counts_lines = f.readlines()
        self.assertEqual(expected_counts_lines, got_counts_lines)
        os.unlink(counts_tsv)

        nextflow_helper.clean_files()
