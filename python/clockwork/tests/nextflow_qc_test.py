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
data_dir = os.path.join(nextflow_helper.data_root_dir, 'nextflow_qc')
modules_dir = os.path.dirname(os.path.abspath(db.__file__))
db_ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')


class TestNextflowQc(unittest.TestCase):
    def test_nextflow_qc_using_database(self):
        '''test nextflow_qc using database'''
        tmp_data_dir = 'tmp.nextflow_qc'
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
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'qc.nf')
        work_dir = 'tmp.nextflow_qc.work'
        dag_file = 'nextflow.qc.dag.db.pdf'
        try:
            os.unlink(dag_file)
        except:
            pass

        command = ' '.join([
            'nextflow run',
            '--dataset_name g1',  # one of the samples is in group2 and should get ignored
            '--ref_id 1',
            '--references_root', os.path.abspath(references_root),
            '--pipeline_root', pipeline_root,
            '--db_config_file', db_ini_file,
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
        got_pipeline_rows = database.get_rows_from_table('Pipeline')
        got_pipeline_rows.sort(key=itemgetter('seqrep_id'))
        expected_pipeline_rows = [
            {'isolate_id': 1, 'seqrep_id': 1, 'seqrep_pool': None, 'version': '0.0.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 1, 'seqrep_id': 1, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'qc', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 2, 'seqrep_pool': None, 'version': '0.0.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 2, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'qc', 'status': 1, 'reference_id': 1},
            {'isolate_id': 3, 'seqrep_id': 3, 'seqrep_pool': None, 'version': '0.0.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 3, 'seqrep_id': 3, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'qc', 'status': -1, 'reference_id': 1},
            {'isolate_id': 4, 'seqrep_id': 4, 'seqrep_pool': None, 'version': '0.0.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
        ]
        self.assertEqual(expected_pipeline_rows, got_pipeline_rows)

        # check QC stats added to database
        got_qc_rows = database.get_rows_from_table('QC')
        got_qc_rows.sort(key=itemgetter('seqrep_id'))
        expected_qc_rows = [{
            'seqrep_id': 1,
            'pipeline_version': clockwork_version,
            'fastqc1_adapter_content': 'pass',
            'fastqc1_basic_statistics': 'pass',
            'fastqc1_gc': 48.0,
            'fastqc1_kmer_content': 'fail',
            'fastqc1_max_sequence_length': 75,
            'fastqc1_min_sequence_length': 75,
            'fastqc1_overrepresented_sequences': 'fail',
            'fastqc1_per_base_n_content': 'pass',
            'fastqc1_per_base_sequence_content': 'fail',
            'fastqc1_per_base_sequence_quality': 'pass',
            'fastqc1_per_sequence_gc_content': 'fail',
            'fastqc1_per_sequence_quality_scores': 'fail',
            'fastqc1_sequence_duplication_levels': 'pass',
            'fastqc1_sequence_length_distribution': 'pass',
            'fastqc1_sequences_flagged_as_poor_quality': 0,
            'fastqc1_total_sequences': 72,
            'fastqc2_adapter_content': 'pass',
            'fastqc2_basic_statistics': 'pass',
            'fastqc2_gc': 48.0,
            'fastqc2_kmer_content': 'fail',
            'fastqc2_max_sequence_length': 75,
            'fastqc2_min_sequence_length': 75,
            'fastqc2_overrepresented_sequences': 'fail',
            'fastqc2_per_base_n_content': 'pass',
            'fastqc2_per_base_sequence_content': 'fail',
            'fastqc2_per_base_sequence_quality': 'pass',
            'fastqc2_per_sequence_gc_content': 'fail',
            'fastqc2_per_sequence_quality_scores': 'fail',
            'fastqc2_sequence_duplication_levels': 'pass',
            'fastqc2_sequence_length_distribution': 'pass',
            'fastqc2_sequences_flagged_as_poor_quality': 0,
            'fastqc2_total_sequences': 72,
            'samtools_average_quality': 40.0,
            'samtools_bases_mapped_cigar': 9900,
            'samtools_bases_trimmed': 0,
            'samtools_error_rate': 0.0,
            'samtools_insert_size_average': 199.6,
            'samtools_insert_size_standard_deviation': 1.0,
            'samtools_inward_oriented_pairs': 66,
            'samtools_outward_oriented_pairs': 0,
            'samtools_pairs_with_other_orientation': 0,
            'samtools_raw_total_sequences': 144,
            'samtools_reads_duplicated': 4,
            'samtools_reads_mapped': 132,
            'het_snp_het_calls': 0,
            'het_snp_positions': 983,
            'het_snp_total_snps': 0,
        },
        {
            'seqrep_id': 2,
            'pipeline_version': clockwork_version,
            'fastqc1_adapter_content': 'pass',
            'fastqc1_basic_statistics': 'pass',
            'fastqc1_gc': 48.0,
            'fastqc1_kmer_content': 'fail',
            'fastqc1_max_sequence_length': 75,
            'fastqc1_min_sequence_length': 75,
            'fastqc1_overrepresented_sequences': 'fail',
            'fastqc1_per_base_n_content': 'pass',
            'fastqc1_per_base_sequence_content': 'fail',
            'fastqc1_per_base_sequence_quality': 'pass',
            'fastqc1_per_sequence_gc_content': 'fail',
            'fastqc1_per_sequence_quality_scores': 'fail',
            'fastqc1_sequence_duplication_levels': 'pass',
            'fastqc1_sequence_length_distribution': 'pass',
            'fastqc1_sequences_flagged_as_poor_quality': 0,
            'fastqc1_total_sequences': 72,
            'fastqc2_adapter_content': 'pass',
            'fastqc2_basic_statistics': 'pass',
            'fastqc2_gc': 49.0,
            'fastqc2_kmer_content': 'fail',
            'fastqc2_max_sequence_length': 75,
            'fastqc2_min_sequence_length': 75,
            'fastqc2_overrepresented_sequences': 'fail',
            'fastqc2_per_base_n_content': 'pass',
            'fastqc2_per_base_sequence_content': 'fail',
            'fastqc2_per_base_sequence_quality': 'pass',
            'fastqc2_per_sequence_gc_content': 'warn',
            'fastqc2_per_sequence_quality_scores': 'fail',
            'fastqc2_sequence_duplication_levels': 'pass',
            'fastqc2_sequence_length_distribution': 'pass',
            'fastqc2_sequences_flagged_as_poor_quality': 0,
            'fastqc2_total_sequences': 72,
            'samtools_average_quality': 40.0,
            'samtools_bases_mapped_cigar': 9900,
            'samtools_bases_trimmed': 0,
            'samtools_error_rate': 0.0,
            'samtools_insert_size_average': 199.7,
            'samtools_insert_size_standard_deviation': 1.1,
            'samtools_inward_oriented_pairs': 66,
            'samtools_outward_oriented_pairs': 0,
            'samtools_pairs_with_other_orientation': 0,
            'samtools_raw_total_sequences': 144,
            'samtools_reads_duplicated': 0,
            'samtools_reads_mapped': 132,
            'het_snp_het_calls': 0,
            'het_snp_positions': 983,
            'het_snp_total_snps': 0,
        }]
        self.assertEqual(expected_qc_rows, got_qc_rows)

        # check QC files got written. No need to check contents, as that is done
        # elsewhere. We're just checking nextflow runs OK here.
        ids = [
            {'sample': 1, 'isolate_id': 1, 'seq_repl': 43},
            {'sample': 2, 'isolate_id': 2, 'seq_repl': 45},
        ]
        for id_dict in ids:
            iso_dir = isolate_dir.IsolateDir(pipeline_root, id_dict['sample'], id_dict['isolate_id'])
            qc_root_dir = iso_dir.pipeline_dir(id_dict['seq_repl'], 'qc', clockwork_version)
            self.assertTrue(os.path.exists(qc_root_dir))
            for method in ['fastqc', 'samtools_qc']:
                this_qc_dir = os.path.join(qc_root_dir, method)
                self.assertTrue(os.path.exists(this_qc_dir))
                self.assertTrue(len(os.listdir(this_qc_dir)) >= 1)

        shutil.rmtree(tmp_data_dir)
        nextflow_helper.clean_files()


    def test_nextflow_qc_using_fastq_input(self):
        '''test nextflow_qc using fastq input'''
        reads1 = os.path.join(data_dir, 'Reads', 'reads.1.1.fq.gz')
        reads2 = os.path.join(data_dir, 'Reads', 'reads.1.2.fq.gz')
        output_dir = 'tmp.test_nextflow_qc_using_fastq_input'
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'qc.nf')
        nextflow_helper.write_config_file()
        work_dir = 'tmp.nextflow_qc.work'
        dag_file = 'nextflow.qc.dag.no_db.pdf'
        try:
            os.unlink(dag_file)
        except:
            pass

        command = ' '.join([
            'nextflow run',
            '--reads_in1', reads1,
            '--reads_in2', reads2,
            '--output_dir', output_dir,
            '--ref_fasta', os.path.join(data_dir, 'Reference', 'ref.fa'),
            '-with-dag', dag_file,
            '-c', nextflow_helper.config_file,
            '-w', work_dir,
            nextflow_file
        ])
        utils.syscall(command)
        os.unlink(nextflow_helper.config_file)
        shutil.rmtree(work_dir)


        self.assertTrue(os.path.exists(output_dir))
        for method in ['fastqc', 'samtools_qc']:
            qc_dir = os.path.join(output_dir, method)
            self.assertTrue(os.path.exists(qc_dir))
            self.assertTrue(len(os.listdir(qc_dir)) >= 1)

        shutil.rmtree(output_dir)
        nextflow_helper.clean_files()
