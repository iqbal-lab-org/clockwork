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
data_dir = os.path.join(nextflow_helper.data_root_dir, 'nextflow_variant_call')
modules_dir = os.path.dirname(os.path.abspath(db.__file__))
db_ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')


def vcf_to_sample(vcf_file):
    with open(vcf_file) as f:
        for line in f:
            if line.startswith('#CHROM\tPOS'):
                return line.rstrip().split()[-1]

    return None


class TestNextflowVarcall(unittest.TestCase):
    def _files_are_present_and_correct(self, pipeline_dir, expected_sample, expect_rmdup_bam=True, expect_ref_check_files=False):
        samtools_dir = os.path.join(pipeline_dir, 'samtools')
        samtools_vcf = os.path.join(samtools_dir, 'samtools.vcf')
        self.assertEqual(expect_rmdup_bam, os.path.exists(os.path.join(samtools_dir, 'rmdup.bam')))
        self.assertTrue(os.path.exists(samtools_vcf))
        self.assertEqual(expected_sample, vcf_to_sample(samtools_vcf))
        cortex_dir = os.path.join(pipeline_dir, 'cortex')
        self.assertTrue(os.path.exists(os.path.join(cortex_dir, 'cortex.in.fofn')))
        self.assertTrue(os.path.exists(os.path.join(cortex_dir, 'cortex.in.index')))
        self.assertTrue(os.path.exists(os.path.join(cortex_dir, 'cortex.in.index_ref.fofn')))
        self.assertTrue(os.path.exists(os.path.join(cortex_dir, 'cortex.log')))
        self.assertTrue(os.path.exists(os.path.join(cortex_dir, 'cortex.out')))
        cortex_vcf_files = [os.path.join(cortex_dir, x) for x in glob.glob(os.path.join(cortex_dir, 'cortex.out', 'vcfs', '**')) if  x.endswith('.vcf')]
        if expect_ref_check_files:
            self.assertEqual(5, len(cortex_vcf_files))
        else:
            self.assertEqual(2, len(cortex_vcf_files))

        for vcf_file in cortex_vcf_files:
            self.assertEqual(expected_sample, vcf_to_sample(vcf_file))
        minos_dir = os.path.join(pipeline_dir, 'minos')
        self.assertTrue(os.path.exists(minos_dir))
        self.assertTrue(os.path.exists(os.path.join(minos_dir, 'final.vcf')))
        simple_merge_dir = os.path.join(pipeline_dir, 'simple_merge')
        self.assertTrue(os.path.exists(simple_merge_dir))
        self.assertTrue(os.path.exists(os.path.join(simple_merge_dir, 'simple_merge.vcf')))

        if expect_ref_check_files:
            self.assertTrue(os.path.exists(samtools_vcf + '.check.stats.tsv'))
            self.assertTrue(os.path.exists(os.path.join(minos_dir, 'final.vcf.check.stats.tsv')))
            self.assertTrue(os.path.exists(os.path.join(simple_merge_dir, 'simple_merge.vcf.check.stats.tsv')))


    def test_nextflow_variant_call_using_database(self):
        '''test nextflow_variant_call using database'''
        tmp_data_dir = 'tmp.nextflow_variant_call_db_input.data'
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
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'variant_call.nf')
        work_dir = 'tmp.nextflow_variant_call_db_input.work'
        dag_file = 'nextflow.variant_call.dag.db.pdf'
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
            '--cortex_mem_height 17',
            '--testing',
            '--truth_ref', os.path.join(tmp_data_dir, 'truth_ref.fa'),
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
        got_rows.sort(key=itemgetter('isolate_id', 'pipeline_name'))
        expected_rows = [
            {'isolate_id': 1, 'seqrep_id': 1, 'seqrep_pool': None, 'version': '0.3.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 1, 'seqrep_id': 2, 'seqrep_pool': None, 'version': '0.3.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 1, 'seqrep_id': None, 'seqrep_pool': '1_2', 'version': clockwork_version, 'pipeline_name': 'variant_call', 'status': 1, 'reference_id': 2},
            {'isolate_id': 2, 'seqrep_id': 3, 'seqrep_pool': None, 'version': '0.3.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 4, 'seqrep_pool': None, 'version': '0.3.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 2, 'seqrep_id': 3, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'variant_call', 'status': 1, 'reference_id': 2},
            {'isolate_id': 2, 'seqrep_id': 4, 'seqrep_pool': None, 'version': clockwork_version, 'pipeline_name': 'variant_call', 'status': 1, 'reference_id': 2},
            {'isolate_id': 3, 'seqrep_id': 5, 'seqrep_pool': None, 'version': '0.3.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
            {'isolate_id': 3, 'seqrep_id': None, 'seqrep_pool': '1', 'version': clockwork_version, 'pipeline_name': 'variant_call', 'status': -1, 'reference_id': 2},
            {'isolate_id': 4, 'seqrep_id': 6, 'seqrep_pool': None, 'version': '0.3.1', 'pipeline_name': 'remove_contam', 'status': 1, 'reference_id': 1},
        ]
        self.assertEqual(expected_rows, got_rows)

        # check VCF files etc got written. No need to check contents, trust the tools
        # We're just checking nextflow runs OK here.
        ids = [
            {'sample': 1, 'seqrep_id': '1_2', 'isolate_id': 1, 'seq_repl': '1_2', 'sample_name': 'site.s1.iso.42.subject.p1.lab_id.l1.seq_reps.1_2'},
            {'sample': 2, 'seqrep_id': 3, 'isolate_id': 2, 'seq_repl': '1', 'sample_name': 'site.s2.iso.43.subject.p2.lab_id.l2.seq_reps.1'},
            {'sample': 2, 'seqrep_id': 4, 'isolate_id': 2, 'seq_repl': '2', 'sample_name': 'site.s2.iso.43.subject.p2.lab_id.l2.seq_reps.2'},
        ]
        for id_dict in ids:
            iso_dir = isolate_dir.IsolateDir(pipeline_root, id_dict['sample'], id_dict['isolate_id'])
            pipeline_dir = iso_dir.pipeline_dir(id_dict['seq_repl'], 'variant_call', clockwork_version, reference_id=2)
            self._files_are_present_and_correct(pipeline_dir, id_dict['sample_name'], expect_ref_check_files=True)

        shutil.rmtree(tmp_data_dir)
        nextflow_helper.clean_files()


    def test_nextflow_variant_call_using_fastq_input(self):
        '''test nextflow_variant_call using fastq input'''
        reads1 = os.path.join(data_dir, 'Reads', 'reads.1.1.fq.gz')
        reads2 = os.path.join(data_dir, 'Reads', 'reads.1.2.fq.gz')
        outdir = os.path.abspath('tmp.test_nextflow_variant_call_fastq_input.out')
        tmp_data_dir = 'tmp.nextflow_variant_call_fastq_input.data'
        if os.path.exists(tmp_data_dir):
            shutil.rmtree(tmp_data_dir)
        shutil.copytree(data_dir, tmp_data_dir)
        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, 'variant_call.nf')
        nextflow_helper.write_config_file()
        work_dir = 'tmp.nextflow_variant_call_fastq_input.work'
        sample_name = 'test_sample_name'
        dag_file = 'nextflow.variant_call.dag.no_db.pdf'
        try:
            os.unlink(dag_file)
        except:
            pass

        command = ' '.join([
            'nextflow run',
            '--reads_in1', reads1,
            '--reads_in2', reads2,
            '--output_dir', outdir,
            '--ref_dir', os.path.join(tmp_data_dir, 'Reference'),
            '--sample_name', sample_name,
            '--cortex_mem_height 17',
            '--testing',
            '-with-dag', dag_file,
            '-c', nextflow_helper.config_file,
            '-w', work_dir,
            nextflow_file
        ])
        utils.syscall(command)
        self._files_are_present_and_correct(outdir, sample_name, expect_rmdup_bam=True, expect_ref_check_files=False)
        os.unlink(nextflow_helper.config_file)
        shutil.rmtree(work_dir)
        shutil.rmtree(tmp_data_dir)
        shutil.rmtree(outdir)
        nextflow_helper.clean_files()

