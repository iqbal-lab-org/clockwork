import unittest
import shutil
import os
from clockwork import fastqc

modules_dir = os.path.dirname(os.path.abspath(fastqc.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'fastqc')


class TestFastqc(unittest.TestCase):
    def test_run_fastqc(self):
        '''test _run_fastqc'''
        reads1 = os.path.join(data_dir, 'reads_1.fq')
        reads2 = os.path.join(data_dir, 'reads_2.fq')
        tmpdir = 'tmp.test.fastqc_run_fastqc'

        fastqc.Fastqc._run_fastqc(tmpdir, [reads1, reads2])

        # We'll trust the contents of fastqc files. Just check they got made.
        for i in '1', '2':
            outdir = os.path.join(tmpdir, 'reads_' + i + '_fastqc')
            self.assertTrue(os.path.exists(outdir))
            self.assertTrue(os.path.isdir(outdir))
            self.assertTrue(os.path.exists(outdir + '.html'))
            self.assertTrue(os.path.exists(outdir + '.zip'))

        shutil.rmtree(tmpdir)



    def test_clean_outdir(self):
        '''test _clean_outdir'''
        tmpdir = 'tmp.test_clean_outdir'
        shutil.copytree(os.path.join(data_dir, 'clean_outdir'), tmpdir)
        fastqc.Fastqc._clean_outdir(tmpdir)

        expected_dirs = {'fastqc_reads_1_fastqc', 'fastqc_reads_2_fastqc'}
        got_dirs = set(os.listdir(tmpdir))
        self.assertEqual(expected_dirs, got_dirs)

        expected_files = {
            'adapter_content.png',
            'duplication_levels.png',
            'fastqc_data.txt',
            'per_base_n_content.png',
            'per_base_quality.png',
            'per_base_sequence_content.png',
            'per_sequence_gc_content.png',
            'per_sequence_quality.png',
            'sequence_length_distribution.png',
        }

        for d in expected_dirs:
            got_files = set(os.listdir(os.path.join(tmpdir, d)))
            self.assertEqual(expected_files, got_files)

        shutil.rmtree(tmpdir)


    def test_seq_length_from_report_string(self):
        '''test _seq_length_from_report_string'''
        self.assertEqual((10, 42), fastqc.Fastqc._seq_length_from_report_string('10-42'))
        self.assertEqual((42, 42), fastqc.Fastqc._seq_length_from_report_string('42'))


    def test_stats_from_report(self):
        '''test _stats_from_report'''
        expected = {
            'filename': 'fastqc_reads_1.fq',
            'total_sequences': 6,
            'min_sequence_length': 42,
            'max_sequence_length': 75,
            'sequences_flagged_as_poor_quality': 0,
            'gc': 49,
            'basic_statistics': 'pass',
            'per_base_sequence_quality': 'pass',
            'per_sequence_quality_scores': 'fail',
            'per_base_sequence_content': 'fail',
            'per_sequence_gc_content': 'fail',
            'per_base_n_content': 'pass',
            'sequence_length_distribution': 'pass',
            'sequence_duplication_levels': 'pass',
            'overrepresented_sequences': 'fail',
            'adapter_content': 'pass',
            'kmer_content': 'pass'
        }
        stats_file = os.path.join(data_dir, 'stats_from_report.txt')
        got = fastqc.Fastqc._stats_from_report(stats_file)
        self.maxDiff = None
        self.assertEqual(expected, got)


    def test_gather_all_stats(self):
        '''test gather_all_stats'''
        indir = os.path.join(data_dir, 'gather_all_stats')
        got = fastqc.Fastqc.gather_all_stats(indir)
        stats1 = fastqc.Fastqc._stats_from_report(os.path.join(indir, 'fastqc_reads_1_fastqc', 'fastqc_data.txt'))
        stats1.pop('filename')
        stats2 = fastqc.Fastqc._stats_from_report(os.path.join(indir, 'fastqc_reads_2_fastqc', 'fastqc_data.txt'))
        stats2.pop('filename')
        expected = {
            'fastqc_reads_1.fq': stats1,
            'fastqc_reads_2.fq': stats2,
        }
        self.assertEqual(expected, got)


    def test_run(self):
        '''test run'''
        reads1 = os.path.join(data_dir, 'reads_1.fq')
        reads2 = os.path.join(data_dir, 'reads_2.fq')
        tmpdir = 'tmp.test.fastqc.run'
        fqc = fastqc.Fastqc(tmpdir, [reads1, reads2])
        fqc.run()
        self.assertTrue(os.path.exists(tmpdir))
        dir1 = os.path.join(tmpdir, 'reads_1_fastqc')
        self.assertTrue(os.path.exists(dir1))
        self.assertTrue(os.path.isdir(dir1))
        dir2 = os.path.join(tmpdir, 'reads_2_fastqc')
        self.assertTrue(os.path.exists(dir2))
        self.assertTrue(os.path.isdir(dir2))
        stats1 = fastqc.Fastqc._stats_from_report(os.path.join(dir1, 'fastqc_data.txt'))
        stats2 = fastqc.Fastqc._stats_from_report(os.path.join(dir2, 'fastqc_data.txt'))
        stats1.pop('filename')
        stats2.pop('filename')
        expected_stats = {
            'reads_1.fq': stats1,
            'reads_2.fq': stats2,
        }
        self.maxDiff = None
        self.assertEqual(expected_stats, fqc.stats)
        shutil.rmtree(tmpdir)
