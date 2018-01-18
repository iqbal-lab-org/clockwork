import unittest
import shutil
import os
from clockwork import samtools_qc

modules_dir = os.path.dirname(os.path.abspath(samtools_qc.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'samtools_qc')

# Made ref files like this:
# fastaq make_random_contigs --seed 42 1 2000 samtools_qc.ref.fa
# bwa index samtools_qc.ref.fa
# art_illumina -ss HS20 -i samtools_qc.ref.fa -o samtools_qc.reads. -l 75 -m 200 -s 5 -f 2 -na
# bwa mem samtools_qc.ref.fa samtools_qc.reads.1.fq samtools_qc.reads.2.fq > samtools_qc.sam

class TestSamtoolsQc(unittest.TestCase):
    def test_map_reads(self):
        '''test _map_reads'''
        # We'll check no errors are thrown and a SAM file gets made.
        # Don't look at contents of SAM. Trust the mapper
        reads1 = os.path.join(data_dir, 'reads.1.fq')
        reads2 = os.path.join(data_dir, 'reads.2.fq')
        ref_fasta = os.path.join(data_dir, 'ref.fa')
        tmp_sam = 'tmp.samtools.qc_map_reads.sam'
        if os.path.exists(tmp_sam):
            os.unlink(tmp_sam)
        samtools_qc.SamtoolsQc._map_reads(ref_fasta, reads1, reads2, tmp_sam)
        self.assertTrue(os.path.exists(tmp_sam))
        os.unlink(tmp_sam)


    def test_make_stats_and_plots(self):
        '''test _make_stats_and_plots'''
        ref_fasta = os.path.join(data_dir, 'ref.fa')
        sam_file = os.path.join(data_dir, 'sam')
        tmp_dir = 'tmp.samtools_qc.make_stats_and_plots'
        os.mkdir(tmp_dir)
        outprefix = os.path.join(tmp_dir, 'test')
        samtools_qc.SamtoolsQc._make_stats_and_plots(sam_file, ref_fasta, outprefix)

        expected_files = [
            'test.plot-acgt-cycles.png',
            'test.plot-coverage.png',
            'test.plot-gc-content.png',
            'test.plot-insert-size.png',
            'test.plot-mism-per-cycle.png',
            'test.plot-quals2.png',
            'test.plot-quals3.png',
            'test.plot-quals-hm.png',
            'test.plot-quals.png',
            'test.stats',
        ]
        expected_files.sort()

        got_files = sorted(list(os.listdir(tmp_dir)))
        self.assertEqual(expected_files, got_files)
        shutil.rmtree(tmp_dir)


    def test_stats_from_report(self):
        '''test stats_from_report'''
        stats_file = os.path.join(data_dir, 'stats_from_report.txt')
        got = samtools_qc.SamtoolsQc.stats_from_report(stats_file)

        expected = {
            'raw_total_sequences': 54,
            'reads_mapped': 54,
            'reads_duplicated': 2,
            'bases_mapped_cigar': 4050,
            'bases_trimmed': 0,
            'error_rate': 7.407407e-03,
            'average_quality': 35.3,
            'insert_size_average': 198.6,
            'insert_size_standard_deviation': 4.1,
            'inward_oriented_pairs': 27,
            'outward_oriented_pairs': 0,
            'pairs_with_other_orientation': 0,
        }

        self.assertEqual(expected, got)


    def test_het_snp_stats_from_summary_file(self):
        '''test het_snp_stats_from_summary_file'''
        infile = os.path.join(data_dir, 'het_snps_report.tsv')
        got = samtools_qc.SamtoolsQc.het_snp_stats_from_summary_file(infile)
        expected = {
            'Total_length': 2000,
            'Positions_used': 1673,
            'Total_SNPs': 200,
            'Het_SNPs': 84,
            'Percent_SNPs_are_het': 42.0,
        }
        self.assertEqual(expected, got)


    def test_run(self):
        '''test run'''
        reads1 = os.path.join(data_dir, 'reads.1.fq')
        reads2 = os.path.join(data_dir, 'reads.2.fq')
        ref_fasta = os.path.join(data_dir, 'ref.fa')
        tmp_dir = 'tmp.test.samtools_qc.run'
        samqc = samtools_qc.SamtoolsQc(ref_fasta, reads1, reads2, tmp_dir)
        samqc.run()

        expected_files = [
            'het_snps.het_calls.vcf',
            'het_snps.per_contig.tsv',
            'het_snps.summary.tsv',
            'samtools_qc.plot-acgt-cycles.png',
            'samtools_qc.plot-coverage.png',
            'samtools_qc.plot-gc-content.png',
            'samtools_qc.plot-insert-size.png',
            'samtools_qc.plot-mism-per-cycle.png',
            'samtools_qc.plot-quals2.png',
            'samtools_qc.plot-quals3.png',
            'samtools_qc.plot-quals-hm.png',
            'samtools_qc.plot-quals.png',
            'samtools_qc.stats',
        ]
        expected_files.sort()

        got_files = sorted(list(os.listdir(tmp_dir)))
        self.assertEqual(expected_files, got_files)
        expected_stats = samtools_qc.SamtoolsQc.stats_from_report(os.path.join(tmp_dir, 'samtools_qc.stats'))
        self.assertEqual(expected_stats, samqc.stats)
        shutil.rmtree(tmp_dir)
