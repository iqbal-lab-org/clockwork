import unittest
import filecmp
import os
from clockwork import simple_vcf_merger

modules_dir = os.path.dirname(os.path.abspath(simple_vcf_merger.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'simple_vcf_merger')


# NOTE: you could argable not write tests for _merge_lists_of_vcf_records and
# _merge_vcf_records and just do one for run()

class TestSimpleVcfMerger(unittest.TestCase):
    # No longer run this test. cluster_vcf_records was updated to check that
    # the REF strings in the input VCF files match the reference genome. That
    # is not true for the test data for this sample. And the simple merger is
    # not being used.
    def _test_run(self):
        '''test run'''
        samtools_vcf = os.path.join(data_dir, 'samtools.vcf')
        cortex_vcf = os.path.join(data_dir, 'cortex.vcf')
        ref_fasta = os.path.join(data_dir, 'NC_000962.3.fa')
        expected_vcf = os.path.join(data_dir, 'expect.vcf')
        tmp_vcf = 'tmp.simple_vcf_merger.out.vcf'
        # ... run the merging ...
        merger = simple_vcf_merger.SimpleVcfMerger(samtools_vcf,cortex_vcf, tmp_vcf, ref_fasta, homozygous_only=True, min_SNP_qual=30.0, min_dp4=5.0, min_GT_conf=5.0)
        merger.run()
        self.assertTrue(filecmp.cmp(expected_vcf, tmp_vcf, shallow=False))
        os.unlink(tmp_vcf)
