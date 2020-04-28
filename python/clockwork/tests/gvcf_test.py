import filecmp
import os
import subprocess
import unittest
from clockwork import gvcf

from cluster_vcf_records import vcf_record

modules_dir = os.path.dirname(os.path.abspath(gvcf.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "gvcf")


def lines_from_vcf_ignore_file_date(vcf):
    with open(vcf) as f:
        return [x for x in f if not x.startswith("##fileDate=")]


class TestGvcf(unittest.TestCase):
    def test_move_info_fields_to_format(self):
        record = vcf_record.VcfRecord(
            "ref\t1\t.\tC\tG\t.\t.\tfoo=bar;spam=eggs\tcleese\tchapman"
        )
        gvcf._move_info_fields_to_format(record)
        assert record.INFO == {}
        assert record.FORMAT == {"foo": "bar", "spam": "eggs", "cleese": "chapman"}

    def test_gvcf_from_minos_vcf_and_samtools_gvcf(self):
        ref_fasta = os.path.join(
            data_dir, "gvcf_from_minos_vcf_and_samtools_gvcf.ref.fa"
        )
        minos_vcf = os.path.join(
            data_dir, "gvcf_from_minos_vcf_and_samtools_gvcf.minos.vcf"
        )
        samtools_vcf = os.path.join(
            data_dir, "gvcf_from_minos_vcf_and_samtools_gvcf.samtools.vcf"
        )
        tmp_out = "tmp.gvcf_from_minos_vcf_and_samtools_gvcf.out.vcf"
        subprocess.check_output(f"rm -f {tmp_out}", shell=True)
        gvcf.gvcf_from_minos_vcf_and_samtools_gvcf(
            ref_fasta, minos_vcf, samtools_vcf, tmp_out
        )
        expect_lines = lines_from_vcf_ignore_file_date(
            os.path.join(data_dir, "gvcf_from_minos_vcf_and_samtools_gvcf.out.vcf")
        )
        got_lines = lines_from_vcf_ignore_file_date(tmp_out)
        self.assertEqual(expect_lines, got_lines)
        os.unlink(tmp_out)

    def test_samtools_vcf_record_to_frs(self):
        record = vcf_record.VcfRecord(
            "ref\t1\t.\tC\tG\t.\t.\tCALLER=samtools\tDP4\t1,2,14,13"
        )
        self.assertEqual(gvcf._samtools_vcf_record_to_frs(record, 0), 0.1)
        self.assertEqual(gvcf._samtools_vcf_record_to_frs(record, 1), 0.9)

    def test_vcf_record_pass_index(self):
        record = vcf_record.VcfRecord(
            "ref\t1\t.\tC\tG\t.\t.\tCALLER=samtools\tGT:DP:DP4\t1/1:20:1,2,14,13"
        )
        self.assertEqual(1, gvcf._vcf_record_pass_index(record, min_frs=0.9, min_dp=5))
        self.assertEqual(
            None, gvcf._vcf_record_pass_index(record, min_frs=0.9, min_dp=21)
        )
        self.assertEqual(
            None, gvcf._vcf_record_pass_index(record, min_frs=0.99, min_dp=5)
        )

        record = vcf_record.VcfRecord(
            "ref\t1\t.\tC\tG\t.\tPASS\tCALLER=minos\tGT:DP:FRS\t1/1:20:0.95"
        )
        self.assertEqual(1, gvcf._vcf_record_pass_index(record))
        self.assertEqual(
            1, gvcf._vcf_record_pass_index(record, require_minos_pass=False)
        )
        self.assertEqual(None, gvcf._vcf_record_pass_index(record, min_frs=0.96))
        self.assertEqual(None, gvcf._vcf_record_pass_index(record, min_dp=21))
        record = vcf_record.VcfRecord(
            "ref\t1\t.\tC\tG\t.\tFAIL\tCALLER=minos\tGT:DP:FRS\t1/1:20:0.95"
        )
        self.assertEqual(None, gvcf._vcf_record_pass_index(record))
        self.assertEqual(
            1, gvcf._vcf_record_pass_index(record, require_minos_pass=False)
        )
        self.assertEqual(
            None,
            gvcf._vcf_record_pass_index(record, require_minos_pass=False, min_frs=0.96),
        )
        self.assertEqual(
            None,
            gvcf._vcf_record_pass_index(record, require_minos_pass=False, min_dp=21),
        )

    def test_gvcf_to_fasta(self):
        vcf = os.path.join(data_dir, "gvcf_to_fasta.vcf")
        tmp_out = "tmp.gvcf_to_fasta.fa"
        subprocess.check_output(f"rm -f {tmp_out}", shell=True)
        gvcf.gvcf_to_fasta(vcf, tmp_out)
        expect_fasta = os.path.join(data_dir, "gvcf_to_fasta.fa")
        self.assertTrue(filecmp.cmp(tmp_out, expect_fasta, shallow=False))
        os.unlink(tmp_out)
