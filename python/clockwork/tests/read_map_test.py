import unittest
import shutil
import filecmp
import os
import pysam
from clockwork import read_map, utils

modules_dir = os.path.dirname(os.path.abspath(read_map.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "read_map")

# note: reads and ref are copies of samtoosl_qc files.
# Then manually copied the last read of each file, to make
# it a duplicate


class TestReadMap(unittest.TestCase):
    def test_map_reads(self):
        """test map_reads"""
        reads1 = os.path.join(data_dir, "reads.1.fq")
        reads2 = os.path.join(data_dir, "reads.2.fq")
        ref_fasta = os.path.join(data_dir, "ref.fa")
        tmp_sam = "tmp.test_map_reads.sam"
        if os.path.exists(tmp_sam):
            os.unlink(tmp_sam)
        read_map.map_reads(
            ref_fasta, reads1, reads2, tmp_sam, read_group=("1", "GROUP_NAME")
        )
        self.assertTrue(os.path.exists(tmp_sam))
        tmp_stats = tmp_sam + ".stats"
        expected_stats = os.path.join(data_dir, "flagstat")
        utils.syscall("samtools flagstat " + tmp_sam + " > " + tmp_stats)
        self.assertTrue(filecmp.cmp(expected_stats, tmp_stats, shallow=False))
        found_rg_line = False
        with open(tmp_sam) as f:
            for line in f:
                if line == "@RG\tLB:LIB\tID:1\tSM:GROUP_NAME\n":
                    found_rg_line = True
                    break

        self.assertTrue(found_rg_line)
        os.unlink(tmp_sam)
        os.unlink(tmp_stats)

    def test_map_reads_secondary_hits_removed(self):
        """test map_reads secondary hits get removed"""
        reads1 = os.path.join(data_dir, "secondary_hits_removed.reads_1.fq")
        reads2 = os.path.join(data_dir, "secondary_hits_removed.reads_2.fq")
        ref_fasta = os.path.join(data_dir, "secondary_hits_removed.ref.fa")
        tmp_sam = "tmp.test_map_reads.sam"
        if os.path.exists(tmp_sam):
            os.unlink(tmp_sam)
        read_map.map_reads(ref_fasta, reads1, reads2, tmp_sam)
        # bwa mem reports one secondary alignment, so the 1 read pair makes
        # 3 SAM records. So should have 2 records after removing secondary match.
        self.assertEqual(2, utils.sam_record_count(tmp_sam))
        os.unlink(tmp_sam)

    def test_map_reads_rmdup(self):
        """test map_reads rmdup"""
        reads1 = os.path.join(data_dir, "reads.1.fq")
        reads2 = os.path.join(data_dir, "reads.2.fq")
        ref_fasta = os.path.join(data_dir, "ref.fa")
        tmp_sam = "tmp.test_map_reads.sam"
        if os.path.exists(tmp_sam):
            os.unlink(tmp_sam)
        read_map.map_reads(ref_fasta, reads1, reads2, tmp_sam, rmdup=True)
        self.assertTrue(os.path.exists(tmp_sam))
        tmp_stats = tmp_sam + ".stats"
        expected_stats = os.path.join(data_dir, "rmdup.flagstat")
        utils.syscall("samtools flagstat " + tmp_sam + " > " + tmp_stats)
        self.assertTrue(filecmp.cmp(expected_stats, tmp_stats, shallow=False))
        os.unlink(tmp_sam)
        os.unlink(tmp_stats)

    def test_map_reads_markdup(self):
        """test map_reads markdup"""
        reads1 = os.path.join(data_dir, "reads.1.fq")
        reads2 = os.path.join(data_dir, "reads.2.fq")
        ref_fasta = os.path.join(data_dir, "ref.fa")
        tmp_sam = "tmp.test_map_reads.sam"
        if os.path.exists(tmp_sam):
            os.unlink(tmp_sam)
        read_map.map_reads(ref_fasta, reads1, reads2, tmp_sam, markdup=True)
        self.assertTrue(os.path.exists(tmp_sam))
        tmp_stats = tmp_sam + ".stats"
        expected_stats = os.path.join(data_dir, "markdup.flagstat")
        utils.syscall("samtools flagstat " + tmp_sam + " > " + tmp_stats)
        self.assertTrue(filecmp.cmp(expected_stats, tmp_stats, shallow=False))
        os.unlink(tmp_sam)
        os.unlink(tmp_stats)

    def test_map_reads_markdup_and_rmdup(self):
        """test map_reads rmdup and markdup"""
        with self.assertRaises(read_map.Error):
            read_map.map_reads(
                "ref_fasta", "reads1", "reads2", "sam", rmdup=True, markdup=True
            )

    def test_map_reads_set(self):
        """test map_reads_set"""
        reads1 = [
            os.path.join(data_dir, "reads_set_" + str(i) + ".1.fq")
            for i in range(1, 4, 1)
        ]
        reads2 = [
            os.path.join(data_dir, "reads_set_" + str(i) + ".2.fq")
            for i in range(1, 4, 1)
        ]
        ref_fasta = os.path.join(data_dir, "ref.fa")
        try:
            os.unlink(tmp_bam)
        except:
            pass
        tmp_bam = "tmp.read_map.map_reads_set.bam"

        def reads_in_bam(bam_file):
            return len(list(pysam.Samfile(bam_file, "rb").fetch(until_eof=True)))

        # check works when list length is 1. Should have 12 read pairs in the bam
        read_map.map_reads_set(ref_fasta, [reads1[0]], [reads2[0]], tmp_bam, rmdup=True)
        self.assertTrue(os.path.exists(tmp_bam))
        self.assertEqual(24, reads_in_bam(tmp_bam))
        os.unlink(tmp_bam)

        # check works when list length > 1, Should have 35 read pairs in the bam
        read_map.map_reads_set(
            ref_fasta,
            reads1,
            reads2,
            tmp_bam,
            rmdup=True,
            read_group=("1", "GROUP_NAME"),
        )
        self.assertEqual(70, reads_in_bam(tmp_bam))
        os.unlink(tmp_bam)
