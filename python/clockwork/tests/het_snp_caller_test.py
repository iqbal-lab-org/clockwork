# This is pretty much a rewrite of:
# https://github.com/sanger-pathogens/vr-codebase/blob/master/t/Pathogens/QC/HetSNPCalculator.t

import unittest
import filecmp
import os
import shutil

from clockwork import het_snp_caller

modules_dir = os.path.dirname(os.path.abspath(het_snp_caller.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "het_snp_caller")


class TestHetSnpCaller(unittest.TestCase):
    def test_run_mpileup(self):
        """test _run_mpileup"""
        # We'll check no errors are thrown and a VCF file gets made.
        # Don't look at contents of VCF, trust samtools
        pass

    def test_vcf_line_is_snp_and_or_het(self):
        """test _vcf_line_is_snp_and_or_het"""
        self.assertEqual(
            ("chrom1", True, True),
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,7,0;ADR=17,8,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            ),
        )

        self.assertEqual(
            ("chrom2", True, True),
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom2\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,7,0;ADR=17,8,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            ),
        )

        self.assertEqual(
            ("chrom1", False, False),
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,1,0;ADR=17,10,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            ),
        )

        self.assertEqual(
            ("chrom1", False, False),
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,10,0;ADR=17,1,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            ),
        )

        self.assertEqual(
            ("chrom1", False, False),
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,1,0;ADR=17,1,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            ),
        )

        self.assertEqual(
            ("chrom1", True, False),
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=1,22,0;ADR=1,17,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            ),
        )

        with self.assertRaises(Exception):
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=1,22,0;ADR=1,17,2,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            )

        with self.assertRaises(Exception):
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADR=1,2,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            )

        with self.assertRaises(Exception):
            het_snp_caller.HetSnpCaller._vcf_line_is_snp_and_or_het(
                "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=1,2,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255",
                4,
                2,
                0.9,
            )

    def test_filter_vcf_and_count_snps(self):
        """test _filter_vcf_and_count_snps"""
        infile = os.path.join(data_dir, "filter_vcf_and_count_snps.in.vcf")
        outfile = "tmp.filter_vcf_and_count_snps.out.vcf"
        got_counts = het_snp_caller.HetSnpCaller._filter_vcf_and_count_snps(
            infile, outfile, 4, 2, 0.9
        )
        expected_counts = {
            "chrom1": {"hets": 1, "positions": 5, "snps": 2},
            "chrom2": {"hets": 1, "positions": 1, "snps": 1},
        }
        self.assertEqual(expected_counts, got_counts)
        expected_vcf = os.path.join(data_dir, "filter_vcf_and_count_snps.expect.vcf")
        self.assertTrue(filecmp.cmp(expected_vcf, outfile, shallow=False))
        os.unlink(outfile)

    def test_write_reports(self):
        """test _write_reports"""
        contig_lengths = {"ctg1": 100, "ctg2": 200, "ctg3": 300}
        snp_data = {
            "ctg1": {"positions": 50, "snps": 2, "hets": 1},
            "ctg2": {"positions": 100, "snps": 0, "hets": 0},
        }

        got_per_contig_file = "tmp.het_snp_caller.write_reports.per_contig"
        got_summary_file = "tmp.het_snp_caller.write_reports.summary"
        het_snp_caller.HetSnpCaller._write_reports(
            snp_data, contig_lengths, got_summary_file, got_per_contig_file
        )
        expected_per_contig_file = os.path.join(data_dir, "write_reports.per_contig")
        expected_summary_file = os.path.join(data_dir, "write_reports.summary")
        self.assertTrue(
            filecmp.cmp(expected_per_contig_file, got_per_contig_file, shallow=False)
        )
        self.assertTrue(
            filecmp.cmp(expected_summary_file, got_summary_file, shallow=False)
        )
        os.unlink(got_per_contig_file)
        os.unlink(got_summary_file)

    def test_run(self):
        """test run"""
        sorted_bam = os.path.join(data_dir, "run.bam")
        ref_fasta = os.path.join(data_dir, "run.ref.fa")
        outprefix = "tmp.het_snp_caller.run"
        hsc = het_snp_caller.HetSnpCaller(sorted_bam, ref_fasta, outprefix)
        hsc.run()

        expected_per_contig = os.path.join(data_dir, "run.expect.per_contig.tsv")
        expected_summary = os.path.join(data_dir, "run.expect.summary.tsv")
        got_per_contig = outprefix + ".per_contig.tsv"
        got_summary = outprefix + ".summary.tsv"
        got_vcf = outprefix + ".het_calls.vcf"

        self.assertTrue(filecmp.cmp(expected_per_contig, got_per_contig, shallow=False))
        self.assertTrue(filecmp.cmp(expected_summary, got_summary, shallow=False))
        expected_lines = [
            "ref1\t251\t.\tC\tA,<*>\t0\t.\tDP=31;ADF=9,5,0;ADR=9,8,0;AD=18,13,0;I16=9,9,5,8,720,28800,520,20800,900,45000,572,25168,279,5591,203,4283;QS=0.580645,0.419355,0;VDB=0.78464;SGB=-0.683931;RPBZ=0.621123;MQBZ=-5.47723;MQSBZ=-0.626654;BQBZ=0;SCBZ=0;FS=0;MQ0F=0\tPL\t255,0,255,255,255,255\n"
        ]
        with open(got_vcf) as f:
            got_lines = [x for x in f if not x.startswith("#")]
        self.assertEqual(expected_lines, got_lines)
        os.unlink(got_per_contig)
        os.unlink(got_summary)
        os.unlink(got_vcf)
