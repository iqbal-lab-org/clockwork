import filecmp
import os
import random
import re
import unittest
import pyfastaq
from clockwork import utils, reference_dir, var_call_one_sample_pipeline


class TestVarCallOneSamplePipeline(unittest.TestCase):
    def test_run(self):
        root_outdir = "tmp.var_call_one_sample"
        utils.syscall(f"rm -rf {root_outdir}")
        os.mkdir(root_outdir)

        ref_fa = os.path.join(root_outdir, "ref.fa")
        ref_fa_mutated = f"{ref_fa}.mut.fa"
        random.seed(42)
        ref_seq = random.choices(["A", "C", "G", "T"], k=1000)
        ref_seq[499] = "A"
        with open(ref_fa, "w") as f:
            print(">ref", "".join(ref_seq), sep="\n", file=f)
        ref_fa_mutated = f"{ref_fa}.mut.fa"
        ref_seq[499] = "T"
        with open(ref_fa_mutated, "w") as f:
            print(">ref_mutated", "".join(ref_seq), sep="\n", file=f)

        reads1 = os.path.join(root_outdir, "reads1.fq")
        reads2 = os.path.join(root_outdir, "reads2.fq")
        utils.syscall(
            f"fastaq to_perfect_reads {ref_fa_mutated} - 200 1 20 75 | fastaq deinterleave - {reads1} {reads2}"
        )

        ref_dir = reference_dir.ReferenceDir(
            directory=os.path.join(root_outdir, "ref_dir")
        )
        ref_dir.make_index_files(ref_fa, False, True, cortex_mem_height=21)
        var_call_out = os.path.join(root_outdir, "varcall")
        var_call_one_sample_pipeline.run(
            [reads1],
            [reads2],
            ref_dir.directory,
            var_call_out,
            sample_name="test_sample",
            debug=False,
            keep_bam=True,
            cortex_mem_height=21,
        )

        got_files = sorted(list(os.listdir(var_call_out)))
        expect_files = [
            "cortex.vcf",
            "final.vcf",
            "map.bam",
            "map.bam.bai",
            "samtools.vcf",
        ]
        self.assertEqual(got_files, expect_files)

        with open(os.path.join(var_call_out, "final.vcf")) as f:
            calls = [x for x in f if not x.startswith("#")]
        self.assertEqual(len(calls), 1)
        fields = calls[0].split("\t")
        self.assertEqual(fields[1], "500")
        self.assertEqual(fields[3], "A")
        self.assertEqual(fields[4], "T")
        utils.syscall(f"rm -r {root_outdir}")
