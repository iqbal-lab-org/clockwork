import unittest
import os
import filecmp
from clockwork import read_trim

modules_dir = os.path.dirname(os.path.abspath(read_trim.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "read_trim")


# We need trimmomatic directory to run trimmomatic tests.
# Assume we're running the tests in the Singularity container, so
# look for trimmomatic directory in there, unless the env variable
# $CLOCKWORK_TRIMMO_DIR is set, in which case use that.
trimmo_root = os.environ.get("CLOCKWORK_TRIMMO_DIR", "/bioinf-tools/Trimmomatic-0.36")
# and if the direcotry is there, then we can run the test, otherwise skip
can_test_trimmo = os.path.isdir(trimmo_root)


class TestReadTrim(unittest.TestCase):
    @unittest.skipUnless(can_test_trimmo, "Trimmoatic dir not found")
    def test_run_trimmomatic(self):
        """test run_trimmomatic"""
        tmp_out1 = "test.read_trim.trimmo.1.fq"
        tmp_out2 = "test.read_trim.trimmo.2.fq"

        read_trim.run_trimmomatic(
            os.path.join(data_dir, "trimmomatic_reads_1.fq"),
            os.path.join(data_dir, "trimmomatic_reads_2.fq.gz"),
            tmp_out1,
            tmp_out2,
            trimmo_root,
            adapters="TruSeq3-PE-2.fa",
            minlen=20,
            verbose=0,
            threads=1,
            qual_trim="SLIDINGWINDOW:4:15",
            adapters_included=True,
        )

        expected_1 = os.path.join(data_dir, "trimmomatic_reads_expect_1.fq")
        expected_2 = os.path.join(data_dir, "trimmomatic_reads_expect_2.fq")

        self.assertTrue(filecmp.cmp(tmp_out1, expected_1, shallow=False))
        self.assertTrue(filecmp.cmp(tmp_out2, expected_2, shallow=False))
        os.unlink(tmp_out1)
        os.unlink(tmp_out2)
