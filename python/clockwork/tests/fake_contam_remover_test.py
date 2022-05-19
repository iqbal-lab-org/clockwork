import unittest
import filecmp
import os
import pathlib
import shutil
from clockwork import fake_contam_remover

modules_dir = os.path.dirname(os.path.abspath(fake_contam_remover.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "fake_contam_remover")


class TestFakeContamRemover(unittest.TestCase):
    def test_write_counts_tsv(self):
        """test _write_counts_tsv"""
        reads1 = os.path.join(data_dir, "reads_1.fq.gz")
        reads2 = os.path.join(data_dir, "reads_2.fq.gz")
        tmp_out = "tmp.fake_contam_remover.write_counts_tsv.out.tsv"
        expected = os.path.join(data_dir, "write_counts_tsv.expect.tsv")
        fake_contam_remover.FakeContamRemover._write_counts_tsv(reads1, reads2, tmp_out)
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)

    def test_symlink_reads_file(self):
        """test _symlink_reads_file"""
        tmpdir = "tmp.symlink_reads_file" ""
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)

        file_to_symlink = os.path.join(tmpdir, "foo")
        symlink_name = os.path.join(tmpdir, "bar")

        # foo doesn't exist
        with self.assertRaises(Exception):
            fake_contam_remover.FakeContamRemover._symlink_reads_file(
                file_to_symlink, symlink_name
            )

        with open(file_to_symlink, "w") as f:
            print("shrubbery", file=f)

        fake_contam_remover.FakeContamRemover._symlink_reads_file(
            file_to_symlink, symlink_name
        )
        self.assertTrue(os.path.exists(symlink_name))
        self.assertTrue(pathlib.Path(symlink_name).is_symlink())
        self.assertTrue(pathlib.Path(symlink_name).exists())
        self.assertEqual("foo", os.readlink(symlink_name))

        # bar is not in the same directory as foo
        with self.assertRaises(Exception):
            fake_contam_remover.FakeContamRemover._symlink_reads_file(
                file_to_symlink, "bar"
            )

        shutil.rmtree(tmpdir)

    def test_run(self):
        """test run"""
        tmpdir = "tmp.fake_contam_remover.run"
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        reads1 = os.path.join(tmpdir, "reads_1.fq.gz")
        reads2 = os.path.join(tmpdir, "reads_2.fq.gz")
        shutil.copyfile(os.path.join(data_dir, "reads_1.fq.gz"), reads1)
        shutil.copyfile(os.path.join(data_dir, "reads_2.fq.gz"), reads2)
        link1 = os.path.join(tmpdir, "remove_contam_1.fq.gz")
        link2 = os.path.join(tmpdir, "remove_contam_2.fq.gz")
        counts_file = os.path.join(tmpdir, "counts.tsv")

        cremover = fake_contam_remover.FakeContamRemover(
            reads1, reads2, link1, link2, counts_file
        )
        cremover.run()

        self.assertTrue(os.path.exists(link1))
        self.assertTrue(os.path.exists(link2))
        self.assertTrue(pathlib.Path(link1).is_symlink())
        self.assertTrue(pathlib.Path(link2).is_symlink())
        self.assertEqual("reads_1.fq.gz", os.readlink(link1))
        self.assertEqual("reads_2.fq.gz", os.readlink(link2))
        expect_counts = os.path.join(data_dir, "write_counts_tsv.expect.tsv")
        self.assertTrue(filecmp.cmp(expect_counts, counts_file))
        shutil.rmtree(tmpdir)
