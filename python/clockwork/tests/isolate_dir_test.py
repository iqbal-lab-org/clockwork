import unittest
import os
import shutil
from clockwork import isolate_dir

modules_dir = os.path.dirname(os.path.abspath(isolate_dir.__file__))
data_dir = os.path.join(modules_dir, "tests", "data")


class TestIsolateDir(unittest.TestCase):
    def test_sample_id_to_dir_numbers(self):
        """test sample_id_to_dir_numbers"""
        self.assertEqual(
            ["00", "00", "00", "00"], isolate_dir.IsolateDir.sample_id_to_dir_numbers(0)
        )
        self.assertEqual(
            ["00", "00", "00", "01"], isolate_dir.IsolateDir.sample_id_to_dir_numbers(1)
        )
        self.assertEqual(
            ["00", "00", "00", "99"],
            isolate_dir.IsolateDir.sample_id_to_dir_numbers(99),
        )
        self.assertEqual(
            ["00", "00", "01", "00"],
            isolate_dir.IsolateDir.sample_id_to_dir_numbers(100),
        )
        self.assertEqual(
            ["42", "00", "43", "00"],
            isolate_dir.IsolateDir.sample_id_to_dir_numbers(42004300),
        )
        self.assertEqual(
            ["12", "34", "56", "78"],
            isolate_dir.IsolateDir.sample_id_to_dir_numbers(12345678),
        )

    def test_make_sample_dir_path(self):
        """test _make_sample_dir_path"""
        tmp_root_dir = "tmp.isolate_dir.make_sample_dir_path"
        expected_path = os.path.join(tmp_root_dir, "42", "43", "44", "45")
        got_path = isolate_dir.IsolateDir._make_sample_dir_path(tmp_root_dir, 42434445)
        self.assertEqual(expected_path, got_path)

    def test_init(self):
        """test __init__"""
        tmp_root_dir = os.path.abspath("tmp.isolate_dir.make_sample_dir_path")
        sample_id = 12345678
        isolate_id = 42
        iso_dir = isolate_dir.IsolateDir(tmp_root_dir, sample_id, isolate_id)

        expected_sample_dir = isolate_dir.IsolateDir._make_sample_dir_path(
            tmp_root_dir, sample_id
        )
        self.assertEqual(expected_sample_dir, iso_dir.sample_dir)

        expected_isolate_dir = os.path.join(expected_sample_dir, str(isolate_id))
        self.assertEqual(expected_isolate_dir, iso_dir.isolate_dir)

        expected_reads_dir = os.path.join(expected_isolate_dir, "Reads")
        self.assertEqual(expected_reads_dir, iso_dir.reads_dir)

    def test_make_essential_dirs(self):
        """test make_essential_dirs"""
        tmp_root_dir = os.path.abspath("tmp.isolate_dir.make_sample_dir_path")
        if os.path.exists(tmp_root_dir):
            shutil.rmtree(tmp_root_dir)

        sample_id = 12345678
        isolate_id = 42
        iso_dir = isolate_dir.IsolateDir(tmp_root_dir, sample_id, isolate_id)
        iso_dir.make_essential_dirs()

        self.assertTrue(os.path.exists(iso_dir.sample_dir))
        self.assertTrue(os.path.exists(iso_dir.isolate_dir))
        self.assertTrue(os.path.exists(iso_dir.reads_dir))
        self.assertTrue(os.path.exists(iso_dir.pipelines_dir))
        shutil.rmtree(tmp_root_dir)

    def test_contamination_counts_filename(self):
        """test contamination_counts_filename"""
        root = os.path.abspath(os.path.join("knights", "ni"))
        sample_id = 42
        isolate_id = 11
        iso_dir = isolate_dir.IsolateDir(root, sample_id, isolate_id)
        expected = os.path.join(
            iso_dir.isolate_dir, "Reads", "reads.remove_contam.1.counts.tsv"
        )
        self.assertEqual(expected, iso_dir.contamination_counts_filename(1))
        expected = os.path.join(
            iso_dir.isolate_dir, "Reads", "reads.remove_contam.2.counts.tsv"
        )
        self.assertEqual(expected, iso_dir.contamination_counts_filename(2))

    def test_reads_filename(self):
        """test reads_filename"""
        root = os.path.abspath(os.path.join("spam", "shrubbery"))
        sample_id = 123456
        isolate_id = 42
        iso_dir = isolate_dir.IsolateDir(root, sample_id, isolate_id)
        self.assertEqual(
            os.path.join(iso_dir.isolate_dir, "Reads", "reads.original.11.1.fq.gz"),
            iso_dir.reads_filename("original", 11, 1),
        )
        self.assertEqual(
            os.path.join(
                iso_dir.isolate_dir, "Reads", "reads.remove_contam.12.2.fq.gz"
            ),
            iso_dir.reads_filename("remove_contam", 12, 2),
        )
        self.assertEqual(
            os.path.join(iso_dir.isolate_dir, "Reads", "reads.contam.11.1.fq.gz"),
            iso_dir.reads_filename("contam", 11, 1),
        )

        with self.assertRaises(Exception):
            iso_dir.reads_filename("original", 11, 42)
        with self.assertRaises(Exception):
            iso_dir.reads_filename("oops_wrong_type", 11, 1)

    def test_pipeline_dir(self):
        """test pipeline_dir"""
        root = os.path.abspath(os.path.join("dave", "lister"))
        sample_id = 42
        isolate_id = 1
        iso_dir = isolate_dir.IsolateDir(root, sample_id, isolate_id)
        self.assertEqual(
            os.path.join(iso_dir.isolate_dir, "Pipelines", "2", "name", "1.0.0"),
            iso_dir.pipeline_dir(2, "name", "1.0.0"),
        )
        self.assertEqual(
            os.path.join(iso_dir.isolate_dir, "Pipelines", "2", "name", "1.0.1.ref.42"),
            iso_dir.pipeline_dir(2, "name", "1.0.1", reference_id=42),
        )

    def test_xml_submission_file(self):
        """test xml_submission_file"""
        root = os.path.abspath(os.path.join("papa", "lazarou"))
        sample_id = 64738
        isolate_id = 42
        iso_dir = isolate_dir.IsolateDir(root, sample_id, isolate_id)
        self.assertEqual(
            os.path.join(iso_dir.sample_dir, "ena_sample_submission.xml"),
            iso_dir.xml_submission_file("sample"),
        )
        self.assertEqual(
            os.path.join(iso_dir.isolate_dir, "ena_experiment_submission.xml"),
            iso_dir.xml_submission_file("experiment"),
        )
        with self.assertRaises(Exception):
            iso_dir.xml_submission_file("run")
        self.assertEqual(
            os.path.join(
                iso_dir.reads_dir, "reads.remove_contam.42.ena_run_submission.xml"
            ),
            iso_dir.xml_submission_file("run", sequence_replicate=42),
        )
