import unittest
import copy
import datetime
import imp
import os
import shutil
import pyfastaq
from operator import itemgetter
from clockwork import (
    db,
    db_connection,
    db_maker,
    db_schema,
    isolate_dir,
    mykrobe,
    reference_dir,
    utils,
)

modules_dir = os.path.dirname(os.path.abspath(db.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "db")
from clockwork import __version__ as clockwork_version

ini_file = os.path.join(modules_dir, "tests", "data", "db.ini")


class TestDb(unittest.TestCase):
    def check_tsv(self, expected_lines, filename):
        with open(filename) as f:
            got = [x.rstrip() for x in f.readlines()]

        self.maxDiff = None
        self.assertEqual(expected_lines, got)

    def setUp(self):
        imp.reload(db_schema)  # because some tests change the schema
        try:
            db_connection.DbConnection(ini_file, destroy=True)
        except:
            pass

        dbm = db_maker.DbMaker(ini_file)
        dbm.run()
        self.db = db.Db(ini_file)

    def tearDown(self):
        self.db.commit_and_close()
        db_connection.DbConnection(ini_file, destroy=True, must_exist=True)

    def test_add_row_to_table(self):
        """test add_row_to_table"""
        db_schema.tables["Test"] = [
            ("col1", "integer"),
            ("col2", "integer"),
            ("col3", "text"),
        ]
        self.db.execute("""CREATE TABLE Test (col1 integer, col2 integer, col3 text)""")
        got = self.db.get_rows_from_table("Test")
        self.assertEqual([], got)
        data_dict1 = {"col1": 1, "col2": 0, "col3": "foo"}
        data_dict2 = {"col1": 2, "col2": 1, "col3": "bar"}
        self.db.add_row_to_table("Test", data_dict1)
        got = self.db.get_rows_from_table("Test")
        self.assertEqual([data_dict1], got)
        self.db.add_row_to_table("Test", data_dict2)
        got = self.db.get_rows_from_table("Test")
        self.assertEqual([data_dict1, data_dict2], got)

    def test_value_to_string(self):
        """test _value_to_string"""
        self.assertEqual("42", db.Db._value_to_string(42, "integer"))
        self.assertEqual('"42"', db.Db._value_to_string(42, "text"))
        self.assertEqual("NULL", db.Db._value_to_string(None, "text"))
        self.assertEqual("NULL", db.Db._value_to_string(None, "integer"))
        self.assertEqual('"2017-03-25"', db.Db._value_to_string("2017-03-25", "date"))

    def test_update_row(self):
        """test update_row"""
        db_schema.tables["Test"] = [
            ("col1", "integer"),
            ("col2", "integer"),
            ("col3", "text"),
        ]
        self.db.execute("""CREATE TABLE Test (col1 integer, col2 integer, col3 text)""")
        got = self.db.get_rows_from_table("Test")
        self.assertEqual([], got)
        rows = [(1, 0, 0), (2, 1, 0)]
        rows = [
            {"col1": 1, "col2": 3, "col3": "foo"},
            {"col1": 2, "col2": 4, "col3": "ni!"},
        ]
        self.db.add_row_to_table("Test", rows[0])
        self.db.add_row_to_table("Test", rows[1])
        got = self.db.get_rows_from_table("Test", order_by="col1")
        self.assertEqual(rows, got)

        self.db.update_row("Test", {"col1": 2}, {"col2": 3})
        rows[1]["col2"] = 3
        got = self.db.get_rows_from_table("Test", order_by="col1")
        self.assertEqual(rows, got)

        rows.append({"col1": 44, "col2": 5, "col3": "ni!"})
        self.db.add_row_to_table("Test", rows[-1])
        got = self.db.get_rows_from_table("Test", order_by="col1")
        self.assertEqual(rows, got)

        self.db.update_row("Test", {"col2": 3}, {"col3": "ping!"})
        rows[0]["col3"] = "ping!"
        rows[1]["col3"] = "ping!"
        got = self.db.get_rows_from_table("Test", order_by="col1")
        self.assertEqual(rows, got)

        self.db.update_row(
            "Test", {"col1": 2, "col2": 3}, {"col2": 100, "col3": "shrubbery"}
        )
        rows[1]["col2"] = 100
        rows[1]["col3"] = "shrubbery"
        got = self.db.get_rows_from_table("Test", order_by="col1")
        self.assertEqual(rows, got)

    def test_query_to_dict(self):
        """test query_to_dict"""
        db_schema.tables["Test"] = [("col1", "integer"), ("col2", "text")]
        self.db.execute("""CREATE TABLE Test (col1 integer, col2 text)""")
        got = self.db.query_to_dict("SELECT * FROM Test")
        self.assertEqual([], got)

        dict1 = {"col1": 1, "col2": "spam"}
        self.db.add_row_to_table("Test", dict1)
        got = self.db.query_to_dict("SELECT * FROM Test")
        self.assertEqual([dict1], got)

        dict2 = {"col1": 2, "col2": "eggs"}
        self.db.add_row_to_table("Test", dict2)
        got = self.db.query_to_dict("SELECT * FROM Test ORDER BY col1")
        self.assertEqual([dict1, dict2], got)

    def test_has_sample_isolate_seq_replicate(self):
        """test has_sample_isolate_seq_replicate"""
        self.assertFalse(self.db.has_sample_isolate_seq_replicate(42, "1", "2"))
        isolate_row = {
            "isolate_id": None,
            "sample_id": 42,
            "isolate_number_from_lab": "lab_isolate_id",
            "pool_sequence_replicates": 0,
            "ena_experiment_accession": None,
        }
        self.db.add_row_to_table("Isolate", isolate_row)
        seqrep_row = {
            "seqrep_id": None,
            "isolate_id": 1,
            "sequence_replicate_number": 3,
            "original_reads_file_1_md5": "abcde",
            "original_reads_file_2_md5": "abcdf",
            "remove_contam_reads_file_1_md5": "12345",
            "remove_contam_reads_file_2_md5": "12346",
            "withdrawn": 0,
            "import_status": 0,
            "submission_date": "20170101",
            "submit_to_ena": 1,
            "instrument_model": "Illumina HiSeq 2500",
            "ena_run_accession": None,
            "ena_on_hold": None,
        }
        self.db.add_row_to_table("Seqrep", seqrep_row)
        self.assertFalse(
            self.db.has_sample_isolate_seq_replicate(41, "lab_isolate_id", 2)
        )
        self.assertFalse(
            self.db.has_sample_isolate_seq_replicate(42, "lab_isolate_id", 2)
        )
        self.assertFalse(
            self.db.has_sample_isolate_seq_replicate(42, "wrong_lab_isolate_id", 3)
        )
        self.assertTrue(
            self.db.has_sample_isolate_seq_replicate(42, "lab_isolate_id", 3)
        )

    def test_get_sample_and_replicate_uniqueness(self):
        """Test _get_sample_and_replicate_uniqueness"""
        data = {
            "subject_id": "p1",
            "site_id": "s1",
            "lab_id": "l1",
            "isolate_number": "42",
            "sequence_replicate_number": "43",
            "submission_date": datetime.date(2017, 12, 25),
            "reads_file_1": "reads_1_1.fq",
            "reads_file_1_md5": "abcdefghijklmnopqrstuvwyx123456",
            "reads_file_2": "reads_1_2.fq",
            "reads_file_2_md5": "abcdefghijklmnopqrstuvwyx123457",
            "dataset_name": "g1",
            "submit_to_ena": "0",
            "instrument_model": "Illumina HiSeq 2500",
            "ena_center_name": "Centre 1",
            "ena_on_hold": "0",
            "ena_run_accession": "ERR123456",
            "ena_sample_accession": "ERS123456",
        }

        self.assertEqual(
            (True, False, None), self.db._get_sample_and_replicate_uniqueness(data)
        )
        self.db.add_one_seqrep(data)
        self.assertEqual(
            (True, True, 1), self.db._get_sample_and_replicate_uniqueness(data)
        )

        data["isolate_number"] = "41"
        self.assertEqual(
            (True, False, 1), self.db._get_sample_and_replicate_uniqueness(data)
        )
        data["isolate_number"] = "42"

        data["sequence_replicate_number"] = "44"
        self.assertEqual(
            (True, False, 1), self.db._get_sample_and_replicate_uniqueness(data)
        )
        data["sequence_replicate_number"] = "43"

        data["subject_id"] = "p2"
        self.db.add_one_seqrep(data)
        self.assertEqual(
            (True, True, 2), self.db._get_sample_and_replicate_uniqueness(data)
        )

    def test_add_one_seqrep(self):
        """test add_one_seqrep"""
        got = self.db.get_rows_from_table("Sample")
        self.assertEqual([], got)

        sample_dict = {
            "subject_id": "p1",
            "site_id": "s1",
            "lab_id": "l1",
            "isolate_number": "42",
            "sequence_replicate_number": 43,
            "submission_date": datetime.date(2017, 12, 25),
            "reads_file_1": "reads_43_1.fq",
            "reads_file_1_md5": "md5_43_1",
            "reads_file_2_md5": "md5_43_2",
            "reads_file_2": "reads_43_2.fq",
            "dataset_name": "g1",
            "submit_to_ena": "0",
            "instrument_model": "Illumina HiSeq 2500",
            "ena_center_name": "Centre 1",
            "ena_on_hold": "0",
            "ena_run_accession": "ERR123456",
            "ena_sample_accession": "ERS123456",
        }

        self.assertEqual((1, 1, 1), self.db.add_one_seqrep(sample_dict))
        expected_sample = [
            {
                "sample_id": 1,
                "subject_id": sample_dict["subject_id"],
                "site_id": sample_dict["site_id"],
                "sample_id_from_lab": sample_dict["lab_id"],
                "dataset_name": sample_dict["dataset_name"],
                "ena_center_name": sample_dict["ena_center_name"],
                "ena_sample_accession": sample_dict["ena_sample_accession"],
                "ena_study_accession": None,
            }
        ]
        got_sample = self.db.get_rows_from_table("Sample", order_by="sample_id")
        self.assertEqual(expected_sample, got_sample)

        expected_isolate = [
            {
                "isolate_id": 1,
                "sample_id": 1,
                "isolate_number_from_lab": sample_dict["isolate_number"],
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
            }
        ]
        got_isolate = self.db.get_rows_from_table("Isolate", order_by="isolate_id")
        self.assertEqual(expected_isolate, got_isolate)

        expected_seqrep = [
            {
                "seqrep_id": 1,
                "isolate_id": 1,
                "sequence_replicate_number": sample_dict["sequence_replicate_number"],
                "original_reads_file_1_md5": sample_dict["reads_file_1_md5"],
                "original_reads_file_2_md5": sample_dict["reads_file_2_md5"],
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "withdrawn": 0,
                "import_status": 0,
                "submission_date": sample_dict["submission_date"],
                "instrument_model": sample_dict["instrument_model"],
                "submit_to_ena": int(sample_dict["submit_to_ena"]),
                "ena_run_accession": sample_dict["ena_run_accession"],
                "ena_on_hold": int(sample_dict["ena_on_hold"]),
            }
        ]
        got_seqrep = self.db.get_rows_from_table("Seqrep", order_by="seqrep_id")
        self.assertEqual(expected_seqrep, got_seqrep)

        sample_dict = {
            "subject_id": "p2",
            "site_id": "s2",
            "lab_id": "l2",
            "isolate_number": "45",
            "sequence_replicate_number": 46,
            "submission_date": datetime.date(2017, 12, 25),
            "reads_file_1": "reads_46_1.fq",
            "reads_file_1_md5": "md5_46_1",
            "reads_file_2_md5": "md5_46_2",
            "reads_file_2": "reads_46_2.fq",
            "dataset_name": "g2",
            "instrument_model": "Illumina HiSeq 2500",
            "submit_to_ena": "1",
            "ena_center_name": "Centre 11",
            "ena_on_hold": "1",
            "ena_run_accession": None,
            "ena_sample_accession": None,
        }

        self.assertEqual((2, 2, 2), self.db.add_one_seqrep(sample_dict))
        expected_sample.append(
            {
                "sample_id": 2,
                "subject_id": sample_dict["subject_id"],
                "site_id": sample_dict["site_id"],
                "sample_id_from_lab": sample_dict["lab_id"],
                "dataset_name": sample_dict["dataset_name"],
                "ena_center_name": sample_dict["ena_center_name"],
                "ena_sample_accession": sample_dict["ena_sample_accession"],
                "ena_study_accession": None,
            }
        )
        got_sample = self.db.get_rows_from_table("Sample", order_by="sample_id")
        self.assertEqual(expected_sample, got_sample)

        expected_isolate.append(
            {
                "isolate_id": 2,
                "sample_id": 2,
                "isolate_number_from_lab": sample_dict["isolate_number"],
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
            }
        )
        got_isolate = self.db.get_rows_from_table("Isolate", order_by="isolate_id")
        self.assertEqual(expected_isolate, got_isolate)

        expected_seqrep.append(
            {
                "seqrep_id": 2,
                "isolate_id": 2,
                "sequence_replicate_number": sample_dict["sequence_replicate_number"],
                "original_reads_file_1_md5": sample_dict["reads_file_1_md5"],
                "original_reads_file_2_md5": sample_dict["reads_file_2_md5"],
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "withdrawn": 0,
                "import_status": 0,
                "submission_date": sample_dict["submission_date"],
                "instrument_model": sample_dict["instrument_model"],
                "submit_to_ena": int(sample_dict["submit_to_ena"]),
                "ena_run_accession": sample_dict["ena_run_accession"],
                "ena_on_hold": int(sample_dict["ena_on_hold"]),
            }
        )
        got_seqrep = self.db.get_rows_from_table("Seqrep", order_by="seqrep_id")
        self.assertEqual(expected_seqrep, got_seqrep)

    def test_make_remove_contam_jobs_tsv(self):
        """test make_remove_contam_jobs_tsv"""
        tmp_out = "tmp.db_test.isolates_to_remove_contam_jobs_tsv.tsv"

        def make_data_dict(i):
            return {
                "subject_id": "subj" + str(i),
                "site_id": "site" + str(i),
                "lab_id": "lab" + str(i),
                "isolate_number": "il" + str(i),
                "sequence_replicate_number": i,
                "submission_date": datetime.date(2017, 1, 1),
                "reads_file_1": "reads1",
                "reads_file_1_md5": "md5." + str(i) + ".1",
                "reads_file_2": "reads2",
                "reads_file_2_md5": "md5." + str(i) + ".2",
                "dataset_name": "group1",
                "instrument_model": "Illumina HiSeq 2500",
                "submit_to_ena": 0,
                "ena_center_name": "Centre " + str(i),
                "ena_on_hold": 0,
                "ena_run_accession": "ERR" + str(i),
                "ena_sample_accession": "ERS" + str(i),
            }

        # import_status = 0, so should ignore
        data_dict1 = make_data_dict(1)
        self.db.add_one_seqrep(data_dict1)

        # withdrawn, so should ignore
        data_dict2 = make_data_dict(2)
        self.db.add_one_seqrep(data_dict2)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.2.1"},
            {"import_status": 1, "withdrawn": 1},
        )

        # these rows pass filter (but different values in Pipeline table)
        data_dict3 = make_data_dict(3)
        self.db.add_one_seqrep(data_dict3)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.3.1"}, {"import_status": 1}
        )
        pipeline_row3 = {
            "isolate_id": 3,
            "seqrep_id": 3,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 0,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row3)

        data_dict4 = make_data_dict(4)
        self.db.add_one_seqrep(data_dict4)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.4.1"}, {"import_status": 1}
        )
        pipeline_row4 = {
            "isolate_id": 4,
            "seqrep_id": 4,
            "seqrep_pool": None,
            "version": "1.0.1",
            "pipeline_name": "remove_contam",
            "status": 0,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row4)

        # Should always get this one because there's no entry in the Pipeline table
        data_dict5 = make_data_dict(5)
        data_dict5["dataset_name"] = "group1"
        self.db.add_one_seqrep(data_dict5)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.5.1"}, {"import_status": 1}
        )

        # Should always get this one because there's no entry in the Pipeline table
        # (but it's from group2 instead of group1)
        data_dict5 = make_data_dict(6)
        data_dict5["dataset_name"] = "group2"
        self.db.add_one_seqrep(data_dict5)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.6.1"}, {"import_status": 1}
        )

        # Should not get this because status=0
        data_dict6 = make_data_dict(7)
        data_dict6["dataset_name"] = "group2"
        self.db.add_one_seqrep(data_dict6)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.7.1"}, {"import_status": 1}
        )
        pipeline_row6 = {
            "isolate_id": 7,
            "seqrep_id": 7,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row6)

        expected_rows = self.db.get_rows_from_table("Pipeline")
        pipeline_root = "/root/"
        self.db.add_reference("ref_name")
        refs_root = "/refs/"
        refdir = reference_dir.ReferenceDir(
            pipeline_references_root_dir=refs_root, reference_id=1
        )

        expected_tsv_header = "\t".join(
            [
                "reads_in1",
                "reads_in2",
                "counts_tsv",
                "reads_contam1",
                "reads_contam2",
                "reads_remove_contam1",
                "reads_remove_contam2",
                "sample_id",
                "seqrep_id",
                "isolate_id",
                "sequence_replicate_number",
                "reference_id",
                "ref_fasta",
                "contam_tsv",
            ]
        )

        def make_expected_tsv_line(i):
            iso_dir = isolate_dir.IsolateDir(pipeline_root, i, i)
            return "\t".join(
                [
                    iso_dir.reads_filename("original", i, 1),
                    iso_dir.reads_filename("original", i, 2),
                    iso_dir.contamination_counts_filename(i),
                    iso_dir.reads_filename("contam", i, 1),
                    iso_dir.reads_filename("contam", i, 2),
                    iso_dir.reads_filename("remove_contam", i, 1),
                    iso_dir.reads_filename("remove_contam", i, 2),
                    str(i),
                    str(i),
                    str(i),
                    str(i),
                    "1",
                    refdir.ref_fasta,
                    refdir.remove_contam_metadata_tsv,
                ]
            )

        expected_tsv_2 = make_expected_tsv_line(2)
        expected_tsv_3 = make_expected_tsv_line(3)
        expected_tsv_4 = make_expected_tsv_line(4)
        expected_tsv_5 = make_expected_tsv_line(5)
        expected_tsv_6 = make_expected_tsv_line(6)

        tmp_out = "tmp.db_test.make_remove_contam_jobs_tsv"
        self.db.make_remove_contam_jobs_tsv(
            tmp_out, pipeline_root, 1, refs_root, dataset_name="group1"
        )
        expected_lines = [expected_tsv_header, expected_tsv_5]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.append(
            {
                "isolate_id": 5,
                "seqrep_id": 5,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "remove_contam",
                "status": 0,
                "reference_id": 1,
            }
        )
        expected_rows.sort(key=itemgetter("seqrep_id"))
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        self.assertEqual(expected_rows, got_rows)

        # getting jobs for the same pipeline again should return nothing
        self.db.make_remove_contam_jobs_tsv(
            tmp_out, pipeline_root, 1, refs_root, dataset_name="group1"
        )
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        self.assertEqual(expected_rows, got_rows)
        self.check_tsv([expected_tsv_header], tmp_out)
        os.unlink(tmp_out)

        # different dataset_name should get more jobs
        self.db.make_remove_contam_jobs_tsv(
            tmp_out, pipeline_root, 1, refs_root, dataset_name="group2"
        )
        expected_lines = [expected_tsv_header, expected_tsv_6]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.append(
            {
                "isolate_id": 6,
                "seqrep_id": 6,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "remove_contam",
                "status": 0,
                "reference_id": 1,
            }
        )
        expected_rows.sort(key=itemgetter("seqrep_id"))
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        self.assertEqual(expected_rows, got_rows)

    def test_make_qc_jobs_tsv(self):
        """test make_qc_jobs_tsv"""
        tmp_out = "tmp.db_test.make_qc_jobs_tsv.tsv"

        def make_data_dict(i):
            return {
                "subject_id": "subj" + str(i),
                "site_id": "site" + str(i),
                "lab_id": "lab" + str(i),
                "isolate_number": "il" + str(i),
                "sequence_replicate_number": i,
                "submission_date": datetime.date(2017, 1, 1),
                "reads_file_1": "reads1",
                "reads_file_1_md5": "md5." + str(i) + ".1",
                "reads_file_2": "reads2",
                "reads_file_2_md5": "md5." + str(i) + ".2",
                "dataset_name": "group1",
                "instrument_model": "Illumina HiSeq 2500",
                "submit_to_ena": 0,
                "ena_center_name": "Centre " + str(i),
                "ena_on_hold": 0,
                "ena_run_accession": "ERR" + str(i),
                "ena_sample_accession": "ERS" + str(i),
            }

        # import_status = 0, so should ignore
        data_dict1 = make_data_dict(1)
        self.db.add_one_seqrep(data_dict1)

        # withdrawn, so should ignore even though remove_contam has been run
        data_dict2 = make_data_dict(2)
        self.db.add_one_seqrep(data_dict2)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.2.1"},
            {"import_status": 1, "withdrawn": 1},
        )
        pipeline_row2 = {
            "isolate_id": 2,
            "seqrep_id": 2,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row2)

        # Have not had remove_contam run, so should get ignored
        data_dict3 = make_data_dict(3)
        self.db.add_one_seqrep(data_dict3)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.3.1"}, {"import_status": 1}
        )
        data_dict4 = make_data_dict(4)
        self.db.add_one_seqrep(data_dict4)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.4.1"}, {"import_status": 1}
        )
        pipeline_row4 = {
            "isolate_id": 4,
            "seqrep_id": 4,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 0,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row4)

        # Should not get used for version 1.0.0. Has had remove_contam run, but already has qc run as well
        data_dict5 = make_data_dict(5)
        self.db.add_one_seqrep(data_dict5)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.5.1"}, {"import_status": 1}
        )
        pipeline_row5 = {
            "isolate_id": 5,
            "seqrep_id": 5,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row5)
        pipeline_row5["pipeline_name"] = "qc"
        self.db.add_row_to_table("Pipeline", pipeline_row5)

        #  Should get picked up to run - has had remove_contam run
        data_dict6 = make_data_dict(6)
        data_dict6["dataset_name"] = "group2"
        self.db.add_one_seqrep(data_dict6)
        self.db.update_row(
            "Seqrep", {"original_reads_file_1_md5": "md5.6.1"}, {"import_status": 1}
        )
        pipeline_row6 = {
            "isolate_id": 6,
            "seqrep_id": 6,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row6)

        expected_rows = self.db.get_rows_from_table("Pipeline")
        pipeline_root = "tmp.db_test.make_qc_jobs.pipeline_root"
        refs_root = "/refs/"
        self.db.add_reference("ref_name")
        refdir = reference_dir.ReferenceDir(
            pipeline_references_root_dir=refs_root, reference_id=1
        )
        expected_tsv_header = "\t".join(
            [
                "reads_in1",
                "reads_in2",
                "output_dir",
                "sample_id",
                "seqrep_id",
                "isolate_id",
                "sequence_replicate_number",
                "reference_id",
                "ref_fasta",
            ]
        )

        def make_expected_tsv_line(pipeline_root, i, pipeline_version):
            iso_dir = isolate_dir.IsolateDir(pipeline_root, i, i)
            return "\t".join(
                [
                    iso_dir.reads_filename("remove_contam", i, 1),
                    iso_dir.reads_filename("remove_contam", i, 2),
                    iso_dir.pipeline_dir(i, "qc", pipeline_version),
                    str(i),
                    str(i),
                    str(i),
                    str(i),
                    "1",
                    refdir.ref_fasta,
                ]
            )

        expected_tsv_4 = make_expected_tsv_line(pipeline_root, 4, "1.0.1")
        expected_tsv_5_0 = make_expected_tsv_line(pipeline_root, 5, "1.0.0")
        expected_tsv_5_1 = make_expected_tsv_line(pipeline_root, 5, "1.0.1")
        expected_tsv_6_0 = make_expected_tsv_line(pipeline_root, 6, "1.0.0")
        expected_tsv_6_1 = make_expected_tsv_line(pipeline_root, 6, "1.0.1")

        tmp_out = "tmp.db_test.make_qc_jobs_tsv"
        self.db.make_qc_jobs_tsv(
            tmp_out, pipeline_root, 1, refs_root, pipeline_version="1.0.0"
        )
        expected_lines = [expected_tsv_header, expected_tsv_6_0]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.append(
            {
                "isolate_id": 6,
                "seqrep_id": 6,
                "seqrep_pool": None,
                "version": "1.0.0",
                "pipeline_name": "qc",
                "status": 0,
                "reference_id": 1,
            }
        )
        expected_rows.sort(key=itemgetter("seqrep_id"))
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        self.assertEqual(expected_rows, got_rows)

        # getting jobs for the same pipeline again should return nothing
        self.db.make_qc_jobs_tsv(
            tmp_out, pipeline_root, 1, refs_root, pipeline_version="1.0.0"
        )
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        self.assertEqual(expected_rows, got_rows)
        self.check_tsv([expected_tsv_header], tmp_out)
        os.unlink(tmp_out)

        # different pipeline version should get more jobs, but test get only group2
        self.db.make_qc_jobs_tsv(
            tmp_out,
            pipeline_root,
            1,
            refs_root,
            pipeline_version="1.0.1",
            dataset_name="group2",
        )
        expected_lines = [expected_tsv_header, expected_tsv_6_1]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.append(
            {
                "isolate_id": 6,
                "seqrep_id": 6,
                "seqrep_pool": None,
                "version": "1.0.1",
                "pipeline_name": "qc",
                "status": 0,
                "reference_id": 1,
            }
        )
        expected_rows.sort(key=itemgetter("seqrep_id"))
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        self.assertEqual(expected_rows, got_rows)

        shutil.rmtree(pipeline_root)

    def test_make_variant_call_or_mykrobe_jobs_tsv(self):
        """test make_variant_call_or_mykrobe_jobs_tsv"""
        tmp_out = "tmp.db_test.make_variant_call_or_mykrobe_jobs_tsv.tsv"

        def make_data_dict(i):
            return {
                "subject_id": "subj" + str(i),
                "site_id": "site" + str(i),
                "lab_id": "lab" + str(i),
                "isolate_number": "il" + str(i),
                "sequence_replicate_number": i,
                "submission_date": datetime.date(2017, 1, 1),
                "reads_file_1": "reads1",
                "reads_file_1_md5": "md5." + str(i) + ".1",
                "reads_file_2": "reads2",
                "reads_file_2_md5": "md5." + str(i) + ".2",
                "dataset_name": "group1",
                "instrument_model": "Illumina HiSeq 2500",
                "submit_to_ena": 0,
                "ena_center_name": "Centre 42",
                "ena_on_hold": 0,
                "ena_run_accession": "ERR" + str(i),
                "ena_sample_accession": "ERS" + str(i),
            }

        # import_status = 0, so should ignore
        data_dict1 = make_data_dict(1)
        self.db.add_one_seqrep(data_dict1)

        # withdrawn, so should ignore
        data_dict2 = make_data_dict(2)
        self.db.add_one_seqrep(data_dict2)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.2.1"},
            {"import_status": 1, "withdrawn": 1},
        )

        # Have not had remove_contam run, so should get ignored
        data_dict3 = make_data_dict(3)
        self.db.add_one_seqrep(data_dict3)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.3.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        data_dict4 = make_data_dict(4)
        self.db.add_one_seqrep(data_dict4)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.4.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row4 = {
            "isolate_id": 4,
            "seqrep_id": 4,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 0,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row4)

        # Should not get used for version 1.0.0. Has had remove_contam run, but already has variant_call run as well
        data_dict5 = make_data_dict(5)
        data_dict5["sequence_replicate_number"] = 1
        self.db.add_one_seqrep(data_dict5)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.5.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row5 = {
            "isolate_id": 5,
            "seqrep_id": 5,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row5)
        pipeline_row5.update(
            {"pipeline_name": "variant_call", "seqrep_id": None, "seqrep_pool": "1"}
        )
        self.db.add_row_to_table("Pipeline", pipeline_row5)

        # Should not get used for version 1.0.1. Has had remove_contam run, but already has variant_call run as well
        # This one is not pooled
        data_dict6 = make_data_dict(6)
        data_dict6["sequence_replicate_number"] = 1
        self.db.add_one_seqrep(data_dict6)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.6.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row6 = {
            "isolate_id": 6,
            "seqrep_id": 6,
            "seqrep_pool": None,
            "version": "1.0.1",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row6)
        pipeline_row6["pipeline_name"] = "variant_call"
        self.db.add_row_to_table("Pipeline", pipeline_row6)
        self.db.update_row(
            "Isolate", {"isolate_id": 6}, {"pool_sequence_replicates": 0}
        )

        #  Should get picked up to run - have remove_contam run. Two sequencing replicates from same isolate that should not be pooled
        data_dict7 = make_data_dict(7)
        data_dict7["sequence_replicate_number"] = 1
        self.db.add_one_seqrep(data_dict7)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.7.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row7 = {
            "isolate_id": 7,
            "seqrep_id": 7,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row7)
        data_dict8 = copy.copy(data_dict7)
        data_dict8["sequence_replicate_number"] = 2
        data_dict8["reads_file_1_md5"] = "md5.8.1"
        data_dict8["reads_file_2_md5"] = "md5.8.1"
        self.db.add_one_seqrep(data_dict8)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.8.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row8 = {
            "isolate_id": 7,
            "seqrep_id": 8,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row8)
        self.db.update_row(
            "Isolate", {"isolate_id": 7}, {"pool_sequence_replicates": 0}
        )

        #  Should get picked up to run - have remove_contam run. Two sequencing replicates from same isolate that need to be pooled
        data_dict9 = make_data_dict(9)
        data_dict9["sequence_replicate_number"] = 1
        data_dict9["dataset_name"] = "group2"
        self.db.add_one_seqrep(data_dict9)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.9.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row9 = {
            "isolate_id": 9,
            "seqrep_id": 9,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row9)
        data_dict10 = copy.copy(data_dict9)
        data_dict10["sequence_replicate_number"] = 2
        data_dict10["reads_file_1_md5"] = "md5.10.1"
        data_dict10["reads_file_2_md5"] = "md5.10.1"
        self.db.add_one_seqrep(data_dict10)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.10.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row10 = {
            "isolate_id": 9,
            "seqrep_id": 10,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row10)

        expected_rows = self.db.get_rows_from_table("Pipeline")
        pipeline_root = "tmp.db_test.make_variant_call_jobs.pipeline_root"
        self.db.add_reference("ref_name1")
        self.db.add_reference("ref_name2")
        refs_root = "/refs/"
        refdir1 = reference_dir.ReferenceDir(
            pipeline_references_root_dir=refs_root, reference_id=1
        )
        refdir2 = reference_dir.ReferenceDir(
            pipeline_references_root_dir=refs_root, reference_id=2
        )
        expected_tsv_header = "\t".join(
            [
                "reads_in1",
                "reads_in2",
                "output_dir",
                "sample_id",
                "pool",
                "isolate_id",
                "seqrep_id",
                "sequence_replicate_number",
                "reference_id",
                "reference_dir",
            ]
        )

        def make_expected_tsv_line(
            pipeline_root,
            i,
            pool,
            seqrep_ids,
            sequence_replicate_numbers,
            pipeline_version,
            reference_id,
            reference_dir,
            sample,
        ):
            iso_dir = isolate_dir.IsolateDir(pipeline_root, i, i)
            sequence_replicate_numbers_string = "_".join(
                [str(x) for x in sequence_replicate_numbers]
            )
            return "\t".join(
                [
                    " ".join(
                        [
                            iso_dir.reads_filename("remove_contam", j, 1)
                            for j in sequence_replicate_numbers
                        ]
                    ),
                    " ".join(
                        [
                            iso_dir.reads_filename("remove_contam", j, 2)
                            for j in sequence_replicate_numbers
                        ]
                    ),
                    iso_dir.pipeline_dir(
                        sequence_replicate_numbers_string,
                        "variant_call",
                        pipeline_version,
                        reference_id=reference_id,
                    ),
                    sample,
                    pool,
                    str(i),
                    "_".join([str(x) for x in seqrep_ids]),
                    "_".join([str(x) for x in sequence_replicate_numbers]),
                    str(reference_id),
                    reference_dir.directory,
                ]
            )

        expected_tsv_5 = make_expected_tsv_line(
            pipeline_root,
            5,
            "1",
            [5],
            [1],
            "1.0.1",
            1,
            refdir1,
            "site.site5.iso.il5.subject.subj5.lab_id.lab5.seq_reps.1",
        )
        expected_tsv_6 = make_expected_tsv_line(
            pipeline_root,
            6,
            "0",
            [6],
            [1],
            "1.0.0",
            1,
            refdir1,
            "site.site6.iso.il6.subject.subj6.lab_id.lab6.seq_reps.1",
        )
        expected_tsv_7_1_0 = make_expected_tsv_line(
            pipeline_root,
            7,
            "0",
            [7],
            [1],
            "1.0.0",
            1,
            refdir1,
            "site.site7.iso.il7.subject.subj7.lab_id.lab7.seq_reps.1",
        )
        expected_tsv_7_2_0 = make_expected_tsv_line(
            pipeline_root,
            7,
            "0",
            [8],
            [2],
            "1.0.0",
            1,
            refdir1,
            "site.site7.iso.il7.subject.subj7.lab_id.lab7.seq_reps.2",
        )
        expected_tsv_7_1_1 = make_expected_tsv_line(
            pipeline_root,
            7,
            "0",
            [7],
            [1],
            "1.0.1",
            1,
            refdir1,
            "site.site7.iso.il7.subject.subj7.lab_id.lab7.seq_reps.1",
        )
        expected_tsv_7_2_1 = make_expected_tsv_line(
            pipeline_root,
            7,
            "0",
            [8],
            [2],
            "1.0.1",
            1,
            refdir1,
            "site.site7.iso.il7.subject.subj7.lab_id.lab7.seq_reps.2",
        )
        expected_tsv_9_pool_0 = make_expected_tsv_line(
            pipeline_root,
            8,
            "1",
            [9, 10],
            [1, 2],
            "1.0.0",
            1,
            refdir1,
            "site.site9.iso.il9.subject.subj9.lab_id.lab9.seq_reps.1_2",
        )
        expected_tsv_9_pool_1 = make_expected_tsv_line(
            pipeline_root,
            8,
            "1",
            [9, 10],
            [1, 2],
            "1.0.1",
            1,
            refdir1,
            "site.site9.iso.il9.subject.subj9.lab_id.lab9.seq_reps.1_2",
        )

        tmp_out = "tmp.db_test.make_varian_call_jobs_tsv"
        self.db.make_variant_call_or_mykrobe_jobs_tsv(
            "variant_call",
            tmp_out,
            pipeline_root,
            1,
            refs_root,
            pipeline_version="1.0.0",
        )
        expected_lines = [
            expected_tsv_header,
            expected_tsv_6,
            expected_tsv_7_1_0,
            expected_tsv_7_2_0,
            expected_tsv_9_pool_0,
        ]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.extend(
            [
                {
                    "isolate_id": 6,
                    "seqrep_id": 6,
                    "seqrep_pool": None,
                    "version": "1.0.0",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
                {
                    "isolate_id": 7,
                    "seqrep_id": 7,
                    "seqrep_pool": None,
                    "version": "1.0.0",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
                {
                    "isolate_id": 7,
                    "seqrep_id": 8,
                    "seqrep_pool": None,
                    "version": "1.0.0",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
                {
                    "isolate_id": 8,
                    "seqrep_id": None,
                    "seqrep_pool": "1_2",
                    "version": "1.0.0",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
            ]
        )

        def sort_pipeline_rows(rows):
            rows.sort(
                key=lambda x: (
                    x["isolate_id"],
                    -1 if x["seqrep_id"] is None else x["seqrep_id"],
                    -1 if x["seqrep_pool"] is None else x["seqrep_pool"],
                )
            )

        sort_pipeline_rows(expected_rows)
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)

        # getting jobs for the same pipeline again should return nothing
        self.db.make_variant_call_or_mykrobe_jobs_tsv(
            "variant_call",
            tmp_out,
            pipeline_root,
            1,
            refs_root,
            pipeline_version="1.0.0",
        )
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)
        self.check_tsv([expected_tsv_header], tmp_out)
        os.unlink(tmp_out)

        # different pipeline version should get more jobs
        self.db.make_variant_call_or_mykrobe_jobs_tsv(
            "variant_call",
            tmp_out,
            pipeline_root,
            1,
            refs_root,
            pipeline_version="1.0.1",
        )
        expected_lines = [
            expected_tsv_header,
            expected_tsv_5,
            expected_tsv_7_1_1,
            expected_tsv_7_2_1,
            expected_tsv_9_pool_1,
        ]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.extend(
            [
                {
                    "isolate_id": 5,
                    "seqrep_id": None,
                    "seqrep_pool": "1",
                    "version": "1.0.1",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
                {
                    "isolate_id": 7,
                    "seqrep_id": 7,
                    "seqrep_pool": None,
                    "version": "1.0.1",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
                {
                    "isolate_id": 7,
                    "seqrep_id": 8,
                    "seqrep_pool": None,
                    "version": "1.0.1",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
                {
                    "isolate_id": 8,
                    "seqrep_id": None,
                    "seqrep_pool": "1_2",
                    "version": "1.0.1",
                    "pipeline_name": "variant_call",
                    "status": 0,
                    "reference_id": 1,
                },
            ]
        )

        sort_pipeline_rows(expected_rows)
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)

        # difference reference version (but same pipeline version) should get more jobs, but limit to group2
        self.db.make_variant_call_or_mykrobe_jobs_tsv(
            "variant_call",
            tmp_out,
            pipeline_root,
            2,
            refs_root,
            pipeline_version="1.0.1",
            dataset_name="group2",
        )
        expected_lines = [expected_tsv_header, expected_tsv_9_pool_1]
        expected_lines = [
            x.replace("\t1\t/refs/1", "\t2\t/refs/2").replace("1.ref.1", "1.ref.2")
            for x in expected_lines
        ]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.append(
            {
                "isolate_id": 8,
                "seqrep_id": None,
                "seqrep_pool": "1_2",
                "version": "1.0.1",
                "pipeline_name": "variant_call",
                "status": 0,
                "reference_id": 2,
            }
        )
        sort_pipeline_rows(expected_rows)
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)

    def test_make_generic_pipeline_jobs_tsv(self):
        """test make_generic_pipeline_jobs_tsv"""
        pipeline_name = "test_pipeline"
        pipeline_root = "tmp.db_test.make_generic_pipeline_jobs.pipeline_root"

        def make_data_dict(i):
            return {
                "subject_id": "subj" + str(i),
                "site_id": "site" + str(i),
                "lab_id": "lab" + str(i),
                "isolate_number": "il" + str(i),
                "sequence_replicate_number": i,
                "submission_date": datetime.date(2017, 1, 1),
                "reads_file_1": "reads1",
                "reads_file_1_md5": "md5." + str(i) + ".1",
                "reads_file_2": "reads2",
                "reads_file_2_md5": "md5." + str(i) + ".2",
                "dataset_name": "group1",
                "instrument_model": "Illumina HiSeq 2500",
                "submit_to_ena": 0,
                "ena_center_name": "Centre 42",
                "ena_on_hold": 0,
                "ena_run_accession": "ERR" + str(i),
                "ena_sample_accession": "ERS" + str(i),
            }

        # import_status = 0, so should ignore
        data_dict1 = make_data_dict(1)
        self.db.add_one_seqrep(data_dict1)

        # withdrawn, so should ignore
        data_dict2 = make_data_dict(2)
        self.db.add_one_seqrep(data_dict2)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.2.1"},
            {"import_status": 1, "withdrawn": 1},
        )

        # Have not had remove_contam run, so should get ignored
        data_dict3 = make_data_dict(3)
        self.db.add_one_seqrep(data_dict3)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.3.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        data_dict4 = make_data_dict(4)
        self.db.add_one_seqrep(data_dict4)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.4.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        # pipeline_row4 = {'isolate_id': 4, 'seqrep_id': 4, 'seqrep_pool': None, 'version': '1.0.0',
        #    'pipeline_name': pipeline_name, 'status': 0, 'reference_id': 1}
        # self.db.add_row_to_table('Pipeline', pipeline_row4)

        #  Have been run on version 1.0.0 already so should get ignored
        data_dict5 = make_data_dict(5)
        data_dict5["sequence_replicate_number"] = 1
        self.db.add_one_seqrep(data_dict5)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.5.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row5 = {
            "isolate_id": 5,
            "seqrep_id": 5,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row5)
        pipeline_row5["pipeline_name"] = pipeline_name
        self.db.add_row_to_table("Pipeline", pipeline_row5)

        # Has not been run, so should get picked up. Two sequencing replicates from the
        # same isolate that should not be pooled
        data_dict6 = make_data_dict(6)
        data_dict6["sequence_replicate_number"] = 1
        self.db.add_one_seqrep(data_dict6)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.6.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row6 = {
            "isolate_id": 6,
            "seqrep_id": 6,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row6)
        data_dict7 = copy.copy(data_dict6)
        data_dict7["sequence_replicate_number"] = 2
        data_dict7["reads_file_1_md5"] = "md5.7.1"
        data_dict7["reads_file_2_md5"] = "md5.7.1"
        self.db.add_one_seqrep(data_dict7)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.7.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row7 = {
            "isolate_id": 6,
            "seqrep_id": 7,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row7)
        self.db.update_row(
            "Isolate", {"isolate_id": 6}, {"pool_sequence_replicates": 0}
        )

        # Has not been run, so should get picked up. Two sequencing replicates from the
        # same isolate that should be pooled
        data_dict8 = make_data_dict(8)
        data_dict8["sequence_replicate_number"] = 1
        data_dict8["dataset_name"] = "group2"
        self.db.add_one_seqrep(data_dict8)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.8.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row8 = {
            "isolate_id": 8,
            "seqrep_id": 8,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row8)
        data_dict9 = copy.copy(data_dict8)
        data_dict9["sequence_replicate_number"] = 2
        data_dict9["reads_file_1_md5"] = "md5.9.1"
        data_dict9["reads_file_2_md5"] = "md5.9.1"
        self.db.add_one_seqrep(data_dict9)
        self.db.update_row(
            "Seqrep",
            {"original_reads_file_1_md5": "md5.9.1"},
            {"import_status": 1, "withdrawn": 0},
        )
        pipeline_row9 = {
            "isolate_id": 8,
            "seqrep_id": 9,
            "seqrep_pool": None,
            "version": "1.0.0",
            "pipeline_name": "remove_contam",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row9)

        expected_rows = self.db.get_rows_from_table("Pipeline")
        expected_tsv_header = "\t".join(
            [
                "reads_in1",
                "reads_in2",
                "output_dir",
                "sample_id",
                "pool",
                "isolate_id",
                "seqrep_id",
                "sequence_replicate_number",
            ]
        )

        def make_expected_tsv_line(
            pipeline_root,
            i,
            pool,
            seqrep_ids,
            sequence_replicate_numbers,
            pipeline_version,
        ):
            iso_dir = isolate_dir.IsolateDir(pipeline_root, i, i)
            sequence_replicate_numbers_string = "_".join(
                [str(x) for x in sequence_replicate_numbers]
            )
            return "\t".join(
                [
                    " ".join(
                        [
                            iso_dir.reads_filename("remove_contam", j, 1)
                            for j in sequence_replicate_numbers
                        ]
                    ),
                    " ".join(
                        [
                            iso_dir.reads_filename("remove_contam", j, 2)
                            for j in sequence_replicate_numbers
                        ]
                    ),
                    iso_dir.pipeline_dir(
                        sequence_replicate_numbers_string,
                        pipeline_name,
                        pipeline_version,
                    ),
                    str(i),
                    pool,
                    str(i),
                    "_".join([str(x) for x in seqrep_ids]),
                    "_".join([str(x) for x in sequence_replicate_numbers]),
                ]
            )

        expected_tsv_5 = make_expected_tsv_line(
            pipeline_root, 5, "1", [5], [1], "1.0.0"
        )
        expected_tsv_6_1 = make_expected_tsv_line(
            pipeline_root, 6, "0", [6], [1], "1.0.0"
        )
        expected_tsv_6_2 = make_expected_tsv_line(
            pipeline_root, 6, "0", [7], [2], "1.0.0"
        )
        expected_tsv_7_pool_1_0_0 = make_expected_tsv_line(
            pipeline_root, 7, "1", [8, 9], [1, 2], "1.0.0"
        )
        expected_tsv_7_pool_1_0_1 = make_expected_tsv_line(
            pipeline_root, 7, "1", [8, 9], [1, 2], "1.0.1"
        )

        tmp_out = "tmp.db_test.make_generic_pipeline_jobs_tsv.tsv"
        self.db.make_generic_pipeline_jobs_tsv(
            tmp_out, pipeline_root, pipeline_name, pipeline_version="1.0.0"
        )
        expected_lines = [
            expected_tsv_header,
            expected_tsv_5,
            expected_tsv_6_1,
            expected_tsv_6_2,
            expected_tsv_7_pool_1_0_0,
        ]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)

        #  Check the new pipline rows got added (with status=0 because they
        # are now "running")
        expected_rows.extend(
            [
                {
                    "isolate_id": 5,
                    "seqrep_id": None,
                    "seqrep_pool": "1",
                    "version": "1.0.0",
                    "pipeline_name": pipeline_name,
                    "status": 0,
                    "reference_id": None,
                },
                {
                    "isolate_id": 6,
                    "seqrep_id": 6,
                    "seqrep_pool": None,
                    "version": "1.0.0",
                    "pipeline_name": pipeline_name,
                    "status": 0,
                    "reference_id": None,
                },
                {
                    "isolate_id": 6,
                    "seqrep_id": 7,
                    "seqrep_pool": None,
                    "version": "1.0.0",
                    "pipeline_name": pipeline_name,
                    "status": 0,
                    "reference_id": None,
                },
                {
                    "isolate_id": 7,
                    "seqrep_id": None,
                    "seqrep_pool": "1_2",
                    "version": "1.0.0",
                    "pipeline_name": pipeline_name,
                    "status": 0,
                    "reference_id": None,
                },
            ]
        )

        def sort_pipeline_rows(rows):
            rows.sort(
                key=lambda x: (
                    x["isolate_id"],
                    -1 if x["seqrep_id"] is None else x["seqrep_id"],
                    "-1" if x["seqrep_pool"] is None else x["seqrep_pool"],
                )
            )

        sort_pipeline_rows(expected_rows)
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)

        # getting jobs for the same pipeline again should return nothing
        self.db.make_generic_pipeline_jobs_tsv(
            tmp_out, pipeline_root, pipeline_name, pipeline_version="1.0.0"
        )
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)
        self.check_tsv([expected_tsv_header], tmp_out)
        os.unlink(tmp_out)

        # different pipeline version should get more jobs, but only from group2
        self.db.make_generic_pipeline_jobs_tsv(
            tmp_out,
            pipeline_root,
            pipeline_name,
            pipeline_version="1.0.1",
            dataset_name="group2",
        )
        expected_lines = [expected_tsv_header, expected_tsv_7_pool_1_0_1]
        self.check_tsv(expected_lines, tmp_out)
        os.unlink(tmp_out)
        expected_rows.extend(
            [
                {
                    "isolate_id": 7,
                    "seqrep_id": None,
                    "seqrep_pool": "1_2",
                    "version": "1.0.1",
                    "pipeline_name": pipeline_name,
                    "status": 0,
                    "reference_id": None,
                }
            ]
        )

        sort_pipeline_rows(expected_rows)
        got_rows = self.db.get_rows_from_table("Pipeline")
        sort_pipeline_rows(got_rows)
        self.assertEqual(expected_rows, got_rows)

    def test_get_rows_from_table(self):
        """test get_rows_from_table"""
        db_schema.tables["Test"] = [
            ("col1", "integer"),
            ("col2", "integer"),
            ("col3", "text"),
        ]
        self.db.execute("""CREATE TABLE Test (col1 integer, col2 integer, col3 text)""")
        self.assertEqual([], self.db.get_rows_from_table("Test"))
        rows = [{"col1": 0, "col2": 1, "col3": "abcdef"}]
        self.db.add_row_to_table("Test", rows[0])
        self.assertEqual(rows, self.db.get_rows_from_table("Test"))
        self.assertEqual(
            [{"col1": 0, "col2": 1}],
            self.db.get_rows_from_table("Test", columns="col1, col2"),
        )

        rows.append({"col1": 1, "col2": 0, "col3": "covfefe"})
        self.db.add_row_to_table("Test", rows[-1])
        rows_cols_1_2 = [{"col1": 0, "col2": 1}, {"col1": 1, "col2": 0}]
        self.assertEqual(
            rows_cols_1_2, self.db.get_rows_from_table("Test", columns="col1, col2")
        )
        rows_cols_1_2.reverse()
        self.assertEqual(
            rows_cols_1_2,
            self.db.get_rows_from_table("Test", columns="col1, col2", order_by="col2"),
        )
        self.assertEqual(
            [rows[1]], self.db.get_rows_from_table("Test", where="col1 = 1")
        )

    def test_isolate_id_to_sample_id(self):
        """test isolate_id_to_sample_id"""
        self.assertEqual(None, self.db.isolate_id_to_sample_id(0))
        new_row = {
            "isolate_id": 42,
            "sample_id": 43,
            "isolate_number_from_lab": "foo",
            "pool_sequence_replicates": 0,
            "ena_experiment_accession": "ENA123",
        }
        self.db.add_row_to_table("Isolate", new_row)
        self.assertEqual(43, self.db.isolate_id_to_sample_id(42))

    def test_seqrep_id_to_sequence_replicate_number(self):
        """test seqrep_id_to_sequence_replicate_number"""
        self.assertEqual(None, self.db.seqrep_id_to_sequence_replicate_number(42))
        new_row = {
            "seqrep_id": 42,
            "isolate_id": 43,
            "sequence_replicate_number": 100,
            "original_reads_file_1_md5": "foo",
            "original_reads_file_2_md5": "bar",
            "remove_contam_reads_file_1_md5": "spam",
            "remove_contam_reads_file_2_md5": "eggs",
            "withdrawn": 0,
            "import_status": 1,
            "submission_date": datetime.date(2017, 12, 25),
            "instrument_model": "Illumina HiSeq 2500",
            "submit_to_ena": 0,
            "ena_run_accession": "ERR1",
            "ena_on_hold": 0,
        }
        self.db.add_row_to_table("Seqrep", new_row)
        self.assertEqual(100, self.db.seqrep_id_to_sequence_replicate_number(42))

    def test_get_vcfs_and_reads_files_for_minos_multi_sample_calling(self):
        """test get_vcfs_and_reads_files_for_minos_multi_sample_calling"""
        ref_id = 42
        dataset_name = "set1"
        pipeline_root = os.path.abspath(os.path.join("foo", "bar"))

        self.db.add_row_to_table(
            "Isolate",
            {
                "isolate_id": 1,
                "sample_id": 1,
                "isolate_number_from_lab": "iso1",
                "pool_sequence_replicates": 0,
                "ena_experiment_accession": "ENA1",
            },
        )
        self.db.add_row_to_table(
            "Pipeline",
            {
                "isolate_id": 1,
                "seqrep_id": None,
                "seqrep_pool": "1",
                "version": "1.2.3",
                "pipeline_name": "variant_call",
                "status": 1,
                "reference_id": 42,
            },
        )

        self.db.add_row_to_table(
            "Isolate",
            {
                "isolate_id": 2,
                "sample_id": 2,
                "isolate_number_from_lab": "iso2",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": "ENA2",
            },
        )
        self.db.add_row_to_table(
            "Pipeline",
            {
                "isolate_id": 2,
                "seqrep_id": None,
                "seqrep_pool": "1_2",
                "version": "1.2.3",
                "pipeline_name": "variant_call",
                "status": 1,
                "reference_id": 42,
            },
        )

        # should be ignored because wrong reference_id
        self.db.add_row_to_table(
            "Isolate",
            {
                "isolate_id": 3,
                "sample_id": 3,
                "isolate_number_from_lab": "iso3",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": "ENA3",
            },
        )
        self.db.add_row_to_table(
            "Pipeline",
            {
                "isolate_id": 3,
                "seqrep_id": None,
                "seqrep_pool": "1",
                "version": "1.2.3",
                "pipeline_name": "variant_call",
                "status": 1,
                "reference_id": 43,
            },
        )

        # should be ignored because wrong pipeline version
        self.db.add_row_to_table(
            "Isolate",
            {
                "isolate_id": 4,
                "sample_id": 4,
                "isolate_number_from_lab": "iso4",
                "pool_sequence_replicates": 0,
                "ena_experiment_accession": "ENA4",
            },
        )
        self.db.add_row_to_table(
            "Pipeline",
            {
                "isolate_id": 4,
                "seqrep_id": None,
                "seqrep_pool": "1",
                "version": "1.2.4",
                "pipeline_name": "variant_call",
                "status": 1,
                "reference_id": 42,
            },
        )
        got = self.db.get_vcfs_and_reads_files_for_minos_multi_sample_calling(
            dataset_name, pipeline_root, ref_id, pipeline_version="1.2.3"
        )
        vcf_iso1 = os.path.join(
            pipeline_root,
            "00",
            "00",
            "00",
            "01",
            "1",
            "Pipelines",
            "1",
            "variant_call",
            "1.2.3.ref.42",
            "minos",
            "final.vcf",
        )
        vcf_iso2 = os.path.join(
            pipeline_root,
            "00",
            "00",
            "00",
            "02",
            "2",
            "Pipelines",
            "1_2",
            "variant_call",
            "1.2.3.ref.42",
            "minos",
            "final.vcf",
        )
        bam1 = os.path.join(
            pipeline_root,
            "00",
            "00",
            "00",
            "01",
            "1",
            "Pipelines",
            "1",
            "variant_call",
            "1.2.3.ref.42",
            "samtools",
            "rmdup.bam",
        )
        bam2 = os.path.join(
            pipeline_root,
            "00",
            "00",
            "00",
            "02",
            "2",
            "Pipelines",
            "1_2",
            "variant_call",
            "1.2.3.ref.42",
            "samtools",
            "rmdup.bam",
        )
        expected = [vcf_iso1 + "\t" + bam1, vcf_iso2 + "\t" + bam2]
        self.assertEqual(expected, got)

    def test_update_remove_contam_stats(self):
        """test _update_remove_contam_stats"""
        counts_file = os.path.join(data_dir, "update_remove_contam_stats.counts.tsv")
        with self.assertRaises(db.Error):
            self.db._update_remove_contam_stats(1, counts_file)

        got_rows = self.db.get_rows_from_table("Read_counts")
        self.assertEqual(0, len(got_rows))
        seqrep_row = {
            "seqrep_id": 1,
            "isolate_id": 1,
            "sequence_replicate_number": 3,
            "original_reads_file_1_md5": "md51",
            "original_reads_file_2_md5": "md52",
            "remove_contam_reads_file_1_md5": "md53",
            "remove_contam_reads_file_2_md5": "md54",
            "withdrawn": 0,
            "import_status": 1,
            "instrument_model": "Illumina HiSeq 2500",
            "submission_date": "20170101",
            "submit_to_ena": 0,
            "ena_run_accession": 0,
            "ena_on_hold": 0,
        }
        self.db.add_row_to_table("Seqrep", seqrep_row)
        self.db._update_remove_contam_stats(1, counts_file)
        got_rows = self.db.get_rows_from_table("Read_counts")
        self.assertEqual(1, len(got_rows))
        expected = {
            "seqrep_id": 1,
            "original_total": 58,
            "contamination": 16,
            "not_contamination": 24,
            "unmapped": 18,
            "total_after_remove_contam": 44,
        }
        self.assertEqual(expected, got_rows[0])

        # should raise error because row already exists for this seqrep
        with self.assertRaises(db.Error):
            self.db._update_remove_contam_stats(1, counts_file)

    def test_update_qc_stats(self):
        """test _update_qc_stats"""
        seqrep_row = {
            "seqrep_id": 1,
            "isolate_id": 2,
            "sequence_replicate_number": 3,
            "original_reads_file_1_md5": "md51",
            "original_reads_file_2_md5": "md52",
            "remove_contam_reads_file_1_md5": "md53",
            "remove_contam_reads_file_2_md5": "md54",
            "withdrawn": 0,
            "import_status": 1,
            "instrument_model": "Illumina HiSeq 2500",
            "submission_date": "20170101",
            "submit_to_ena": 0,
            "ena_run_accession": 0,
            "ena_on_hold": 0,
        }
        # seqrep_row = (0, 1, 2, 3, 'md51', 'md52', 'mds53', 'md54', 0, 0, 1, '20170101', 0, 0, 0)
        self.db.add_row_to_table("Seqrep", seqrep_row)
        isolate_row = {
            "isolate_id": 2,
            "sample_id": 1,
            "isolate_number_from_lab": "i1",
            "pool_sequence_replicates": 1,
            "ena_experiment_accession": None,
        }
        self.db.add_row_to_table("Isolate", isolate_row)
        pipeline_root = os.path.join(data_dir, "update_qc_stats_pipeline_root")
        self.db._update_qc_stats(1, "1.0.0", pipeline_root)
        got_rows = self.db.get_rows_from_table("QC")
        expected_rows = [
            {
                "pipeline_version": "1.0.0",
                "seqrep_id": 1,
                "fastqc1_adapter_content": "pass",
                "fastqc1_basic_statistics": "pass",
                "fastqc1_gc": 49.0,
                "fastqc1_kmer_content": "pass",
                "fastqc1_max_sequence_length": 75,
                "fastqc1_min_sequence_length": 75,
                "fastqc1_overrepresented_sequences": "fail",
                "fastqc1_per_base_n_content": "pass",
                "fastqc1_per_base_sequence_content": "fail",
                "fastqc1_per_base_sequence_quality": "pass",
                "fastqc1_per_sequence_gc_content": "fail",
                "fastqc1_per_sequence_quality_scores": "fail",
                "fastqc1_sequence_duplication_levels": "pass",
                "fastqc1_sequence_length_distribution": "pass",
                "fastqc1_sequences_flagged_as_poor_quality": 0,
                "fastqc1_total_sequences": 6,
                "fastqc2_adapter_content": "pass",
                "fastqc2_basic_statistics": "pass",
                "fastqc2_gc": 50.0,
                "fastqc2_kmer_content": "pass",
                "fastqc2_max_sequence_length": 75,
                "fastqc2_min_sequence_length": 75,
                "fastqc2_overrepresented_sequences": "fail",
                "fastqc2_per_base_n_content": "pass",
                "fastqc2_per_base_sequence_content": "fail",
                "fastqc2_per_base_sequence_quality": "pass",
                "fastqc2_per_sequence_gc_content": "fail",
                "fastqc2_per_sequence_quality_scores": "fail",
                "fastqc2_sequence_duplication_levels": "pass",
                "fastqc2_sequence_length_distribution": "pass",
                "fastqc2_sequences_flagged_as_poor_quality": 0,
                "fastqc2_total_sequences": 6,
                "samtools_average_quality": 35.3,
                "samtools_bases_mapped_cigar": 4050,
                "samtools_bases_trimmed": 0,
                "samtools_error_rate": 0.00740741,
                "samtools_insert_size_average": 198.6,
                "samtools_insert_size_standard_deviation": 4.1,
                "samtools_inward_oriented_pairs": 27,
                "samtools_outward_oriented_pairs": 0,
                "samtools_pairs_with_other_orientation": 0,
                "samtools_positions_with_depth_of_0": 2,
                "samtools_positions_with_depth_atleast_2": 7,
                "samtools_positions_with_depth_atleast_5": 5,
                "samtools_positions_with_depth_atleast_10": 2,
                "samtools_positions_with_depth_atleast_20": 1,
                "samtools_positions_with_depth_atleast_100": 1,
                "samtools_raw_total_sequences": 54,
                "samtools_reads_duplicated": 2,
                "samtools_reads_mapped": 54,
                "het_snp_het_calls": 84,
                "het_snp_positions": 1673,
                "het_snp_total_snps": 200,
            }
        ]
        self.maxDiff = None
        self.assertEqual(expected_rows, got_rows)

    def test_update_finished_pipeline_run(self):
        """test update_finished_pipeline_run"""

        # Error because nothing found in Pipeline table
        with self.assertRaises(db.Error):
            self.db.update_finished_pipeline_run(1, 1, None, "pipeline_name", 1)

        # Error because status is 1, which means has already been run and updated
        pipeline_row = {
            "isolate_id": 1,
            "seqrep_id": 1,
            "seqrep_pool": None,
            "version": clockwork_version,
            "pipeline_name": "pipeline_name",
            "status": 1,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row)
        with self.assertRaises(db.Error):
            self.db.update_finished_pipeline_run(
                1, 1, None, "pipeline_name", 1, pipeline_version=clockwork_version
            )

        # Error because two rows found
        pipeline_row["seqrep_id"] = 2
        self.db.add_row_to_table("Pipeline", pipeline_row)
        self.db.add_row_to_table("Pipeline", pipeline_row)
        with self.assertRaises(db.Error):
            self.db.update_finished_pipeline_run(
                1, 2, None, "pipeline_name", 1, pipeline_version=clockwork_version
            )

        # Should run with no errors
        pipeline_row["seqrep_id"] = 3
        pipeline_row["status"] = 0
        self.db.add_row_to_table("Pipeline", pipeline_row)
        self.db.update_finished_pipeline_run(
            1, 3, None, "pipeline_name", 1, pipeline_version=clockwork_version
        )
        expected_rows = [
            {
                "isolate_id": 1,
                "seqrep_id": 3,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 1,
                "reference_id": 1,
            }
        ]
        got_rows = self.db.get_rows_from_table("Pipeline", where="seqrep_id = 3")
        self.assertEqual(expected_rows, got_rows)

        # test new_pipeline status -1
        pipeline_row["seqrep_id"] = 4
        pipeline_row["reference_id"] = 2
        self.db.add_row_to_table("Pipeline", pipeline_row)
        self.db.update_finished_pipeline_run(
            1,
            4,
            None,
            "pipeline_name",
            -1,
            pipeline_version=clockwork_version,
            reference_id=2,
        )
        expected_rows = [
            {
                "isolate_id": 1,
                "seqrep_id": 4,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": -1,
                "reference_id": 2,
            }
        ]
        got_rows = self.db.get_rows_from_table("Pipeline", where="seqrep_id = 4")
        self.assertEqual(expected_rows, got_rows)

        # Test when pooled
        pipeline_row = {
            "isolate_id": 2,
            "seqrep_id": None,
            "seqrep_pool": "42_43",
            "version": clockwork_version,
            "pipeline_name": "pipeline_name",
            "status": 0,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row)
        self.db.update_finished_pipeline_run(
            2,
            None,
            "42_43",
            "pipeline_name",
            1,
            pipeline_version=clockwork_version,
            reference_id=1,
        )
        expected_rows = [
            {
                "isolate_id": 2,
                "seqrep_id": None,
                "seqrep_pool": "42_43",
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 1,
                "reference_id": 1,
            }
        ]
        got_rows = self.db.get_rows_from_table("Pipeline", where="isolate_id = 2")
        self.assertEqual(expected_rows, got_rows)

    def test_update_finished_pipeline_run_remove_contam(self):
        """test update_finished_pipeline_run remove_contam"""
        # Test remove_contam pipeline. Should add md5s to Seqrep table
        pipeline_root = "tmp.test_update_finished_pipeline_run_remove_contam.root"
        if os.path.exists(pipeline_root):
            shutil.rmtree(pipeline_root)

        os.mkdir(pipeline_root)
        sequence_replicate_number = 2
        sample_dict = {
            "subject_id": "p1",
            "site_id": "s1",
            "lab_id": "l1",
            "isolate_number": "42",
            "sequence_replicate_number": sequence_replicate_number,
            "submission_date": datetime.date(2017, 12, 25),
            "reads_file_1": "reads_43_1.fq",
            "reads_file_1_md5": "md5_43_1",
            "reads_file_2_md5": "md5_43_2",
            "reads_file_2": "reads_43_2.fq",
            "dataset_name": "g1",
            "submit_to_ena": "0",
            "instrument_model": "Illumina HiSeq 2500",
            "ena_center_name": "Centre 42",
            "ena_on_hold": "0",
            "ena_run_accession": "ERR123456",
            "ena_sample_accession": "ERS123456",
        }
        (seqrep_id, isolate_id, sample_id) = self.db.add_one_seqrep(sample_dict)
        iso_dir = isolate_dir.IsolateDir(pipeline_root, sample_id, isolate_id)
        iso_dir.make_essential_dirs()
        reads1 = iso_dir.reads_filename("remove_contam", sequence_replicate_number, 1)
        reads2 = iso_dir.reads_filename("remove_contam", sequence_replicate_number, 2)
        f = pyfastaq.utils.open_file_write(reads1)
        print("foo", file=f)
        pyfastaq.utils.close(f)
        f = pyfastaq.utils.open_file_write(reads2)
        print("bar", file=f)
        pyfastaq.utils.close(f)
        expect_md5_1 = utils.md5(reads1)
        expect_md5_2 = utils.md5(reads2)

        pipeline_row = {
            "isolate_id": isolate_id,
            "seqrep_id": seqrep_id,
            "seqrep_pool": None,
            "version": clockwork_version,
            "pipeline_name": "remove_contam",
            "status": 0,
            "reference_id": 1,
        }
        self.db.add_row_to_table("Pipeline", pipeline_row)

        to_copy_counts_file = os.path.join(
            data_dir, "update_finished_pipeline_run_remove_contam", "counts.tsv"
        )
        counts_file = iso_dir.contamination_counts_filename(sequence_replicate_number)
        shutil.copyfile(to_copy_counts_file, counts_file)

        self.db.update_finished_pipeline_run(
            isolate_id,
            seqrep_id,
            None,
            "remove_contam",
            1,
            reference_id=1,
            pipeline_version=clockwork_version,
            pipeline_root=pipeline_root,
        )
        got_rows = self.db.get_rows_from_table(
            "Seqrep", where="seqrep_id =" + str(seqrep_id)
        )
        self.assertEqual(1, len(got_rows))
        self.assertEqual(expect_md5_1, got_rows[0]["remove_contam_reads_file_1_md5"])
        self.assertEqual(expect_md5_2, got_rows[0]["remove_contam_reads_file_2_md5"])
        shutil.rmtree(pipeline_root)

    def test_load_success_jobs_file(self):
        """test _load_success_jobs_file"""
        expected = {(1, "3", "1", False), (2, "4", "42", False)}
        got = db.Db._load_success_jobs_file(
            os.path.join(data_dir, "load_success_jobs_file.no_pool_column.tsv")
        )
        self.assertEqual(expected, got)

        expected = {
            (1, "3", "1", False),
            (2, "4_5", "42_43", True),
            (3, "6", "2", True),
        }
        got = db.Db._load_success_jobs_file(
            os.path.join(data_dir, "load_success_jobs_file.with_pool_column.tsv")
        )
        self.assertEqual(expected, got)

    def test_update_finished_pipeline_run_failed_jobs_not_pooled(self):
        """test update_finished_pipeline_run_failed_jobs not pooled"""
        for i in [1, 2, 3]:
            pipeline_row = {
                "isolate_id": i,
                "seqrep_id": i,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 0,
                "reference_id": 1,
            }
            self.db.add_row_to_table("Pipeline", pipeline_row)

        success_tsv = os.path.join(
            data_dir, "update_finished_pipeline_run_failed_jobs.success.tsv"
        )
        jobs_tsv = os.path.join(
            data_dir, "update_finished_pipeline_run_failed_jobs.jobs.tsv"
        )
        self.db.update_finished_pipeline_run_failed_jobs(
            jobs_tsv, success_tsv, "pipeline_name"
        )
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("seqrep_id"))
        expected_rows = [
            {
                "isolate_id": 1,
                "seqrep_id": 1,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 0,
                "reference_id": 1,
            },
            {
                "isolate_id": 2,
                "seqrep_id": 2,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": -1,
                "reference_id": 1,
            },
            {
                "isolate_id": 3,
                "seqrep_id": 3,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": -1,
                "reference_id": 1,
            },
        ]
        self.assertEqual(expected_rows, got_rows)

    def test_update_finished_pipeline_run_failed_jobs_pooled(self):
        """test update_finished_pipeline_run_failed_jobs pooled"""
        pipeline_rows = [
            {
                "isolate_id": 1,
                "seqrep_id": 1,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 0,
                "reference_id": 1,
            },
            {
                "isolate_id": 1,
                "seqrep_id": 2,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 0,
                "reference_id": 1,
            },
            {
                "isolate_id": 2,
                "seqrep_id": None,
                "seqrep_pool": "1",
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 0,
                "reference_id": 1,
            },
            {
                "isolate_id": 3,
                "seqrep_id": None,
                "seqrep_pool": "42_43",
                "version": clockwork_version,
                "pipeline_name": "pipeline_name",
                "status": 0,
                "reference_id": 1,
            },
        ]
        for row in pipeline_rows:
            self.db.add_row_to_table("Pipeline", row)

        got = self.db.get_rows_from_table("Pipeline")

        success_tsv = os.path.join(
            data_dir, "update_finished_pipeline_run_failed_jobs_pooled.success.tsv"
        )
        jobs_tsv = os.path.join(
            data_dir, "update_finished_pipeline_run_failed_jobs_pooled.jobs.tsv"
        )
        self.db.update_finished_pipeline_run_failed_jobs(
            jobs_tsv, success_tsv, "pipeline_name", reference_id=1
        )
        got_rows = self.db.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("isolate_id"))
        for i in [0, 1, 3]:
            pipeline_rows[i]["status"] = -1
        self.assertEqual(pipeline_rows, got_rows)

    def test_add_reference(self):
        """test add_reference"""
        got_rows = self.db.get_rows_from_table("Reference")
        self.assertEqual([], got_rows)
        got_id = self.db.add_reference("ref_name")
        self.assertEqual(1, got_id)
        got_rows = self.db.get_rows_from_table("Reference")
        expected_rows = [{"reference_id": 1, "name": "ref_name"}]
        self.assertEqual(expected_rows, got_rows)
        with self.assertRaises(db.Error):
            self.db.add_reference("ref_name")
        got_id = self.db.add_reference("ref_name2")
        self.assertEqual(2, got_id)
        expected_rows.append({"reference_id": 2, "name": "ref_name2"})
        got_rows = self.db.get_rows_from_table("Reference")
        got_rows.sort(key=itemgetter("reference_id"))
        self.assertEqual(expected_rows, got_rows)

    def test_add_mykrobe_custom_panel(self):
        """test add_mykrobe_custom_panel"""
        tmp_probes = "tmp.test_add_mykrobe_custom_panel.probes.fa"
        tmp_json = "tmp.test_add_mykrobe_custom_panel.json"
        with open(tmp_probes, "w"):
            pass
        with open(tmp_json, "w"):
            pass
        ref_root = "tmp.test_add_mykrobe_custom_panel.refs"
        os.mkdir(ref_root)
        name1 = "mykrobe_test"
        species = "tb"
        ref_id = self.db.add_mykrobe_custom_panel(
            species, name1, ref_root, probes_fasta=tmp_probes, var_to_res_json=tmp_json
        )
        self.assertEqual(1, ref_id)
        got_rows = self.db.get_rows_from_table("Reference")
        expected_rows = [{"reference_id": 1, "name": name1}]
        panel_dir = self.db.get_reference_dir(ref_id, ref_root)
        panel = mykrobe.Panel(panel_dir.directory)
        self.assertTrue(os.path.exists(panel.probes_fasta))
        self.assertTrue(os.path.exists(panel.var_to_res_json))
        self.assertEqual(species, panel.metadata["species"])
        self.assertEqual(name1, panel.metadata["name"])
        self.assertFalse(panel.metadata["is_built_in"])
        with self.assertRaises(db.Error):
            self.db.add_mykrobe_custom_panel(
                species, name1, ref_root, tmp_probes, tmp_json
            )
        os.unlink(tmp_probes)
        os.unlink(tmp_json)

        # Test adding a panel that is built-in to mykrobe
        name2 = "walker-2015"
        ref_id = self.db.add_mykrobe_custom_panel(species, name2, ref_root)
        self.assertEqual(2, ref_id)
        got_rows = self.db.get_rows_from_table("Reference")
        expected_rows = [
            {"reference_id": 1, "name": name1},
            {"reference_id": 2, "name": name2},
        ]
        panel_dir = self.db.get_reference_dir(ref_id, ref_root)
        panel = mykrobe.Panel(panel_dir.directory)
        self.assertFalse(os.path.exists(panel.probes_fasta))
        self.assertFalse(os.path.exists(panel.var_to_res_json))
        self.assertEqual(species, panel.metadata["species"])
        self.assertEqual(name2, panel.metadata["name"])
        self.assertTrue(panel.metadata["is_built_in"])
        shutil.rmtree(ref_root)

    def test_has_reference(self):
        """test has_reference"""
        self.assertFalse(self.db.has_reference(1))
        self.db.add_reference("kryten")
        self.assertTrue(self.db.has_reference(1))

    def test_get_reference_dir(self):
        """test get_reference_dir"""
        root = "/root/"
        with self.assertRaises(db.Error):
            self.db.get_reference_dir(1, root)
        self.db.add_reference("refname")
        got = self.db.get_reference_dir(1, root)
        expected = reference_dir.ReferenceDir(
            pipeline_references_root_dir=root, reference_id=1
        )
        self.assertEqual(expected, got)

    def test_backup(self):
        """test backup"""
        tmpfile = "tmp.db.backup"
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)
        self.db.backup(outfile=tmpfile)
        # trust the mysqldump command, just cgeck backup file got made
        self.assertTrue(os.path.exists(tmpfile))
        with self.assertRaises(db.Error):
            self.db.backup(outfile=tmpfile)
        os.unlink(tmpfile)

        tmpdir = "tmp.db.backup.dir"
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)

        os.mkdir(tmpdir)
        self.db.backup(backup_dir=tmpdir)
        backup_files = os.listdir(tmpdir)
        self.assertEqual(1, len(backup_files))
        backup_file = os.path.join(tmpdir, backup_files[0])
        shutil.rmtree(tmpdir)
