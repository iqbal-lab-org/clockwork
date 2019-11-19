import unittest
import copy
import datetime
import shutil
import filecmp
import os
import re
from clockwork import (
    spreadsheet_importer,
    db,
    db_connection,
    db_maker,
    lock_file,
    utils,
)

modules_dir = os.path.dirname(os.path.abspath(spreadsheet_importer.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "spreadsheet_importer")
db_ini_file = os.path.join(modules_dir, "tests", "data", "db.ini")


class TestSpreadsheetImporter(unittest.TestCase):
    def setUp(self):
        try:
            db_connection.DbConnection(db_ini_file, destroy=True)
        except:
            pass

        dbm = db_maker.DbMaker(db_ini_file)
        dbm.run()
        self.db = db.Db(db_ini_file)

    def tearDown(self):
        self.db.commit_and_close()
        db_connection.DbConnection(db_ini_file, destroy=True, must_exist=True)

    def test_create_template_spreadsheet(self):
        """Test create_template_spreadsheet"""
        filename = "tmp.create_template_spreadsheet.xlsx"
        if os.path.exists(filename):
            os.unlink(filename)
        spreadsheet_importer.create_template_spreadsheet(filename)
        self.assertTrue(os.path.exists(filename))
        os.unlink(filename)

    def test_validate_data(self):
        """Test validate data"""
        dropbox_dir = "tmp.spreadsheet_importer.validate_data.dropbox"
        os.mkdir(dropbox_dir)

        data = [
            {
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
                "instrument_model": "Illumina HiSeq 2000",
                "ena_center_name": "Center 1",
                "submit_to_ena": "0",
                "ena_on_hold": "0",
                "ena_run_accession": "ERR123456",
                "ena_sample_accession": "ERS123456",
            }
        ]

        expected = []

        utils.make_empty_file(os.path.join(dropbox_dir, data[0]["reads_file_1"]))
        utils.make_empty_file(os.path.join(dropbox_dir, data[0]["reads_file_2"]))

        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )

        # check error when bas date format
        data[0]["submission_date"] = "oops"
        expected = [
            "Date format error: p1, s1, l1, 42, 43, oops, reads_1_1.fq, abcdefghijklmnopqrstuvwyx123456, reads_1_2.fq, abcdefghijklmnopqrstuvwyx123457, g1, Illumina HiSeq 2000, Center 1, 0, 0, ERR123456, ERS123456"
        ]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        data[0]["submission_date"] = datetime.date(2017, 12, 25)

        # check error when files not found
        os.unlink(os.path.join(dropbox_dir, data[0]["reads_file_1"]))
        os.unlink(os.path.join(dropbox_dir, data[0]["reads_file_2"]))
        expected = [
            "Reads file not found: " + data[0]["reads_file_1"],
            "Reads file not found: " + data[0]["reads_file_2"],
        ]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        utils.make_empty_file(os.path.join(dropbox_dir, data[0]["reads_file_1"]))
        utils.make_empty_file(os.path.join(dropbox_dir, data[0]["reads_file_2"]))

        # check error when non-unique replicate
        new_data = copy.copy(data[0])
        new_data["reads_file_1"] = "reads_2_1.fq"
        new_data["reads_file_1_md5"] = "12345"
        new_data["reads_file_2"] = "reads_2_2.fq"
        new_data["reads_file_2_md5"] = "12346"
        data.append(new_data)
        utils.make_empty_file(os.path.join(dropbox_dir, data[1]["reads_file_1"]))
        utils.make_empty_file(os.path.join(dropbox_dir, data[1]["reads_file_2"]))
        expected = [
            "Replicate subject_id,site_id,lab_id,isolate_number,sequence_replicate_number p1,s1,l1,42,43 found 2 times in spreadsheet"
        ]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        del data[1]

        # check error when duplicate file name
        data[0]["reads_file_1"] = data[0]["reads_file_2"]
        expected = ["Reads file " + data[0]["reads_file_2"] + " found 2 times"]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        data[0]["reads_file_1"] = "reads_1_1.fq"

        # check error when no md5
        data[0]["reads_file_1_md5"] = None
        expected = ["No md5 for reads file " + data[0]["reads_file_1"]]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        data[0]["reads_file_1_md5"] = "abcdefghijklmnopqrstuvwyx123456"

        # check error when md5 mismatch
        md5_file = os.path.join(dropbox_dir, "reads_1_1.fq.md5")
        with open(md5_file, "w") as f:
            print("2ead32a9cbca30ca4b3b9acf59852966  reads_1_1.fq", file=f)
        expected = ["Mismatch in md5 info for reads file " + data[0]["reads_file_1"]]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        data[0]["reads_file_1_md5"] = "abcdefghijklmnopqrstuvwyx123456"
        os.unlink(md5_file)

        # check replicates unique
        self.assertEqual(
            [],
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        self.db.add_one_seqrep(data[0])
        expected = [
            "Replicate already found for subject_id,site_id,lab_id,isolate_number,sequence_replicate_number: p1,s1,l1,42,43"
        ]
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )

        # check patient+site+lab unique
        sample_row = {
            "sample_id": None,
            "subject_id": "p1",
            "site_id": "s1",
            "sample_id_from_lab": "l1",
            "dataset_name": "g1",
            "ena_center_name": "Center 1",
            "ena_sample_accession": "ERS123456",
            "ena_study_accession": "ERP123456",
        }
        self.db.add_row_to_table("Sample", sample_row)
        expected = [
            "Subject(p1) + site(s1) + lab(l1) found more than once in database. Something very wrong!"
        ] + expected
        self.assertEqual(
            expected,
            spreadsheet_importer.SpreadsheetImporter._validate_data(
                self.db, data, dropbox_dir
            ),
        )
        shutil.rmtree(dropbox_dir)

    def test_archive_spreadsheet(self):
        """test _archive_spreadsheet"""
        archive_dir = "tmp.test_archive_spreadsheet_archive"
        os.mkdir(archive_dir)
        xlsx_file = "tmp.test_archive_spreadsheet.xlsx"
        xlsx_file2 = "tmp.test_archive_spreadsheet2.xlsx"
        with open(xlsx_file, "w"):
            pass

        date = utils.date_string_from_file_mtime(xlsx_file)
        date_dir = os.path.join(archive_dir, date)
        got_filename = spreadsheet_importer.SpreadsheetImporter._archive_spreadsheet(
            xlsx_file, archive_dir
        )
        expected_filename = os.path.join(date_dir, xlsx_file)
        self.assertEqual(expected_filename, got_filename)
        self.assertFalse(os.path.exists(xlsx_file))
        self.assertTrue(os.path.exists(expected_filename))

        with open(xlsx_file, "w"):
            pass

        got_filename = spreadsheet_importer.SpreadsheetImporter._archive_spreadsheet(
            xlsx_file, archive_dir
        )
        expected_filename = os.path.join(date_dir, xlsx_file + ".1")
        self.assertEqual(expected_filename, got_filename)
        self.assertFalse(os.path.exists(xlsx_file))
        self.assertTrue(os.path.exists(expected_filename))

        with open(xlsx_file2, "w"):
            pass

        got_filename = spreadsheet_importer.SpreadsheetImporter._archive_spreadsheet(
            xlsx_file2, archive_dir
        )
        expected_filename = os.path.join(date_dir, xlsx_file2)
        self.assertEqual(expected_filename, got_filename)
        self.assertFalse(os.path.exists(xlsx_file2))
        self.assertTrue(os.path.exists(expected_filename))

        shutil.rmtree(archive_dir)

    def test_import_reads_and_update_db(self):
        """test _import_reads_and_update_db"""
        archive_dir = "tmp.test_spreadsheet_importer_run.archive"
        os.mkdir(archive_dir)
        db_backup_dir = "tmp.test_spreadsheet_importer_run.db_backup"
        os.mkdir(db_backup_dir)
        tsv_file = "tmp.test_spreadsheet_importer_run.out.tsv"

        # need to copy the dropbox directory, because run() will move the
        # xlsx file
        original_dropbox_dir = os.path.join(data_dir, "run.dropbox")
        dropbox_dir = "tmp.test_spreadsheet_importer_run.out.dropbox"
        shutil.copytree(original_dropbox_dir, dropbox_dir)
        xlsx_file = os.path.join(dropbox_dir, "import.xlsx")
        date = utils.date_string_from_file_mtime(xlsx_file)

        importer = spreadsheet_importer.SpreadsheetImporter(
            dropbox_dir,
            xlsx_file,
            db_ini_file,
            archive_dir,
            tsv_file,
            db_backup_dir=db_backup_dir,
        )
        importer._import_reads_and_update_db()

        # check xlsx file and import jobs file got archived
        xlsx_archive_file = os.path.join(archive_dir, date, "import.xlsx")
        self.assertTrue(os.path.exists(xlsx_archive_file))
        self.assertTrue(os.path.exists(xlsx_archive_file + ".import_jobs.tsv"))
        self.assertFalse(os.path.exists(xlsx_file))
        shutil.rmtree(archive_dir)

        # check database got backed up
        backup_files = os.listdir(db_backup_dir)
        self.assertEqual(1, len(backup_files))
        shutil.rmtree(db_backup_dir)

        # check tsv file is correct
        reads_prefix = os.path.abspath(os.path.join(dropbox_dir, "reads"))
        expected_tsv_lines = [
            "\t".join(
                [
                    "seqrep_id",
                    "sample_id",
                    "isolate_id",
                    "sequence_replicate_number",
                    "reads1",
                    "reads2",
                    "reads1_md5",
                    "reads2_md5",
                ]
            ),
            "\t".join(
                [
                    "1",
                    "1",
                    "1",
                    "43",
                    reads_prefix + ".1_1.fq",
                    reads_prefix + ".1_2.fq",
                    "abcdefghijklmnopqrstuvwyx123456",
                    "abcdefghijklmnopqrstuvwyx123457",
                ]
            ),
            "\t".join(
                [
                    "2",
                    "2",
                    "2",
                    "45",
                    reads_prefix + ".2_1.fq",
                    reads_prefix + ".2_2.fq",
                    "a73817805eb1d44ca88eb5cb794c7de7",
                    "d468360c689d482b227256d887a05996",
                ]
            ),
        ]
        with open(tsv_file) as f:
            got_tsv_lines = [line.rstrip() for line in f]

        self.assertEqual(expected_tsv_lines, got_tsv_lines)
        os.unlink(tsv_file)

        # check database is correct
        got_sample = self.db.get_rows_from_table("Sample", order_by="sample_id")
        expected_sample = [
            {
                "sample_id": 1,
                "subject_id": "p1",
                "site_id": "s1",
                "sample_id_from_lab": "l1",
                "dataset_name": "g1",
                "ena_center_name": "Center 1",
                "ena_sample_accession": "ERS123456",
                "ena_study_accession": None,
            },
            {
                "sample_id": 2,
                "subject_id": "p2",
                "site_id": "s2",
                "sample_id_from_lab": "l2",
                "dataset_name": "g2",
                "ena_center_name": "Center 1",
                "ena_sample_accession": None,
                "ena_study_accession": None,
            },
        ]
        self.maxDiff = None
        self.assertEqual(expected_sample, got_sample)

        got_seqrep = self.db.get_rows_from_table("Seqrep", order_by="seqrep_id")
        expected_seqrep = [
            {
                "seqrep_id": 1,
                "isolate_id": 1,
                "sequence_replicate_number": 43,
                "original_reads_file_1_md5": "abcdefghijklmnopqrstuvwyx123456",
                "original_reads_file_2_md5": "abcdefghijklmnopqrstuvwyx123457",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "withdrawn": 0,
                "import_status": 0,
                "submission_date": datetime.date(2017, 12, 25),
                "instrument_model": "Illumina HiSeq 2000",
                "submit_to_ena": 0,
                "ena_run_accession": "ERR123456",
                "ena_on_hold": 0,
            },
            {
                "seqrep_id": 2,
                "isolate_id": 2,
                "sequence_replicate_number": 45,
                "original_reads_file_1_md5": "a73817805eb1d44ca88eb5cb794c7de7",
                "original_reads_file_2_md5": "d468360c689d482b227256d887a05996",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "withdrawn": 0,
                "import_status": 0,
                "submission_date": datetime.date(2017, 12, 26),
                "instrument_model": "Illumina HiSeq 2000",
                "submit_to_ena": 1,
                "ena_run_accession": None,
                "ena_on_hold": 1,
            },
        ]

        self.assertEqual(expected_seqrep, got_seqrep)
        shutil.rmtree(dropbox_dir)

    def test_run(self):
        """test run"""
        archive_dir = "tmp.test_spreadsheet_importer_run.archive"
        os.mkdir(archive_dir)
        tsv_file = "tmp.test_spreadsheet_importer_run.out.tsv"

        # need to copy the dropbox directory, because run() will move the
        # xlsx file
        original_dropbox_dir = os.path.join(data_dir, "run.dropbox")
        dropbox_dir = "tmp.test_spreadsheet_importer_run.out.dropbox"
        shutil.copytree(original_dropbox_dir, dropbox_dir)
        xlsx_file = os.path.join(dropbox_dir, "import.xlsx")
        date = utils.date_string_from_file_mtime(xlsx_file)

        importer = spreadsheet_importer.SpreadsheetImporter(
            dropbox_dir, xlsx_file, db_ini_file, archive_dir, tsv_file
        )

        # test lock file stops it running
        utils.make_empty_file(importer.lock_file)

        with self.assertRaises(lock_file.Error):
            importer.run()

        os.unlink(importer.lock_file)

        # we'll just run it - the details are checked in test_import_reads_and_update_db
        importer.run()
        shutil.rmtree(archive_dir)
        shutil.rmtree(dropbox_dir)
        os.unlink(tsv_file)
