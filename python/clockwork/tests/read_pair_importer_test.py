import unittest
import os
import shutil
from clockwork import (
    db,
    db_connection,
    db_maker,
    isolate_dir,
    read_pair_importer,
    utils,
    fqtools,
)

modules_dir = os.path.dirname(os.path.abspath(read_pair_importer.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "read_pair_importer")
ini_file = os.path.join(modules_dir, "tests", "data", "db.ini")


class TestReadPairImporter(unittest.TestCase):
    def setUp(self):
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

    def test_copy_reads_file(self):
        """test _copy_reads_file"""
        original_reads_file = os.path.join(data_dir, "copy_reads_file.fq")
        correct_md5 = "1bf5fb5b641e793ba40285779a059132"
        tmp_reads_file = "tmp.read_pair_importer.copy_reads_file.out"
        if os.path.exists(tmp_reads_file):
            os.unlink(tmp_reads_file)

        read_pair_importer.ReadPairImporter._copy_reads_file(
            original_reads_file, tmp_reads_file, correct_md5
        )
        self.assertTrue(os.path.exists(tmp_reads_file))
        self.assertEqual(correct_md5, utils.md5(tmp_reads_file))

        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._copy_reads_file(
                original_reads_file, tmp_reads_file, "wrong md5"
            )

        os.unlink(tmp_reads_file)

    def test_check_database(self):
        """test _check_database"""
        seqrep_row = {
            "seqrep_id": 1,
            "isolate_id": 1,
            "sequence_replicate_number": 1,
            "original_reads_file_1_md5": "md51",
            "original_reads_file_2_md5": "md52",
            "remove_contam_reads_file_1_md5": None,
            "remove_contam_reads_file_2_md5": None,
            "withdrawn": 0,
            "import_status": 0,
            "instrument_model": "Illumina HiSeq 2000",
            "submission_date": "20170101",
            "submit_to_ena": 0,
            "ena_run_accession": None,
            "ena_on_hold": 0,
        }
        self.db.add_row_to_table("Seqrep", seqrep_row)
        read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 1)

        # wrong isolate_id
        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._check_database(self.db, 1, 2, 1)

        # wrong sequence_replicate_number
        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 2)

        # unknown import status
        self.db.update_row("Seqrep", {"seqrep_id": 1}, {"import_status": 2})
        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 1)

        # already imported (status 1)
        self.db.update_row("Seqrep", {"seqrep_id": 1}, {"import_status": 1})
        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 1)

        # failed import (status -1)
        self.db.update_row("Seqrep", {"seqrep_id": 1}, {"import_status": -1})
        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 1)

        # should be ok
        self.db.update_row("Seqrep", {"seqrep_id": 1}, {"import_status": 0})
        read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 1)

    def test_update_database(self):
        """test _update_database"""
        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._update_database(self.db, 1, 1, 3)

        seqrep_row = {
            "seqrep_id": None,
            "isolate_id": 1,
            "sequence_replicate_number": 3,
            "original_reads_file_1_md5": "md51",
            "original_reads_file_2_md5": "md52",
            "remove_contam_reads_file_1_md5": None,
            "remove_contam_reads_file_2_md5": None,
            "withdrawn": 0,
            "import_status": 0,
            "instrument_model": "Illumina HiSeq 2000",
            "submission_date": "20170101",
            "submit_to_ena": 0,
            "ena_run_accession": None,
            "ena_on_hold": 0,
        }
        self.db.add_row_to_table("Seqrep", seqrep_row)
        read_pair_importer.ReadPairImporter._check_database(self.db, 1, 1, 3)

        with self.assertRaises(read_pair_importer.Error):
            read_pair_importer.ReadPairImporter._update_database(self.db, 1, 1, 4)

        read_pair_importer.ReadPairImporter._update_database(
            self.db, 1, 1, 3, import_status=42
        )
        rows = self.db.get_rows_from_table("Seqrep", where="seqrep_id=1")
        self.assertEqual(1, len(rows))
        self.assertEqual(42, rows[0]["import_status"])

    def _test_run(self, original_reads1, original_reads2, expected_import_status):
        """test run"""
        pipeline_root = "tmp.read_pair_importer.run.root"
        if os.path.exists(pipeline_root):
            shutil.rmtree(pipeline_root)
        os.mkdir(pipeline_root)
        seqrep_id = 1
        sample = 3
        isolate = 2
        sequence_replicate_number = 42

        # copy the reads because the pipeline will delete them later
        reads1 = "tmp.read_pair_importer.reads1.fq"
        reads2 = "tmp.read_pair_importer.reads2.fq"
        md5_1 = utils.rsync_and_md5(original_reads1, reads1)
        md5_2 = utils.rsync_and_md5(original_reads2, reads2)

        # write an md5 file, to check it gets deleted later
        md5_file = reads1 + ".md5"
        utils.syscall("md5sum " + reads1 + " > " + md5_file)

        importer = read_pair_importer.ReadPairImporter(
            ini_file,
            pipeline_root,
            seqrep_id,
            isolate,
            sample,
            sequence_replicate_number,
            reads1,
            reads2,
            md5_1,
            md5_2,
        )

        # no row in Seqrep table
        with self.assertRaises(read_pair_importer.Error):
            importer.run()

        seqrep_row = {
            "seqrep_id": seqrep_id,
            "isolate_id": isolate,
            "sequence_replicate_number": sequence_replicate_number,
            "original_reads_file_1_md5": md5_1,
            "original_reads_file_2_md5": md5_2,
            "remove_contam_reads_file_1_md5": None,
            "remove_contam_reads_file_2_md5": None,
            "withdrawn": 0,
            "import_status": 0,
            "instrument_model": "Illumina HiSeq 2000",
            "submission_date": "20170101",
            "submit_to_ena": 0,
            "ena_run_accession": None,
            "ena_on_hold": 0,
        }
        self.db.add_row_to_table("Seqrep", seqrep_row)
        self.db.commit()

        # need to create a new object so the database changes get picked up.
        # Separate connections do not see changes made by each other.
        importer = read_pair_importer.ReadPairImporter(
            ini_file,
            pipeline_root,
            seqrep_id,
            isolate,
            sample,
            sequence_replicate_number,
            reads1,
            reads2,
            md5_1,
            md5_2,
        )

        # reads file doesn't exit
        importer.reads_file_1 = "oops"
        with self.assertRaises(read_pair_importer.Error):
            importer.run()
        importer.reads_file_1 = reads1

        # check lock file works
        iso_dir = isolate_dir.IsolateDir(pipeline_root, sample, isolate)
        iso_dir.make_essential_dirs()
        lock_file = os.path.join(iso_dir.reads_dir, "import_lock." + str(seqrep_id))
        utils.make_empty_file(lock_file)
        with self.assertRaises(read_pair_importer.Error):
            importer.run()
        os.unlink(lock_file)

        # should run ok
        where_query = " and ".join(
            [
                "sample_id=" + str(sample),
                "isolate_number=" + str(isolate),
                "sequence_replicate_number=" + str(seqrep_id),
            ]
        )
        rows = self.db.get_rows_from_table(
            "Seqrep", where="seqrep_id=" + str(seqrep_id)
        )
        self.assertEqual(1, len(rows))
        self.assertEqual(0, rows[0]["import_status"])

        importer = read_pair_importer.ReadPairImporter(
            ini_file,
            pipeline_root,
            seqrep_id,
            isolate,
            sample,
            sequence_replicate_number,
            reads1,
            reads2,
            md5_1,
            md5_2,
        )
        importer.run()
        # reconnect so that we pick up the changes made by the previous line
        self.db.reconnect()
        reads_out_1 = iso_dir.reads_filename("original", sequence_replicate_number, 1)
        reads_out_2 = iso_dir.reads_filename("original", sequence_replicate_number, 2)
        rows = self.db.get_rows_from_table(
            "Seqrep", where="seqrep_id=" + str(seqrep_id)
        )
        self.assertEqual(1, len(rows))
        self.assertEqual(expected_import_status, rows[0]["import_status"])

        # Files either copied/deleted or not depending on if import was successful
        if expected_import_status == 1:
            self.assertTrue(os.path.exists(reads_out_1))
            self.assertTrue(os.path.exists(reads_out_2))
            self.assertFalse(os.path.exists(reads1))
            self.assertFalse(os.path.exists(reads2))
            self.assertFalse(os.path.exists(md5_file))
        else:
            self.assertFalse(os.path.exists(reads_out_1))
            self.assertFalse(os.path.exists(reads_out_2))
            self.assertTrue(os.path.exists(reads1))
            self.assertTrue(os.path.exists(reads2))
            self.assertTrue(os.path.exists(md5_file))
            os.unlink(reads1)
            os.unlink(reads2)
            os.unlink(md5_file)

        shutil.rmtree(pipeline_root)

    def test_run_good_reads(self):
        reads1 = os.path.join(data_dir, "run.reads.1.fq")
        reads2 = os.path.join(data_dir, "run.reads.2.fq")
        self._test_run(reads1, reads2, 1)

    def test_run_bad_reads(self):
        reads1 = os.path.join(data_dir, "run.bad_reads.1.fq")
        reads2 = os.path.join(data_dir, "run.bad_reads.2.fq")
        self._test_run(reads1, reads2, -1)
