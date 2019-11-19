import unittest
import datetime
import glob
import sys
import shutil
import os
from operator import itemgetter

from clockwork import db, db_connection, db_maker, isolate_dir, utils

sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
import nextflow_helper

data_dir = os.path.join(nextflow_helper.data_root_dir, "nextflow_import")
modules_dir = os.path.dirname(os.path.abspath(db.__file__))
db_ini_file = os.path.join(modules_dir, "tests", "data", "db.ini")


class TestNextflowImport(unittest.TestCase):
    def test_nextflow_import(self):
        """test nextflow_import"""
        nextflow_helper.write_config_file()
        pipeline_root = "tmp.nextflow_import.pipeline_root"
        os.mkdir(pipeline_root)
        try:
            db_connection.DbConnection(db_ini_file, destroy=True)
        except:
            pass

        dbm = db_maker.DbMaker(db_ini_file)
        dbm.run()

        dropbox_dir = "tmp.nextflow_import.dropbox"
        shutil.copytree(os.path.join(data_dir, "dropbox"), dropbox_dir)
        xlsx_archive_dir = "tmp.nextflow_import.xlsx_archive"
        os.mkdir(xlsx_archive_dir)
        expected_xlsx_files = [
            os.path.basename(x) for x in glob.glob(os.path.join(dropbox_dir, "*.xlsx"))
        ]

        nextflow_file = os.path.join(nextflow_helper.nextflow_dir, "import.nf")
        work_dir = "tmp.nextflow_import.work"
        dag_file = "nextflow.import.dag.pdf"
        try:
            os.unlink(dag_file)
        except:
            pass

        command = " ".join(
            [
                "nextflow run",
                "--dropbox_dir",
                dropbox_dir,
                "--pipeline_root",
                pipeline_root,
                "--db_config_file",
                db_ini_file,
                "--xlsx_archive_dir",
                xlsx_archive_dir,
                "-with-dag",
                dag_file,
                "-c",
                nextflow_helper.config_file,
                "-w",
                work_dir,
                nextflow_file,
            ]
        )
        utils.syscall(command)
        os.unlink(nextflow_helper.config_file)
        shutil.rmtree(work_dir)

        # All files should be gone from the dropbox
        self.assertEqual([], os.listdir(dropbox_dir))
        shutil.rmtree(dropbox_dir)

        # The two spreadsheets should have been archived
        got_xlsx_files = [
            os.path.basename(x)
            for x in glob.glob(os.path.join(xlsx_archive_dir, "**", "*.xlsx"))
        ]
        self.assertEqual(expected_xlsx_files, got_xlsx_files)
        shutil.rmtree(xlsx_archive_dir)

        # Check database updated correctly
        database = db.Db(db_ini_file)
        expected_sample_rows = [
            {
                "subject_id": "p1",
                "site_id": "s1",
                "sample_id_from_lab": "l1",
                "dataset_name": "g1",
                "ena_center_name": "Center A",
                "ena_sample_accession": "ERS123456",
                "ena_study_accession": None,
            },
            {
                "subject_id": "p2",
                "site_id": "s2",
                "sample_id_from_lab": "l2",
                "dataset_name": "g2",
                "ena_center_name": "Center A",
                "ena_sample_accession": None,
                "ena_study_accession": None,
            },
            {
                "subject_id": "p1",
                "site_id": "s3",
                "sample_id_from_lab": "l1",
                "dataset_name": "g1",
                "ena_center_name": "Center B",
                "ena_sample_accession": None,
                "ena_study_accession": None,
            },
        ]
        got_sample_rows = sorted(
            database.get_rows_from_table("Sample"), key=itemgetter("site_id")
        )
        # the rows also have the sample_id, which is made by mysql auto increment,
        # We don't know the order in which things are added, so can't check the sample_id.
        for row in got_sample_rows:
            del row["sample_id"]

        self.assertEqual(expected_sample_rows, got_sample_rows)

        expected_rows = [
            {
                "sequence_replicate_number": 1,
                "original_reads_file_1_md5": "edc176f367fe8e5a014c819b9ec9b05c",
                "original_reads_file_2_md5": "0dd551a0d76d90059808f6f7ddbb0e02",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "pool_sequence_replicates": 1,
                "withdrawn": 0,
                "import_status": 1,
                "submission_date": datetime.date(2017, 12, 25),
                "submit_to_ena": 0,
                "ena_run_accession": "ERR123456",
                "ena_on_hold": 0,
                "isolate_number_from_lab": "1",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
                "instrument_model": "Illumina HiSeq 2000",
            },
            {
                "sequence_replicate_number": 1,
                "original_reads_file_1_md5": "fe5cd28cf9394be14794f0a56a2fe845",
                "original_reads_file_2_md5": "d026fd9a439294ed42795bd7f1e7df10",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "pool_sequence_replicates": 1,
                "withdrawn": 0,
                "import_status": 1,
                "submission_date": datetime.date(2017, 12, 26),
                "submit_to_ena": 1,
                "ena_run_accession": None,
                "ena_on_hold": 1,
                "isolate_number_from_lab": "1",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
                "instrument_model": "Illumina HiSeq 2000",
            },
            {
                "sequence_replicate_number": 1,
                "original_reads_file_1_md5": "aa8f077673c158c4f2a19fc3c50e3fa7",
                "original_reads_file_2_md5": "ae6bafef67da3c26576e799c32985ac9",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "pool_sequence_replicates": 1,
                "withdrawn": 0,
                "import_status": 1,
                "submission_date": datetime.date(2017, 12, 26),
                "submit_to_ena": 1,
                "ena_run_accession": None,
                "ena_on_hold": 1,
                "isolate_number_from_lab": "2",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
                "instrument_model": "Illumina HiSeq 2000",
            },
            {
                "sequence_replicate_number": 1,
                "original_reads_file_1_md5": "6b9a34ed492dad739ac03e084f3b2ab9",
                "original_reads_file_2_md5": "7ceffc5314ff7e305b4ab5bd859850c9",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "pool_sequence_replicates": 1,
                "withdrawn": 0,
                "import_status": 1,
                "submission_date": datetime.date(2017, 12, 25),
                "submit_to_ena": 1,
                "ena_run_accession": None,
                "ena_on_hold": 0,
                "isolate_number_from_lab": "1",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
                "instrument_model": "Illumina HiSeq 2500",
            },
            {
                "sequence_replicate_number": 2,
                "original_reads_file_1_md5": "ec0377e321c59c0b1b6392a3c6dfc2dc",
                "original_reads_file_2_md5": "d541ffdb43a0648233ec7408c3626bfd",
                "remove_contam_reads_file_1_md5": None,
                "remove_contam_reads_file_2_md5": None,
                "pool_sequence_replicates": 1,
                "withdrawn": 0,
                "import_status": 1,
                "submission_date": datetime.date(2017, 12, 25),
                "submit_to_ena": 1,
                "ena_run_accession": None,
                "ena_on_hold": 0,
                "isolate_number_from_lab": "1",
                "pool_sequence_replicates": 1,
                "ena_experiment_accession": None,
                "instrument_model": "Illumina HiSeq 2500",
            },
        ]

        expected_rows.sort(key=itemgetter("original_reads_file_1_md5"))
        query = "SELECT * FROM (Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id)"
        got_rows = database.query_to_dict(query)
        got_rows.sort(key=itemgetter("original_reads_file_1_md5"))

        # Check reads files etc written correctly
        for isolate_data in got_rows:
            iso_dir = isolate_dir.IsolateDir(
                pipeline_root, isolate_data["sample_id"], isolate_data["isolate_id"]
            )
            self.assertTrue(os.path.exists(iso_dir.reads_dir))

            for i in [1, 2]:
                self.assertTrue(
                    os.path.exists(
                        iso_dir.reads_filename(
                            "original", isolate_data["sequence_replicate_number"], i
                        )
                    )
                )

        # similar to above, we don't know the sample_id, seqrep_id or isolate_id, which are auto generated.
        for row in got_rows:
            del row["sample_id"]
            del row["seqrep_id"]
            del row["isolate_id"]

        self.assertEqual(expected_rows, got_rows)

        shutil.rmtree(pipeline_root)
        nextflow_helper.clean_files()
        database.commit_and_close()
        db_connection.DbConnection(db_ini_file, destroy=True, must_exist=True)
