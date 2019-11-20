import unittest
import glob
import sys
import shutil
import os
from operator import itemgetter

from clockwork import db, db_connection, isolate_dir, utils
from clockwork import __version__ as clockwork_version

sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
import nextflow_helper

data_dir = os.path.join(nextflow_helper.data_root_dir, "nextflow_generic_pipeline")
modules_dir = os.path.dirname(os.path.abspath(db.__file__))
db_ini_file = os.path.join(modules_dir, "tests", "data", "db.ini")


class TestNextflowGenericPipeline(unittest.TestCase):
    def test_nextflow_generic_pipeline(self):
        """test nextflow generic pipeline using database"""
        tmp_data_dir = "tmp.nextflow_generic_pipeline_db_input.data"
        if os.path.exists(tmp_data_dir):
            shutil.rmtree(tmp_data_dir)
        shutil.copytree(data_dir, tmp_data_dir)
        nextflow_helper.write_config_file()
        mysql_config_file = os.path.join(data_dir, "db.cnf")
        mysql_dump = os.path.join(data_dir, "mysql.dump")
        db_config_data = db_connection.DbConnection._parse_config_file(db_ini_file)
        utils.syscall(
            "mysql --defaults-file="
            + mysql_config_file
            + ' -e "DROP DATABASE IF EXISTS '
            + db_config_data["db"]
            + "; CREATE DATABASE "
            + db_config_data["db"]
            + '"'
        )
        utils.syscall(
            "mysql --defaults-file="
            + mysql_config_file
            + " "
            + db_config_data["db"]
            + " < "
            + mysql_dump
        )
        pipeline_root = os.path.join(tmp_data_dir, "Pipeline_root")
        nextflow_file = os.path.join(
            nextflow_helper.nextflow_dir, "generic_pipeline.nf"
        )
        work_dir = "tmp.nextflow_generic_pipeline.work"
        dag_file = "nextflow.generic_pipeline.dag.pdf"
        pipeline_name = "generic_pipeline"
        script = os.path.join(data_dir, "script.pl")

        try:
            os.unlink(dag_file)
        except:
            pass

        command = " ".join(
            [
                "nextflow run",
                "--dataset_name g1",  # one read pair is from group 2 and should get ignored
                "--pipeline_name",
                pipeline_name,
                "--pipeline_root",
                pipeline_root,
                "--script",
                script,
                "--db_config_file",
                db_ini_file,
                "--max_ram",
                "0.5",
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

        # check database Pipeline table updated as expected
        database = db.Db(db_ini_file)
        got_rows = database.get_rows_from_table("Pipeline")
        got_rows.sort(key=itemgetter("isolate_id", "pipeline_name"))
        expected_rows = [
            {
                "isolate_id": 1,
                "seqrep_id": 1,
                "seqrep_pool": None,
                "version": "0.1.2",
                "pipeline_name": "remove_contam",
                "status": 1,
                "reference_id": 1,
            },
            {
                "isolate_id": 1,
                "seqrep_id": 2,
                "seqrep_pool": None,
                "version": "0.1.2",
                "pipeline_name": "remove_contam",
                "status": 1,
                "reference_id": 1,
            },
            {
                "isolate_id": 1,
                "seqrep_id": None,
                "seqrep_pool": "1_2",
                "version": clockwork_version,
                "pipeline_name": pipeline_name,
                "status": 1,
                "reference_id": None,
            },
            {
                "isolate_id": 2,
                "seqrep_id": 3,
                "seqrep_pool": None,
                "version": "0.1.2",
                "pipeline_name": "remove_contam",
                "status": 1,
                "reference_id": 1,
            },
            {
                "isolate_id": 2,
                "seqrep_id": 4,
                "seqrep_pool": None,
                "version": "0.1.2",
                "pipeline_name": "remove_contam",
                "status": 1,
                "reference_id": 1,
            },
            {
                "isolate_id": 2,
                "seqrep_id": 3,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": pipeline_name,
                "status": 1,
                "reference_id": None,
            },
            {
                "isolate_id": 2,
                "seqrep_id": 4,
                "seqrep_pool": None,
                "version": clockwork_version,
                "pipeline_name": pipeline_name,
                "status": 1,
                "reference_id": None,
            },
            {
                "isolate_id": 3,
                "seqrep_id": 5,
                "seqrep_pool": None,
                "version": "0.1.2",
                "pipeline_name": "remove_contam",
                "status": 1,
                "reference_id": 1,
            },
            {
                "isolate_id": 3,
                "seqrep_id": None,
                "seqrep_pool": "1",
                "version": clockwork_version,
                "pipeline_name": pipeline_name,
                "status": -1,
                "reference_id": None,
            },
            {
                "isolate_id": 4,
                "seqrep_id": 6,
                "seqrep_pool": None,
                "version": "0.1.2",
                "pipeline_name": "remove_contam",
                "status": 1,
                "reference_id": 1,
            },
        ]
        expected_rows.sort(key=itemgetter("isolate_id", "pipeline_name"))
        self.assertEqual(expected_rows, got_rows)

        # check that the expected output file from the script.pl
        # got made (except for the sample that is expected to fail)

        ids = [
            {"sample": 1, "seqrep_id": "1_2", "isolate_id": 1, "seq_repl": "1_2"},
            {"sample": 2, "seqrep_id": 3, "isolate_id": 2, "seq_repl": "1"},
            {"sample": 2, "seqrep_id": 4, "isolate_id": 2, "seq_repl": "2"},
        ]
        for id_dict in ids:
            iso_dir = isolate_dir.IsolateDir(
                pipeline_root, id_dict["sample"], id_dict["isolate_id"]
            )
            pipeline_dir = iso_dir.pipeline_dir(
                id_dict["seq_repl"], pipeline_name, clockwork_version
            )
            counts_file = os.path.join(pipeline_dir, "count.txt")
            self.assertTrue(os.path.exists(counts_file))

        shutil.rmtree(tmp_data_dir)
        nextflow_helper.clean_files()
