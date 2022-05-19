from operator import itemgetter
import os

import pyfastaq

from clockwork import db, isolate_dir


mysql_seqrep_isolate_sample_join = "Seqrep join Isolate on Seqrep.isolate_id = Isolate.isolate_id join Sample on Isolate.sample_id = Sample.sample_id"


class DataFinder:
    def __init__(
        self,
        db_ini_file,
        pipeline_root,
        include_withdrawn=False,
        include_internal_ids=False,
        dataset_name=None,
    ):
        self.db = db.Db(db_ini_file)
        self.pipeline_root = os.path.abspath(pipeline_root)
        if not os.path.exists(self.pipeline_root):
            raise Exception(
                'Pipeline root directory "'
                + self.pipeline_root
                + '" not found. Cannot continue'
            )
        self.dataset_name = dataset_name
        self.include_withdrawn = include_withdrawn
        self.include_internal_ids = include_internal_ids
        self.dataset_name = dataset_name

    def write_seqrep_data_to_file(self, outfile):
        query = "select * from " + mysql_seqrep_isolate_sample_join
        columns = [
            "site_id",
            "subject_id",
            "sample_id_from_lab",  # = lab_id in import spreadsheet
            "isolate_number_from_lab",  # = isolate_number in import spreadhseet
            "sequence_replicate_number",
            "pool_sequence_replicates",
            "submission_date",
            "import_status",
            "dataset_name",
            "submit_to_ena",
            "ena_on_hold",
            "ena_center_name",
            "ena_experiment_accession",
            "ena_run_accession",
            "ena_sample_accession",
            "ena_study_accession",
            "instrument_model",
            "original_reads_file_1_md5",
            "original_reads_file_2_md5",
            "isolate_directory",
            "remove_contam_reads_1",
            "remove_contam_reads_file_1_md5",
            "remove_contam_reads_2",
            "remove_contam_reads_file_2_md5",
        ]

        if self.include_internal_ids:
            columns.extend(["sample_id", "isolate_id", "seqrep_id"])

        where_fields = []

        if self.include_withdrawn:
            columns.append("withdrawn")
        else:
            where_fields.append("withdrawn=0")

        if self.dataset_name is not None:
            where_fields.append('dataset_name="' + self.dataset_name + '"')

        if len(where_fields):
            query += " where " + " AND ".join(where_fields)

        rows = self.db.query_to_dict(query)
        rows.sort(
            key=itemgetter(
                "site_id",
                "subject_id",
                "sample_id_from_lab",
                "isolate_number_from_lab",
                "sequence_replicate_number",
            )
        )
        f = pyfastaq.utils.open_file_write(outfile)
        print(*columns, sep="\t", file=f)

        for row in rows:
            iso_dir_obj = isolate_dir.IsolateDir(
                self.pipeline_root, row["sample_id"], row["isolate_id"]
            )
            row["isolate_directory"] = iso_dir_obj.isolate_dir
            row["remove_contam_reads_1"] = iso_dir_obj.reads_filename(
                "remove_contam", row["sequence_replicate_number"], 1
            )
            row["remove_contam_reads_2"] = iso_dir_obj.reads_filename(
                "remove_contam", row["sequence_replicate_number"], 2
            )
            print(*[row[x] for x in columns], sep="\t", file=f)

        pyfastaq.utils.close(f)

    def write_pipeline_data_to_file(
        self, outfile, pipeline_name, pipeline_version=None, reference_id=None
    ):
        where_fields = ['pipeline_name="' + pipeline_name + '"']
        pooling_seqreps = pipeline_name in {"mykrobe_predict", "variant_call"}
        seqreps_column = (
            "sequence_replicate_numbers"
            if pooling_seqreps
            else "sequence_replicate_number"
        )

        if pipeline_version is not None:
            where_fields.append('version="' + pipeline_version + '"')

        if reference_id is not None:
            where_fields.append("reference_id=" + str(reference_id))

        if self.dataset_name is not None:
            where_fields.append('dataset_name="' + self.dataset_name + '"')

        columns = [
            "pipeline_name",
            "version",
            "dataset_name",
            "reference_id",
            "site_id",
            "subject_id",
            "sample_id_from_lab",  # = lab_id in import spreadsheet
            "isolate_number_from_lab",  # = isolate_number in import spreadhseet
            seqreps_column,
            "pipeline_directory",
        ]

        if pipeline_name == "remove_contam":
            columns.extend(
                [
                    "remove_contam_reads_1",
                    "remove_contam_reads_file_1_md5",
                    "remove_contam_reads_2",
                    "remove_contam_reads_file_2_md5",
                ]
            )

        if self.include_internal_ids:
            columns.extend(["sample_id", "isolate_id", "seqrep_id"])

        if pooling_seqreps:
            pipeline_join = " join Pipeline on Pipeline.isolate_id = Isolate.isolate_id"
        else:
            pipeline_join = " join Pipeline on Pipeline.seqrep_id = Seqrep.seqrep_id"

        query = (
            "select * from "
            + mysql_seqrep_isolate_sample_join
            + " "
            + pipeline_join
            + " where "
            + " AND ".join(where_fields)
        )

        rows = self.db.query_to_dict(query)
        rows.sort(
            key=itemgetter(
                "site_id",
                "subject_id",
                "sample_id_from_lab",
                "isolate_number_from_lab",
                "sequence_replicate_number",
            )
        )
        f = pyfastaq.utils.open_file_write(outfile)
        print(*columns, sep="\t", file=f)

        # The joining of tales means that if the pipelne pools samples, could
        # have duplicate rows. eg pool is 1_2. We'll get a row for replicate
        # 1 and replicate 2. Track what we've seen already and skip the duplicates
        used_pools = set()

        for row in rows:
            iso_dir_obj = isolate_dir.IsolateDir(
                self.pipeline_root, row["sample_id"], row["isolate_id"]
            )
            ref_id = None if pipeline_name == "qc" else row["reference_id"]
            if pooling_seqreps:
                row[seqreps_column] = row["seqrep_pool"]
                pool_tuple = (
                    row["version"],
                    row["isolate_id"],
                    row["seqrep_pool"],
                    row["reference_id"],
                )
                if pool_tuple in used_pools:
                    continue
                used_pools.add(pool_tuple)

            row["pipeline_directory"] = iso_dir_obj.pipeline_dir(
                row[seqreps_column], pipeline_name, row["version"], reference_id=ref_id
            )
            if pipeline_name == "remove_contam":
                row["remove_contam_reads_1"] = iso_dir_obj.reads_filename(
                    "remove_contam", row["sequence_replicate_number"], 1
                )
                row["remove_contam_reads_2"] = iso_dir_obj.reads_filename(
                    "remove_contam", row["sequence_replicate_number"], 2
                )
            print(*[row[x] for x in columns], sep="\t", file=f)

        pyfastaq.utils.close(f)
