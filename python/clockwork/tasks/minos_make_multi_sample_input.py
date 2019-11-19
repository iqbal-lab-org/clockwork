import sys
from clockwork import db


def run(options):
    database = db.Db(options.db_config_file)
    lines = database.get_vcfs_and_reads_files_for_minos_multi_sample_calling(
        options.dataset_name,
        options.pipeline_root,
        options.reference_id,
        pipeline_version=options.pipeline_version,
    )

    print(*lines, sep="\n")
