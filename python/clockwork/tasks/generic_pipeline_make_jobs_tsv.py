import os
from clockwork import db, utils, lock_file


def run(options):
    lock = lock_file.LockFile(
        os.path.join(options.pipeline_root, "generic_pipeline.lock")
    )
    database = db.Db(options.db_config_file)
    database.make_generic_pipeline_jobs_tsv(
        options.outfile,
        options.pipeline_root,
        options.pipeline_name,
        pipeline_version=options.pipeline_version,
        dataset_name=options.dataset_name,
    )
    database.commit_and_close()
    lock.stop()
