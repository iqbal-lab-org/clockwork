import os
from clockwork import db, lock_file


def run(options):
    lock = lock_file.LockFile(os.path.join(options.pipeline_root, "remove_contam.lock"))
    database = db.Db(options.db_config_file)
    database.make_remove_contam_jobs_tsv(
        options.outfile,
        options.pipeline_root,
        0,
        "/fake/path/to/refs/",
        dataset_name=options.dataset_name,
        faking_it=True,
    )
    database.commit_and_close()
    lock.stop()
