import os
from clockwork import db, utils, lock_file

def run(options):
    lock = lock_file.LockFile(os.path.join(options.pipeline_root, 'remove_contam.lock'))
    database = db.Db(options.db_config_file)
    database.make_remove_contam_jobs_tsv(
        options.outfile,
        options.pipeline_root,
        options.reference_id,
        options.reference_root,
        pipeline_version=options.pipeline_version,
        dataset_name=options.dataset_name,
    )
    database.commit_and_close()
    lock.stop()

