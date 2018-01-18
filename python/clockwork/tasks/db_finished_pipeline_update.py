import sys
from clockwork import db

def run(options):
    if options.pool == 1:
        options.seqrep_id = None
    else:
        options.seqrep_pool = None
        options.seqrep_id = int(options.seqrep_id)

    database = db.Db(options.db_config_file)
    database.update_finished_pipeline_run(
        options.isolate_id,
        options.seqrep_id,
        options.seqrep_pool,
        options.pipeline_name,
        options.new_pipeline_status,
        reference_id=options.reference_id,
        pipeline_version=options.pipeline_version,
        pipeline_root=options.pipeline_root,
    )

    database.commit_and_close()

