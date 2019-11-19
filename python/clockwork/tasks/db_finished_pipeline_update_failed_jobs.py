from clockwork import db


def run(options):
    database = db.Db(options.db_config_file)
    database.update_finished_pipeline_run_failed_jobs(
        options.jobs_tsv,
        options.success_jobs_file,
        options.pipeline_name,
        reference_id=options.reference_id,
        pipeline_version=options.pipeline_version,
    )

    database.commit_and_close()
