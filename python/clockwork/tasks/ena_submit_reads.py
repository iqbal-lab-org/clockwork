from clockwork.ena import dataset_submitter

def run(options):
    submitter = dataset_submitter.DatasetSubmitter(
       options.ini_file,
       options.dataset_name,
       options.pipeline_root,
       options.taxon_id,
       fq_upload_threads=options.fq_upload_threads,
       use_test_server=options.use_test_server,
    )
    submitter.run()
