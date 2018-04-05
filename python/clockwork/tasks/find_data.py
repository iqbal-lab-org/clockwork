from clockwork import data_finder

def run(options):
    finder = data_finder.DataFinder(
        options.db_config_file,
        options.pipeline_root,
        include_withdrawn=options.include_withdrawn,
        include_internal_ids=options.include_internal_ids,
        dataset_name=options.dataset_name
    )

    if options.pipeline_name:
        finder.write_pipeline_data_to_file(
            options.outfile,
            options.pipeline_name,
            pipeline_version=options.pipeline_version,
            reference_id=options.reference_id
        )
    else:
        finder.write_seqrep_data_to_file(outfile)

