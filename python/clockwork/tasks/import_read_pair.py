from clockwork import read_pair_importer


def run(options):
    importer = read_pair_importer.ReadPairImporter(
        options.db_config_file,
        options.pipeline_root,
        options.seqrep_id,
        options.isolate_id,
        options.sample_id,
        options.sequence_rep_number,
        options.reads_file_1,
        options.reads_file_2,
        options.reads_file_md5_1,
        options.reads_file_md5_2,
    )
    importer.run()
