import os
import sys
from clockwork import reference_dir, db, lock_file

def run(options):
    using_db = None not in (options.db_config_file, options.pipeline_references_root, options.name)
    if (using_db and options.outdir):
        print('Error! If adding to database, must use --db_config_file,--pipeline_references_root,--name.', file=sys.stderr)
        print('Otherwise, use --outdir.', file=sys.stderr)
        sys.exit(1)

    if using_db:
        lock = lock_file.LockFile(os.path.join(options.pipeline_references_root, 'add_reference.lock'))
        database = db.Db(options.db_config_file)
        ref_id = database.add_reference(options.name)
        database.commit_and_close()
        lock.stop()
    else:
        ref_id = None

    ref_dir = reference_dir.ReferenceDir(
        pipeline_references_root_dir=options.pipeline_references_root,
        reference_id=ref_id,
        directory=options.outdir,
    )

    genome_is_big = options.contam_tsv is not None
    using_cortex = options.contam_tsv is None
    ref_dir.make_index_files(options.fasta_file, genome_is_big, using_cortex, cortex_mem_height=options.cortex_mem_height)

    if options.contam_tsv is not None:
        ref_dir.add_remove_contam_metadata_tsv(options.contam_tsv)

