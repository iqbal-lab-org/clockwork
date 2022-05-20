import os
from clockwork import utils, var_call_one_sample_pipeline


def run(options):
    if len(options.reads_files) % 2 != 0:
        raise Exception(
            f"Must provide even number of reads files. Got these files: {', '.join(options.reads_files)}"
        )
    reads1 = [f for i, f in enumerate(options.reads_files) if i % 2 == 0]
    reads2 = [f for i, f in enumerate(options.reads_files) if i % 2 != 0]
    assert len(reads1) == len(reads2)
    if options.force and os.path.exists(options.outdir):
        utils.syscall(f"rm -r {options.outdir}")

    var_call_one_sample_pipeline.run(
        reads1,
        reads2,
        options.ref_dir,
        options.outdir,
        sample_name=options.sample_name,
        cortex_mem_height=options.mem_height,
        debug=options.debug,
        keep_bam=options.keep_bam,
    )
