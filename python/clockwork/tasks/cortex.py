from clockwork import cortex

def run(options):
    ctx = cortex.CortexRunCalls(
        options.ref_dir,
        options.reads_file,
        options.output_dir,
        options.sample_name,
        mem_height=options.mem_height,
    )
    ctx.run()

