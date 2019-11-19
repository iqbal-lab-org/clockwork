from clockwork import samtools_qc


def run(options):
    sqc = samtools_qc.SamtoolsQc(
        options.ref_fasta, options.reads1, options.reads2, options.output_dir
    )
    sqc.run()
