from clockwork import fastqc


def run(options):
    fqc = fastqc.Fastqc(options.output_dir, options.reads_files)
    fqc.run()
