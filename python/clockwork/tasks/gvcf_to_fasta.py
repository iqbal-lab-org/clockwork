from clockwork import gvcf


def run(options):
    if not 0 <= options.min_frs <= 1:
        raise ValueError(f"--min_frs {options.min_frs} must be between 0 and 1")

    gvcf.gvcf_to_fasta(
        options.gvcf,
        options.outfile,
        require_minos_pass=not options.ignore_minos_pass,
        min_frs=options.min_frs,
        min_dp=options.min_dp,
    )
