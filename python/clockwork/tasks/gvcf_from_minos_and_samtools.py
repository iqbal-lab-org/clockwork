from clockwork import gvcf


def run(options):
    gvcf.gvcf_from_minos_vcf_and_samtools_gvcf(
        options.ref_fasta, options.minos_vcf, options.samtools_vcf, options.outfile,
    )
