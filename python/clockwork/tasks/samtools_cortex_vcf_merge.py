from clockwork import simple_vcf_merger

def run(options):
    merger = simple_vcf_merger.SimpleVcfMerger(
        options.samtools_vcf,
        options.cortex_vcf,
        options.output_vcf,
        options.ref_fasta,
        homozygous_only=not options.include_homozygous,
        max_REF_len=None,
        min_SNP_qual=options.min_SNP_qual,
        min_dp4=options.min_dp4,
        min_GT_conf=options.min_GT_conf,
    )
    merger.run()

