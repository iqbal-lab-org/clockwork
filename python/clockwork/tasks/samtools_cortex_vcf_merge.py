from clockwork import simple_vcf_merger

def run(options):
    merger = simple_vcf_merger.SimpleVcfMerger(
            options.samtools_vcf,
            options.cortex_vcf,
            options.output_vcf,
            options.reference_seqs,
            options.homozygous_only,
            options.max_REF_len,
            options.min_SNP_qual,
            options.min_dp4,
            options.min_GT_conf,
            mem_height=options.mem_height,
            )
    merger.run()


