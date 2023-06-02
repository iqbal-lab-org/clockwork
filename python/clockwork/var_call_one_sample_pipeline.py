import logging
import os

from clockwork import cortex, gvcf, read_map, read_trim, reference_dir, utils


def run(
    reads1_list,
    reads2_list,
    ref_dir,
    outdir,
    sample_name="sample",
    cortex_mem_height=22,
    debug=False,
    keep_bam=False,
    trim_reads=True,
):
    if len(reads1_list) != len(reads2_list):
        raise Exception(
            "Must give same number of forward and reverse reads files. Got:\nForward:{reads1_list}\nReverse:{reads2_list}"
        )

    os.mkdir(outdir)

    trimmed_reads_1 = []
    trimmed_reads_2 = []
    for i in range(len(reads1_list)):
        if not trim_reads:
            trimmed_reads_1.append(reads1_list[i])
            trimmed_reads_2.append(reads2_list[i])
            continue

        trimmed_reads_1.append(os.path.join(outdir, f"trimmed_reads.{i}.1.fq.gz"))
        trimmed_reads_2.append(os.path.join(outdir, f"trimmed_reads.{i}.2.fq.gz"))
        read_trim.run_trimmomatic(
            reads1_list[i],
            reads2_list[i],
            trimmed_reads_1[-1],
            trimmed_reads_2[-1],
        )

    refdir = reference_dir.ReferenceDir(directory=ref_dir)
    rmdup_bam = os.path.join(outdir, "map.bam")
    read_map.map_reads_set(
        refdir.ref_fasta,
        trimmed_reads_1,
        trimmed_reads_2,
        rmdup_bam,
        rmdup=True,
        read_group=("1", sample_name),
    )
    utils.syscall(f"samtools index {rmdup_bam}")
    if trim_reads and not debug:
        for filename in trimmed_reads_1 + trimmed_reads_2:
            os.unlink(filename)

    samtools_vcf = os.path.join(outdir, "samtools.vcf")
    cmd = f"bcftools mpileup --output-type u -f {refdir.ref_fasta} {rmdup_bam} | bcftools call -vm -O v -o {samtools_vcf}"
    utils.syscall(cmd)

    cortex_dir = os.path.join(outdir, "cortex")
    ctx = cortex.CortexRunCalls(
        refdir.directory,
        rmdup_bam,
        cortex_dir,
        sample_name,
        mem_height=cortex_mem_height,
    )
    ctx.run(run_mccortex_view_kmers=False)
    ctx_vcf_dir = os.path.join(cortex_dir, "cortex.out", "vcfs")
    cortex_vcfs = [
        os.path.join(ctx_vcf_dir, x)
        for x in os.listdir(ctx_vcf_dir)
        if x.endswith("raw.vcf")
    ]
    if len(cortex_vcfs) != 1:
        raise Exception("Error running cortex. Could not find output VCF file")
    cortex_vcf = os.path.join(outdir, "cortex.vcf")
    os.rename(cortex_vcfs[0], cortex_vcf)
    if not debug:
        utils.syscall(f"rm -rf {cortex_dir}")

    minos_dir = os.path.join(outdir, "minos")
    cmd = f"minos adjudicate --reads {rmdup_bam} {minos_dir} {refdir.ref_fasta} {samtools_vcf} {cortex_vcf}"
    utils.syscall(cmd)
    final_vcf = os.path.join(outdir, "final.vcf")
    os.rename(os.path.join(minos_dir, "final.vcf"), final_vcf)

    samtools_gvcf = os.path.join(outdir, "samtools.gvcf")
    cmd = f"bcftools mpileup -I --output-type u -f {refdir.ref_fasta} {rmdup_bam} | bcftools call -c -O v -o {samtools_gvcf}"
    utils.syscall(cmd)
    final_gvcf = os.path.join(outdir, "final.gvcf")
    gvcf.gvcf_from_minos_vcf_and_samtools_gvcf(
        refdir.ref_fasta, final_vcf, samtools_gvcf, final_gvcf,
    )
    if not debug:
        os.unlink(samtools_gvcf)
    gvcf.gvcf_to_fasta(final_gvcf, f"{final_gvcf}.fasta")

    if not debug:
        utils.syscall(f"rm -rf {minos_dir}")

    if not (keep_bam or debug):
        os.unlink(rmdup_bam)
        os.unlink(rmdup_bam + ".bai")

    final_vcf = os.path.join(outdir, "final.vcf")
    if not os.path.exists(final_vcf):
        raise Exception(f"Error. Final VCF file not found: {final_vcf}")

    logging.info(f"Finished variant calling. Final VCF file: {final_vcf}")
