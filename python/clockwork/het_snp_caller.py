import os
import re
import pyfastaq
from clockwork import utils


class Error(Exception):
    pass


# This is a rewrite into Python of the Perl package here:
# https://github.com/sanger-pathogens/vr-codebase/blob/master/modules/Pathogens/QC/HetSNPCalculator.pm

adf_regex = re.compile(r"""[\t;]ADF=(?P<adf>[0-9,]+)[\t;]""")
adr_regex = re.compile(r"""[\t;]ADR=(?P<adr>[0-9,]+)[\t;]""")


class HetSnpCaller:
    def __init__(
        self,
        sorted_bam,
        ref_fasta,
        outprefix,
        min_total_depth=4,
        min_second_depth=2,
        max_allele_freq=0.9,
    ):
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.bam = os.path.abspath(sorted_bam)
        self.outprefix = os.path.abspath(outprefix)
        self.min_total_depth = min_total_depth
        self.min_second_depth = min_second_depth
        self.max_allele_freq = max_allele_freq

    @classmethod
    def _run_mpileup(cls, bam, ref, outfile):
        cmd = " ".join(
            [
                "samtools mpileup --skip-indels -d 500 -t INFO/AD,INFO/ADF,INFO/ADR -C50 -uv",
                "-f",
                ref,
                bam,
                ">",
                outfile,
            ]
        )

        utils.syscall(cmd)

    @classmethod
    def _vcf_line_is_snp_and_or_het(
        cls, vcf_line, min_total_depth, min_second_depth, max_allele_freq
    ):
        """returns tuple: (chromosome name, is a SNP, is a het SNP)"""
        fields = vcf_line.rstrip().split("\t")

        if "," not in fields[4]:
            return (fields[0], False, False)

        adf_search = adf_regex.search(fields[7])
        adr_search = adr_regex.search(fields[7])
        if None in [adf_search, adr_search]:
            raise Error("Error! Must have ADF and ADR in info column: " + vcf_line)

        adf_list = [int(x) for x in adf_search.group("adf").split(",")]
        adr_list = [int(x) for x in adr_search.group("adr").split(",")]
        if len(adf_list) != len(adr_list):
            raise Error("Mismatch in lengths of ADF and ADR lists: " + vcf_line)

        adf_sum = sum(adf_list)
        if adf_sum < min_total_depth:
            return (fields[0], False, False)

        adr_sum = sum(adr_list)
        if adr_sum < min_total_depth:
            return (fields[0], False, False)

        adr_max = max(adr_list)
        adf_max = max(adf_list)

        for i, adf in enumerate(adf_list):
            adr = adr_list[i]

            if (
                adf >= min_second_depth
                and adf / adf_sum <= max_allele_freq
                and adr >= min_second_depth
                and adr / adr_sum <= max_allele_freq
            ):
                return (fields[0], True, True)

        if adf_list[0] < adf_max and adr_list[0] < adr_max:
            return (fields[0], True, False)
        else:
            return (fields[0], False, False)

    @classmethod
    def _filter_vcf_and_count_snps(
        cls,
        vcf_file_in,
        vcf_file_out,
        min_total_depth,
        min_second_depth,
        max_allele_freq,
    ):
        results = {}

        with open(vcf_file_in) as f_in, open(vcf_file_out, "w") as f_out:
            for line in f_in:
                if line.startswith("#"):
                    print(line, end="", file=f_out)
                    continue

                chrom, is_snp, is_het = HetSnpCaller._vcf_line_is_snp_and_or_het(
                    line, min_total_depth, min_second_depth, max_allele_freq
                )
                if chrom not in results:
                    results[chrom] = {x: 0 for x in ["positions", "snps", "hets"]}

                results[chrom]["positions"] += 1
                if is_snp:
                    results[chrom]["snps"] += 1
                if is_het:
                    results[chrom]["hets"] += 1
                    print(line, end="", file=f_out)

        return results

    @classmethod
    def _write_reports(cls, snp_data, contig_lengths, summary_out, per_contig_out):
        totals = {x: 0 for x in ["length", "positions", "snps", "hets"]}

        with open(per_contig_out, "w") as f:
            print(
                "Contig",
                "Length",
                "Positions_mapped",
                "SNPs",
                "Het_SNPs",
                sep="\t",
                file=f,
            )

            for contig in sorted(contig_lengths):
                totals["length"] += contig_lengths[contig]

                if contig in snp_data:
                    for key, number in snp_data[contig].items():
                        totals[key] += number

                    print(
                        contig,
                        contig_lengths[contig],
                        snp_data[contig]["positions"],
                        snp_data[contig]["snps"],
                        snp_data[contig]["hets"],
                        sep="\t",
                        file=f,
                    )
                else:
                    print(contig, 0, 0, 0, 0, sep="\t", file=f)

        with open(summary_out, "w") as f:
            if totals["snps"] > 0:
                percent_snps = round(100 * totals["hets"] / totals["snps"], 2)
            else:
                percent_snps = 0
            print(
                "Total_length",
                "Positions_used",
                "Total_SNPs",
                "Het_SNPs",
                "Percent_SNPs_are_het",
                sep="\t",
                file=f,
            )
            print(
                totals["length"],
                totals["positions"],
                totals["snps"],
                totals["hets"],
                percent_snps,
                sep="\t",
                file=f,
            )

    def run(self):
        unfiltered_vcf = self.outprefix + ".unfiltered.vcf"
        filtered_vcf = self.outprefix + ".het_calls.vcf"
        ref_lengths = {}
        pyfastaq.tasks.lengths_from_fai(self.ref_fasta + ".fai", ref_lengths)
        HetSnpCaller._run_mpileup(self.bam, self.ref_fasta, unfiltered_vcf)
        snp_data = HetSnpCaller._filter_vcf_and_count_snps(
            unfiltered_vcf,
            filtered_vcf,
            self.min_total_depth,
            self.min_second_depth,
            self.max_allele_freq,
        )
        os.unlink(unfiltered_vcf)
        summary_tsv = self.outprefix + ".summary.tsv"
        per_contig_tsv = self.outprefix + ".per_contig.tsv"
        HetSnpCaller._write_reports(snp_data, ref_lengths, summary_tsv, per_contig_tsv)
