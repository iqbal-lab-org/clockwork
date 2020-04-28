import datetime

from cluster_vcf_records import vcf_file_read, vcf_record
import pyfastaq

from clockwork import utils


def _combine_minos_and_samtools_header(minos_header, samtools_header, ref_seqs):
    header_start = [
        "##fileformat=VCFv4.2",
        "##source=clockwork merge samtools gvcf and minos vcf",
        "##fileDate=" + str(datetime.date.today()),
    ]
    new_header = [
        '##FILTER=<ID=NO_DATA,Description="No information from minos or samtools">',
        '##INFO=<ID=CALLER,Number=1,Description="Origin of call, one of minos, samtools, or none if there was no depth">',
    ]
    new_header.extend(
        [f"##contig=<ID={k},length={len(v)}>" for k, v in ref_seqs.items()]
    )

    exclude = [
        "##fileformat",
        "##fileDate",
        "##minos_max_read_length",
        "##contig",
        "#CHROM",
    ]

    for l in minos_header, samtools_header:
        for line in l:
            skip = False
            for prefix in exclude:
                if line.startswith(prefix):
                    skip = True
                    break

            if not skip:
                new_header.append(line.replace("##INFO", "##FORMAT"))

    new_header.sort()
    assert minos_header[-1].startswith("#CHROM\t")
    new_header.append(minos_header[-1])
    return header_start + new_header


def _move_info_fields_to_format(record):
    """Changes VCF record in place. Moves all the key/values in INFO column
    into the FORMAT column"""
    for k, v in sorted(record.INFO.items()):
        record.set_format_key_value(k, v)
    record.INFO = {}


def gvcf_from_minos_vcf_and_samtools_gvcf(ref_fasta, minos_vcf, samtools_vcf, out_vcf):
    minos_header, minos_records = vcf_file_read.vcf_file_to_dict(minos_vcf)
    samtools_header, samtools_records = vcf_file_read.vcf_file_to_dict(samtools_vcf)
    ref_seqs = {}
    pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)

    with open(out_vcf, "w") as f:
        print(
            *_combine_minos_and_samtools_header(
                minos_header, samtools_header, ref_seqs
            ),
            sep="\n",
            file=f,
        )
        for ref_name, ref_seq in sorted(ref_seqs.items()):
            if ref_name in samtools_records:
                samtools_records[ref_name] = {
                    x.POS: x for x in samtools_records[ref_name]
                }
            else:
                samtools_records[ref_name] = {}

            if ref_name in minos_records:
                minos_iter = iter(minos_records[ref_name])
                minos_record = next(minos_iter)
            else:
                minos_record = None

            ref_pos = 0

            while ref_pos < len(ref_seq):
                while minos_record is None or ref_pos < minos_record.POS:
                    if ref_pos in samtools_records[ref_name]:
                        record = samtools_records[ref_name][ref_pos]
                        _move_info_fields_to_format(record)
                        record.INFO = {"CALLER": "samtools"}
                        print(record, file=f)
                    else:
                        print(
                            ref_name,
                            ref_pos + 1,
                            ".",
                            ref_seqs[ref_name][ref_pos],
                            ".",
                            ".",
                            "NO_DATA",
                            "CALLER=none",
                            "GT:DP",
                            "./.:0",
                            file=f,
                            sep="\t",
                        )
                    ref_pos += 1
                    if ref_pos >= len(ref_seq):
                        break

                if minos_record is not None:
                    minos_record.INFO["CALLER"] = "minos"
                    print(minos_record, file=f)
                    ref_pos = minos_record.ref_end_pos() + 1
                    try:
                        minos_record = next(minos_iter)
                    except StopIteration:
                        minos_record = None


def _samtools_vcf_record_to_frs(record, geno_index):
    dp4 = [int(x) for x in record.FORMAT["DP4"].split(",")]
    dp = sum(dp4)
    if dp == 0:
        return 0
    elif geno_index == 0:
        return (dp4[0] + dp4[1]) / dp
    else:
        return (dp4[2] + dp4[3]) / dp


def _vcf_record_pass_index(record, require_minos_pass=True, min_frs=0.9, min_dp=5):
    """If the VCF record passes the filters, then returns the genotype index
    of the called allele (0=REF, 1=ALT1, 2=ALT2, ...etc). Returns
    None if the record fails."""
    geno_indexes = record.FORMAT["GT"].split("/")
    if "." in geno_indexes or len(set(geno_indexes)) > 1:
        return None

    geno_index = int(geno_indexes[0])
    alt = record.ALT[geno_index - 1]

    if geno_index > 0 and len(alt) != len(record.REF) and alt[0] != record.REF[0]:
        return None
    elif int(record.FORMAT["DP"]) < min_dp:
        return None
    elif record.INFO["CALLER"] == "minos":
        if (require_minos_pass and record.FILTER != {"PASS"}) or float(
            record.FORMAT["FRS"]
        ) < min_frs:
            return None
        else:
            return geno_index
    else:
        assert record.INFO["CALLER"] == "samtools"
        if _samtools_vcf_record_to_frs(record, geno_index) >= min_frs:
            return geno_index
        else:
            return None


def gvcf_to_fasta(gvcf_file, outfile, require_minos_pass=True, min_frs=0.9, min_dp=5):
    sample = "unknown"
    out_seqs = {}
    expect_lengths = {}

    with open(gvcf_file) as f:
        for line in f:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                sample = line.rstrip().split()[-1]
                continue

            record = vcf_record.VcfRecord(line)
            expect_lengths[record.CHROM] = record.POS + 1
            if record.CHROM not in out_seqs:
                out_seqs[record.CHROM] = []
            out_seq = out_seqs[record.CHROM]

            geno_index = _vcf_record_pass_index(
                record,
                require_minos_pass=require_minos_pass,
                min_frs=min_frs,
                min_dp=min_dp,
            )

            if geno_index is None:
                out_seq.extend("N" * len(record.REF))
            elif geno_index == 0:
                out_seq.extend(record.REF)
            else:
                alt = record.ALT[geno_index - 1]
                # This is an indel or a complex variant. VCF convetion says that
                # the nucleotide before the variant is included, and
                # the first nucleotide of the ref and alt is the same. We can
                # put that first one in the FASTA, then put the rest as Ns.
                if len(alt) == len(record.REF) == 1:
                    out_seq.extend(record.ALT[geno_index - 1])
                elif len(alt) != len(record.REF) and alt[0] == record.REF[0]:
                    out_seq.extend(alt[0] + "N" * (len(record.REF) - 1))
                else:
                    out_seq.extend("N" * len(record.REF))

    with open(outfile, "w") as f:
        for name, nucleotides in sorted(out_seqs.items()):
            seq = pyfastaq.sequences.Fasta(f"{name}.{sample}", "".join(nucleotides))
            assert len(seq) == expect_lengths[name]
            print(seq, file=f)
