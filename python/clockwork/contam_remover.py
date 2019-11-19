import os
import pysam
import pyfastaq
from clockwork import utils


class Error(Exception):
    pass


class ContamRemover:
    """Takes a BAM file in read name order (ie has not been sorted by reference).
       The reference should have sequence(s) that we want, and the rest are assumed to
       be contamination.
       Writes FASTQ files, splitting the read pairs into 3 categories:
       1. Reads we want, where either read matches a wanted sequence
       2. Contaminated reads, where either read is mapped, but not to a wanted sequence.
          Only written if contam_out_{1,2} are not None - these are the output filenames
       3. Unmapped, where both reads unmapped. By default, these reads are included
          in reads_out_{1,2}. But if no_match_out_{1,2} given, then written to these
          files instead"""

    def __init__(
        self,
        metadata_file,
        bam_file,
        counts_outfile,
        reads_out_1,
        reads_out_2,
        contam_out_1=None,
        contam_out_2=None,
        no_match_out_1=None,
        no_match_out_2=None,
        done_file=None,
    ):
        self.metadata_file = os.path.abspath(metadata_file)
        self.bam_file = bam_file
        self.counts_outfile = os.path.abspath(counts_outfile)
        self.files = {
            "wanted1": os.path.abspath(reads_out_1),
            "wanted2": os.path.abspath(reads_out_2),
        }

        assert contam_out_1 == None == contam_out_2 or None not in {
            contam_out_1,
            contam_out_2,
        }
        assert no_match_out_1 == None == no_match_out_2 or None not in {
            no_match_out_1,
            no_match_out_2,
        }

        if contam_out_1 is not None:
            self.files["contam1"] = os.path.abspath(contam_out_1)
            self.files["contam2"] = os.path.abspath(contam_out_2)

        if no_match_out_1 is not None:
            self.files["no_match1"] = os.path.abspath(no_match_out_1)
            self.files["no_match2"] = os.path.abspath(no_match_out_2)

        self.done_file = None if done_file is None else os.path.abspath(done_file)

    @classmethod
    def _load_metadata_file(cls, infile):
        data_by_group = {}
        sequence_is_contam = {}

        with open(infile) as f:
            for line in f:
                try:
                    group, is_contam, *names = line.rstrip().split("\t")
                except:
                    raise Error("Error parsing line:\n" + line)

                is_contam = {"0": False, "1": True}[is_contam]
                if (
                    group in data_by_group
                    and data_by_group[group]["contam"] != is_contam
                ):
                    raise Error(
                        'Error! is_contam for one group must always be the same (column 2). Check group "'
                        + group
                        + '"'
                    )
                elif group not in data_by_group:
                    data_by_group[group] = {"contam": is_contam, "sequences": set()}

                data_by_group[group]["sequences"].update(names)

                for name in names:
                    if name in sequence_is_contam:
                        raise Error('Error! Duplicate sequence name "' + name + '"')
                    sequence_is_contam[name] = is_contam

        return data_by_group, sequence_is_contam

    @staticmethod
    def _sam_to_fastq(sam):
        """Given a pysam alignment, returns the sequence a Fastq object.
           Reverse complements as required and add suffix /1 or /2 as appropriate from the flag"""
        name = sam.qname
        if sam.is_read1:
            name += "/1"
        elif sam.is_read2:
            name += "/2"
        else:
            raise Error(
                "Read "
                + name
                + " must be first or second of pair according to flag. Cannot continue"
            )

        seq = pyfastaq.sequences.Fastq(
            name, utils.decode(sam.seq), utils.decode(sam.qual)
        )

        if sam.is_reverse:
            seq.revcomp()

        return seq

    @staticmethod
    def _read_mapped_and_wanted(sam_record, contam_dict):
        """Takes a pysam.AlignedSegment object. Returns tuple of bools:
        (is mapped, is wanted (ie does not match contamination genome)"""
        if not sam_record.is_unmapped:
            return (True, not contam_dict[sam_record.reference_name])
        else:
            return (False, False)

    @classmethod
    def _write_read_counts_by_group_file(cls, sequences_by_group, read_counts, outfile):
        group_counts = {}
        for group in sequences_by_group:
            group_counts[group] = sum(
                [
                    read_counts.get(sequence, 0)
                    for sequence in sequences_by_group[group]["sequences"]
                ]
            )

        with open(outfile, "w") as f:
            print("Name\tIs_contam\tReads", file=f)
            for group, count in sorted(group_counts.items()):
                is_contam = 1 if sequences_by_group[group]["contam"] else 0
                print(group, is_contam, count, sep="\t", file=f)

            print("Unmapped", 0, read_counts.get("Unmapped", 0), sep="\t", file=f)
            print(
                "Reads_kept_after_remove_contam",
                0,
                read_counts.get("reads_kept_after_remove_contam", 0),
                sep="\t",
                file=f,
            )

    def run(self):
        ref_seqs_by_group, ref_seqs_contam = ContamRemover._load_metadata_file(
            self.metadata_file
        )
        filehandles = {
            f: pyfastaq.utils.open_file_write(self.files[f]) for f in self.files
        }
        open_type = "rb" if self.bam_file.endswith(".bam") else "r"
        sam_reader = pysam.Samfile(self.bam_file, open_type)
        previous_sam = None
        read_counts = {"reads_kept_after_remove_contam": 0}

        for current_sam in sam_reader.fetch(until_eof=True):
            if previous_sam is None:
                previous_sam = current_sam
                continue

            # reads should be in read pair order in the input, but handle if
            # the second read is before the first
            if previous_sam.is_read1:
                assert current_sam.is_read2
            else:
                assert current_sam.is_read1
                previous_sam, current_sam = current_sam, previous_sam

            read_1_mapped, read_1_wanted = ContamRemover._read_mapped_and_wanted(
                previous_sam, ref_seqs_contam
            )
            read_2_mapped, read_2_wanted = ContamRemover._read_mapped_and_wanted(
                current_sam, ref_seqs_contam
            )

            key1 = previous_sam.reference_name if read_1_mapped else "Unmapped"
            key2 = current_sam.reference_name if read_2_mapped else "Unmapped"
            read_counts[key1] = read_counts.get(key1, 0) + 1
            read_counts[key2] = read_counts.get(key2, 0) + 1

            if read_1_wanted or read_2_wanted:
                out_type = "wanted"
            elif read_1_mapped or read_2_mapped:
                out_type = "contam" if "contam1" in self.files else None
            else:
                out_type = "no_match" if "no_match1" in self.files else "wanted"

            if out_type == "wanted":
                read_counts["reads_kept_after_remove_contam"] += 2

            if out_type is not None:
                print(
                    ContamRemover._sam_to_fastq(previous_sam),
                    file=filehandles[out_type + "1"],
                )
                print(
                    ContamRemover._sam_to_fastq(current_sam),
                    file=filehandles[out_type + "2"],
                )

            previous_sam = None

        for filehandle in filehandles.values():
            pyfastaq.utils.close(filehandle)

        if previous_sam is not None:
            raise Error(
                "Odd number of sequences in input BAM. Expected paired reads only."
            )

        ContamRemover._write_read_counts_by_group_file(
            ref_seqs_by_group, read_counts, self.counts_outfile
        )

        if self.done_file is not None:
            utils.make_empty_file(self.done_file)
