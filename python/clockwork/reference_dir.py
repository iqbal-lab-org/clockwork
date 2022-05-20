import os
import pyfastaq
from clockwork import contam_remover, cortex, utils


class ReferenceDir:
    def __init__(
        self, pipeline_references_root_dir=None, reference_id=None, directory=None, minimap2_preset="sr",
    ):
        self.pipeline_references_root_dir = pipeline_references_root_dir
        self.reference_id = reference_id
        self.directory = directory
        self.minimap2_preset = minimap2_preset

        if self.directory is not None:
            self.directory = os.path.abspath(directory)
        elif None not in {pipeline_references_root_dir, reference_id}:
            self.pipeline_references_root_dir = os.path.abspath(
                pipeline_references_root_dir
            )
            self.directory = os.path.join(
                self.pipeline_references_root_dir, str(self.reference_id)
            )
        else:
            raise Exception(
                "Must provide directory, or both of pipeline_references_root_dir,reference_id"
            )

        self.ref_fasta_prefix = os.path.join(self.directory, "ref")
        self.ref_fasta = self.ref_fasta_prefix + ".fa"
        self.minimap2_index = self.ref_fasta + ".minimap2_idx"
        self.ref_fai = self.ref_fasta + ".fai"
        self.remove_contam_metadata_tsv = os.path.join(
            self.directory, "remove_contam_metadata.tsv"
        )

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def make_index_files(
        self, fasta_in, genome_is_big, using_cortex, cortex_mem_height=22
    ):
        # seqtk just hangs if the input file doesn't exist, so check
        # it first instead
        if not os.path.exists(fasta_in):
            raise Exception("File not found: " + fasta_in)

        try:
            os.makedirs(self.directory)
        except:
            raise Exception("Error mkdir " + self.directory)

        # ensure sequence lines are 60 nt long, and remove comments from
        # header lines
        utils.syscall("seqtk seq -C -l 60 " + fasta_in + " > " + self.ref_fasta)
        utils.syscall("samtools faidx " + self.ref_fasta)
        utils.syscall(f"minimap2 -x {self.minimap2_preset} -d {self.minimap2_index} {self.ref_fasta}")

        if using_cortex:
            cortex.make_run_calls_index_files(
                self.ref_fasta, self.ref_fasta_prefix, mem_height=cortex_mem_height
            )

    def add_remove_contam_metadata_tsv(self, infile):
        utils.rsync_and_md5(infile, self.remove_contam_metadata_tsv)
        data_by_group, sequence_is_contam = contam_remover.ContamRemover._load_metadata_file(
            self.remove_contam_metadata_tsv
        )
        names_in_contam = set(sequence_is_contam.keys())
        names_in_fasta = {}
        pyfastaq.tasks.lengths_from_fai(self.ref_fai, names_in_fasta)
        names_in_fasta = set(names_in_fasta.keys())
        if names_in_fasta != names_in_contam:
            raise Exception("Mismtach in names from metadata tsv and fasta files")
