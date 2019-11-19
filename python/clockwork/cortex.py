import glob
import os
import tarfile
import shutil
import sys

import pyfastaq
from clockwork import utils


class Error(Exception):
    pass


def _replace_sample_name_in_vcf(infile, outfile, sample_name):
    changed_name = False

    with open(infile) as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                if len(fields) < 10:
                    raise Error("Not enough columns in header line of VCF: " + line)
                elif len(fields) == 10:
                    fields[9] = sample_name
                    print(*fields, sep="\t", file=f_out)
                    changed_name = True
                else:
                    raise Error(
                        "More than one sample in VCF, from header line: " + line
                    )
            else:
                print(line, end="", file=f_out)

    if not changed_name:
        raise Error("No #CHROM line found in VCF file " + infile)


def _find_binaries(
    cortex_root=None, stampy_script=None, vcftools_dir=None, mccortex=None
):
    if cortex_root is None:
        cortex_root = os.environ.get("CLOCKWORK_CORTEX_DIR", "/bioinf-tools/cortex")
    cortex_root = os.path.abspath(cortex_root)

    if not os.path.isdir(cortex_root):
        raise Error("Cortex directory " + cortex_root + " not found. Cannot continue")

    cortex_run_calls = os.path.join(cortex_root, "scripts", "calling", "run_calls.pl")
    if not os.path.exists(cortex_run_calls):
        raise Error(
            "Cortex run_call.pl script "
            + cortex_run_calls
            + " not found. Cannot continue"
        )

    if stampy_script is None:
        stampy_script = os.environ.get(
            "CLOCKWORK_STAMPY_SCRIPT", "/bioinf-tools/stampy-1.0.32/stampy.py"
        )
    stampy_script = os.path.abspath(stampy_script)

    if not os.path.exists(stampy_script):
        raise Error("Stampy script " + stampy_script + " not found. Cannot contine")

    if vcftools_dir is None:
        vcftools_dir = os.environ.get(
            "CLOCKWORK_VCFTOOLS_DIR", "/bioinf-tools/vcftools-0.1.15/"
        )
    vcftools_dir = os.path.abspath(vcftools_dir)

    if not os.path.isdir(vcftools_dir):
        raise Error(
            "vcftools directory " + vcftools_dir + " not found. Cannot continue"
        )

    if mccortex is None:
        mccortex = os.environ.get("CLOCKWORK_MCCORTEX", "/bioinf-tools/mccortex31")
    mccortex = os.path.abspath(mccortex)

    return cortex_root, cortex_run_calls, stampy_script, vcftools_dir, mccortex


def make_run_calls_index_files(
    ref_fasta,
    outprefix,
    cortex_root=None,
    stampy_script=None,
    vcftools_dir=None,
    mem_height=22,
):
    cortex_root, cortex_run_calls, stampy_script, vcftools_dir, mccortex = _find_binaries(
        cortex_root=cortex_root, stampy_script=stampy_script, vcftools_dir=vcftools_dir
    )
    cortex_var = os.path.join(cortex_root, "bin", "cortex_var_31_c1")
    if not os.path.exists(cortex_var):
        raise Error("Cortex script not found: " + cortex_var)

    fofn = outprefix + ".fofn"
    with open(fofn, "w") as f:
        print(os.path.abspath(ref_fasta), file=f)

    utils.syscall(
        " ".join(
            [
                cortex_var,
                "--kmer_size 31",
                "--mem_height",
                str(mem_height),
                "--mem_width 100",
                "--se_list",
                fofn,
                "--max_read_len 10000",
                "--dump_binary",
                outprefix + ".k31.ctx",
                "--sample_id REF",
            ]
        )
    )

    os.unlink(fofn)

    utils.syscall(" ".join([stampy_script, "-G", outprefix + ".stampy", ref_fasta]))

    utils.syscall(
        " ".join(
            [stampy_script, "-g", outprefix + ".stampy", "-H", outprefix + ".stampy"]
        )
    )


class CortexRunCalls:
    def __init__(
        self,
        ref_dir,
        reads_infile,
        outdir,
        sample_name,
        cortex_root=None,
        stampy_script=None,
        vcftools_dir=None,
        mccortex=None,
        mem_height=22,
    ):
        self.ref_dir = os.path.abspath(ref_dir)
        self.reads_infile = os.path.abspath(reads_infile)
        self.outdir = os.path.abspath(outdir)
        self.sample_name = sample_name
        self.cortex_root = cortex_root
        self.stampy_script = stampy_script
        self.vcftools_dir = vcftools_dir
        self.cortex_log = os.path.join(self.outdir, "cortex.log")

        if os.path.exists(self.outdir):
            raise Error("Error! Output directory already exists " + self.outdir)

        self.cortex_outdir = os.path.join(self.outdir, "cortex.out")
        self.cortex_reads_fofn = os.path.join(self.outdir, "cortex.in.fofn")
        self.cortex_reads_index = os.path.join(self.outdir, "cortex.in.index")
        self.cortex_ref_fofn = os.path.join(self.outdir, "cortex.in.index_ref.fofn")
        self.cortex_root, self.cortex_run_calls, self.stampy_script, self.vcftools_dir, self.mccortex = _find_binaries(
            cortex_root=cortex_root,
            stampy_script=stampy_script,
            vcftools_dir=vcftools_dir,
            mccortex=mccortex,
        )
        self.mem_height = mem_height
        self.kmer_counts_file = os.path.join(self.outdir, "kmer_counts.txt.gz")

        # Cortex uses the sample name in some of its output filenames.
        # If the sample name is long, then the filename can end up too long for
        # unix, giving the error 'File name too long'. We'll use a short name
        # while running cortex, then rename the sample inside the VCF file at the end
        self.tmp_sample_name = "sample"

    def _make_input_files(self):
        try:
            os.mkdir(self.outdir)
        except:
            raise Error("Error making directory " + self.outdir + ". Cannot continue")

        with open(self.cortex_reads_fofn, "w") as f:
            print(self.reads_infile, file=f)

        with open(self.cortex_reads_index, "w") as f:
            print(
                self.tmp_sample_name, self.cortex_reads_fofn, ".", ".", sep="\t", file=f
            )

        with open(self.cortex_ref_fofn, "w") as f:
            print(os.path.join(self.ref_dir, "ref.fa"), file=f)

    def _tidy_files(self):
        shutil.rmtree(os.path.join(self.cortex_outdir, "tmp_filelists"))

        for filename in glob.glob(
            os.path.join(self.cortex_outdir, "binaries", "uncleaned", "**", "**")
        ):
            if not (filename.endswith("log") or filename.endswith(".covg")):
                os.unlink(filename)

        for filename in glob.glob(os.path.join(self.cortex_outdir, "calls", "**")):
            if os.path.isdir(filename):
                for filename2 in os.listdir(filename):
                    if not filename2.endswith("log"):
                        os.unlink(os.path.join(filename, filename2))
            elif not (
                filename.endswith("log") or filename.endswith("callsets.genotyped")
            ):
                os.unlink(filename)

        for filename in glob.glob(os.path.join(self.cortex_outdir, "vcfs", "**")):
            if filename.endswith(".vcf"):
                tmp_vcf = filename + ".tmp"
                _replace_sample_name_in_vcf(filename, tmp_vcf, self.sample_name)
                utils.rsync_and_md5(tmp_vcf, filename)
                os.unlink(tmp_vcf)

            if not (
                (filename.endswith(".vcf") and "FINAL" in filename)
                or filename.endswith("log")
                or filename.endswith("aligned_branches")
            ):
                if os.path.isdir(filename):
                    shutil.rmtree(filename)
                else:
                    os.unlink(filename)

    def _run_mccortex_view_kmers(self):
        # Example filename we're looking for:
        # cortex.out/binaries/cleaned/k31/sample_name.kmer31.q5cleaned_1.ctx
        cleaned_dir = os.path.join(self.cortex_outdir, "binaries", "cleaned")
        if not os.path.exists(cleaned_dir):
            print(
                "Cleaned directory not found "
                + cleaned_dir
                + " ... cannot run mccortex view kmers",
                file=sys.stderr,
            )
            return

        kmer_dirs = [x for x in os.listdir(cleaned_dir) if x.startswith("k")]
        if len(kmer_dirs) != 1:
            print(
                "Error finding kmers directory inside "
                + cleaned_dir
                + " ... cannot run mccortex view kmers",
                file=sys.stderr,
            )

        kmer_dir = os.path.join(cleaned_dir, kmer_dirs[0])
        ctx_files = [x for x in os.listdir(kmer_dir) if x.endswith(".ctx")]

        if len(ctx_files) != 1:
            print(
                "Error finding ctx file inside "
                + kmer_dir
                + " ... cannot run mccortex view kmers",
                file=sys.stderr,
            )

        ctx_file = os.path.join(kmer_dir, ctx_files[0])
        assert os.path.exists(ctx_file)
        command = " ".join(
            [
                self.mccortex,
                "view",
                "--kmers",
                ctx_file,
                r"""| awk '{print $1,$2}' | gzip -9 > """,
                self.kmer_counts_file,
            ]
        )
        utils.syscall(command)

    def run(self):
        self._make_input_files()
        ref_fai = os.path.join(self.ref_dir, "ref.fa.fai")
        genome_size = pyfastaq.tasks.stats_from_fai(ref_fai)["total_length"]

        cmd = " ".join(
            [
                self.cortex_run_calls,
                "--fastaq_index",
                self.cortex_reads_index,
                "--auto_cleaning yes",
                "--first_kmer 31",
                "--bc yes",
                "--pd no",
                "--outdir",
                self.cortex_outdir,
                "--outvcf cortex",
                "--ploidy 2",
                "--stampy_hash",
                os.path.join(self.ref_dir, "ref.stampy"),
                "--stampy_bin",
                self.stampy_script,
                "--list_ref_fasta",
                self.cortex_ref_fofn,
                "--refbindir",
                self.ref_dir,
                "--genome_size",
                str(genome_size),
                "--qthresh 5",
                "--mem_height",
                str(self.mem_height),
                "--mem_width 100",
                "--vcftools_dir",
                self.vcftools_dir,
                "--do_union yes",
                "--ref CoordinatesAndInCalling",
                "--workflow independent",
                "--logfile",
                self.cortex_log,
            ]
        )

        utils.syscall(cmd)
        self._tidy_files()
        self._run_mccortex_view_kmers()
