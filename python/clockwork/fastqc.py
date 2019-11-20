import copy
import shutil
import os
from clockwork import utils


class Error(Exception):
    pass


class Fastqc:
    def __init__(self, outdir, infiles):
        assert isinstance(infiles, list)
        self.infiles = copy.copy(infiles)
        self.outdir = outdir

    @staticmethod
    def _run_fastqc(outdir, infiles):
        assert isinstance(infiles, list)
        try:
            os.mkdir(outdir)
        except:
            raise Error("Error mkdir " + outdir)

        command = "fastqc --threads 1 --extract -o " + outdir + " " + " ".join(infiles)
        utils.syscall(command)

    @classmethod
    def _seq_length_from_report_string(cls, string):
        if "-" in string:
            return tuple([int(x) for x in string.split("-")])
        else:
            return int(string), int(string)

    @staticmethod
    def _stats_from_report(filename):
        stats = {}
        interesting_line_starts = {
            "Total Sequences",
            "Sequences flagged as poor quality",
            "Sequence length",
            "%GC",
            "Filename",
        }

        with open(filename) as f:
            for line in f:
                if (
                    line.startswith(">>") and not line.startswith(">>END_MODULE")
                ) or line.split("\t")[0] in interesting_line_starts:
                    try:
                        stat, value = line.lstrip(">").rstrip().split("\t")
                    except:
                        raise Error("Error splitting line: " + line)

                    stat = "_".join(stat.split()).lower().replace("%", "")
                    if stat == "sequence_length":
                        stats["min_sequence_length"], stats[
                            "max_sequence_length"
                        ] = Fastqc._seq_length_from_report_string(value)
                    else:
                        try:
                            value = int(value)
                        except:
                            pass
                        stats[stat] = value

        return stats

    @classmethod
    def gather_all_stats(cls, input_dir):
        all_stats = {}

        for run_directory in os.listdir(input_dir):
            if not run_directory.endswith("fastqc"):
                continue
            report = os.path.join(input_dir, run_directory, "fastqc_data.txt")
            stats = Fastqc._stats_from_report(report)
            filename = stats.pop("filename")
            all_stats[filename] = stats

        return all_stats

    @staticmethod
    def _clean_outdir(outdir):
        assert os.path.exists(outdir)
        assert os.path.isdir(outdir)

        dirs = [
            os.path.join(outdir, x)
            for x in os.listdir(outdir)
            if os.path.isdir(os.path.join(outdir, x))
        ]

        for directory in dirs:
            for filename in [directory + ".zip", directory + ".html"]:
                try:
                    os.unlink(filename)
                except:
                    raise Error(
                        'Error deleting file "' + filename + '". Cannot continue'
                    )

            images_dir = os.path.join(directory, "Images")
            images = [x for x in os.listdir(images_dir) if x.endswith(".png")]

            for image in images:
                old = os.path.join(images_dir, image)
                new = os.path.join(directory, image)
                try:
                    os.rename(old, new)
                except:
                    raise Error("Error: mv " + old + " " + new)

            shutil.rmtree(images_dir)
            shutil.rmtree(os.path.join(directory, "Icons"))

            for filename in ["fastqc.fo", "fastqc_report.html", "summary.txt"]:
                full_filename = os.path.join(directory, filename)

                try:
                    os.unlink(full_filename)
                except:
                    raise Error("Error rm " + full_filename)

    def run(self):
        Fastqc._run_fastqc(self.outdir, self.infiles)
        Fastqc._clean_outdir(self.outdir)
        self.stats = Fastqc.gather_all_stats(self.outdir)
