import os
import glob
from clockwork import het_snp_caller, read_map, utils


class Error(Exception):
    pass


class SamtoolsQc:
    def __init__(self, ref_fasta, reads1, reads2, outdir):
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.reads1 = os.path.abspath(reads1)
        self.reads2 = os.path.abspath(reads2)
        self.outdir = os.path.abspath(outdir)

        if os.path.exists(self.outdir):
            raise Error("Error! Output directory already exists " + self.outdir)

    @classmethod
    def _map_reads(cls, ref_fasta, reads1, reads2, outfile):
        read_map.map_reads(ref_fasta, reads1, reads2, outfile, markdup=True)

    @classmethod
    def _make_depth_stats(cls, samfile, outprefix):
        depth_file = outprefix + ".depths"
        cmd = " ".join(["samtools depth", "-a", samfile, ">", depth_file])
        utils.syscall(cmd)

    @classmethod
    def depth_stats(cls, filename):
        depths = {
            "eq_0": 0,
            "atleast_2": 0,
            "atleast_5": 0,
            "atleast_10": 0,
            "atleast_20": 0,
            "atleast_100": 0,
        }

        with open(filename) as f:
            for line in f:
                _, _, depth = line.rstrip().split("\t")
                depth = int(depth)
                if depth == 0:
                    depths["eq_0"] += 1
                if depth >= 2:
                    depths["atleast_2"] += 1
                if depth >= 5:
                    depths["atleast_5"] += 1
                if depth >= 10:
                    depths["atleast_10"] += 1
                if depth >= 20:
                    depths["atleast_20"] += 1
                if depth >= 100:
                    depths["atleast_100"] += 1
        return depths

    @classmethod
    def _make_stats_and_plots(cls, samfile, ref_fasta, outprefix):
        stats_file = outprefix + ".stats"

        cmd = " ".join(["samtools stats", "-r", ref_fasta, samfile, ">", stats_file])
        utils.syscall(cmd)

        cmd = " ".join(["plot-bamstats", "-p", outprefix + ".plot", stats_file])
        utils.syscall(cmd)

        for filename in glob.glob(outprefix + "*"):
            if filename.endswith(".gp") or filename.endswith(".html"):
                os.unlink(filename)

    @classmethod
    def stats_from_report(cls, filename):
        stats = {}
        wanted_keys = {
            "raw total sequences:": False,
            "reads mapped:": False,
            "reads duplicated:": False,
            "bases mapped (cigar):": False,
            "bases trimmed:": False,
            "error rate:": True,
            "average quality:": True,
            "insert size average:": True,
            "insert size standard deviation:": True,
            "inward oriented pairs:": False,
            "outward oriented pairs:": False,
            "pairs with other orientation:": False,
        }

        with open(filename) as f:
            for line in f:
                if len(stats) > 0 and not line.startswith("SN"):
                    break

                if line.startswith("SN"):
                    fields = line.rstrip().split("\t")
                    if fields[1] in wanted_keys:
                        key = (
                            fields[1]
                            .replace(" ", "_")
                            .rstrip(":")
                            .replace("(", "")
                            .replace(")", "")
                        )
                        value = (
                            float(fields[2])
                            if wanted_keys[fields[1]]
                            else int(fields[2])
                        )

                        stats[key] = value

        return stats

    @classmethod
    def het_snp_stats_from_summary_file(cls, filename):
        with open(filename) as f:
            lines = [x.rstrip().split("\t") for x in f.readlines()]
            assert len(lines) == 2
            stats = {}
            for i, key in enumerate(lines[0]):
                if key == "Percent_SNPs_are_het":
                    stats[key] = float(lines[1][i])
                else:
                    stats[key] = int(lines[1][i])

        return stats

    def run(self):
        try:
            os.mkdir(self.outdir)
        except:
            raise Error("Error making output directory " + self.outdir)

        outprefix = os.path.join(self.outdir, "samtools_qc")
        samfile = os.path.join(self.outdir, "tmp.sam")
        SamtoolsQc._map_reads(self.ref_fasta, self.reads1, self.reads2, samfile)
        SamtoolsQc._make_depth_stats(samfile, outprefix)
        SamtoolsQc._make_stats_and_plots(samfile, self.ref_fasta, outprefix)
        hsc = het_snp_caller.HetSnpCaller(
            samfile, self.ref_fasta, os.path.join(self.outdir, "het_snps")
        )
        hsc.run()
        self.stats = SamtoolsQc.stats_from_report(outprefix + ".stats")
        os.unlink(samfile)
