import itertools
import multiprocessing
import os
import re
import requests
import shutil
import sys
from clockwork import spreadsheet_helper, utils


# https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html
run_pattern = re.compile("^[EDS]RR[0-9]{6,}$")


def _download(download_cmd):
    print(download_cmd)
    try:
        utils.syscall(download_cmd)
    except Exception:
        print("Error running: ", download_cmd)


def _download_run(run_id, outdir):
    expected_1 = os.path.join(outdir, run_id + "_1.fastq.gz")
    expected_2 = os.path.join(outdir, run_id + "_2.fastq.gz")
    if os.path.exists(expected_1) and os.path.exists(expected_2):
        print("Already got", outdir, " so do not download again")
        return

    cmd = " ".join(["enaDataGet", "-f fastq", "-d", outdir, run_id,])
    _download(cmd)


def _download_sample(sample_id, outdir):
    cmd = " ".join(["enaGroupGet", "-f fastq", "-d", outdir, sample_id,])
    _download(cmd)


def _download_run_or_sample(accession, outdir):
    print("Start downloading", accession)
    if run_pattern.match(accession):
        _download_run(accession, outdir)
    else:
        _download_sample(accession, outdir)
    print("End downloading", accession)


def _write_md5_file(infile):
    print("Start calculate md5 of", infile)
    with open(infile + ".md5", "w") as f:
        print(utils.md5(infile), os.path.basename(infile), sep="  ", file=f)
    print("End calculate md5 of", infile)


class EnaDownloader:
    def __init__(
        self,
        data_infile,
        outdir,
        site_id,
        lab_id,
        submission_date,
        dataset_name,
        download_threads=1,
        md5_threads=1,
    ):
        self.data_infile = os.path.abspath(data_infile)
        self.outdir = outdir
        self.site_id = site_id
        self.lab_id = lab_id
        self.submission_date = submission_date
        self.dataset_name = dataset_name
        self.download_threads = download_threads
        self.md5_threads = md5_threads

    @classmethod
    def _load_data_infile(cls, infile):
        data = {}
        seen_accessions = set()

        with open(infile) as f:
            for line in f:
                fields = line.rstrip().split("\t")
                if len(fields) != 2:
                    raise Exception(
                        "Expect two columns in each line. Error at this line:\n" + line
                    )
                if fields[0] in data:
                    raise Exception('Error. Non-unique sample name "' + fields[0] + '"')
                accessions = set(fields[1].split(","))
                accessions_seen_already = seen_accessions.intersection(accessions)
                if len(accessions_seen_already) > 0:
                    accessions_seen_already = sorted(list(accessions_seen_already))
                    raise Exception(
                        "Error. Non-unique acccession(s): "
                        + ",".join(accessions_seen_already)
                    )

                seen_accessions.update(accessions)
                data[fields[0]] = accessions

        return data

    @classmethod
    def _rename_files(cls, input_data, outdir):
        filename_data = {}

        for sample, accession_set in input_data.items():
            filename_data[sample] = set()

            for accession in sorted(list(accession_set)):
                # If we're rerunning, then the files might already be there, so
                # no files to rename. But do still need to update the filename_data dict
                expected_1 = os.path.join(outdir, accession + "_1.fastq.gz")
                expected_2 = os.path.join(outdir, accession + "_2.fastq.gz")
                print(expected_1, expected_2)
                if os.path.exists(expected_1) and os.path.exists(expected_2):
                    filename_data[sample].add(accession)
                    continue

                accession_dir = os.path.join(outdir, accession)
                if run_pattern.match(accession):
                    # new_run_accessions = [accession]
                    new_files = [
                        os.path.join(
                            accession_dir, accession + "_" + str(i) + ".fastq.gz"
                        )
                        for i in (1, 2)
                    ]
                    new_files = [x for x in new_files if os.path.exists(x)]
                    if len(new_files) > 0:
                        filename_data[sample].add(accession)
                    delete_dirs = set()
                else:
                    new_files = []
                    try:
                        new_run_accessions = os.listdir(os.path.join(outdir, accession))
                    except:
                        print("No files found for sample", accession)
                        new_run_accessions = []
                    delete_dirs = set(new_run_accessions)

                    for run_accession in new_run_accessions:
                        new_files.extend(
                            [
                                os.path.join(
                                    accession_dir,
                                    run_accession,
                                    run_accession + "_" + str(i) + ".fastq.gz",
                                )
                                for i in (1, 2)
                            ]
                        )

                    new_files = [x for x in new_files if os.path.exists(x)]
                    if len(new_files) > 0:
                        filename_data[sample].update(new_run_accessions)

                if len(new_files) == 0:
                    print("No paired files for accession", accession, file=sys.stderr)

                for filename in new_files:
                    basename = os.path.basename(filename)
                    new_name = os.path.join(outdir, basename)
                    if not os.path.exists(new_name):
                        print("rename", filename, new_name)
                        os.rename(filename, new_name)
                    else:
                        print("skip rename, new file found already", filename, new_name)

                for d in delete_dirs:
                    shutil.rmtree(os.path.join(accession_dir, d))

            for accession in sorted(list(accession_set)):
                # If there are unpaired reads, then those files will be left here.
                # We don't want them. Hence rmtree instead rmdir.
                # There may be no dir to delete, hence try/except pass
                try:
                    shutil.rmtree(os.path.join(outdir, accession))
                except:
                    pass

        filename_data = {
            x: filename_data[x] for x in filename_data if len(filename_data[x]) > 0
        }
        return filename_data

    @classmethod
    def _ena_run_to_sample_and_instrument_model(cls, run_id):
        wanted_fields = ["sample_accession", "instrument_model", "center_name"]
        url = "http://www.ebi.ac.uk/ena/portal/api/filereport?"
        data = {
            "accession": run_id,
            "result": "read_run",
            "fields": ",".join(wanted_fields),
        }

        try:
            r = requests.get(url, data)
        except:
            raise Exception(
                "Error querying ENA to get sample from run " + run_id + "\n" + r.url
            )

        if r.status_code != requests.codes.ok:
            raise Exception(
                "Error requesting data. Error code: "
                + str(r.status_code)
                + "\n"
                + r.url
            )

        # rstrip explicitly on newline because the center name can be empty, which
        # means the line ends in a tab character, which we want to keep so that
        # later when splitting on tab we end up with 3 fields, not 2.
        #
        # Also, even though we request fields:
        #    sample_accession,instrument_model,center_name
        # it also returns the run_accession. This change must have happened
        # in the last few months (today is 12/11/2020). So let's just look
        # for the fields we actually want and ignore everything else, so there's
        # a chance that if they change it again this code will still work.
        lines = r.text.rstrip("\n").split("\n")
        error_message = f"Unexpected format from ENA request {r.url}\nGot:\n{lines}"
        if len(lines) != 2:
            raise Exception(error_message)

        field_names = lines[0].split("\t")
        field_vals = lines[1].split("\t")
        if len(field_names) != len(field_vals):
            raise Exception(error_message)

        try:
            results = dict(zip(field_names, field_vals))
        except:
            raise Exception(error_message)

        for field in wanted_fields:
            if field not in results:
                raise Exepction(error_message)

        return [results[k] for k in wanted_fields]

    @classmethod
    def _write_import_tsv(
        cls,
        outfile,
        input_data,
        filename_data,
        site_id,
        lab_id,
        submission_date,
        dataset_name,
        test_ena_dict=None,
    ):
        with open(outfile, "w") as f:
            print(*spreadsheet_helper.columns, sep="\t", file=f)

            for sample in sorted(input_data):
                if sample not in filename_data:
                    print("Skipping", sample, "from final TSV file", outfile)
                    continue

                run_list = sorted(list(filename_data[sample]))
                for i, run in enumerate(run_list):
                    if test_ena_dict is not None:
                        ena_sample_id = test_ena_dict[run]
                        instrument_model = "Illumina HiSeq 2500"
                        center_name = "CENTER 42"
                    else:
                        metadata = EnaDownloader._ena_run_to_sample_and_instrument_model(
                            run
                        )
                        try:
                            ena_sample_id, instrument_model, center_name = metadata
                        except:
                            print(
                                "Error getting ENA metadata for run ",
                                run,
                                ". Got:",
                                metadata,
                                ". Skipping",
                                sep="",
                            )
                            continue

                    if center_name == "":
                        print("WARNING: empty center_name for run", run)
                        center_name = "UNKNOWN"

                    print(
                        sample,
                        site_id,
                        lab_id,
                        1,
                        i + 1,
                        submission_date,
                        run + "_1.fastq.gz",
                        0,
                        run + "_2.fastq.gz",
                        0,
                        dataset_name,
                        instrument_model,
                        center_name,
                        0,
                        0,
                        run,
                        ena_sample_id,
                        sep="\t",
                        file=f,
                    )

    def run(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        input_data = EnaDownloader._load_data_infile(self.data_infile)

        # Gather list of all accession to download, then run in parallel
        all_accessions = []
        for accession_set in input_data.values():
            all_accessions.extend(list(accession_set))
        all_accessions.sort()
        pool = multiprocessing.Pool(self.download_threads)
        pool.starmap(
            _download_run_or_sample, zip(all_accessions, itertools.repeat(self.outdir))
        )

        # Rename files made by ENA download script
        filename_data = EnaDownloader._rename_files(input_data, self.outdir)

        # Make md5 file of each fastq.gz file, in parallel
        fastq_files = []
        for run_set in filename_data.values():
            for run in run_set:
                for i in ("1", "2"):
                    filename = os.path.join(self.outdir, run + "_" + i + ".fastq.gz")
                    assert os.path.exists(filename)
                    fastq_files.append(filename)

        pool = multiprocessing.Pool(self.md5_threads)
        pool.starmap(_write_md5_file, zip(fastq_files))

        # write import TSV file
        tsv = os.path.join(self.outdir, "import.tsv")
        EnaDownloader._write_import_tsv(
            tsv,
            input_data,
            filename_data,
            self.site_id,
            self.lab_id,
            self.submission_date,
            self.dataset_name,
        )
