import itertools
import multiprocessing
import os
import re
import requests
from clockwork import spreadsheet_importer, utils

class Error (Exception): pass

run_pattern = re.compile('^[EDS]RR[0-9]{6,7}$')


def _download_run(run_id, outdir):
    cmd = ' '.join([
        'enaDataGet',
        '-f fastq',
        '-d', outdir,
        run_id,
    ])
    print(cmd)
    utils.syscall(cmd)


def _download_sample(sample_id, outdir):
    cmd = ' '.join([
        'enaGroupGet',
        '-f fastq',
        '-d', outdir,
        sample_id,
    ])
    print(cmd)
    utils.syscall(cmd)


def _download_run_or_sample(accession, outdir):
    print('Start downloading', accession)
    if run_pattern.match(accession):
        _download_run(accession, outdir)
    else:
        _download_sample(accession, outdir)
    print('End downloading', accession)


def _write_md5_file(infile):
    print('Start calculate md5 of', infile)
    with open(infile + '.md5', 'w') as f:
        print(utils.md5(infile), os.path.basename(infile), sep='  ', file=f)
    print('End calculate md5 of', infile)


class EnaDownloader:
    def __init__(self, data_infile, outdir, site_id, lab_id, submission_date, dataset_name, download_threads=1, md5_threads=1):
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
                fields = line.rstrip().split('\t')
                if len(fields) != 2:
                    raise Error('Expect two columns in each line. Error at this line:\n' + line)
                if fields[0] in data:
                    raise Error('Error. Non-unique sample name "' + fields[0] + '"')
                accessions = set(fields[1].split(','))
                accessions_seen_already = seen_accessions.intersection(accessions)
                if len(accessions_seen_already) > 0:
                    accessions_seen_already = sorted(list(accessions_seen_already))
                    raise Error('Error. Non-unique acccession(s): ' + ','.join(accessions_seen_already))

                seen_accessions.update(accessions)
                data[fields[0]] = accessions

        return data


    @classmethod
    def _rename_files(cls, input_data, outdir):
        filename_data = {}

        for sample, accession_set in input_data.items():
            filename_data[sample] = set()

            for accession in sorted(list(accession_set)):
                accession_dir = os.path.join(outdir, accession)
                if run_pattern.match(accession):
                    new_run_accessions = [accession]
                    new_files = [os.path.join(accession_dir, accession + '_' + str(i) + '.fastq.gz') for i in (1,2)]
                    filename_data[sample].add(accession)
                    delete_dirs = set()
                else:
                    new_files = []
                    new_run_accessions = os.listdir(os.path.join(outdir, accession))
                    filename_data[sample].update(new_run_accessions)
                    delete_dirs = set(new_run_accessions)

                    for run_accession in new_run_accessions:
                        new_files.extend([os.path.join(accession_dir, run_accession, run_accession + '_' + str(i) + '.fastq.gz') for i in (1,2)])

                for filename in new_files:
                    basename = os.path.basename(filename)
                    new_name = os.path.join(outdir, basename)
                    assert not os.path.exists(new_name)
                    print('rename', filename, new_name)
                    os.rename(filename, new_name)

                for d in delete_dirs:
                    os.rmdir(os.path.join(accession_dir, d))

            for accession in sorted(list(accession_set)):
                os.rmdir(os.path.join(outdir, accession))

        assert input_data.keys() == filename_data.keys()
        return filename_data


    @classmethod
    def _ena_run_to_sample_and_instrument_model(cls, run_id):
        url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?'
        data = {'accession': run_id, 'result': 'read_run', 'download': 'txt', 'fields': 'sample_accession,instrument_model,center_name'}

        try:
            r = requests.get(url, data)
        except:
            raise Error('Error querying ENA to get sample from run ' + run_id + '\n' + r.url)

        if r.status_code != requests.codes.ok:
            raise Error('Error requesting data. Error code: ' + str(r.status_code) + '\n' + r.url)

        lines = r.text.rstrip().split('\n')
        if len(lines) != 2 or lines[0] != 'sample_accession\tinstrument_model\tcenter_name':
            raise Error('Unexpected format from ENA request ' + r.url + '\nGot:' + str(lines))

        return lines[1].split('\t')


    @classmethod
    def _write_import_tsv(cls, outfile, input_data, filename_data, site_id, lab_id, submission_date, dataset_name, test_ena_dict=None):
        with open(outfile, 'w') as f:
            print(*spreadsheet_importer.columns, sep='\t', file=f)

            for sample in sorted(input_data):
                assert sample in filename_data
                run_list = sorted(list(filename_data[sample]))
                for i, run in enumerate(run_list):
                    if test_ena_dict is not None:
                        ena_sample_id = test_ena_dict[run]
                        instrument_model = 'Illumina HiSeq 2500'
                        center_name = 'CENTER 42'
                    else:
                        ena_sample_id, instrument_model, center_name = EnaDownloader._ena_run_to_sample_and_instrument_model(run)

                    print(sample, site_id, lab_id, 1, i + 1, submission_date,
                        run + '_1.fastq.gz', 0, run + '_2.fastq.gz', 0,
                        dataset_name, instrument_model, center_name, 0, 0, run, ena_sample_id, sep='\t', file=f)


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
        pool.starmap(_download_run_or_sample, zip(all_accessions, itertools.repeat(self.outdir)))

        # Rename files made by ENA download script
        filename_data = EnaDownloader._rename_files(input_data, self.outdir)

        # Make md5 file of each fastq.gz file, in parallel
        fastq_files = []
        for run_set in filename_data.values():
            for run in run_set:
                for i in ('1', '2'):
                    filename = os.path.join(self.outdir, run + '_' + i + '.fastq.gz')
                    assert os.path.exists(filename)
                    fastq_files.append(filename)

        pool = multiprocessing.Pool(self.md5_threads)
        pool.starmap(_write_md5_file, zip(fastq_files))

        # write import TSV file
        tsv = os.path.join(self.outdir, 'import.tsv')
        EnaDownloader._write_import_tsv(tsv, input_data, filename_data, self.site_id,
            self.lab_id, self.submission_date, self.dataset_name)

