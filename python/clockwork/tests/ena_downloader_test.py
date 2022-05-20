import unittest
import copy
import shutil
import filecmp
import os
import re
from clockwork import ena_downloader

modules_dir = os.path.dirname(os.path.abspath(ena_downloader.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "ena_downloader")


class TestEnaDownloader(unittest.TestCase):
    def test_load_data_infile(self):
        """Test _load_data_infile"""
        infile_bad_columns = os.path.join(data_dir, "load_data_infile.bad_columns.tsv")
        with self.assertRaises(Exception):
            ena_downloader.EnaDownloader._load_data_infile(infile_bad_columns)

        infile_repeat_sample_name = os.path.join(
            data_dir, "load_data_infile.repeat_sample_name.tsv"
        )
        with self.assertRaises(Exception):
            ena_downloader.EnaDownloader._load_data_infile(infile_repeat_sample_name)

        infile_repeat_accession = os.path.join(
            data_dir, "load_data_infile.repeat_accesion.tsv"
        )
        with self.assertRaises(Exception):
            ena_downloader.EnaDownloader._load_data_infile(infile_repeat_accession)

        infile_good = os.path.join(data_dir, "load_data_infile.good.tsv")
        expected = {
            "sample1": {"ERR42"},
            "sample2": {"ERS43"},
            "sample3": {"ERR44", "ERR45"},
        }
        got = ena_downloader.EnaDownloader._load_data_infile(infile_good)
        self.assertEqual(expected, got)

    def test_rename_files(self):
        """Test _rename_files"""
        indir = "tmp.ena_doanloader.test_rename_files"
        shutil.copytree(os.path.join(data_dir, "rename_files"), indir)
        input_data = ena_downloader.EnaDownloader._load_data_infile(
            os.path.join(data_dir, "rename_files.data.tsv")
        )
        got_dict = ena_downloader.EnaDownloader._rename_files(input_data, indir)
        expected_dict = {
            "id1": {"ERR0000012", "ERR0000013", "ERR0000014"},
            "id2": {"ERR0000050"},
            "id3": {"ERR0001000", "ERR0001001"},
        }
        self.assertEqual(expected_dict, got_dict)

        runs = []
        for id_set in expected_dict.values():
            runs.extend(list(id_set))

        expected_file_list = sorted(
            [x + "_" + str(i) + ".fastq.gz" for i in (1, 2) for x in runs]
        )
        got_file_list = sorted(os.listdir(indir))
        self.assertEqual(expected_file_list, got_file_list)

        # Rerun and check no errors. The code can't handle when a rerun includes
        # samples. They end up getting ignored if they were already downloaded.
        # The documentation of the main script says to use run accessions anyway, not
        # sample accessions. Test rerunning here, but have to remove the data that
        # arose from sample IDs in the input file.
        got_dict = ena_downloader.EnaDownloader._rename_files(input_data, indir)
        del expected_dict["id3"]
        expected_dict["id1"].remove("ERR0000014")
        self.assertEqual(expected_dict, got_dict)
        shutil.rmtree(indir)

    def test_ena_run_to_sample_and_instrument_model(self):
        """test _ena_run_to_sample_and_instrument_model"""
        run_accession = "ERR550803"
        (
            got_sample,
            got_instrument,
            got_center_name,
        ) = ena_downloader.EnaDownloader._ena_run_to_sample_and_instrument_model(
            run_accession
        )
        self.assertEqual("SAMEA2533482", got_sample)
        self.assertEqual("Illumina HiSeq 2500", got_instrument)
        self.assertEqual("FZB", got_center_name)

    def test_write_import_tsv(self):
        """Test _write_import_tsv"""
        input_data = {
            "id1": {"ERS123"},
            "id2": {"ERR1", "ERR2"},
        }

        runs = {
            "id1": {"ERR10", "ERR11"},
            "id2": {"ERR1", "ERR2"},
        }

        outfile = "tmp.ena_downloader.write_import_tsv.out"
        test_ena_dict = {x: x + "_sample" for x in ["ERR1", "ERR2", "ERR10", "ERR11"]}
        ena_downloader.EnaDownloader._write_import_tsv(
            outfile,
            input_data,
            runs,
            "site1",
            "lab1",
            "20170101",
            "group1",
            test_ena_dict=test_ena_dict,
        )
        expected_file = os.path.join(data_dir, "write_import_tsv.expected.tsv")
        self.assertTrue(filecmp.cmp(expected_file, outfile, shallow=False))
        os.unlink(outfile)
