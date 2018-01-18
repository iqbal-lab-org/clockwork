import unittest
import filecmp
import os
import re
import pysam
import pyfastaq
from clockwork import contam_remover

modules_dir = os.path.dirname(os.path.abspath(contam_remover.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'contam_remover')


class TestContamRemover(unittest.TestCase):
    def test_load_metadata_file(self):
        '''test _load_metadata_file'''
        infile = os.path.join(data_dir, 'load_metadata_file.tsv')
        expected_by_group = {
            'Mtb': {'contam': False, 'sequences': {'Mtb_seq'}},
            'human': {'contam': True, 'sequences': {'human.chr1', 'human.chr2', 'human.chr3'}},
            'virus': {'contam': True, 'sequences': {'influenza', 'HIV', 'EBV'}},
            'other': {'contam': True, 'sequences': {'other.1', 'other.2', 'other.3'}},
        }

        expected_sequence_is_contam = {x: True for x in ('Mtb_seq', 'human.chr1', 'human.chr2', 'human.chr3', 'influenza', 'HIV', 'EBV', 'other.1', 'other.2', 'other.3')}
        expected_sequence_is_contam['Mtb_seq'] = False

        got_by_group, got_is_contam = contam_remover.ContamRemover._load_metadata_file(infile)
        self.assertEqual(expected_by_group, got_by_group)
        self.assertEqual(expected_sequence_is_contam, got_is_contam)


    def test_sam_to_fastq(self):
        '''test sam_to_fastq'''
        expected = [
        pyfastaq.sequences.Fastq('read1/1', 'GAACACGCGTCTCATCGCCGTATCACGTACATTAAAGGTAAATCGCTCCAGGTAATGAAT', 'IIIIIIIIIIIIIIIIIIIIIIIIIIGGGGHHHHHFFFFFEEEEDDDDDDDCCCCCCCCB'),
        pyfastaq.sequences.Fastq('read1/2', 'ACCCCTCTCTATGCGATAGACGTATTGGATGTTATATGACCATTGAGTCTAGGTCCTTCT', 'ABCDEFGHIIIIIIIIIIIIIHHHHIHIHIHIHIIIIIIIIIIIIIIIHHHHHGGGGGGG')
                                                                        ]
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'sam_to_fastq.sam'), 'r')
        i = 0
        for s in sam_reader.fetch(until_eof=True):
            self.assertEqual(expected[i], contam_remover.ContamRemover._sam_to_fastq(s))
            i += 1


    def test_read_mapped_and_wanted(self):
        '''test _read_mapped_and_wanted'''
        contam_dict = {'ref1': False, 'foo2': True}
        expected = [
            (True, True),
            (True, False),
            (False, False)
        ]

        sam_reader = pysam.Samfile(os.path.join(data_dir, 'read_mapped_and_wanted.sam'), 'r')

        i = 0
        for s in sam_reader.fetch(until_eof=True):
            self.assertEqual(expected[i], contam_remover.ContamRemover._read_mapped_and_wanted(s, contam_dict))
            i += 1


    def test_write_read_counts_by_group_file(self):
        '''test_write_read_counts_by_group_file'''
        sequences_by_group = {
            'group1': {'contam': False, 'sequences': {'seq1', 'seq2'}},
            'group2': {'contam': True, 'sequences': {'seq3'}},
            'group3': {'contam': True, 'sequences':  {'seq4', 'seq5'}},
        }

        read_counts = {
            'seq1': 10,
            'seq2': 20,
            'seq4': 42,
            'reads_kept_after_remove_contam': 32,
        }

        outfile = 'tmp.write_read_counts_by_group_file.out'
        contam_remover.ContamRemover._write_read_counts_by_group_file(sequences_by_group, read_counts, outfile)
        expected_file = os.path.join(data_dir, 'write_read_counts_by_group_file.expected.tsv')
        self.assertTrue(filecmp.cmp(expected_file, outfile, shallow=False))
        os.unlink(outfile)


    def test_run_no_match_not_wanted(self):
        '''test run when no_match_are_wanted=False'''
        metadata_file = os.path.join(data_dir, 'run.metadata.tsv')
        samfile = os.path.join(data_dir, 'run.sam')
        outprefix = 'tmp.contam_remover.test_run'
        counts_file = outprefix + '.counts.tsv'
        reads1 = outprefix + '.wanted_1.fq'
        reads2 = outprefix + '.wanted_2.fq'
        no_match1 = outprefix  + '.no_match_1.fq'
        no_match2 = outprefix  + '.no_match_2.fq'
        contam_1 = outprefix  + '.contam_1.fq'
        contam_2 = outprefix  + '.contam_2.fq'
        done_file = outprefix + '.done'
        if os.path.exists(done_file):
            os.unlink(done_file)
        cremover = contam_remover.ContamRemover(metadata_file, samfile, counts_file, reads1, reads2, contam_out_1=contam_1, contam_out_2=contam_2, no_match_out_1=no_match1, no_match_out_2=no_match2, done_file=done_file)

        cremover.run()
        for read_type in 'contam', 'no_match', 'wanted':
            for i in 1,2:
                got_file = outprefix + '.' + read_type + '_' + str(i) + '.fq'
                self.assertTrue(os.path.exists(got_file))
                expected_file = os.path.join(data_dir, 'run.no_match_not_wanted.' + read_type + '_' + str(i) + '.fq')
                self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
                os.unlink(got_file)

        expected_counts_file = os.path.join(data_dir, 'run.expected_counts.tsv')
        self.assertTrue(filecmp.cmp(expected_counts_file, counts_file, shallow=False))
        os.unlink(counts_file)
        self.assertTrue(os.path.exists(done_file))
        os.unlink(done_file)


    def test_run_no_match_not_wanted_discard_contam(self):
        '''test run when no_match_are_wanted=False and discard contam reads'''
        metadata_file = os.path.join(data_dir, 'run.metadata.tsv')
        samfile = os.path.join(data_dir, 'run.sam')
        outprefix = 'tmp.contam_remover.test_run'
        counts_file = outprefix + '.counts.tsv'
        reads1 = outprefix + '.wanted_1.fq'
        reads2 = outprefix + '.wanted_2.fq'
        no_match1 = outprefix  + '.no_match_1.fq'
        no_match2 = outprefix  + '.no_match_2.fq'
        done_file = outprefix + '.done'
        if os.path.exists(done_file):
            os.unlink(done_file)
        cremover = contam_remover.ContamRemover(metadata_file, samfile, counts_file, reads1, reads2, no_match_out_1=no_match1, no_match_out_2=no_match2, done_file=done_file)
        cremover.run()
        for read_type in 'no_match', 'wanted':
            for i in 1,2:
                got_file = outprefix + '.' + read_type + '_' + str(i) + '.fq'
                self.assertTrue(os.path.exists(got_file))
                expected_file = os.path.join(data_dir, 'run.no_match_not_wanted.' + read_type + '_' + str(i) + '.fq')
                self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
                os.unlink(got_file)

        expected_counts_file = os.path.join(data_dir, 'run.expected_counts.tsv')
        self.assertTrue(filecmp.cmp(expected_counts_file, counts_file, shallow=False))
        os.unlink(counts_file)
        self.assertTrue(os.path.exists(done_file))
        os.unlink(done_file)


    def test_run_no_match_are_wanted(self):
        '''test run when no_match_are_wanted=True'''
        metadata_file = os.path.join(data_dir, 'run.metadata.tsv')
        samfile = os.path.join(data_dir, 'run.sam')
        outprefix = 'tmp.contam_remover.test_run'
        counts_file = outprefix + '.counts.tsv'
        reads1 = outprefix + '.wanted_1.fq'
        reads2 = outprefix + '.wanted_2.fq'
        contam_1 = outprefix  + '.contam_1.fq'
        contam_2 = outprefix  + '.contam_2.fq'
        done_file = outprefix + '.done'
        if os.path.exists(done_file):
            os.unlink(done_file)
        cremover = contam_remover.ContamRemover(metadata_file, samfile, counts_file, reads1, reads2, contam_out_1=contam_1, contam_out_2=contam_2, done_file=done_file)
        cremover.run()
        for read_type in 'contam', 'wanted':
            for i in 1,2:
                got_file = outprefix + '.' + read_type + '_' + str(i) + '.fq'
                self.assertTrue(os.path.exists(got_file))
                expected_file = os.path.join(data_dir, 'run.no_match_are_wanted.' + read_type + '_' + str(i) + '.fq')
                self.assertTrue(filecmp.cmp(expected_file, got_file, shallow=False))
                os.unlink(got_file)

        expected_counts_file = os.path.join(data_dir, 'run.no_match_are_wanted.counts.tsv')
        self.assertTrue(filecmp.cmp(expected_counts_file, counts_file, shallow=False))
        os.unlink(counts_file)
        self.assertTrue(os.path.exists(done_file))
        os.unlink(done_file)

