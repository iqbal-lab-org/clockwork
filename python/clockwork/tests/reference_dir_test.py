import filecmp
import unittest
import os
import shutil
from clockwork import reference_dir

modules_dir = os.path.dirname(os.path.abspath(reference_dir.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data', 'reference_dir')


class TestReferenceDir(unittest.TestCase):
    def test_init(self):
        '''test __init__'''
        with self.assertRaises(reference_dir.Error):
            refdir = reference_dir.ReferenceDir(pipeline_references_root_dir='foo')
        with self.assertRaises(reference_dir.Error):
            refdir = reference_dir.ReferenceDir(reference_id=42)

        refdir = reference_dir.ReferenceDir(pipeline_references_root_dir='foo', reference_id=42)
        self.assertEqual(refdir.directory, os.path.join(os.getcwd(), 'foo', '42'))

        refdir = reference_dir.ReferenceDir(directory='bar')
        self.assertEqual(refdir.directory, os.path.join(os.getcwd(), 'bar'))


    def test_make_index_files(self):
        '''test make_index_files'''
        tmp_root_dir = 'tmp.reference_dir.make_index_files'
        if os.path.exists(tmp_root_dir):
            shutil.rmtree(tmp_root_dir)
        fasta_in = os.path.join(data_dir, 'make_index_files.ref.in.fa.gz')
        expected_ref = os.path.join(data_dir, 'make_index_files.ref.expected.fa')
        ref_dir = reference_dir.ReferenceDir(pipeline_references_root_dir=tmp_root_dir, reference_id=42)
        with self.assertRaises(reference_dir.Error):
            ref_dir.make_index_files('file_does_not_exist', False, True, cortex_mem_height=17)

        ref_dir.make_index_files(fasta_in, False, True, cortex_mem_height=17)
        self.assertTrue(os.path.exists(ref_dir.directory))
        self.assertTrue(os.path.exists(ref_dir.ref_fasta))
        self.assertTrue(filecmp.cmp(ref_dir.ref_fasta, expected_ref, shallow=False))
        self.assertTrue(os.path.exists(ref_dir.ref_fai))
        self.assertTrue(os.path.exists(ref_dir.ref_fasta + '.bwt'))
        self.assertTrue(os.path.exists(ref_dir.ref_fasta_prefix + '.stampy.sthash'))
        self.assertTrue(os.path.exists(ref_dir.ref_fasta_prefix + '.stampy.stidx'))
        self.assertTrue(os.path.exists(ref_dir.ref_fasta_prefix + '.k31.ctx'))
        shutil.rmtree(tmp_root_dir)


    def test_add_remove_contam_metadata_tsv(self):
        '''test add_remove_contam_metadata_tsv'''
        tmp_root_dir = 'tmp.reference_dir.add_remove_contam_metadata_tsv'
        if os.path.exists(tmp_root_dir):
            shutil.rmtree(tmp_root_dir)
        fasta_in =  os.path.join(data_dir, 'add_remove_contam_metadata_tsv.ref.fa')
        ref_dir = reference_dir.ReferenceDir(pipeline_references_root_dir=tmp_root_dir, reference_id=42)
        ref_dir.make_index_files(fasta_in, False, False, cortex_mem_height=17)

        bad_files = [
            os.path.join(data_dir, 'add_remove_contam_metadata_tsv.ref.missing_from_tsv.tsv'),
            os.path.join(data_dir, 'add_remove_contam_metadata_tsv.ref.extra_in_tsv.tsv'),
        ]

        for bad_file in bad_files:
            with self.assertRaises(reference_dir.Error):
                ref_dir.add_remove_contam_metadata_tsv(bad_file)

        ref_dir.add_remove_contam_metadata_tsv(os.path.join(data_dir, 'add_remove_contam_metadata_tsv.ref.good.tsv'))
        shutil.rmtree(tmp_root_dir)
