import unittest
import os
import shutil
from clockwork.ena import object_creator, xml_create

modules_dir = os.path.join(os.path.dirname(os.path.abspath(object_creator.__file__)), os.pardir)
data_dir = os.path.normpath(os.path.join(modules_dir, 'tests', 'data', 'ena', 'object_creator'))


class TestObjectCreator(unittest.TestCase):
    def test_init_bad_parameters(self):
        '''test init with bad parameters'''
        with self.assertRaises(object_creator.Error):
            object_creator.ObjectCreator('ini_file', 'not_a_project', 'obj.xml', 'obj_alias', 'sub_alias', 'center 42', 'title')

        with self.assertRaises(object_creator.Error):
            # missing project_description
            object_creator.ObjectCreator('ini_file', 'project', 'obj.xml', 'obj_alias', 'sub_alias', 'center 42', 'title')

        with self.assertRaises(object_creator.Error):
            # missing taxon_id
            object_creator.ObjectCreator('ini_file', 'sample', 'obj.xml', 'obj_alias', 'sub_alias', 'center 42', 'title')

        with self.assertRaises(object_creator.Error):
            # missing study_accession, sample_accession, library_name, platform, instrument
            object_creator.ObjectCreator('ini_file', 'experiment', 'obj.xml', 'obj_alias', 'sub_alias', 'center 42', 'title')

        with self.assertRaises(object_creator.Error):
            # missing experiment_accession, reads_1, md5_1, reads_2, md5_2
            object_creator.ObjectCreator('ini_file', 'run', 'obj.xml', 'obj_alias', 'sub_alias', 'center 42', 'title')


    def test_run_project(self):
        '''test run making project'''
        obj_xml = 'tmp.object_creator.project.obj.xml'
        ini_file = os.path.join(data_dir, 'conf.ini')
        obj = object_creator.ObjectCreator(ini_file, 'project', obj_xml, 'objct alias', 'sub alias', 'center 42', 'title', project_description='project description', unit_test='success')
        obj.run()
        self.assertTrue(obj.submission_receipt.successful)
        os.unlink(obj_xml)
        self.assertTrue(os.path.exists(obj.submission_xml))
        os.unlink(obj.submission_xml)
        self.assertTrue(os.path.exists(obj.receipt_xml))
        os.unlink(obj.receipt_xml)


    def test_run_sample(self):
        '''test run making sample'''
        obj_xml = 'tmp.object_creator.sample.obj.xml'
        ini_file = os.path.join(data_dir, 'conf.ini')
        obj = object_creator.ObjectCreator(ini_file, 'sample', obj_xml, 'objct alias', 'sub alias', 'center 42', 'title', taxon_id=42, unit_test='success')
        obj.run()
        self.assertTrue(obj.submission_receipt.successful)
        os.unlink(obj_xml)
        self.assertTrue(os.path.exists(obj.submission_xml))
        os.unlink(obj.submission_xml)
        self.assertTrue(os.path.exists(obj.receipt_xml))
        os.unlink(obj.receipt_xml)


    def test_run_experiment(self):
        '''test run making experiment'''
        obj_xml = 'tmp.object_creator.experiment.obj.xml'
        ini_file = os.path.join(data_dir, 'conf.ini')
        obj = object_creator.ObjectCreator(ini_file, 'experiment', obj_xml, 'objct alias', 'sub alias', 'center 42', 'title', study_accession='ERP123', sample_accession='ERS42', library_name='lib name', platform='ILLUMINA', instrument='HISEQ', unit_test='success')
        obj.run()
        self.assertTrue(obj.submission_receipt.successful)
        os.unlink(obj_xml)
        self.assertTrue(os.path.exists(obj.submission_xml))
        os.unlink(obj.submission_xml)
        self.assertTrue(os.path.exists(obj.receipt_xml))
        os.unlink(obj.receipt_xml)


    def test_run_run(self):
        '''test run making run'''
        obj_xml = 'tmp.object_creator.run.obj.xml'
        ini_file = os.path.join(data_dir, 'conf.ini')
        obj = object_creator.ObjectCreator(ini_file, 'run', obj_xml, 'object alias', 'sub alias', 'center 42', 'title', experiment_accession='ERX123', reads_1='reads1.fq', md5_1='md51', reads_2='reads2.fq', md5_2='md52', unit_test='success')
        obj.run()
        self.assertTrue(obj.submission_receipt.successful)
        os.unlink(obj_xml)
        self.assertTrue(os.path.exists(obj.submission_xml))
        os.unlink(obj.submission_xml)
        self.assertTrue(os.path.exists(obj.receipt_xml))
        os.unlink(obj.receipt_xml)

