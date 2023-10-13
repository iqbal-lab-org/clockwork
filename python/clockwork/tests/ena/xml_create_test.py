import unittest
import filecmp
import os
import xmltodict
import xml.etree.ElementTree as ET
from clockwork.ena import xml_create

modules_dir = os.path.join(os.path.dirname(os.path.abspath(xml_create.__file__)), os.pardir)
data_dir = os.path.normpath(os.path.join(modules_dir, 'tests', 'data', 'ena', 'xml_create'))


def xml_files_the_same(file1, file2):
    with open(file1, "rb") as f:
        contents1 = xmltodict.parse(f)
    with open(file2, "rb") as f:
        contents2 = xmltodict.parse(f)
    return contents1 == contents2

class TestXmlCreate(unittest.TestCase):
    def test_element_tree_to_file(self):
        '''test element_tree_to_file'''
        tmp_file = 'tmp.ena.xml_create.element_tree_to_file.xml'
        root = ET.Element('ROOT')
        sub_element_1 = ET.SubElement(root, 'SUB_ELEMENT1', {'foo': 'bar'})
        sub_element_2 = ET.SubElement(root, 'SUB_ELEMENT2', {'knight': 'ni!'})
        xml_create.element_tree_to_file(root, tmp_file)
        expected_file = os.path.join(data_dir, 'element_tree_to_file.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)


    def test_make_add_submission_xml(self):
        '''test make_add_submission_xml'''
        tmp_file = 'tmp.ena.xml_create.make_add_submission_xml.xml'
        xml_create.make_add_submission_xml(tmp_file, 'alias 42', 'Royston Vasey', 'source.xml', 'schema_name', hold_for_two_years=False)
        expected_file = os.path.join(data_dir, 'make_add_submission_xml.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)

        xml_create.make_add_submission_xml(tmp_file, 'alias 42', 'Royston Vasey', 'source.xml', 'schema_name', broker_name='Broker 11', hold_for_two_years=False)
        expected_file = os.path.join(data_dir, 'make_add_submission_xml.with_broker.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)


    def test_make_project_xml(self):
        '''test make_project_xml'''
        tmp_file = 'tmp.ena.xml_create.make_project_xml.xml'
        xml_create.make_project_xml(tmp_file, 'project_alias', 'center name', 'title text', 'description text')
        expected_file = os.path.join(data_dir, 'make_project_xml.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)


    def test_make_sample_xml(self):
        '''test make_sample_xml'''
        tmp_file = 'tmp.ena.xml_create.make_sample_xml.xml'
        xml_create.make_sample_xml(tmp_file, 'sample alias', 'center name', 'title', 42, {'foo': 'bar'})
        expected_file = os.path.join(data_dir, 'make_sample_xml.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)


    def test_make_experiment_xml(self):
        '''test make_experiment_xml'''
        tmp_file = 'tmp.ena.xml_create.make_experiment_xml.xml'
        xml_create.make_experiment_xml(tmp_file, 'experiment alias', 'center name', 'title', 'study_acc', 'sample_acc', 'library name', 'platform', 'instrument', {'spam': 'eggs'})
        expected_file = os.path.join(data_dir, 'make_experiment_xml.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)


    def test_make_paired_fastq_run_xml(self):
        '''test make_paired_fastq_run_xml'''
        tmp_file = 'tmp.ena_xml_create.make_paired_fastq_run_xml.xml'
        xml_create.make_paired_fastq_run_xml(tmp_file, 'run alias', 'center name', 'experiment_acc', 'reads1.fq', 'md51', 'reads2.gq', 'md52')
        expected_file = os.path.join(data_dir, 'make_paired_fastq_run_xml.xml')
        self.assertTrue(xml_files_the_same(expected_file, tmp_file))
        os.unlink(tmp_file)
