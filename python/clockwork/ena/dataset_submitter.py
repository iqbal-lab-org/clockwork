import configparser
import itertools
import multiprocessing
import os
from clockwork import db, isolate_dir
from clockwork.ena import object_creator, submission_receipt, submit_files, xml_create

class Error (Exception): pass

def _upload_fastq_file_pair(names_tuple, ini_file, unit_test=None):
    seqrep_id, to_upload_1, uploaded_1, to_upload_2, uploaded_2 = names_tuple
    if unit_test is not None:
        return (seqrep_id, unit_test == 'success')

    try:
        submit_files.upload_file_to_ena_ftp(ini_file, to_upload_1, uploaded_1)
        submit_files.upload_file_to_ena_ftp(ini_file, to_upload_2, uploaded_2)
    except:
        return (seqrep_id, False)

    return (seqrep_id, True)


class DatasetSubmitter:
    def __init__(self, ini_file, dataset_name, pipeline_root, taxon_id, fq_upload_threads=1, use_test_server=False, unit_test=None):
        self.ini_file = os.path.abspath(ini_file)
        self.db = db.Db(self.ini_file)
        self.dataset_name = dataset_name
        self.pipeline_root = os.path.abspath(pipeline_root)
        self.taxon_id = taxon_id
        self.fq_upload_threads = fq_upload_threads
        self.use_test_server = use_test_server
        self.unit_test = unit_test
        self.project_xml_dir = DatasetSubmitter.dataset_xml_dir(self.pipeline_root)
        self.project_xml = DatasetSubmitter.dataset_xml_file(self.pipeline_root, self.dataset_name)
        self.centre_number_to_name = DatasetSubmitter._get_centres_from_ini_file(self.ini_file)
        self.broker_name = DatasetSubmitter._get_broker_name_from_ini_file(self.ini_file)
        self.study_prefix = DatasetSubmitter._get_key_from_ini_file(self.ini_file, 'ena_login', 'study_prefix')
        if self.study_prefix is None:
            raise Error('Error! Must provide study_prefix in [ena_login] section of ini file ' + self.ini_file)


    @classmethod
    def dataset_xml_dir(cls, pipeline_root):
        return os.path.join(os.path.abspath(pipeline_root), 'Project_xmls')


    @classmethod
    def dataset_xml_file(cls, pipeline_root, dataset_name):
        return os.path.join(DatasetSubmitter.dataset_xml_dir(pipeline_root), dataset_name + '.dataset_submission.xml')


    def _get_data_from_db(self):
        columns = ','.join([
            'Seqrep.seqrep_id',
            'sequence_replicate_number',
            'remove_contam_reads_file_1_md5',
            'remove_contam_reads_file_2_md5',
            'ena_run_accession',
            'Isolate.isolate_id',
            'isolate_number_from_lab',
            'ena_experiment_accession',
            'Sample.sample_id',
            'subject_id',
            'ena_center_name',
            'site_id',
            'sample_id_from_lab',
            'instrument_model',
            'ena_sample_accession',
            'ena_study_accession',
        ])
        join = 'Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id JOIN Sample ON Isolate.sample_id = Sample.sample_id'
        where = 'submit_to_ena=1 AND import_status=1 AND ena_run_accession IS NULL AND dataset_name="' + self.dataset_name + '"'
        query = ' '.join(['SELECT', columns, 'FROM', '(' + join + ')', 'WHERE', '(' + where + ')'])
        return self.db.query_to_dict(query)


    @classmethod
    def _get_centres_from_ini_file(cls, ini_file):
        config = configparser.ConfigParser()
        try:
            config.read(ini_file)
        except:
            raise Error('Error reading config file ' + ini_file)

        if 'sequencing_centres' not in config:
            return {}
        else:
            return config['sequencing_centres']


    @classmethod
    def _get_key_from_ini_file(cls, ini_file, section, key):
        config = configparser.ConfigParser()
        try:
            config.read(ini_file)
        except:
            raise Error('Error reading config file ' + ini_file)

        if section not in config:
            return None
        else:
            return config[section].get(key, None)


    @classmethod
    def _get_broker_name_from_ini_file(cls, ini_file):
        return DatasetSubmitter._get_key_from_ini_file(ini_file, 'ena_login', 'broker_name')


    @classmethod
    def _ena_center_name_from_db_data(cls, data_in, number_to_name_dict=None):
        if number_to_name_dict is None:
            number_to_name_dict = {}

        center_names = {x['ena_center_name'] for x in data_in}
        if len(center_names) == 1:
            center_name = center_names.pop()
        else:
            raise Error('Error getting unique ena_center_name from: ' + str(data_in))

        return number_to_name_dict.get(center_name, center_name)


    def _submit_study_object(self, data_in):
        if not os.path.exists(self.project_xml_dir):
            os.mkdir(self.project_xml_dir)

        study_accessions_from_db = {x['ena_study_accession'] for x in data_in}
        if len(study_accessions_from_db) > 1:
            raise Error('Error! More than one study ID found for dataset ' + self.dataset_name + '. Got: ' + str(study_accessions_from_db))

        if not os.path.exists(self.project_xml):
            assert study_accessions_from_db == {None}
            project_alias = 'project.' + self.dataset_name
            submit_alias = 'submit.' + project_alias
            center_name = DatasetSubmitter._ena_center_name_from_db_data(data_in, number_to_name_dict=self.centre_number_to_name)
            title = self.study_prefix + '. ' + center_name + '. ' + self.dataset_name
            project_description = title
            project_creator = object_creator.ObjectCreator(self.ini_file, 'project',
              self.project_xml, project_alias, submit_alias, center_name, title,
              project_description, use_test_server=self.use_test_server,
              unit_test=self.unit_test, broker_name=self.broker_name
            )
            project_creator.run()
            if not project_creator.submission_receipt.successful:
                raise Error('Error submitting project to ena. XML file: ' + self.project_xml)

            ena_study_accession = project_creator.submission_receipt.accessions.get('PROJECT', None)
            if ena_study_accession is None:
                raise Error('Error getting proejct accession from ' + project_creator.receipt_xml)

            for row in data_in:
                row['ena_study_accession'] = ena_study_accession
                self.db.update_row('Sample', {'sample_id': row['sample_id']}, {'ena_study_accession': ena_study_accession})
            self.db.commit()
        else:
            assert len(study_accessions_from_db) == 1


    def _submit_sample_objects(self, data_in):
        submitted_samples = {} # sample id -> ena accession

        for row in data_in:
            if row['ena_sample_accession'] is not None:
                continue
            elif row['sample_id'] in submitted_samples:
                row['ena_sample_accession'] = submitted_samples[row['sample_id']]
            else:
                assert row['ena_sample_accession'] is None
                iso_dir = isolate_dir.IsolateDir(self.pipeline_root, row['sample_id'], row['isolate_id'])
                object_xml = iso_dir.xml_submission_file('sample')
                object_alias = 'sample.' + str(row['sample_id'])
                submit_alias = 'submit.' + object_alias
                center_name = DatasetSubmitter._ena_center_name_from_db_data(data_in, number_to_name_dict=self.centre_number_to_name)
                title = row['subject_id'] + '. ' + center_name + '. ' + row['sample_id_from_lab']
                obj_creator = object_creator.ObjectCreator(self.ini_file, 'sample', object_xml, object_alias, submit_alias, center_name, title, taxon_id=self.taxon_id, use_test_server=self.use_test_server, unit_test=self.unit_test, broker_name=self.broker_name)
                obj_creator.run()
                if obj_creator.submission_receipt.successful:
                    try:
                        sample_accession = obj_creator.submission_receipt.accessions['SAMPLE']
                    except:
                        sample_accession = 'FAIL'
                else:
                    sample_accession = 'FAIL'

                row['ena_sample_accession'] = sample_accession
                self.db.update_row('Sample', {'sample_id': row['sample_id']}, {'ena_sample_accession': sample_accession})
                self.db.commit()
                submitted_samples[row['sample_id']] = sample_accession


    def _submit_experiment_objects(self, data_in):
        submitted_isolates = {} # isolate id -> ena accession

        for row in data_in:
            if row['ena_experiment_accession'] is not None or row['ena_sample_accession'] == 'FAIL':
                continue
            elif row['isolate_id'] in submitted_isolates:
                row['ena_experiment_accession'] = submitted_isolates[row['isolate_id']]
            else:
                assert row['ena_experiment_accession'] is None
                iso_dir = isolate_dir.IsolateDir(self.pipeline_root, row['sample_id'], row['isolate_id'])
                object_xml = iso_dir.xml_submission_file('experiment')
                object_alias = 'experiment.' + str(row['isolate_id'])
                submit_alias = 'submit.' + object_alias
                center_name = DatasetSubmitter._ena_center_name_from_db_data(data_in, number_to_name_dict=self.centre_number_to_name)
                title = row['subject_id'] + '. ' + center_name + '. ' + row['sample_id_from_lab'] + '. ' + row['isolate_number_from_lab']
                library_name = title
                obj_creator = object_creator.ObjectCreator(self.ini_file, 'experiment', object_xml, object_alias,
                    submit_alias, center_name, title,
                    study_accession=row['ena_study_accession'],
                    sample_accession=row['ena_sample_accession'],
                    library_name=library_name,
                    platform='ILLUMINA',
                    instrument=row['instrument_model'],
                    use_test_server=self.use_test_server, unit_test=self.unit_test,
                    broker_name=self.broker_name,
                )
                obj_creator.run()
                if obj_creator.submission_receipt.successful:
                    try:
                        experiment_accession = obj_creator.submission_receipt.accessions['EXPERIMENT']
                    except:
                        experiment_accession = 'FAIL'
                else:
                    experiment_accession = 'FAIL'

                row['ena_experiment_accession'] = experiment_accession
                self.db.update_row('Isolate', {'isolate_id': row['isolate_id']}, {'ena_experiment_accession': experiment_accession})
                self.db.commit()
                submitted_isolates[row['isolate_id']] = experiment_accession


    def _submit_runs(self, data_in):
        # Note: reads have to be in the dropbox before submitting the Run object.
        # Upload all the reads first, in parallel, then submit the runs.
        fq_pairs_to_upload = []  # (seqrep, full path on disk1, dropbox name1, full path on disk2, dropbox name2)
        for row in data_in:
            iso_dir = isolate_dir.IsolateDir(self.pipeline_root, row['sample_id'], row['isolate_id'])
            fq_pairs_to_upload.append((
                row['seqrep_id'],
                iso_dir.reads_filename('remove_contam', row['sequence_replicate_number'], 1),
                str(row['seqrep_id']) + '.1.' + row['remove_contam_reads_file_1_md5'] + '.fq.gz',
                iso_dir.reads_filename('remove_contam', row['sequence_replicate_number'], 2),
                str(row['seqrep_id']) + '.2.' + row['remove_contam_reads_file_2_md5'] + '.fq.gz',
            ))


        self.pool = multiprocessing.Pool(self.fq_upload_threads)
        upload_return_values = self.pool.starmap(_upload_fastq_file_pair, zip(fq_pairs_to_upload, itertools.repeat(self.ini_file), itertools.repeat(self.unit_test)))
        upload_success = {x[0]: x[1] for x in upload_return_values}
        fq_pairs_to_upload = {x[0]: x for x in fq_pairs_to_upload}

        # Fastqs are uploaded, now submit the xmls and update the database
        for row in data_in:
            assert row['seqrep_id'] in fq_pairs_to_upload
            assert row['seqrep_id'] in upload_success
            assert row['ena_run_accession'] is None
            assert row['ena_experiment_accession'] is not None

            iso_dir = isolate_dir.IsolateDir(self.pipeline_root, row['sample_id'], row['isolate_id'])
            object_xml = iso_dir.xml_submission_file('run', sequence_replicate = row['sequence_replicate_number'])
            object_alias = 'run.' + str(row['isolate_id'])
            submit_alias = 'submit.' + object_alias
            center_name = DatasetSubmitter._ena_center_name_from_db_data(data_in, number_to_name_dict=self.centre_number_to_name)
            title = None # not needed for a run

            obj_creator = object_creator.ObjectCreator(self.ini_file, 'run', object_xml, object_alias,
                submit_alias, center_name, title,
                experiment_accession=row['ena_experiment_accession'],
                reads_1=fq_pairs_to_upload[row['seqrep_id']][2],
                md5_1=row['remove_contam_reads_file_1_md5'],
                reads_2=fq_pairs_to_upload[row['seqrep_id']][4],
                md5_2=row['remove_contam_reads_file_2_md5'],
                use_test_server=self.use_test_server,
                unit_test=self.unit_test,
                broker_name=self.broker_name,
            )
            obj_creator.run()

            if obj_creator.submission_receipt.successful:
                try:
                    run_accession = obj_creator.submission_receipt.accessions['RUN']
                except:
                    run_accession = 'FAIL'
            else:
                run_accession = 'FAIL'

            row['ena_run_accession'] = run_accession
            self.db.update_row('Seqrep', {'seqrep_id': row['seqrep_id']}, {'ena_run_accession': run_accession})
            self.db.commit()


    def run(self):
        db_data = self._get_data_from_db()
        self._submit_study_object(db_data)
        self._submit_sample_objects(db_data)
        self._submit_experiment_objects(db_data)
        self._submit_runs(db_data)
