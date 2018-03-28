import csv
import copy
import datetime
import os
import re
import sys
import tempfile
from operator import attrgetter, itemgetter
from clockwork import db_connection, db_schema, isolate_dir, fastqc, lock_file, mykrobe, reference_dir, samtools_qc, utils
from clockwork import __version__ as clockwork_version

class Error (Exception): pass

class Db:
    def __init__(self, ini_file):
        '''ini_file must be the filename of a config file, and needs the following:

        [db_login]
        user = username
        password = my_password
        host = hostname
        db = db_name
        port = 1234

        The port line is optional. If not there, the default will be used.'''
        self.connection = db_connection.DbConnection(ini_file)
        self.cursor = self.connection.connection.cursor()


    def reconnect(self):
        self.connection.reconnect()
        self.cursor = self.connection.connection.cursor()


    def execute(self, cmd, row=None):
        '''Runs the given command. Raises an error if the command to be run
        is not successful'''
        if row is None:
            try:
                self.cursor.execute(cmd)
            except:
                raise Error('Error running the following mysql command: ' + cmd)
        else:
            try:
                self.cursor.execute(cmd, row)
            except:
                raise Error('Error running the following mysql command: ' + cmd + ', ' + str(row))


    def add_row_to_table(self, table, data):
        '''Adds a row to the given table.
        Data should be a dictionary of column name => value'''
        cmd = ['INSERT INTO', 'table']
        columns = []
        values = []
        for column, col_type in db_schema.tables[table]:
            columns.append(column)
            values.append(Db._value_to_string(data[column], col_type))

        if not( len(data) == len(columns) == len(values)):
            print('Mismatch in number of keys/values when trying to add row to ' + table + ' table.', file=sys.stderr)
            print('Columns should be:', *columns, file=sys.stderr)
            print('Data dict is:', data, file=sys.stderr)
            raise Error('Error adding row to table')

        cmd = ' '.join([
            'INSERT INTO', table,
            '(' + ', '.join(columns) + ')',
            'VALUES (' + ', '.join(values) + ')'
        ])

        self.execute(cmd)


    @classmethod
    def _value_to_string(cls, value, var_type):
        if value is None:
            return 'NULL'
        elif var_type in {'date', 'text'}:
            return '"' + str(value) + '"'
        else:
            return str(value)


    def update_row(self, table, where_dict, update_dict):
        '''Updates a row in the given table. row(s) updated are
        determined by the dictionary where_dict, by taking all
        rows where col1==value1, col2==value2 .. .etc for all
        key/values  colX=>valueX in the dict.
        Row(s) are updated by being given the values in the update_dict
        dictionary, where key=column name and value = new value.'''
        update_keys, update_values = zip(*update_dict.items())
        where_keys, where_values = zip(*where_dict.items())

        set_values = []
        where_values = []

        for column, col_type in db_schema.tables[table]:
            if column in update_dict:
                set_values.append(column + ' = ' + Db._value_to_string(update_dict[column], col_type))
            if column in where_dict:
                value_string = Db._value_to_string(where_dict[column], col_type)
                if value_string == 'NULL':
                    where_values.append(column + ' IS NULL')
                else:
                    where_values.append(column + ' = ' + value_string)

        assert len(set_values) == len(update_dict)
        assert len(where_values) == len(where_dict)

        cmd = ' '.join([
            'UPDATE', table,
            'SET', ', '.join(set_values),
            'WHERE', ' AND '.join(where_values)
        ])
        self.execute(cmd)


    def query_to_dict(self, query):
        '''Runs the given query and returns a list of results.
        Each element of the list is a dictionary of column name => value.'''
        self.execute(query)
        columns = [description[0] for description in self.cursor.description]
        return [dict(zip(columns, row)) for row in self.cursor.fetchall()]


    def get_rows_from_table(self, table, columns='*', where=None, order_by=None):
        '''Gets rows from the given table, as a list of dictionaryies.
        All rows and columns retrieved unordered by default. Filter them by using
        the strings 'columns' and 'where'.'''
        cmd = ['SELECT', columns, 'FROM', table]
        if where is not None:
            cmd.extend(['WHERE', where])

        if order_by is not None:
            cmd.extend(['ORDER BY', order_by])

        return self.query_to_dict(' '.join(cmd))


    def isolate_id_to_sample_id(self, isolate_id):
        got_rows = self.get_rows_from_table('Isolate', where='isolate_id = ' + str(isolate_id))
        if len(got_rows) == 0:
            return None
        elif len(got_rows) == 1:
            return got_rows[0]['sample_id']
        else:
            raise Error('Error! More than one row found in Isolate table for isolate_id ' + str(isolate_id))


    def seqrep_id_to_sequence_replicate_number(self, seqrep_id):
        got_rows = self.get_rows_from_table('Seqrep', where='seqrep_id = ' + str(seqrep_id))
        if len(got_rows) == 0:
            return None
        elif len(got_rows) == 1:
            return got_rows[0]['sequence_replicate_number']
        else:
            raise Error('Error! More than one row found in Seqrep table for seqrep_id ' + str(seqrep_id))


    def has_sample_isolate_seq_replicate(self, sample_id, isolate_number_from_lab, sequence_replicate):
        '''Returns true iff the give sample_id, isolate_number_from_lab and sequence_replicate
        are in the database. Throws error if if more than one isolate is found with the given
        sample_id and isolate_number_from_lab.'''
        where_query = ' '.join([
            'sample_id =', str(sample_id),
            'AND isolate_number_from_lab = "' +  isolate_number_from_lab + '"',
        ])
        isolate_matches = self.get_rows_from_table('Isolate', where=where_query)

        if len(isolate_matches) == 0:
            return False
        elif len(isolate_matches) > 1:
            raise Error('Error! Found more than one Isolate with sample_id ' + str(sample_id) + ' and isolate_number_from_lab ' + isolate_number_from_lab + '.')

        where_query = ' '.join([
            'isolate_id =', str(isolate_matches[0]['isolate_id']),
            'AND sequence_replicate_number = "' + str(sequence_replicate) + '"',
        ])
        seqrep_matches = self.get_rows_from_table('Seqrep', where=where_query)
        return len(seqrep_matches) > 0


    def _get_sample_and_replicate_uniqueness(self, sample_dict):
        '''Returns tuple: patient_site_lab_unique (bool), replicates_exist (bool), sample_id (int).
        sample_dict should be an element of the list
        made by SpreadsheetImporter.load_data_from_spreadsheet'''
        where_query = ''.join([
            'subject_id="', sample_dict['subject_id'], '"'
            ' and site_id="', sample_dict['site_id'], '"'
            ' and sample_id_from_lab="', sample_dict['lab_id'], '"'
        ])
        sample_rows = self.get_rows_from_table('Sample', where=where_query)
        sample_id = None
        replicates_exist = False
        patient_site_lab_unique = len(sample_rows) <= 1
        sample_id = sample_rows[0]['sample_id'] if len(sample_rows) == 1 else None

        if len(sample_rows) > 0:
            for sample_row in sample_rows:
                if self.has_sample_isolate_seq_replicate(sample_row['sample_id'], sample_dict['isolate_number'], sample_dict['sequence_replicate_number']):
                    replicates_exist = True
                    break

        return patient_site_lab_unique, replicates_exist, sample_id


    def add_one_seqrep(self, sample_dict):
        '''Adds one sequencing replicate to the database using data
        in sample_dict, which is an element of list made by
        SpreadsheetImporter.load_data_from_spreadsheet.
        Returns tuple (seqrep_id, isolate_id, sample_id) - row IDs of this seqrep
        (makes new row in Seqrep table, but other ids may not be new)'''
        # check if we have seen this sample before
        patient_site_lab_unique, replicates_exist, sample_id = self._get_sample_and_replicate_uniqueness(sample_dict)

        if not patient_site_lab_unique:
            raise Error('Error! Should not find >1 row for one sample. Cannot continue. Data was:\n' + sample_dict)

        if replicates_exist:
            raise Error('Error! isolate_number and sequence_replicate_number already in database for data:\n' + sample_dict)

        if sample_id is None:
            sample_row = {
                'sample_id': sample_id,
                'subject_id': sample_dict['subject_id'],
                'site_id': sample_dict['site_id'],
                'sample_id_from_lab': sample_dict['lab_id'],
                'dataset_name': sample_dict['dataset_name'],
                'ena_center_name': sample_dict['ena_center_name'],
                'ena_sample_accession': sample_dict['ena_sample_accession'],
                'ena_study_accession': None,
            }
            self.add_row_to_table('Sample', sample_row)
            sample_id = self.cursor.lastrowid

        assert sample_id is not None

        isolate_id = None
        isolate_where = 'sample_id = ' + str(sample_id) + ' AND isolate_number_from_lab = "' + sample_dict['isolate_number'] + '"'
        isolate_rows = self.get_rows_from_table('Isolate', where=isolate_where)
        if len(isolate_rows) == 0:
            isolate_row = {
                'isolate_id': None,
                'sample_id': sample_id,
                'isolate_number_from_lab': sample_dict['isolate_number'],
                'pool_sequence_replicates': 1,
                'ena_experiment_accession': None,
            }
            self.add_row_to_table('Isolate', isolate_row)
            isolate_id = self.cursor.lastrowid
        elif len(isolate_rows) == 1:
            isolate_id = isolate_rows[0]['isolate_id']
        else:
            raise Error('Error! More than one row found from query: ' + isolate_where)

        assert isolate_id is not None


        seqrep_row = {
            'seqrep_id': None,
            'isolate_id': isolate_id,
            'sequence_replicate_number': sample_dict['sequence_replicate_number'],
            'original_reads_file_1_md5': sample_dict['reads_file_1_md5'],
            'original_reads_file_2_md5': sample_dict['reads_file_2_md5'],
            'remove_contam_reads_file_1_md5': None,
            'remove_contam_reads_file_2_md5': None,
            'withdrawn': 0,
            'import_status': 0,
            'submission_date': sample_dict['submission_date'],
            'instrument_model': sample_dict['instrument_model'],
            'submit_to_ena': sample_dict['submit_to_ena'],
            'ena_run_accession': sample_dict['ena_run_accession'],
            'ena_on_hold': sample_dict['ena_on_hold'],
        }
        self.add_row_to_table('Seqrep', seqrep_row)
        seqrep_id = self.cursor.lastrowid
        return seqrep_id, isolate_id, sample_id


    def commit(self):
        try:
            self.connection.commit()
        except:
            raise Error('Error committing changes to database. Cannot continue')


    def make_remove_contam_jobs_tsv(self, outfile, pipeline_root, reference_id, references_root, dataset_name=None):
        '''Writes TSV file of job info for isolates that need to have remove_contam
           run on them. Used by nextflow remove_contam pipeline.
           remove_comtan can only be run on a sequence replicate once,
           regardless of the version of the pipeline that was run'''
        refdir = self.get_reference_dir(reference_id, os.path.abspath(references_root))
        pipeline_root = os.path.abspath(pipeline_root)

        pipeline_select = 'SELECT seqrep_id FROM Pipeline where Pipeline.pipeline_name = "remove_contam"'
        from_query = 'FROM (Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id JOIN Sample ON Isolate.sample_id = Sample.sample_id) WHERE'
        if dataset_name is None:
            dataset_query = ''
        else:
            dataset_query = 'AND dataset_name = "' + dataset_name + '"'

        query = ' '.join([
            'SELECT Sample.sample_id, Isolate.isolate_id, seqrep_id, sequence_replicate_number',
            from_query, '(',
              'seqrep_id NOT IN (' + pipeline_select + ')',
              dataset_query,
              'AND Seqrep.import_status = 1 AND Seqrep.withdrawn != 1)'
        ])

        rows = sorted(self.query_to_dict(query), key=itemgetter('sample_id', 'isolate_id', 'seqrep_id', 'sequence_replicate_number'))

        with open(outfile, 'w') as f:
            print('reads_in1', 'reads_in2', 'counts_tsv', 'reads_contam1', 'reads_contam2', 'reads_remove_contam1', 'reads_remove_contam2',
                  'sample_id', 'seqrep_id', 'isolate_id', 'sequence_replicate_number', 'reference_id', 'ref_fasta', 'contam_tsv',
                  sep='\t', file=f)

            for row in rows:
                iso_dir = isolate_dir.IsolateDir(pipeline_root, row['sample_id'], row['isolate_id'])
                reads_in1 = iso_dir.reads_filename('original', row['sequence_replicate_number'], 1)
                reads_in2 = iso_dir.reads_filename('original', row['sequence_replicate_number'], 2)
                counts_tsv = iso_dir.contamination_counts_filename(row['sequence_replicate_number'])
                reads_contam1 = iso_dir.reads_filename('contam', row['sequence_replicate_number'], 1)
                reads_contam2 = iso_dir.reads_filename('contam', row['sequence_replicate_number'], 2)
                reads_remove_contam1 = iso_dir.reads_filename('remove_contam', row['sequence_replicate_number'], 1)
                reads_remove_contam2 = iso_dir.reads_filename('remove_contam', row['sequence_replicate_number'], 2)
                print(reads_in1, reads_in2, counts_tsv, reads_contam1, reads_contam2, reads_remove_contam1, reads_remove_contam2,
                  row['sample_id'], row['seqrep_id'], row['isolate_id'], row['sequence_replicate_number'], reference_id,
                  refdir.ref_fasta, refdir.remove_contam_metadata_tsv,
                  sep='\t', file=f)

        for row in rows:
            db_row = {'isolate_id': row['isolate_id'], 'seqrep_id': row['seqrep_id'], 'seqrep_pool': None,
                'version': clockwork_version, 'pipeline_name': 'remove_contam', 'status': 0, 'reference_id': reference_id}
            self.add_row_to_table('Pipeline', db_row)

        self.commit()


    def make_qc_jobs_tsv(self, outfile, pipeline_root, reference_id, references_root, pipeline_version=None, dataset_name=None):
        '''Writes TSV file of job info for isolates that need to have
        QC pipeline run on them. Used by nextflow qc pipeline.'''
        pipeline_version = clockwork_version if pipeline_version is None else pipeline_version
        refdir = self.get_reference_dir(reference_id, os.path.abspath(references_root))
        pipeline_root = os.path.abspath(pipeline_root)

        pipeline_remove_contam_select = 'SELECT seqrep_id FROM Pipeline as Pipeline2 where Pipeline2.pipeline_name = "remove_contam" and Pipeline2.status = 1'
        pipeline_qc_select = 'SELECT seqrep_id FROM Pipeline where Pipeline.version = "' + pipeline_version + '" and Pipeline.pipeline_name = "qc"'
        from_query = 'FROM (Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id JOIN Sample ON Isolate.sample_id = Sample.sample_id) WHERE'
        if dataset_name is None:
            dataset_query = ''
        else:
            dataset_query = 'AND dataset_name = "' + dataset_name + '"'

        query = ' '.join([
            'SELECT Sample.sample_id, Isolate.isolate_id, seqrep_id, sequence_replicate_number',
            from_query, '(',
              'seqrep_id NOT IN (' + pipeline_qc_select + ') ',
              'AND seqrep_id IN (' + pipeline_remove_contam_select + ')',
              dataset_query,
              'AND Seqrep.import_status = 1 AND Seqrep.withdrawn != 1)'
        ])
        rows = sorted(self.query_to_dict(query), key=itemgetter('sample_id', 'isolate_id', 'seqrep_id', 'sequence_replicate_number'))

        with open(outfile, 'w') as f:
            print('reads_in1', 'reads_in2', 'output_dir', 'sample_id', 'seqrep_id',
                'isolate_id', 'sequence_replicate_number', 'reference_id', 'ref_fasta', sep='\t', file=f)

            for row in rows:
                iso_dir = isolate_dir.IsolateDir(pipeline_root, row['sample_id'], row['isolate_id'])
                reads_in1 = iso_dir.reads_filename('remove_contam', row['sequence_replicate_number'], 1)
                reads_in2 = iso_dir.reads_filename('remove_contam', row['sequence_replicate_number'], 2)
                output_dir = iso_dir.pipeline_dir(row['sequence_replicate_number'], 'qc', pipeline_version)
                assert not os.path.exists(output_dir)
                try:
                    os.makedirs(output_dir)
                except:
                    raise Error('Error making directory ' + output_dir + ' -- cannot continue')
                print(reads_in1, reads_in2, output_dir, row['sample_id'], row['seqrep_id'],
                    row['isolate_id'], row['sequence_replicate_number'], reference_id, refdir.ref_fasta, sep='\t', file=f)

        for row in rows:
            db_row = {'isolate_id': row['isolate_id'], 'seqrep_id': row['seqrep_id'], 'seqrep_pool': None,
                'version': pipeline_version, 'pipeline_name': 'qc', 'status': 0, 'reference_id': reference_id}
            self.add_row_to_table('Pipeline', db_row)

        self.commit()



    def make_variant_call_or_mykrobe_jobs_tsv(self, pipeline_name, outfile, pipeline_root, reference_id, references_root, pipeline_version=None, dataset_name=None):
        '''Writes TSV file of job info for isolates that need to have
        variant_call or mykrobe pipeline run on them.
        Used by nextflow variant_call or mykrobe pipeline.
        pipeline_name must be 'variant_call', or 'mykrobe'.'''
        assert pipeline_name in {'mykrobe', 'variant_call'}
        pipeline_version = clockwork_version if pipeline_version is None else pipeline_version
        refdir = self.get_reference_dir(reference_id, os.path.abspath(references_root))
        pipeline_root = os.path.abspath(pipeline_root)

        # For pooled, need to get the isolates that have been run already, then
        # for each run: check if it was run using all sequence replicates. If not,
        # then need to rerun it. Also, must check that each of the sequence
        # replicates have had remove_contam run already, and that there is
        # not a pending variant_call/mykrobe run (ie pipeline status = 0).
        # This could probably be done with some awesome mysql-fu. I'll
        # use python dicts instead (maybe less efficient?)

        # Get all rows that have had remove contam run
        isolates_and_seqreps_remove_contam_run = self.get_rows_from_table('Pipeline', columns='isolate_id,seqrep_id,seqrep_pool', where='pipeline_name = "remove_contam" AND status = 1')
        seqreps_have_had_remove_contam_run = {x['seqrep_id'] for x in isolates_and_seqreps_remove_contam_run}

        # Get all variant_call/mykrobe jobs for this pipeline version
        # that have been run or are pending to run
        jobs_run_already = self.query_to_dict('SELECT isolate_id, seqrep_id, seqrep_pool FROM Pipeline where pipeline_name = "' + pipeline_name + '" AND reference_id = ' + str(reference_id) + ' AND version = "' + pipeline_version + '"')
        pools_have_had_pipeline_run = {(x['isolate_id'], x['seqrep_pool']) for x in jobs_run_already if x['seqrep_pool'] is not None}
        unpooled_seqreps_have_had_pipeline_run = {x['seqrep_id'] for x in jobs_run_already if x['seqrep_id'] is not None}

        # Get isolates/seqreps, for filtering later for those that could have
        # the pipeline run on them
        seqrep_isolate_sample_join = 'Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id JOIN Sample ON Isolate.sample_id = Sample.sample_id'
        where_query = 'WHERE Seqrep.import_status = 1 AND Seqrep.withdrawn != 1'
        dataset_query = '' if dataset_name is None else ' AND dataset_name = "' + dataset_name + '"'
        query = 'SELECT Sample.site_id, Sample.subject_id, Sample.sample_id_from_lab, Isolate.isolate_id, Isolate.isolate_number_from_lab, Isolate.sample_id, Isolate.pool_sequence_replicates, Seqrep.seqrep_id, Seqrep.sequence_replicate_number FROM (' + seqrep_isolate_sample_join + ')' + where_query + dataset_query
        all_possible_rows = self.query_to_dict(query)
        sample_data = {}
        isolate_id_to_number_from_lab = {}

        # pooled: loop over each isolate. Get list of all seqreps for that isolate.
        #   1. Check if all seqreps have had remove_contam run
        #   2. Check if pipeline has been run for that pool of seqreps
        pooled_isolates_to_check = {}

        for seqrep_data in all_possible_rows:
            isolate_id_to_number_from_lab[seqrep_data['isolate_id']] = seqrep_data['isolate_number_from_lab']
            sample_data[seqrep_data['sample_id']] = {
                'subject_id': seqrep_data['subject_id'],
                'site_id': seqrep_data['site_id'],
                'lab_id': seqrep_data['sample_id_from_lab'],
            }
            if seqrep_data['pool_sequence_replicates'] == 1:
                if seqrep_data['isolate_id'] not in pooled_isolates_to_check:
                    pooled_isolates_to_check[seqrep_data['isolate_id']] = {
                        'sample_id': seqrep_data['sample_id'], 'seqrep_ids': [], 'sequence_replicate_numbers': []}
                assert pooled_isolates_to_check[seqrep_data['isolate_id']]['sample_id'] == seqrep_data['sample_id']
                pooled_isolates_to_check[seqrep_data['isolate_id']]['seqrep_ids'].append(seqrep_data['seqrep_id'])
                pooled_isolates_to_check[seqrep_data['isolate_id']]['sequence_replicate_numbers'].append(seqrep_data['sequence_replicate_number'])

        data_to_print = []

        for isolate_id, isolate_data in pooled_isolates_to_check.items():
            seqrep_set = set(isolate_data['seqrep_ids'])
            if not seqrep_set.issubset(seqreps_have_had_remove_contam_run):
                continue

            seqrep_string = '_'.join([str(x) for x in sorted(isolate_data['sequence_replicate_numbers'])])
            if (isolate_id, seqrep_string) not in pools_have_had_pipeline_run:
                data_to_print.append({'sample': isolate_data['sample_id'], 'isolate': isolate_id, 'seqrep_ids': isolate_data['seqrep_ids'], 'sequence_replicate_numbers': isolate_data['sequence_replicate_numbers'], 'pool': True})



        # not pooled: loop over each seqrep and check if it's:
        #   1. had remove_contam run
        #   2. not had pipeline run (or pending to run)
        for seqrep_data in all_possible_rows:
            if seqrep_data['pool_sequence_replicates'] == 0 and seqrep_data['seqrep_id'] in seqreps_have_had_remove_contam_run and seqrep_data['seqrep_id'] not in unpooled_seqreps_have_had_pipeline_run:
                data_to_print.append({'sample': seqrep_data['sample_id'], 'isolate': seqrep_data['isolate_id'], 'seqrep_ids': [seqrep_data['seqrep_id']], 'sequence_replicate_numbers': [seqrep_data['sequence_replicate_number']], 'pool': False})


        data_to_print.sort(key=itemgetter('sample', 'isolate'))
        database_rows = []

        with open(outfile, 'w') as f:
            print('reads_in1', 'reads_in2', 'output_dir', 'sample_id', 'pool',
            'isolate_id', 'seqrep_id', 'sequence_replicate_number', 'reference_id', 'reference_dir', sep='\t', file=f)

            for data in data_to_print:
                seqrep_numbers_string = '_'.join([str(x) for x in sorted(data['sequence_replicate_numbers'])])
                sample_name = '.'.join([
                    'site', sample_data[data['sample']]['site_id'],
                    'iso', isolate_id_to_number_from_lab[data['isolate']],
                    'subject', sample_data[data['sample']]['subject_id'],
                    'lab_id', sample_data[data['sample']]['lab_id'],
                    'seq_reps', seqrep_numbers_string,
                ])

                iso_dir = isolate_dir.IsolateDir(pipeline_root, data['sample'], data['isolate'])
                reads_in1 = []
                reads_in2 = []
                for seqrep_number in data['sequence_replicate_numbers']:
                    reads_in1.append(iso_dir.reads_filename('remove_contam', seqrep_number, 1))
                    reads_in2.append(iso_dir.reads_filename('remove_contam', seqrep_number, 2))

                seqrep_ids_string = '_'.join([str(x) for x in sorted(data['seqrep_ids'])])

                print(
                    ' '.join(reads_in1),
                    ' '.join(reads_in2),
                    iso_dir.pipeline_dir(seqrep_numbers_string, pipeline_name, pipeline_version, reference_id=reference_id),
                    sample_name,
                    1 if data['pool'] else 0,
                    data['isolate'],
                    seqrep_ids_string,
                    seqrep_numbers_string,
                    reference_id,
                    refdir.directory,
                    sep='\t', file=f,
                )


                database_rows.append({
                    'isolate_id': data['isolate'],
                    'seqrep_id': None if data['pool'] else data['seqrep_ids'][0],
                    'version': pipeline_version,
                    'pipeline_name': pipeline_name,
                    'status': 0,
                    'seqrep_pool': seqrep_numbers_string if data['pool'] else None,
                    'reference_id': reference_id,
                })

        for row in database_rows:
            self.add_row_to_table('Pipeline', row)

        self.commit()


    def make_generic_pipeline_jobs_tsv(self, outfile, pipeline_root, pipeline_name, pipeline_version=None, dataset_name=None):
        pipeline_version = clockwork_version if pipeline_version is None else pipeline_version
        pipeline_root = os.path.abspath(pipeline_root)

        # For pooled, need to get the isolates that have been run already, then
        # for each run: check if it was run using all sequence replicates. If not,
        # then need to rerun it. Also, must check that each of the sequence
        # replicates have had remove_contam run already, and that there is
        # not a pending pipeline run (ie pipeline status = 0).
        # This could probably be done with some awesome mysql-fu. I'll
        # use python dicts instead (maybe less efficient?)

        # Get all rows that have had remove contam run
        isolates_and_seqreps_remove_contam_run = self.get_rows_from_table('Pipeline', columns='isolate_id,seqrep_id,seqrep_pool', where='pipeline_name = "remove_contam" AND status = 1')
        seqreps_have_had_remove_contam_run = {x['seqrep_id'] for x in isolates_and_seqreps_remove_contam_run}

        # Get all jobs for this pipeline version
        # that have been run or are pending to run
        jobs_run_already = self.query_to_dict('SELECT isolate_id, seqrep_id, seqrep_pool FROM Pipeline where pipeline_name = "' + pipeline_name + '" AND version = "' + pipeline_version + '"')
        pools_have_had_pipeline_run = {(x['isolate_id'], x['seqrep_pool']) for x in jobs_run_already if x['seqrep_pool'] is not None}
        unpooled_seqreps_have_had_pipeline_run = {x['seqrep_id'] for x in jobs_run_already if x['seqrep_id'] is not None}

        # Get isolates/seqreps, for filtering later for those that could have
        # pipeline run on them
        seqrep_isolate_sample_join = 'Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id JOIN Sample ON Isolate.sample_id = Sample.sample_id'
        where_query = 'WHERE Seqrep.import_status = 1 AND Seqrep.withdrawn != 1'
        dataset_query = '' if dataset_name is None else ' AND dataset_name = "' + dataset_name + '"'
        query = 'SELECT Isolate.isolate_id, Isolate.sample_id, Isolate.pool_sequence_replicates, Seqrep.seqrep_id, Seqrep.sequence_replicate_number FROM (' + seqrep_isolate_sample_join + ')' + where_query + dataset_query
        all_possible_rows = self.query_to_dict(query)

        # pooled: loop over each isolate. Get list of all seqreps for that isolate.
        #   1. Check if all seqreps have had remove_contam run
        #   2. Check if pipeline has been run for that pool of seqreps
        pooled_isolates_to_check = {}

        for seqrep_data in all_possible_rows:
            if seqrep_data['pool_sequence_replicates'] == 1:
                if seqrep_data['isolate_id'] not in pooled_isolates_to_check:
                    pooled_isolates_to_check[seqrep_data['isolate_id']] = {
                        'sample_id': seqrep_data['sample_id'], 'seqrep_ids': [], 'sequence_replicate_numbers': []}
                assert pooled_isolates_to_check[seqrep_data['isolate_id']]['sample_id'] == seqrep_data['sample_id']
                pooled_isolates_to_check[seqrep_data['isolate_id']]['seqrep_ids'].append(seqrep_data['seqrep_id'])
                pooled_isolates_to_check[seqrep_data['isolate_id']]['sequence_replicate_numbers'].append(seqrep_data['sequence_replicate_number'])

        data_to_print = []

        for isolate_id, isolate_data in pooled_isolates_to_check.items():
            seqrep_set = set(isolate_data['seqrep_ids'])
            if not seqrep_set.issubset(seqreps_have_had_remove_contam_run):
                continue

            seqrep_string = '_'.join([str(x) for x in sorted(isolate_data['sequence_replicate_numbers'])])
            if (isolate_id, seqrep_string) not in pools_have_had_pipeline_run:
                data_to_print.append({'sample': isolate_data['sample_id'], 'isolate': isolate_id, 'seqrep_ids': isolate_data['seqrep_ids'], 'sequence_replicate_numbers': isolate_data['sequence_replicate_numbers'], 'pool': True})

        # not pooled: loop over each seqrep and check if it's:
        #   1. had remove_contam run
        #   2. not had pipeline run (or pending to run)
        for seqrep_data in all_possible_rows:
            if seqrep_data['pool_sequence_replicates'] == 0 and seqrep_data['seqrep_id'] in seqreps_have_had_remove_contam_run and seqrep_data['seqrep_id'] not in unpooled_seqreps_have_had_pipeline_run:
                data_to_print.append({'sample': seqrep_data['sample_id'], 'isolate': seqrep_data['isolate_id'], 'seqrep_ids': [seqrep_data['seqrep_id']], 'sequence_replicate_numbers': [seqrep_data['sequence_replicate_number']], 'pool': False})


        data_to_print.sort(key=itemgetter('sample', 'isolate'))
        database_rows = []

        with open(outfile, 'w') as f:
            print('reads_in1', 'reads_in2', 'output_dir', 'sample_id', 'pool',
            'isolate_id', 'seqrep_id', 'sequence_replicate_number', sep='\t', file=f)

            for data in data_to_print:
                iso_dir = isolate_dir.IsolateDir(pipeline_root, data['sample'], data['isolate'])
                reads_in1 = []
                reads_in2 = []
                for seqrep_number in data['sequence_replicate_numbers']:
                    reads_in1.append(iso_dir.reads_filename('remove_contam', seqrep_number, 1))
                    reads_in2.append(iso_dir.reads_filename('remove_contam', seqrep_number, 2))

                seqrep_numbers_string = '_'.join([str(x) for x in sorted(data['sequence_replicate_numbers'])])
                seqrep_ids_string = '_'.join([str(x) for x in sorted(data['seqrep_ids'])])

                print(
                    ' '.join(reads_in1),
                    ' '.join(reads_in2),
                    iso_dir.pipeline_dir(seqrep_numbers_string, pipeline_name, pipeline_version),
                    data['sample'],
                    1 if data['pool'] else 0,
                    data['isolate'],
                    seqrep_ids_string,
                    seqrep_numbers_string,
                    sep='\t', file=f,
                )


                database_rows.append({
                    'isolate_id': data['isolate'],
                    'seqrep_id': None if data['pool'] else data['seqrep_ids'][0],
                    'version': pipeline_version,
                    'pipeline_name': pipeline_name,
                    'status': 0,
                    'seqrep_pool': seqrep_numbers_string if data['pool'] else None,
                    'reference_id': None,
                })

        for row in database_rows:
            self.add_row_to_table('Pipeline', row)

        self.commit()


    def get_vcfs_and_reads_files_for_minos_multi_sample_calling(self, dataset_name, pipeline_root, reference_id, pipeline_version=None):
        '''Returns a list of lines in the right format for
        minos's multi_sample_pipeline TDV data input file'''
        # We need the isolates that have had the variant calling pipeline run, from
        # there get the sample id, so we can get paths to VCFs and BAMs
        pipeline_version = clockwork_version if pipeline_version is None else pipeline_version
        query = 'SELECT Isolate.isolate_id, Isolate.sample_id, Pipeline.seqrep_pool FROM Isolate JOIN Pipeline on Isolate.isolate_id = Pipeline.isolate_id where Pipeline.version = "' + pipeline_version + '" and Pipeline.pipeline_name = "variant_call" and Pipeline.status = 1 and Pipeline.reference_id = ' + str(reference_id)
        all_wanted_db_rows = self.query_to_dict(query)
        files_list = []

        for row in all_wanted_db_rows:
            iso_dir = isolate_dir.IsolateDir(pipeline_root, row['sample_id'], row['isolate_id'])
            pipeline_dir = iso_dir.pipeline_dir(row['seqrep_pool'], 'variant_call', pipeline_version, reference_id=reference_id)
            vcf = os.path.join(pipeline_dir, 'minos', 'final.vcf')
            bam = os.path.join(pipeline_dir, 'samtools', 'rmdup.bam')
            files_list.append(vcf + '\t' + bam)

        return files_list


    def commit_and_close(self):
        self.commit()

        try:
            self.connection.close()
        except:
            raise Error('Error closing database connection. Cannot continue')


    def _update_remove_contam_stats(self, seqrep_id, read_counts_file):
        '''Updates the Read counts table for the given sequence replicate, using
        numbers in the file read_counts_file'''
        read_counts = {x: 0 for x in ['contamination', 'not_contamination', 'unmapped']}
        total_after_remove_contam = 0

        with open(read_counts_file) as f:
            for line in f:
                if line == 'Name\tIs_contam\tReads\n':
                    continue
                name, is_contam, number = line.rstrip().split('\t')
                if name == 'Unmapped':
                    read_counts['unmapped'] = int(number)
                elif name == 'Reads_kept_after_remove_contam':
                    total_after_remove_contam = int(number)
                elif is_contam == '1':
                    read_counts['contamination'] += int(number)
                else:
                    assert is_contam == '0'
                    read_counts['not_contamination'] += int(number)

        seqrep_rows = self.get_rows_from_table('Seqrep', where='seqrep_id =' + str(seqrep_id))
        if len(seqrep_rows) != 1:
            raise Error('Error. seqrep_id ' + str(seqrep_id) + ' not found in Seqrep table. Cannot continue')
        seqrep_rows = self.get_rows_from_table('Read_counts', where='seqrep_id =' + str(seqrep_id))
        if len(seqrep_rows) != 0:
            raise Error('Error. Found entry for seqrep_id ' + str(seqrep_id) + ' in Read_counts table.')

        total_reads = sum(read_counts.values())

        new_row = {
            'seqrep_id': seqrep_id,
            'original_total': total_reads,
            'contamination': read_counts['contamination'],
            'not_contamination': read_counts['not_contamination'],
            'unmapped': read_counts['unmapped'],
            'total_after_remove_contam': total_after_remove_contam,
        }
        self.add_row_to_table('Read_counts', new_row)


    def _update_qc_stats(self, seqrep_id, pipeline_version, pipeline_root):
        '''Updates QC table for the give nsequence replicate. Assumes that
        FASTQC and samtools QC have been run.'''
        query = 'SELECT Isolate.sample_id, Isolate.isolate_id, seqrep_id, sequence_replicate_number FROM (Seqrep JOIN Isolate ON Seqrep.isolate_id = Isolate.isolate_id) WHERE seqrep_id = ' + str(seqrep_id)
        got_rows = self.query_to_dict(query)
        if len(got_rows) != 1:
            raise Error('Error, expceted exactly 1 row from this query: ' + query + '\nbut got:' + str(got_rows))
        seqrep_data = got_rows[0]
        iso_dir = isolate_dir.IsolateDir(pipeline_root, seqrep_data['sample_id'], seqrep_data['isolate_id'])
        qc_dir = iso_dir.pipeline_dir(seqrep_data['sequence_replicate_number'], 'qc', pipeline_version)
        samtools_stats = samtools_qc.SamtoolsQc.stats_from_report(os.path.join(qc_dir, 'samtools_qc', 'samtools_qc.stats'))
        fastqc_stats = fastqc.Fastqc.gather_all_stats(os.path.join(qc_dir, 'fastqc'))
        assert len(fastqc_stats) == 2
        new_row = {'seqrep_id' : seqrep_id, 'pipeline_version' : pipeline_version}
        new_row.update({'samtools_' + x: samtools_stats[x] for x in samtools_stats})

        # Depending on the reads, fastqc info reported can vary. Eg per tile sequence quality
        # may or may not be there. So only keep keys that we expect to add to the database.
        wanted_keys = set([x[0] for x in db_schema.tables['QC']])
        for i in ('1', '2'):
            fastqc_filename_regex = re.compile(r'''.*[^0-9]([?P=<one_or_two>''' + i + r'''])\.[^0-9]+$''')
            new_key_prefix = 'fastqc' + i
            fastqc_stats_keys = [x for x in fastqc_stats if fastqc_filename_regex.match(x) is not None]
            assert len(fastqc_stats_keys) == 1
            fastqc_stats_key = fastqc_stats_keys[0]
            stats_dict = fastqc_stats[fastqc_stats_key]
            new_row.update({new_key_prefix + '_' + x: stats_dict[x] for x in stats_dict if new_key_prefix + '_' + x in wanted_keys})

        het_stats = samtools_qc.SamtoolsQc.het_snp_stats_from_summary_file(os.path.join(qc_dir, 'samtools_qc', 'het_snps.summary.tsv'))
        new_row['het_snp_positions'] = het_stats['Positions_used']
        new_row['het_snp_total_snps'] = het_stats['Total_SNPs']
        new_row['het_snp_het_calls'] = het_stats['Het_SNPs']

        self.add_row_to_table('QC', new_row)


    def update_finished_pipeline_run(self, isolate_id, seqrep_id, seqrep_pool, pipeline_name, new_pipeline_status, reference_id=None, pipeline_version=None, pipeline_root=None):
        assert (seqrep_id is None) != (seqrep_pool is None)
        if new_pipeline_status == 1:
            if pipeline_name == 'qc' and pipeline_root is None:
                raise Error('Error! Must supply pipeline_root when pipeline_name is qc')

            if pipeline_name == 'remove_contam' and pipeline_root is None:
                raise Error('Error! must supply contam_reads_count_file when pipeline_name is remove_contam')

        if pipeline_version is None:
            pipeline_version = clockwork_version

        pool = 'IS NULL' if seqrep_pool is None else '= "' + seqrep_pool + '"'

        where_query = [
            'isolate_id =', str(isolate_id),
            'AND seqrep_id ', 'IS NULL' if seqrep_id is None else '= ' + str(seqrep_id),
            'AND version = "' + pipeline_version + '"',
            'AND pipeline_name = "' + pipeline_name + '"',
            'AND seqrep_pool ', pool,
        ]

        if reference_id is not None:
            where_query.append('AND reference_id = ' + str(reference_id))

        pipeline_rows = self.get_rows_from_table('Pipeline', where=' '.join(where_query))

        if len(pipeline_rows) == 0:
            raise Error('Error! No match found in database for isolate_id=' + str(isolate_id) + ', seqrep_id=' + str(seqrep_id) + ', seqrep_pool=' + str(seqrep_pool) + ', pipeline_version=' + pipeline_version + ', pipeline_name=' + pipeline_name + ',reference_id=' + str(reference_id))
        elif len(pipeline_rows) > 1:
            raise Error('Error! More than one row found in database for isolate_id=' + str(isolate_id) + ', seqrep_id=' + str(seqrep_id) + ', seqrep_pool=' + str(seqrep_pool) + ', pipeline_version=' + pipeline_version + ', pipeline_name=' + pipeline_name + ',reference_id=' + str(reference_id) + '\n' + '\n'.join([str(row) for row in pipeline_rows]))

        row = pipeline_rows[0]
        if row['status'] != 0:
            raise Error('Error! Expected status 0 but got ' + str(row['status']) + ' for seqrep_id=' + str(seqrep_id) + ', pipeline_version=' + pipeline_version + ', pipeline_name=' + pipeline_name)

        where_dict = {
            'isolate_id': isolate_id,
            'seqrep_id': seqrep_id,
            'seqrep_pool': seqrep_pool,
            'version': pipeline_version,
            'pipeline_name': pipeline_name,
            'status': 0,
        }
        if reference_id is not None:
            where_dict['reference_id'] = reference_id
        update_dict = copy.copy(where_dict)
        update_dict['status'] = new_pipeline_status

        if new_pipeline_status == 1:
            if pipeline_name == 'qc':
                assert pipeline_root is not None
                self._update_qc_stats(seqrep_id, pipeline_version, pipeline_root)
            elif pipeline_name == 'remove_contam':
                sample_id = self.isolate_id_to_sample_id(isolate_id)
                if sample_id is None:
                    raise Error('Error getting sample_id for isolate_id ' + str(isolate_id))
                sequence_replicate_number = self.seqrep_id_to_sequence_replicate_number(seqrep_id)
                if sequence_replicate_number is None:
                    raise Error('Error getting sequence_replicate_number for seqrep_id ' + str(seqrep_id))
                iso_dir = isolate_dir.IsolateDir(pipeline_root, sample_id, isolate_id)
                contamination_counts_filename = iso_dir.contamination_counts_filename(sequence_replicate_number)
                self._update_remove_contam_stats(seqrep_id, contamination_counts_filename)
                md5_dict = {
                    'remove_contam_reads_file_1_md5': utils.md5(iso_dir.reads_filename('remove_contam', sequence_replicate_number, 1)),
                    'remove_contam_reads_file_2_md5': utils.md5(iso_dir.reads_filename('remove_contam', sequence_replicate_number, 2)),
                }
                self.update_row('Seqrep', {'seqrep_id': seqrep_id}, md5_dict)


        self.update_row('Pipeline', where_dict, update_dict)
        self.commit()


    @classmethod
    def _load_success_jobs_file(cls, infile):
        successful_jobs = set()
        with open(infile) as f:
            csv_reader = csv.DictReader(f, delimiter='\t')
            for column in ['isolate_id', 'seqrep_id', 'sequence_replicate_number']:
                if column not in csv_reader.fieldnames:
                    raise Error('Error! column "' + column + '" not found in file ' + infile)

            for row in csv_reader:
                pool = row.get('pool', '0') == '1'
                successful_jobs.add((int(row['isolate_id']), row['seqrep_id'], row['sequence_replicate_number'], pool))

        return successful_jobs


    def update_finished_pipeline_run_failed_jobs(self, jobs_tsv, success_jobs_file, pipeline_name, pipeline_version=None, reference_id=None):
        if pipeline_version is None:
            pipeline_version = clockwork_version

        successful_jobs = Db._load_success_jobs_file(success_jobs_file)

        with open(jobs_tsv) as f:
            csv_reader = csv.DictReader(f, delimiter='\t')
            for column in ['isolate_id', 'seqrep_id', 'sequence_replicate_number']:
                if column not in csv_reader.fieldnames:
                    raise Error('Error! column "' + column + '" not found in file ' + jobs_tsv)

            using_pools = 'pool' in csv_reader.fieldnames

            for row in csv_reader:
                pool = row.get('pool', '0') == '1'
                key = (int(row['isolate_id']), row['seqrep_id'], row['sequence_replicate_number'], pool)

                if key not in successful_jobs:
                    seqrep = None if pool else row['seqrep_id']
                    seqrep_pool = row['sequence_replicate_number'] if pool else None
                    self.update_finished_pipeline_run(row['isolate_id'], seqrep, seqrep_pool, pipeline_name, -1, pipeline_version=pipeline_version, reference_id=reference_id)


    def add_reference(self, name):
        got_rows = self.get_rows_from_table('Reference', columns='*', where='name ="' + name + '"')
        if len(got_rows):
            raise Error('Error! Reference with name "' + name + '" already in database. Cannot add it again')
        self.add_row_to_table('Reference', {'reference_id': None, 'name': name})
        reference_id = self.cursor.lastrowid
        self.commit()
        return reference_id


    def add_mykrobe_custom_panel(self, species, name, pipeline_references_root, probes_fasta, var_to_res_json):
        lock = lock_file.LockFile(os.path.join(pipeline_references_root, 'add_mykrobe_panel.lock'))
        reference_id = self.add_reference(name)
        self.commit()
        lock.stop()
        panel_dir = reference_dir.ReferenceDir(
            pipeline_references_root_dir=pipeline_references_root,
            reference_id=reference_id,
        )
        panel = mykrobe.CustomPanel(panel_dir.directory)
        panel.setup_files(species, probes_fasta, var_to_res_json)
        return reference_id


    def has_reference(self, reference_id):
        got_rows = self.get_rows_from_table('Reference', columns='*', where='reference_id =' + str(reference_id))
        return len(got_rows) > 0


    def get_reference_dir(self, reference_id, reference_root):
        if not self.has_reference(reference_id):
            raise Error('Reference with ID ' + str(reference_id) + ' not found.')
        return reference_dir.ReferenceDir(pipeline_references_root_dir=reference_root, reference_id=reference_id)


    def backup(self, backup_dir=None, outfile=None):
        if backup_dir is not None:
            backup_dir = os.path.abspath(backup_dir)
            assert outfile is None
            if not os.path.isdir(backup_dir):
                raise Error('Error! mysql backup directory not found: ' + backup_dir)
            now = datetime.datetime.now()
            now_string = '-'.join([str(x).zfill(2) for x in now.timetuple()[:6]])
            outfile = os.path.join(backup_dir, 'backup.' + now_string)
            if os.path.exists(outfile):
                i = 2
                while os.path.exists(outfile + '.' + str(i)):
                    i += 1
                outfile = outfile + '.' + str(i)
        elif outfile is not None:
            assert backup_dir is None
        else:
            raise Error('Error! Must provide back_dir or outfile. Cannot continue')

        outfile = os.path.abspath(outfile)

        if os.path.exists(outfile):
            raise Error('Error! mysql backup file already exists: ' + outfile)

        # to use the mysql dump command, need to make a config file
        # that has the login info in it
        tmp_fh = tempfile.NamedTemporaryFile(mode='w', prefix=outfile + '.tmp.cfg', delete=False)
        print('[client]',
            'host = ' + self.connection.host,
            'user = ' + self.connection.user,
            'password = ' + self.connection.password,
            sep='\n', file=tmp_fh)

        if self.connection.port is not None:
            print('port =', self.connection.port, file=tmp_fh)

        tmp_fh.close()

        backup_command = ' '.join([
            'mysqldump',
            '--defaults-file=' + tmp_fh.name,
            self.connection.db,
            '>', outfile,
        ])
        utils.syscall(backup_command)
        os.unlink(tmp_fh.name)

