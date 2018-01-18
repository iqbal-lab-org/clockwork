import copy
import datetime
import dateutil.parser
import os
import re
import string
import xlsxwriter
import openpyxl
from clockwork import db, lock_file, utils
from clockwork.common_data import allowed_sequencing_instruments


class Error (Exception): pass


columns = [
    'subject_id',
    'site_id',
    'lab_id',
    'isolate_number',
    'sequence_replicate_number',
    'submission_date',
    'reads_file_1',
    'reads_file_1_md5',
    'reads_file_2',
    'reads_file_2_md5',
    'dataset_name',
    'instrument_model',
    'ena_center_name',
    'submit_to_ena',
    'ena_on_hold',
    'ena_run_accession',
    'ena_sample_accession',
]

date_iso_regex = re.compile(r'^[0-9]{4}-[0-9]{2}-[0-9]{2}$')

def create_template_spreadsheet(filename):
    workbook = xlsxwriter.Workbook(filename)
    main_sheet = workbook.add_worksheet('data')

    # Want to lock the column headings (added later( of main_sheet.
    # See https://stackoverflow.com/questions/40885097/xlsxwriter-lock-only-specific-cells
    unlocked = workbook.add_format({'locked': False})
    bold_locked = workbook.add_format({'bold': True, 'locked': True})
    main_sheet.protect()

    admin_sheet = workbook.add_worksheet('admin')
    admin_sheet.protect()
    sequencing_instruments = sorted(list(allowed_sequencing_instruments))
    admin_sheet.write(0, 0, 'Allowed sequencing instruments')
    for i, instrument in enumerate(sequencing_instruments):
        admin_sheet.write(i + 1, 0, instrument)

    check_sequencing_instrument = {'validate': 'list', 'source': '=admin!$A$2:$A$' + str(len(sequencing_instruments) + 1)}
    integer_greater_than_zero = {'validate': 'integer', 'criteria': '>', 'value': 0, 'error_message': 'Must be an integer greater than zero'}
    zero_or_one = {'validate': 'integer', 'criteria': 'between', 'minimum': 0, 'maximum': 1, 'error_message': 'Must be 0 or 1'}
    date_later_than_2000 = {'validate': 'date', 'criteria': '>', 'value': datetime.date(2000, 1, 1), 'error_message': 'Must be a date from this century'}
    is_alphanumeric = lambda column: 'ISNUMBER(SUMPRODUCT(FIND(MID(' + column + '2,ROW(INDIRECT("1:"&LEN(' + column + '2))),1),"abcdefghijklmnopqrstuvwxyz0123456789")))'
    looks_like_md5 = lambda column: '=AND(32=LEN(' + column + '2), ' + is_alphanumeric(column) + ')'
    validations = {
        'isolate_number': copy.copy(integer_greater_than_zero),
        'sequence_replicate_number': copy.copy(integer_greater_than_zero),
        'instrument_model': copy.copy(check_sequencing_instrument),
        'submission_date': copy.copy(date_later_than_2000),
        'submit_to_ena': copy.copy(zero_or_one),
        'ena_on_hold': copy.copy(zero_or_one),
        'reads_file_1_md5': None,
        'reads_file_2_md5': None,
    }
    for i, column in enumerate(columns):
        main_sheet.set_column(i, i, len(column), unlocked)
        if column in validations:
            column_letter = string.ascii_uppercase[i]
            if column.endswith('md5'):
                validation = {'validate': 'custom', 'value': looks_like_md5(column_letter), 'error_message': 'MD5 sum must be 32 characters long and only contain numbers and lower case letters'}
            else:
                validation = validations[column]
            main_sheet.data_validation(column_letter + '2:' + column_letter + '1048576', validation)

    main_sheet.write_row('A1', columns, bold_locked)
    workbook.close()


class SpreadsheetImporter:
    def __init__(self, dropbox_dir, xlsx_file, db_ini_file, xlsx_archive_dir, jobs_outfile, db_backup_dir=None):
        self.dropbox_dir = os.path.abspath(dropbox_dir)
        self.xlsx_file = os.path.join(self.dropbox_dir, os.path.basename(xlsx_file))
        self.db_ini_file = os.path.abspath(db_ini_file)
        self.xlsx_archive_dir = os.path.abspath(xlsx_archive_dir)
        self.jobs_outfile = os.path.abspath(jobs_outfile)
        self.lock_file = os.path.join(self.xlsx_archive_dir, 'import_spreadsheet.lock')
        self.db_backup_dir = None if db_backup_dir is None else os.path.abspath(db_backup_dir)


    @classmethod
    def load_data_from_spreadsheet(cls, infile):
        if infile.endswith('.xlsx'):
            workbook = openpyxl.load_workbook(infile)
            sheet = workbook.active
            fields_list = []
            for row in sheet:
                fields_list.append([str(x.value) for x in row if x.value is not None])
        elif infile.endswith('.tsv'):
            with open(infile) as f:
                fields_list = []
                for line in f:
                    fields_list.append([x.strip() for x in line.rstrip().split('\t')])

        data = []


        for line_counter, fields in enumerate(fields_list):
            if line_counter == 0:
                if fields != columns:
                    fields_str = ', '.join([str(x) for x in fields])
                    expected_str = ', '.join(columns)
                    raise Error('Error in first line of spreadsheet "' + infile +
                            '". Column names not correct.\nExpected: ' + expected_str +
                            '\nGot:      ' + fields_str)
            elif len(fields) != len(columns):
                raise Error('Wrong number of fields in line ' + str(line_counter) + ' of spreadsheet: ' +
                        ', '.join([str(x) for x in fields]))
            else:
                new_data = dict(zip(columns, fields))

                for key in ['ena_run_accession', 'ena_sample_accession', 'reads_file_1_md5', 'reads_file_2_md5']:
                    if new_data[key] in {0, '0'}:
                        new_data[key] = None

                try:
                    date = dateutil.parser.parse(new_data['submission_date']).date()
                except:
                    raise Error('Date format error: ' + SpreadsheetImporter._row_data_dict_to_string(new_data))

                assert SpreadsheetImporter._looks_like_date(str(date))
                new_data['submission_date'] = date
                data.append(new_data)

        return data


    @classmethod
    def _row_data_dict_to_string(cls, row_data):
        return ', '.join([str(row_data[columns[i]]) for i in range(len(columns))])


    @classmethod
    def _looks_like_date(cls, date_string):
        '''Returns true iff the date_string looks like a date of the form YYYY-MM-DD.
           We're always getting ISO date strings from the datetime module, but
           check anyway'''
        if not date_iso_regex.match(date_string):
            return False

        year, month, day = date_string[:4], date_string[5:7], date_string[8:]

        try:
            datetime.date(int(year), int(month), int(day))
        except:
            return False

        return True


    @classmethod
    def _validate_data(cls, database, data, dropbox_dir):
        '''Input should be data, made by load_data_from_spreadsheet().
           Sanity checks that it is ok, and returns a list of error messages.
           If the list has length zero, then all is OK.'''
        errors = []
        all_filenames = {}
        all_replicates = {}
        replicate_keys = ('subject_id', 'site_id', 'lab_id', 'isolate_number', 'sequence_replicate_number')

        for data_dict in data:
            if type(data_dict['submission_date']) is not datetime.date:
                errors.append('Date format error: ' + SpreadsheetImporter._row_data_dict_to_string(data_dict))

            for i in [1, 2]:
                read_file_key = 'reads_file_' + str(i)
                filename = data_dict[read_file_key]
                md5_key = read_file_key + '_md5'

                if not os.path.exists(os.path.join(dropbox_dir, filename)):
                    errors.append('Reads file not found: ' + filename)

                all_filenames[filename] = all_filenames.get(filename, 0) + 1

                md5_file = os.path.join(dropbox_dir, filename + '.md5')
                if os.path.exists(md5_file):
                    md5sum_from_file = utils.load_md5_from_file(md5_file)
                else:
                    md5sum_from_file = None

                if md5sum_from_file is None and data_dict[md5_key] is None:
                    errors.append('No md5 for reads file ' + filename)
                elif md5sum_from_file is not None and data_dict[md5_key] is not None and md5sum_from_file != data_dict[md5_key]:
                    errors.append('Mismatch in md5 info for reads file ' + filename)
                elif data_dict[md5_key] is None and md5sum_from_file is not None:
                    data_dict[md5_key] = md5sum_from_file

            replicate = tuple([data_dict[x] for x in replicate_keys])
            all_replicates[replicate] = all_replicates.get(replicate, 0) + 1

            patient_site_lab_unique, replicates_exist, sample_id = database._get_sample_and_replicate_uniqueness(data_dict)

            if not patient_site_lab_unique:
                errors.append('Subject(' + data_dict['subject_id'] + ') + site(' + data_dict['site_id'] + ') + lab(' + data_dict['lab_id'] + ') found more than once in database. Something very wrong!')

            if replicates_exist:
                errors.append('Replicate already found for ' + ','.join(replicate_keys) + ': ' + ','.join([data_dict[x] for x in replicate_keys]))

        for filename, count in sorted(all_filenames.items()):
            if count > 1:
                errors.append('Reads file ' + filename + ' found ' + str(count) + ' times')


        for replicate, count in sorted(all_replicates.items()):
            if count > 1:
                errors.append('Replicate ' + ','.join(replicate_keys) + ' ' + ','.join(replicate) + ' found ' + str(count) + ' times in spreadsheet')

        return errors


    @classmethod
    def _archive_spreadsheet(cls, xlsx_file, archive_dir):
        date_string = utils.date_string_from_file_mtime(xlsx_file)
        date_directory = os.path.join(archive_dir, date_string)
        if not os.path.exists(date_directory):
            try:
                os.mkdir(date_directory)
            except:
                raise Error('Error mkdir ' + date_directory)

        xlsx_basename = os.path.basename(xlsx_file)
        existing_xlsx_files = set(os.listdir(date_directory))

        # this is unlikely to happen, but if name of archived xlsx
        # file already exists, append .N on the end (where N is smallest
        # int such that file does not already exist)
        if xlsx_basename in existing_xlsx_files:
            i = 1
            while xlsx_basename + '.' + str(i) in existing_xlsx_files:
                i += 1
            xlsx_basename += '.' + str(i)

        new_name = os.path.join(date_directory, xlsx_basename)
        utils.rsync_and_md5(xlsx_file, new_name)
        os.unlink(xlsx_file)
        return new_name


    def _import_reads_and_update_db(self):
        database = db.Db(self.db_ini_file)
        data = SpreadsheetImporter.load_data_from_spreadsheet(self.xlsx_file)
        xlsx_dir = os.path.dirname(self.xlsx_file)
        data_errors = SpreadsheetImporter._validate_data(database, data, self.dropbox_dir)

        if len(data_errors) > 0:
            raise Error('Error(s) importing spreadsheet:\n' + '\n'.join(data_errors))

        try:
            f_out = open(self.jobs_outfile, 'w')
        except:
            raise Error('Error opening file "' + self.jobs_outfile + '". Cannot continue')

        print('seqrep_id', 'sample_id', 'isolate_id', 'sequence_replicate_number', 'reads1', 'reads2', 'reads1_md5', 'reads2_md5', sep='\t', file=f_out)

        for data_dict in data:
            reads1 = os.path.join(xlsx_dir, data_dict['reads_file_1'])
            reads2 = os.path.join(xlsx_dir, data_dict['reads_file_2'])
            assert os.path.exists(reads1) and os.path.exists(reads2)
            seqrep_id, isolate_id, sample_id = database. add_one_seqrep(data_dict)
            print(seqrep_id, sample_id, isolate_id, data_dict['sequence_replicate_number'], reads1, reads2, data_dict['reads_file_1_md5'], data_dict['reads_file_2_md5'], sep='\t', file=f_out)

        f_out.close()
        xlsx_backup_file = SpreadsheetImporter._archive_spreadsheet(self.xlsx_file, self.xlsx_archive_dir)
        jobs_backup_file = xlsx_backup_file + ".import_jobs.tsv"
        assert not os.path.exists(jobs_backup_file)
        utils.rsync_and_md5(self.jobs_outfile, jobs_backup_file)
        database.commit_and_close()

        if self.db_backup_dir is not None:
            database.backup(self.db_backup_dir)


    def run(self):
        lock = lock_file.LockFile(self.lock_file)
        try:
            self._import_reads_and_update_db()
        except:
            lock.stop()
            raise Error('Error immport reads or updating database')

        lock.stop()


