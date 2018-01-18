import multiprocessing
import os

from clockwork import spreadsheet_importer, utils
from clockwork.common_data import allowed_sequencing_instruments

class Error (Exception): pass


# This is here for multiprocessing
def _check_md5(filename, expect_md5):
    got_md5 = utils.md5(filename)
    if expect_md5 == got_md5:
        return (filename, None, None)
    else:
        return (filename, expect_md5, got_md5)


class SpreadsheetValidator:
    def __init__(self, spreasheet_xlsx, data_root_dir, outfile, md5_threads=1):
        self.spreasheet_xlsx = spreasheet_xlsx
        self.data_root_dir = data_root_dir
        self.outfile = outfile
        self.md5_threads = md5_threads


    @classmethod
    def _check_no_blank_values(cls, list_of_dicts):
        errors = []
        skip_keys = {'ena_run_accession', 'ena_sample_accession'}
        for i, d in enumerate(list_of_dicts):
            for key, value in sorted(d.items()):
                if key in skip_keys:
                    continue

                if value in {None, ''}:
                    errors.append('Empty_field\t' + key + '\t' + str(i+2))
        return errors


    @classmethod
    def _get_non_unique_dict_list_values(cls, list_of_dicts, key=None, key_list=None):
        '''Given a list of dictionaries, returns a dict of values corresponding
        to the key or key_list and the indexes in the list where they appear, for values
        that are found more than once'''
        assert (key_list is None or key is None) and key_list != key
        if key_list is None:
            key_list = [key]
        counts = {}

        for i, d in enumerate(list_of_dicts):
            value_tuple = tuple(d[k] for k in key_list)
            if value_tuple not in counts:
                counts[value_tuple] = []
            counts[value_tuple].append(i)

        if key is None:
            return {x : counts[x] for x in counts if len(counts[x]) > 1}
        else:
            return {x[0]: counts[x] for x in counts if len(counts[x]) > 1}


    @classmethod
    def _check_uniqueness_of_values(cls, all_data):
        errors = []
        # Check filenames and md5s are unique
        keys_should_be_unique = ['reads_file_1', 'reads_file_1_md5', 'reads_file_2', 'reads_file_2_md5']
        for key in keys_should_be_unique:
            duplicates = SpreadsheetValidator._get_non_unique_dict_list_values(all_data, key=key)
            for dup_key, value in sorted(duplicates.items()):
                # We have a list of data, indexes start at one. But
                # input file had header line. So need to add 2 to get
                # the row number from input file
                lines_string = ','.join([str(x+2) for x in value])
                errors.append('Non_unique\t' + key + '\t' + dup_key + '\t' + lines_string)

        # Check no repeated combination of subject_id/site_id/lab_id/isolate_number/sequence_replicate_number
        ids_combo = ['subject_id', 'site_id', 'lab_id', 'isolate_number', 'sequence_replicate_number']
        duplicates = SpreadsheetValidator._get_non_unique_dict_list_values(all_data, key_list=ids_combo)
        keys_string = ','.join(ids_combo)
        for dup_key, value in sorted(duplicates.items()):
            values_string = ','.join(dup_key)
            lines_string = ','.join([str(x+2) for x in value])
            errors.append('Non_unique\t' + keys_string + '\t' + values_string + '\t' + lines_string)

        return errors


    @classmethod
    def _check_global_file_and_md5_column_intersection(cls, all_data):
        files_1 = set([d['reads_file_1'] for d in all_data])
        files_2 = set([d['reads_file_2'] for d in all_data])
        file_intersection = sorted(list(files_1.intersection(files_2)))
        md5_1 = set([d['reads_file_1_md5'] for d in all_data])
        md5_2 = set([d['reads_file_2_md5'] for d in all_data])
        md5_intersection = sorted(list(md5_1.intersection(md5_2)))
        errors = ['Filename_in_both_columns\t' + x for x in file_intersection]
        errors.extend(['md5_in_both_columns\t' +  x for x in md5_intersection])
        return errors


    @classmethod
    def _check_files_exist_and_md5(cls, all_data, data_root_dir, md5_threads=1):
        original_dir = os.getcwd()
        os.chdir(data_root_dir)
        errors = []
        files_that_exist = []
        for i, d in enumerate(all_data):
            for key in 'reads_file_1', 'reads_file_2':
                filename = os.path.join(d[key])
                if os.path.exists(filename):
                    files_that_exist.append((filename, d[key + '_md5']))
                else:
                    errors.append('File_not_found\t' + d[key] + '\t' + str(i+2))

        if md5_threads == 0:
            os.chdir(data_root_dir)
            return errors
        elif md5_threads > 1:
            pool = multiprocessing.Pool(md5_threads)
            md5_check = pool.starmap(_check_md5, files_that_exist)
        else:
            md5_check = [_check_md5(t[0], t[1]) for t in files_that_exist]

        md5_check.sort()
        errors.extend(['md5_error\t' + '\t'.join(t) for t in md5_check if t[1] is not None])
        os.chdir(original_dir)
        return errors


    @classmethod
    def _check_integers(cls, all_data, key, min_value=None, max_value=None):
        errors = []
        for i, d in enumerate(all_data):
            try:
                test_value = int(d[key])
            except:
                errors.append('Not_an_integer\t' + key + '\t' + str(d[key]) + '\t' + str(i+2))
                continue

            if (min_value is not None and test_value < min_value) or (max_value is not None and test_value > max_value):
                errors.append('Integer_out_of_range\t' + key + '\t' + str(test_value) + '\t' + str(i+2))

        return errors


    @classmethod
    def _check_instrument_model(cls, all_data):
        errors = []
        for i, d in enumerate(all_data):
            if d['instrument_model'] not in allowed_sequencing_instruments:
                errors.append('Unknown_instrument_model\t' + d['instrument_model'] + '\t' + str(i+2))
        return errors


    def run(self):
        all_data = spreadsheet_importer.SpreadsheetImporter.load_data_from_spreadsheet(self.spreasheet_xlsx)
        errors = SpreadsheetValidator._check_no_blank_values(all_data)
        errors.extend(SpreadsheetValidator._check_uniqueness_of_values(all_data))
        errors.extend(SpreadsheetValidator._check_global_file_and_md5_column_intersection(all_data))
        errors.extend(SpreadsheetValidator._check_files_exist_and_md5(all_data, self.data_root_dir, self.md5_threads))
        errors.extend(SpreadsheetValidator._check_integers(all_data, 'isolate_number', min_value=1))
        errors.extend(SpreadsheetValidator._check_integers(all_data, 'sequence_replicate_number', min_value=1))
        errors.extend(SpreadsheetValidator._check_integers(all_data, 'submit_to_ena', min_value=0, max_value=1))
        errors.extend(SpreadsheetValidator._check_integers(all_data, 'ena_on_hold', min_value=0, max_value=1))
        errors.extend(SpreadsheetValidator._check_instrument_model(all_data))

        with open(self.outfile, 'w') as f:
            for line in errors:
                print(line, file=f)

