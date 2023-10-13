import datetime
import dateutil.parser
import re

import openpyxl

from clockwork import spreadsheet_importer


date_iso_regex = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")

columns = [
    "subject_id",
    "site_id",
    "lab_id",
    "isolate_number",
    "sequence_replicate_number",
    "submission_date",
    "reads_file_1",
    "reads_file_1_md5",
    "reads_file_2",
    "reads_file_2_md5",
    "dataset_name",
    "instrument_model",
    "ena_center_name",
    "submit_to_ena",
    "ena_on_hold",
    "ena_run_accession",
    "ena_sample_accession",
]


def looks_like_date(date_string):
    """Returns true iff the date_string looks like a date of the form YYYY-MM-DD.
       We're always getting ISO date strings from the datetime module, but
       check anyway"""
    if not date_iso_regex.match(date_string):
        return False

    year, month, day = date_string[:4], date_string[5:7], date_string[8:]

    try:
        datetime.date(int(year), int(month), int(day))
    except:
        return False

    return True


def row_data_dict_to_string(row_data):
    return ", ".join([str(row_data[columns[i]]) for i in range(len(columns))])


def load_data_from_spreadsheet(infile):
    if infile.endswith(".xlsx"):
        workbook = openpyxl.load_workbook(infile)
        sheet = workbook.active
        fields_list = []
        for row in sheet:
            fields_list.append([str(x.value) for x in row if x.value is not None])
    elif infile.endswith(".tsv"):
        with open(infile) as f:
            fields_list = []
            for line in f:
                fields_list.append([x.strip() for x in line.rstrip().split("\t")])

    data = []

    for line_counter, fields in enumerate(fields_list):
        if line_counter == 0:
            if fields != columns:
                fields_str = ", ".join([str(x) for x in fields])
                expected_str = ", ".join(columns)
                raise Exception(
                    'Error in first line of spreadsheet "'
                    + infile
                    + '". Column names not correct.\nExpected: '
                    + expected_str
                    + "\nGot:      "
                    + fields_str
                )
        elif len(fields) != len(columns):
            raise Exception(
                "Wrong number of fields in line "
                + str(line_counter)
                + " of spreadsheet: "
                + ", ".join([str(x) for x in fields])
            )
        else:
            new_data = dict(zip(columns, fields))

            for key in [
                "ena_run_accession",
                "ena_sample_accession",
                "reads_file_1_md5",
                "reads_file_2_md5",
            ]:
                if new_data[key] in {0, "0"}:
                    new_data[key] = None

            try:
                date = dateutil.parser.parse(new_data["submission_date"]).date()
            except:
                raise Exception(
                    "Date format error: "
                    + spreadsheet_importer.SpreadsheetImporter.row_data_dict_to_string(
                        new_data
                    )
                )

            assert looks_like_date(str(date))
            new_data["submission_date"] = date
            data.append(new_data)

    return data
