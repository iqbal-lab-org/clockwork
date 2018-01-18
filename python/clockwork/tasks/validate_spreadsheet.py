from clockwork import spreadsheet_validator

def run(options):
    validator = spreadsheet_validator.SpreadsheetValidator(
        options.infile,
        options.data_dir,
        options.outfile,
        md5_threads=options.threads,
    )
    validator.run()

