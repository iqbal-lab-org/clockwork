from clockwork import spreadsheet_importer

def run(options):
    importer = spreadsheet_importer.SpreadsheetImporter(
        options.dropbox_dir,
        options.xlsx_file,
        options.db_config_file,
        options.xls_archive_dir,
        options.jobs_outfile,
        db_backup_dir=options.db_backup_dir,
    )
    importer.run()

