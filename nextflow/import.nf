params.help = false
params.dropbox_dir = ""
params.pipeline_root = ""
params.db_config_file = ""
params.xlsx_archive_dir = ""


if (params.help){
    log.info"""
        Clockwork import pipeline. Imports data from dropbox directory.

        Usage: nextflow run import.nf <arguments>

        Mandatory arguments:
          --dropbox dir        Directory with import spreadsheet and fastq files
          --pipeline_root      Root directory of pipeline files
          --xlsx_archive_dir   Directory where spreadsheets are archived
          --db_config_file     Name of database config file
    """.stripIndent()
    exit 0
}


def check_parameter(params, parameter_name){
   if ( !params[parameter_name]){
      error "You must specifiy a " + parameter_name
   } else {
      variable = params[parameter_name]
      return variable
   }
}


dropbox_dir = file(check_parameter(params, "dropbox_dir"))
pipeline_root = file(check_parameter(params, "pipeline_root"))
db_config_file = file(check_parameter(params, "db_config_file"))
xlsx_archive_dir = file(check_parameter(params, "xlsx_archive_dir"))

if (!dropbox_dir.isDirectory()) {
    exit 1, "Dropbox directory not found: ${params.dropbox_dir} -- aborting"
}

if (!pipeline_root.isDirectory()) {
    exit 1, "Samples directory not found: ${params.pipeline_root} -- aborting"
}

if (!xlsx_archive_dir.isDirectory() ) {
    exit 1, "Excel sheets directory not found: ${params.xlsx_archive_dir} -- aborting"
}

if (!db_config_file.exists()) {
    exit 1, "Sqlite file not found: ${params.db_config_file} -- aborting"
}

Channel
    .fromPath("${dropbox_dir}/*.{tsv,xlsx}")
    .set{parse_xlsx_channel}


process parse_xlxs {
    maxForks 1
    memory '1000 MB'

    input:
    file xlsx_file from parse_xlsx_channel

    output:
    file jobs_tsv into import_files

    """
    clockwork import_spreadsheet ${dropbox_dir} ${db_config_file} ${xlsx_file} ${xlsx_archive_dir} jobs_tsv
    """
}


import_files.splitCsv(header: true, sep:'\t').set{tsv_lines}


process import_sample {
    maxForks 20
    memory '500 MB'

    input:
    val fields from tsv_lines

    """
    clockwork import_read_pair ${db_config_file} ${pipeline_root} ${fields.seqrep_id} ${fields.isolate_id} ${fields.sample_id} ${fields.sequence_replicate_number} ${fields.reads1} ${fields.reads2} ${fields.reads1_md5} ${fields.reads2_md5}
    """
}

