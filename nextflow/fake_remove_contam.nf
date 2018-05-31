params.help = false
params.pipeline_root = ""
params.db_config_file = ""
params.dataset_name = ""
params.max_forks = 20


if (params.help){
    log.info"""
        Clockwork fake_remove_contam pipeline. Fakes a run of remove_contam pipeline,
        updating database as if remove_Contam had really been run. Use this if
        you do not want to run remove_contam pipeline. It results in "decontaminated"
        reads that are the same as the original reads.

        Usage: nextflow run fake_remove_contam.nf <arguments>

        Required arguments:
          --pipeline_root DIRECTORY
                                Root directory of pipeline files
          --db_config_file FILENAME
                                Name of database config file

        Optional:
          --dataset_name DATASET_NAME
                                Limit to all samples in the given dataset_name
          --max_forks INT       Limit number of concurrent jobs [${params.max_forks}]
    """.stripIndent()

    exit 0
}


pipeline_root = file(params.pipeline_root).toAbsolutePath()
db_config_file = file(params.db_config_file)

if (!pipeline_root.isDirectory()) {
    exit 1, "Pipeline root directory not found: ${params.pipeline_root} -- aborting"
}

if (!db_config_file.exists()) {
    exit 1, "Sqlite file not found: ${params.db_config_file} -- aborting"
}

if (params.dataset_name != ""){
    dataset_name_opt = "--dataset_name " + params.dataset_name
}
else {
    dataset_name_opt = ""
}


process make_jobs_tsv {
    maxForks 1
    memory '1 GB'

    input:
    file db_config_file

    output:
    file jobs_tsv

    script:
    """
    clockwork fake_remove_contam_make_jobs_tsv ${dataset_name_opt} ${pipeline_root} ${db_config_file} jobs_tsv
    """
}


process fake_remove_contam {
    maxForks params.max_forks
    memory '500 MB'
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    val tsv_fields from jobs_tsv.splitCsv(header: true, sep:'\t')

    output:
    val "${tsv_fields.isolate_id}\t${tsv_fields.seqrep_id}\t${tsv_fields.sequence_replicate_number}" into update_database_worked

    script:
    """
    clockwork fake_remove_contam ${tsv_fields.reads_in1} ${tsv_fields.reads_in2} ${tsv_fields.reads_remove_contam1} ${tsv_fields.reads_remove_contam2} ${tsv_fields.counts_tsv}
    clockwork db_finished_pipeline_update --pipeline_root ${pipeline_root} ${db_config_file} 0 ${tsv_fields.isolate_id} ${tsv_fields.seqrep_id} ${tsv_fields.sequence_replicate_number} remove_contam
    """
}


process update_database_failed_jobs {
    maxForks 1
    memory '1 GB'

    input:
    val seqrep_worked from update_database_worked.collect()
    file jobs_tsv

    """
    cat <<"EOF" > success.seqrep_ids.txt
isolate_id	seqrep_id	sequence_replicate_number
${seqrep_worked.join('\n')}
EOF
    clockwork db_finished_pipeline_update_failed_jobs ${db_config_file} jobs_tsv success.seqrep_ids.txt remove_contam
    """
}

