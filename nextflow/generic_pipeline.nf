#!/usr/bin/env nextflow

if (params.help){
    log.info"""
        Clockwork generic pipeline. Runs a user-provided
        script on every sample. The arguments provided to the
        script are:
          outdir
          sample_id
          pool
          isolate_id
          seqrep_id
          seqrep_number
          reads_fwd1
          reads_fwd2
        and optionally more pairs of fwd/rev reads filenames at the end.

        Usage: nextflow run generic_pipeline.nf <arguments>

        Required arguments:
          --pipeline_name STRING
                                Name of pipeline to be run
          --script FILENAME     Name of script to be run
          --pipeline_root DIRECTORY
                                Root directory of pipeline files
          --db_config_file FILENAME
                                Name of database config file

        Optional:
          --dataset_name DATASET_NAME
                                Limit to all samples in the given dataset_name

          --max_forks INT       Max number of jobs to run at the same time [${params.max_forks}]
                                Default is 100.
          --max_ram FLOAT       RAM limit in GB for each job. Default is [${params.max_ram}]
    """.stripIndent()

    exit 0
}



script = file(params.script).toAbsolutePath()
pipeline_root = file(params.pipeline_root).toAbsolutePath()
db_config_file = file(params.db_config_file)

if (params.pipeline_name == "") {
    exit 1, "Must use --pipeline_name -- aborting"
}

if (!script.exists()) {
    exit 1, "Script from option --script not found: ${params.script} -- aborting"
}

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


// This writes a tsv file, one line per job to be run (plus a header line).
process make_jobs_tsv {
    maxForks 1
    memory '1 GB'

    input:
    file db_config_file

    output:
    file jobs_tsv

    script:
    """
    clockwork generic_pipeline_make_jobs_tsv ${dataset_name_opt} ${params.pipeline_name} ${pipeline_root} ${db_config_file} jobs_tsv
    """
}


process run_script {
    maxForks params.max_forks
    memory params.max_ram + ' GB'
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    val tsv_fields from jobs_tsv.splitCsv(header: true, sep:'\t')

    output:
    //file done_file
    set file(done_file), val(tsv_fields) into run_script_out_channel

    """
    rm -fr ${tsv_fields.output_dir}
    mkdir -p ${tsv_fields.output_dir}
    reads_string=\$(echo ${tsv_fields.reads_in1} ${tsv_fields.reads_in2} | awk '{ORS=" "; for (i=1; i<=NF/2; i++) {print \$i, \$(i+NF/2)}; print "\\n"}')
    echo "read_string: \${reads_string}"
    ${script} ${tsv_fields.output_dir} ${tsv_fields.sample_id} ${tsv_fields.pool} ${tsv_fields.isolate_id} ${tsv_fields.seqrep_id} ${tsv_fields.sequence_replicate_number} \${reads_string}
    touch done_file
    """
}

process update_database {
    maxForks 10
    time '3m'
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    set file(done_file), val(tsv_fields) from run_script_out_channel

    output:
    val "${tsv_fields.isolate_id}\t${tsv_fields.seqrep_id}\t${tsv_fields.sequence_replicate_number}\t${tsv_fields.pool}" into update_database_worked

    """
    clockwork db_finished_pipeline_update ${db_config_file} ${tsv_fields.pool} ${tsv_fields.isolate_id} ${tsv_fields.seqrep_id} ${tsv_fields.sequence_replicate_number} ${params.pipeline_name}
    """
}

process update_database_failed_jobs {
    maxForks 1

    input:
    val seqrep_worked from update_database_worked.collect()
    file jobs_tsv

    script:
    """
    cat <<"EOF" > success.seqrep_ids.txt
isolate_id	seqrep_id	sequence_replicate_number	pool
${seqrep_worked.join('\n')}
EOF
    clockwork db_finished_pipeline_update_failed_jobs ${db_config_file} jobs_tsv success.seqrep_ids.txt ${params.pipeline_name}
    """
}

