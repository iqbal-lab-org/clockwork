params.help = false
params.ref_id = ""
params.references_root = ""
params.pipeline_root = ""
params.db_config_file = ""
params.dataset_name = ""
params.testing = false

if (params.testing) {
    test_opt_string = '--testing'
}
else {
    test_opt_string = ''
}

if (params.help){
    log.info"""
        Clockwork mykrobe pipeline.
        Runs mykrobe predict on all samples.

        Usage: nextflow run mykrobe.nf <arguments>

        Required arguments:

          --ref_id INT          Reference ID in database
          --references_root DIRECTORY
                                Root directory of pipeline reference files
          --pipeline_root DIRECTORY
                                Root directory of pipeline files
          --db_config_file FILENAME
                                Name of database config file

        Optional:
          --dataset_name DATASET_NAME
                                Limit to all samples in the given dataset_name

    """.stripIndent()

    exit 0
}



references_root = file(params.references_root).toAbsolutePath()
pipeline_root = file(params.pipeline_root).toAbsolutePath()
db_config_file = file(params.db_config_file)

if (!references_root.isDirectory()) {
    exit 1, "References root directory not found: ${params.pipeline_root} -- aborting"
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



/* This writes a tsv file, one line per mykrobe job to be run
   (plus a header line).
*/
process make_jobs_tsv {
    maxForks 1
    memory '1 GB'

    input:
    file db_config_file

    output:
    file jobs_tsv

    """
    clockwork mykrobe_predict_make_jobs_tsv ${dataset_name_opt} ${pipeline_root} ${db_config_file} jobs_tsv ${params.ref_id} ${references_root}
    """
}


process run_mykrobe {
    memory {params.testing ? '1 GB' : '3 GB'}
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    val tsv_fields from jobs_tsv.splitCsv(header: true, sep:'\t')

    output:
    file(done_file)
    val(tsv_fields) into run_mykrobe_out

    """
    rm -rf ${tsv_fields.output_dir} mykrobe.out
    mkdir -p ${tsv_fields.output_dir}
    reads_string=\$(echo ${tsv_fields.reads_in1} ${tsv_fields.reads_in2} | awk '{ORS=" "; for (i=1; i<=NF/2; i++) {print \$i, \$(i+NF/2)}; print "\\n"}')
    clockwork mykrobe_predict ${test_opt_string} ${tsv_fields.sample_id} mykrobe.out ${tsv_fields.reference_dir} \${reads_string}
    rsync -av mykrobe.out/ ${tsv_fields.output_dir}/
    rsync -av mykrobe.out/ ${tsv_fields.output_dir}/
    rm -r mykrobe.out
    touch done_file
    """
}


process update_database {
    maxForks 10
    time '3m'
    errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    val tsv_fields from run_mykrobe_out

    output:
    val "${tsv_fields.isolate_id}\t${tsv_fields.seqrep_id}\t${tsv_fields.sequence_replicate_number}\t${tsv_fields.pool}" into update_database_worked

    script:
    """
    clockwork db_finished_pipeline_update --reference_id ${params.ref_id} ${db_config_file} ${tsv_fields.pool} ${tsv_fields.isolate_id} ${tsv_fields.seqrep_id} ${tsv_fields.sequence_replicate_number} mykrobe_predict
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
    clockwork db_finished_pipeline_update_failed_jobs --reference_id ${params.ref_id} ${db_config_file} jobs_tsv success.seqrep_ids.txt mykrobe_predict
    """
}

