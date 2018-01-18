params.help = false
params.ref_id = ""
params.references_root = ""
params.pipeline_root = ""
params.db_config_file = ""
params.dataset_name = ""
params.ref_fasta = ""
params.reads_in1 = ""
params.reads_in2 = ""
params.output_dir = ""
params.max_forks_samtools_qc = 100
params.max_forks_fastqc = 100


if (params.help){
    log.info"""
        Clockwork QC pipeline.
        Can either be run on one pair of fastq files,
        or on multiple samples using a database.

        Usage: nextflow run qc.nf <arguments>

        Required if using a database:
          --ref_id INT          Reference ID in database
          --references_root DIRECTORY
                                Root directory of pipeline reference files
          --pipeline_root DIRECTORY
                                Root directory of pipeline files
          --db_config_file FILENAME
                                Name of database config file

        Optional if using a database:
          --dataset_name DATASET_NAME
                                Limit to all samples in the given dataset

        Required if running on one pair of FASTQ files:
          --ref_fasta FILENAME  Name of reference fasta file
          --reads_in1 FILENAME  Name of forwards reads file
          --reads_in2 FILENAME  Name of reverse reads file
          --output_dir DIRECTORY
                                Name of output directory (must not
                                already exist)

        Other options:
          --max_forks_fastqc INT
                                Limit number of concurrent fastqc jobs [100]
          --max_forks_samtools_qc INT
                                Limit number of concurrent samtools_qc jobs [100]
    """.stripIndent()

    exit 0
}


using_db_input = params.pipeline_root != "" && params.references_root != "" && params.db_config_file != "" && params.ref_id != ""
using_reads_input = params.ref_fasta != "" && params.reads_in1 != "" && params.reads_in2 != "" && params.output_dir != ""

if (using_db_input == using_reads_input) {
    exit 1, "Must use either all of ref_fasta,reads_in1,reads_in2,output_dir, or all of ref_id,references_root,pipeline_root,db_config_file -- aborting"
}

if (using_db_input) {
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
}
else if (using_reads_input) {
    ref_fasta = file(params.ref_fasta).toAbsolutePath()
    reads_in1 = file(params.reads_in1).toAbsolutePath()
    reads_in2 = file(params.reads_in2).toAbsolutePath()
    output_dir = file(params.output_dir).toAbsolutePath()

    if (!ref_fasta.exists()) {
        exit 1, "Reference fasta file not found: ${params.ref_fasta} -- aborting"
    }

    if (!reads_in1.exists()) {
        exit 1, "Reads file 1 not found: ${params.reads_in1} -- aborting"
    }

    if (!reads_in2.exists()) {
        exit 1, "Reads file 2 not found: ${params.reads_in2} -- aborting"
    }

    if (output_dir.exists()) {
        exit 1, "Output directory already exists: ${params.output_dir} -- aborting"
    }
    if (!output_dir.mkdirs()) {
        exit 1, "Cannot create output directory: ${params.output_dir} -- check file system permissions"
    }
}
else {
    exit 1, "Must use either all of reads_in1,reads_in2,output_dir, or all of pipeline_root,db_config_file -- aborting"
}


/* This writes a tsv file, one line per QC job to be run
   (plus a header line).
   If we're using the database, then we need to gather the jobs from it.
   If we're not using the database, then manually write the tsv file, which
   is then simply the header line, plus a line that defines the single
   remove contamination job to be run. Put dots for the database-related
   fields, which get ignored.
*/
process make_jobs_tsv {
    maxForks 1
    memory '1 GB'

    input:
    if (using_db_input) {
        file db_config_file
    }

    output:
    file jobs_tsv
    file jobs_tsv into jobs_tsv_channel

    script:
    if (using_db_input)
        """
        clockwork qc_make_jobs_tsv ${dataset_name_opt} ${pipeline_root} ${db_config_file} jobs_tsv ${params.ref_id} ${references_root}
        """
    else if (using_reads_input)
        """
        echo "reads_in1 reads_in2 output_dir sample_id seqrep_id isolate_id sequence_replicate_number reference_id ref_fasta" > tmp.tsv
        echo "${reads_in1} ${reads_in2} ${output_dir} . . . . . ${ref_fasta}" >> tmp.tsv
        sed 's/ /\t/g' tmp.tsv > jobs_tsv
        rm tmp.tsv
        """
    else
        error "Must get input from database or provide reads files -- aborting"
}


jobs_tsv_channel.splitCsv(header: true, sep:'\t').into{samtools_qc_tsv_fields; fastqc_tsv_fields}

process samtools_qc {
    maxForks params.max_forks_samtools_qc
    memory '4 GB'
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    val tsv_fields from samtools_qc_tsv_fields

    output:
    val tsv_fields into update_database_channel_from_samtools_qc

    """
    echo ${tsv_fields}
    rm -rf ${tsv_fields.output_dir}/samtools_qc
    clockwork samtools_qc ${tsv_fields.ref_fasta} ${tsv_fields.reads_in1} ${tsv_fields.reads_in2} samtools_qc
    rsync -a samtools_qc ${tsv_fields.output_dir}/
    """
}


process fastqc {
    maxForks params.max_forks_fastqc
    memory '2 GB'
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    val tsv_fields from fastqc_tsv_fields

    output:
    val tsv_fields into update_database_channel_from_fastqc

    """
    echo ${tsv_fields}
    rm -rf ${tsv_fields.output_dir}/fastqc
    clockwork fastqc fastqc ${tsv_fields.reads_in1} ${tsv_fields.reads_in2}
    rsync -a fastqc ${tsv_fields.output_dir}/
    """
}


if (using_db_input) {
    /*
      Ensure the same sample is input from each of samtools_qc and fastqc into
      the update_database process by phasing the two channels, by taking
      pairs where the seqrep_id is the same.
    */
    phased_tsv_fields = update_database_channel_from_samtools_qc.phase(update_database_channel_from_fastqc){ it -> it.seqrep_id }


    process update_database {
        maxForks 10
        memory '1 GB'

        when:
        using_db_input

        input:
        val tsv_fields from phased_tsv_fields

        output:
        val "${tsv_fields[0].isolate_id}\t${tsv_fields[0].seqrep_id}\t${tsv_fields[0].sequence_replicate_number}" into update_database_worked

        script:
        """
        clockwork db_finished_pipeline_update --pipeline_root ${pipeline_root} ${db_config_file} 0 ${tsv_fields[0].isolate_id} ${tsv_fields[0].seqrep_id} ${tsv_fields[0].sequence_replicate_number} qc
        """
    }


    process update_database_failed_jobs {
        maxForks 1
        memory '1 GB'

        input:
        val seqrep_worked from update_database_worked.collect()
        file jobs_tsv

        script:
        """
        cat <<"EOF" > success.seqrep_ids.txt
isolate_id	seqrep_id	sequence_replicate_number
${seqrep_worked.join('\n')}
EOF
        clockwork db_finished_pipeline_update_failed_jobs ${db_config_file} jobs_tsv success.seqrep_ids.txt qc
        """
    }
}
