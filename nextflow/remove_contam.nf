params.help = false
params.reads_in1 = ""
params.reads_in2 = ""
params.outprefix = ""
params.ref_metadata_tsv = ""
params.pipeline_root = ""
params.references_root = ""
params.db_config_file = ""
params.dataset_name = ""
params.ref_fasta = ""
params.ref_id = ""
params.testing = false
params.max_forks_map_reads = 100
params.max_forks_sam_to_fastq_files = 100
params.mapping_threads = 1


if (params.help){
    log.info"""
        Clockwork remove_contam pipeline. Removes reads that are contaminated.
        Can be run on one pair of fastq files, or on multiple samples using a database.

        Usage: nextflow run remove_contam.nf <arguments>

        Required arguments:

        Required if using a database:
          --ref_id INT          Reference ID in database
          --references_root DIRECTORY
                                Root directory of pipeline reference files
          --pipeline_root DIRECTORY
                                Root directory of pipeline files
          --db_config_file FILENAME
                                Name of database config file

        Optional if using a database:
          --dataset_name DATASET_NAME    Limit to all samples in the given dataset_name

        Required if running on one pair of FASTQ files:
          --ref_fasta FILENAME  Name of reference fasta file
          --ref_metadata_tsv FILENAME
                                Name of metadata TSV file, each sequence in
                                reference fasta must be in this file
          --reads_in1 FILENAME  Name of forwards reads file
          --reads_in2 FILENAME  Name of reverse reads file
          --outprefix PATH      Prefix of name of output files

        Other options:
          --mapping_threads INT
                                Number of threads used by each read mapping process.
                                This option is only recommended if running on one
                                pair of FASTQ files [${params.mapping_threads}]
          --max_forks_map_reads INT
                                Limit number of concurrent map_reads jobs [${params.max_forks_map_reads}]
          --max_forks_sam_to_fastq_files INT
                                Limit number of concurrent
                                sam_to_fastq_files jobs [${params.max_forks_sam_to_fastq_files}]
    """.stripIndent()

    exit 0
}


using_db_input = params.pipeline_root != "" && params.references_root != "" && params.db_config_file != "" && params.ref_id != ""
using_reads_input = params.ref_fasta != "" && params.ref_metadata_tsv != "" && params.reads_in1 != "" && params.reads_in2 != "" && params.outprefix != ""

if (using_db_input == using_reads_input) {
    exit 1, "Must use either all of ref_fasta,ref_metadata_tsv,reads_in1,reads_in2,outprefix, or all of ref_id,references_root,pipeline_root,db_config_file -- aborting"
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
    ref_metadata_tsv = file(params.ref_metadata_tsv).toAbsolutePath()
    reads_in1 = file(params.reads_in1).toAbsolutePath()
    reads_in2 = file(params.reads_in2).toAbsolutePath()
    outprefix = file(params.outprefix).toAbsolutePath()

    if (!ref_fasta.exists()) {
        exit 1, "Reference fasta file not found: ${params.ref_fasta} -- aborting"
    }

    if (!ref_metadata_tsv.exists()) {
        exit 1, "Reference metadata file not found: ${params.ref_metadata_tsv} -- aborting"
    }

    if (!reads_in1.exists()) {
        exit 1, "Reads file 1 not found: ${params.reads_in1} -- aborting"
    }

    if (!reads_in2.exists()) {
        exit 1, "Reads file 2 not found: ${params.reads_in2} -- aborting"
    }

}
else {
    exit 1, "Must use either all of reads_in1,reads_in2,outprefix, or all of pipeline_root,db_config_file -- aborting"
}




/* This writes a tsv file, one line per remove contamination job to be run
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

    script:
    if (using_db_input)
        """
        clockwork remove_contam_make_jobs_tsv ${dataset_name_opt} ${pipeline_root} ${db_config_file} jobs_tsv ${params.ref_id} ${references_root}
        """
    else if (using_reads_input)
        """
        echo "reads_in1 reads_in2 counts_tsv reads_contam1 reads_contam2 reads_remove_contam1 reads_remove_contam2 sample_id seqrep_id isolate_id sequence_replicate_number reference_id ref_fasta contam_tsv" > tmp.tsv
        echo "${reads_in1} ${reads_in2} ${outprefix}.counts.tsv ${outprefix}.contam.1.fq.gz ${outprefix}.contam.2.fq.gz ${outprefix}.remove_contam.1.fq.gz ${outprefix}.remove_contam.2.fq.gz . . . . . ${ref_fasta} ${ref_metadata_tsv}" >> tmp.tsv
        sed 's/ /\t/g' tmp.tsv > jobs_tsv
        rm tmp.tsv
        """
    else
        error "Must get input from database or provide reads files -- aborting"
}


process map_reads {
    maxForks params.max_forks_map_reads
    memory {params.testing ? '3 GB' : '9 GB'}
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    val tsv_fields from jobs_tsv.splitCsv(header: true, sep:'\t')

    output:
    set file(contam_sam), val(tsv_fields) into sam_to_fastq_files_in

    script:
    """
    clockwork map_reads --threads ${params.mapping_threads} --unsorted_sam sample_name ${tsv_fields.ref_fasta} contam_sam ${tsv_fields.reads_in1} ${tsv_fields.reads_in2}
    """
}


process sam_to_fastq_files {
    maxForks params.max_forks_sam_to_fastq_files
    memory '2 GB'
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    set file(contam_sam), val(tsv_fields) from sam_to_fastq_files_in

    output:
    val tsv_fields into update_database_tsv_fields_channel

    script:
    """
    clockwork remove_contam --contam_out_1 ${tsv_fields.reads_contam1} --contam_out_2 ${tsv_fields.reads_contam2} ${tsv_fields.contam_tsv} contam_sam ${tsv_fields.counts_tsv} ${tsv_fields.reads_remove_contam1} ${tsv_fields.reads_remove_contam2}
    """
}


if (using_db_input) {
    process update_database {
        maxForks 10
        memory '1 GB'

        input:
        val tsv_fields from update_database_tsv_fields_channel

        output:
        val "${tsv_fields.isolate_id}\t${tsv_fields.seqrep_id}\t${tsv_fields.sequence_replicate_number}" into update_database_worked

        """
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
}
