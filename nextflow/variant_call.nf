params.help = false
params.ref_id = ""
params.references_root = ""
params.pipeline_root = ""
params.db_config_file = ""
params.dataset_name = ""
params.ref_dir = ""
params.reads_in1 = ""
params.reads_in2 = ""
params.output_dir = ""
params.sample_name = ""
params.testing = false
params.cortex_mem_height = 22
params.max_forks_trim_reads = 100
params.max_forks_map_reads = 100
params.max_forks_samtools = 100
params.max_forks_cortex = 100
params.max_forks_combine_variant_calls = 100
params.minos_max_read_length = 200
params.truth_ref = ""


if (params.help){
    log.info"""
        Clockwork variant_call pipeline.
        Calls variants using Cortex and BWA MEM + Samtools.
        Can either be run on one pair of fastq files,
        or on multiple samples using a database.

        Usage: nextflow run variant_call.nf <arguments>

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
          --dataset_name DATASET_NAME
                                Limit to all samples in the given dataset_name

        Required if running on one pair of FASTQ files:
          --ref_dir DIRECTORY   Name of reference data directory, made by:
                                clockwork prepare_reference
          --reads_in1 FILENAME  Name of forwards reads file
          --reads_in2 FILENAME  Name of reverse reads file
          --output_dir DIRECTORY
                                Name of output directory (must not
                                already exist)
          --sample_name STRING  Sample name - this will be put into
                                output VCF files

        Other options:
          --max_forks_trim_reads INT
                                Limit number of concurrent trim_reads jobs [${params.max_forks_trim_reads}]
          --max_forks_map_reads INT
                                Limit number of concurrent map_reads jobs [${params.max_forks_map_reads}]
          --max_forks_samtools INT
                                Limit number of concurrent samtools jobs [${params.max_forks_samtools}]
          --max_forks_cortex INT
                                Limit number of concurrent cortex jobs [${params.max_forks_cortex}]
          --minos_max_read_length INT
                                --max_read_length passed to 'minos adjudicate'
                                when combining variant calls [${params.minos_max_read_length}]
          --truth_ref FILENAME
                                If this option is used, checks the called variants using this
                                provided FASTA file as the "truth" genome
    """.stripIndent()

    exit 0
}


using_db_input = params.pipeline_root != "" && params.references_root != "" && params.db_config_file != "" && params.ref_id != ""
using_reads_input = params.ref_dir != "" && params.reads_in1 != "" && params.reads_in2 != "" && params.output_dir != "" && params.sample_name != ""

if (using_db_input == using_reads_input) {
    exit 1, "Must use either all of ref_dir,reads_in1,reads_in2,output_dir,sample_name or all of ref_id,references_root,pipeline_root,db_config_file -- aborting"
}

if (params.truth_ref) {
    truth_ref = file(params.truth_ref).toAbsolutePath()
    if (!truth_ref.exists()) {
        exit 1, "Truth reference file not found: ${params.truth_ref} -- aborting"
    }
}
else {
    truth_ref = ""
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
    reference_dir = file(params.ref_dir).toAbsolutePath()
    reads_in1 = file(params.reads_in1).toAbsolutePath()
    reads_in2 = file(params.reads_in2).toAbsolutePath()
    output_dir = file(params.output_dir).toAbsolutePath()

    if (!reference_dir.exists()) {
        exit 1, "Reference files directory not found: ${params.ref_dir} -- aborting"
    }

    if (!reads_in1.exists()) {
        exit 1, "Reads file 1 not found: ${params.reads_in1} -- aborting"
    }

    if (!reads_in2.exists()) {
        exit 1, "Reads file 2 not found: ${params.reads_in2} -- aborting"
    }

    if (params.sample_name == "") {
        exit 1, "Must supply sample name with --sample_name foo -- aborting"
    }
}
else {
    exit 1, "Must use either all of reads_in1,reads_in2,output_dir, or all of pipeline_root,db_config_file -- aborting"
}



/* This writes a tsv file, one line per variant calling contamination job to be run
   (plus a header line).
   If we're using the database, then we need to gather the jobs from it.
   If we're not using the database, then manually write the tsv file, which
   is then simply the header line, plus a line that defines the single
   variant calling job to be run. Put dots for the database-related
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
        clockwork variant_call_make_jobs_tsv ${dataset_name_opt} ${pipeline_root} ${db_config_file} jobs_tsv ${params.ref_id} ${references_root}
        """
    else if (using_reads_input)
        """
        echo "reads_in1 reads_in2 output_dir sample_id pool isolate_id seqrep_id sequence_replicate_number reference_id reference_dir" > tmp.tsv
        echo "${reads_in1} ${reads_in2} ${output_dir} ${params.sample_name} . . . . . ${reference_dir}" >> tmp.tsv
        sed 's/ /\t/g' tmp.tsv > jobs_tsv
        rm tmp.tsv
        """
    else
        error "Must get input from database or provide reads files -- aborting"
}


process trim_reads {
    maxForks params.max_forks_trim_reads
    memory {params.testing ? '3 GB' : '4 GB'}
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    val tsv_fields from jobs_tsv.splitCsv(header: true, sep:'\t')

    output:
    set file(trimmed_reads_dir), val(tsv_fields) into map_reads_in

    """
    rm -fr trimmed_reads_dir
    mkdir trimmed_reads_dir
    reads_string=\$(echo ${tsv_fields.reads_in1} ${tsv_fields.reads_in2} | awk '{ORS=" "; for (i=1; i<=NF/2; i++) {print \$i, \$(i+NF/2)}; print "\\n"}')
    echo "read_string: \${reads_string}"
    clockwork trim_reads trimmed_reads_dir/reads \${reads_string}
    """
}


process map_reads {
    maxForks params.max_forks_map_reads
    memory {params.testing ? '3 GB' : '4 GB'}
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    set file(trimmed_reads_dir), val(tsv_fields) from map_reads_in

    output:
    set(file("rmdup.bam"), val(tsv_fields)) into call_vars_samtools_in
    set(file("rmdup.bam"), val(tsv_fields)) into call_vars_cortex_in

    script:
    if (using_db_input)
        """
        reads_string=\$(ls trimmed_reads_dir/reads*fq | sort)
        clockwork map_reads ${tsv_fields.sample_id} ${tsv_fields.reference_dir}/ref.fa rmdup.bam \${reads_string}
        """
    else if (using_reads_input)
        """
        reads_string=\$(ls trimmed_reads_dir/reads*fq | sort)
        clockwork map_reads ${params.sample_name} ${tsv_fields.reference_dir}/ref.fa rmdup.bam \${reads_string}
        """
    else
        error "Must get input from database or provide reads files -- aborting"
}


process call_vars_samtools {
    maxForks params.max_forks_samtools
    memory {params.testing ? '3 GB' : '4 GB'}
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    set(file("rmdup.bam"), val(tsv_fields)) from call_vars_samtools_in

    output:
    set(val(tsv_fields), file("rmdup.bam")) into combine_variant_calls_channel_from_samtools

    """
    samtools mpileup -ugf ${tsv_fields.reference_dir}/ref.fa "rmdup.bam" | bcftools call -vm -O v -o samtools.vcf
    rm -rf ${tsv_fields.output_dir}/samtools/
    mkdir -p ${tsv_fields.output_dir}/samtools/
    rsync --copy-links samtools.vcf rmdup.bam ${tsv_fields.output_dir}/samtools/
    """
}


process call_vars_cortex {
    maxForks params.max_forks_cortex
    memory {params.testing ? '3 GB' : '12 GB'}
    if (using_db_input) {
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3
    }

    input:
    set(file("rmdup.bam"), val(tsv_fields)) from call_vars_cortex_in

    output:
    val(tsv_fields) into combine_variant_calls_channel_from_cortex

    script:
    if (using_db_input)
        """
        clockwork cortex --mem_height ${params.cortex_mem_height} ${tsv_fields.reference_dir} "rmdup.bam" cortex ${tsv_fields.sample_id}
        rm -fr ${tsv_fields.output_dir}/cortex/
        mkdir -p  ${tsv_fields.output_dir}
        rsync -a cortex ${tsv_fields.output_dir}/
        """
    else if (using_reads_input)
        """
        clockwork cortex --mem_height ${params.cortex_mem_height} ${tsv_fields.reference_dir} "rmdup.bam" cortex ${params.sample_name}
        rm -fr ${tsv_fields.output_dir}/cortex/
        mkdir -p  ${tsv_fields.output_dir}
        rsync -a cortex ${tsv_fields.output_dir}/
        """
    else
        error "Must get input from database or provide reads files -- aborting"
}


/*
  Ensure the same sample is input from each of call_vars_samtools and call_vars_cortex into
  the combine_variant_calls process by phasing the two channels, by taking
  pairs where the tsv fields are all the same
*/
call_vars_joined_channel = combine_variant_calls_channel_from_samtools.join(combine_variant_calls_channel_from_cortex, by:0)


process combine_variant_calls_minos {
    maxForks params.max_forks_combine_variant_calls
    memory { 5.GB * task.attempt }
    if (using_db_input) {
       errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
       maxRetries 3
    }

    input:
    set (val(tsv_fields), file("rmdup.bam")) from call_vars_joined_channel

    output:
    val tsv_fields into combine_variant_calls_minos_out

    script:
    """
    rm -rf ${tsv_fields.output_dir}/minos
    cortex_vcf=\$(find ${tsv_fields.output_dir}/cortex/cortex.out/vcfs/ -name "*FINAL*raw.vcf")
    minos adjudicate --max_read_length ${params.minos_max_read_length} --force --reads rmdup.bam minos ${tsv_fields.reference_dir}/ref.fa ${tsv_fields.output_dir}/samtools/samtools.vcf \$cortex_vcf
    rsync -av minos/ ${tsv_fields.output_dir}/minos
    rm -r minos
    """
}



if (truth_ref) {
    process check_calls_using_truth_ref {
        memory '2 GB'
        errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
        maxRetries 3

        input:
        val tsv_fields from combine_variant_calls_minos_out

        output:
        val tsv_fields into check_calls_using_truth_ref_out

        script:
        """
        samtools_vcf=\$(find ${tsv_fields.output_dir}/samtools/ -name samtools.vcf)
        if [ \$samtools_vcf -a -f \$samtools_vcf ]; then
            minos check_with_ref \$samtools_vcf ${tsv_fields.reference_dir}/ref.fa ${truth_ref} \$samtools_vcf.check
        fi

        cortex_vcf=\$(find ${tsv_fields.output_dir}/cortex/cortex.out/vcfs/ -name "*FINAL*raw.vcf")
        if [ \$cortex_vcf -a -f \$cortex_vcf ]; then
            minos check_with_ref \$cortex_vcf ${tsv_fields.reference_dir}/ref.fa ${truth_ref} \$cortex_vcf.check
        fi

        minos_vcf=\$(find ${tsv_fields.output_dir}/minos/ -name final.vcf)
        if [ \$minos_vcf -a -f \$minos_vcf ]; then
            minos check_with_ref \$minos_vcf ${tsv_fields.reference_dir}/ref.fa ${truth_ref} \$minos_vcf.check
        fi
        """
    }
}


if (using_db_input) {
    process update_database {
        maxForks 10
        time '3m'
        if (using_db_input) {
            errorStrategy {task.attempt < 3 ? 'retry' : 'ignore'}
            maxRetries 3
        }

        if (truth_ref) {
            input_channel = check_calls_using_truth_ref_out
        }
        else{
           input_channel = combine_variant_calls_minos_out
        }

        when:
        using_db_input


        input:
        val tsv_fields from input_channel

        output:
        val "${tsv_fields.isolate_id}\t${tsv_fields.seqrep_id}\t${tsv_fields.sequence_replicate_number}\t${tsv_fields.pool}" into update_database_worked

        script:
        """
        clockwork db_finished_pipeline_update --reference_id ${params.ref_id} ${db_config_file} ${tsv_fields.pool} ${tsv_fields.isolate_id} ${tsv_fields.seqrep_id} ${tsv_fields.sequence_replicate_number} variant_call
        """
    }


    process update_database_failed_jobs {
        maxForks 1

        when:
        using_db_input

        input:
        val seqrep_worked from update_database_worked.collect()
        file jobs_tsv

        script:
        """
        cat <<"EOF" > success.seqrep_ids.txt
isolate_id	seqrep_id	sequence_replicate_number	pool
${seqrep_worked.join('\n')}
EOF
        clockwork db_finished_pipeline_update_failed_jobs --reference_id ${params.ref_id} ${db_config_file} jobs_tsv success.seqrep_ids.txt variant_call
        """
    }
}
