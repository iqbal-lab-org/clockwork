version = 1

tables = {
    'Isolate': [
        ('isolate_id', 'integer'),
        ('sample_id', 'integer'),
        ('isolate_number_from_lab', 'text'),
        ('pool_sequence_replicates', 'integer'),
        ('ena_experiment_accession', 'text'),
    ],

    'Pipeline': [
        ('isolate_id', 'integer'),
        ('seqrep_id', 'integer'),
        ('seqrep_pool', 'text'),
        ('version', 'text'),
        ('pipeline_name', 'text'),
        ('status', 'integer'),
        ('reference_id', 'integer'),
    ],

    'QC': [
        ('seqrep_id', 'integer'),
        ('pipeline_version', 'text'),
        ('fastqc1_gc', 'float'),
        ('fastqc1_adapter_content', 'text'),
        ('fastqc1_basic_statistics', 'text'),
        ('fastqc1_kmer_content', 'text'),
        ('fastqc1_max_sequence_length', 'integer'),
        ('fastqc1_min_sequence_length', 'integer'),
        ('fastqc1_overrepresented_sequences', 'text'),
        ('fastqc1_per_base_n_content', 'text'),
        ('fastqc1_per_base_sequence_content', 'text'),
        ('fastqc1_per_base_sequence_quality', 'text'),
        ('fastqc1_per_sequence_gc_content', 'text'),
        ('fastqc1_per_sequence_quality_scores', 'text'),
        ('fastqc1_sequence_duplication_levels', 'text'),
        ('fastqc1_sequence_length_distribution', 'text'),
        ('fastqc1_sequences_flagged_as_poor_quality', 'integer'),
        ('fastqc1_total_sequences', 'integer'),
        ('fastqc2_gc', 'float'),
        ('fastqc2_adapter_content', 'text'),
        ('fastqc2_basic_statistics', 'text'),
        ('fastqc2_kmer_content', 'text'),
        ('fastqc2_max_sequence_length', 'integer'),
        ('fastqc2_min_sequence_length', 'integer'),
        ('fastqc2_overrepresented_sequences', 'text'),
        ('fastqc2_per_base_n_content', 'text'),
        ('fastqc2_per_base_sequence_content', 'text'),
        ('fastqc2_per_base_sequence_quality', 'text'),
        ('fastqc2_per_sequence_gc_content', 'text'),
        ('fastqc2_per_sequence_quality_scores', 'text'),
        ('fastqc2_sequence_duplication_levels', 'text'),
        ('fastqc2_sequence_length_distribution', 'text'),
        ('fastqc2_sequences_flagged_as_poor_quality', 'integer'),
        ('fastqc2_total_sequences', 'integer'),
        ('samtools_raw_total_sequences', 'integer'),
        ('samtools_reads_mapped', 'integer'),
        ('samtools_reads_duplicated', 'integer'),
        ('samtools_bases_mapped_cigar', 'integer'),
        ('samtools_bases_trimmed', 'integer'),
        ('samtools_error_rate', 'float'),
        ('samtools_average_quality', 'float'),
        ('samtools_insert_size_average', 'float'),
        ('samtools_insert_size_standard_deviation', 'float'),
        ('samtools_inward_oriented_pairs', 'integer'),
        ('samtools_outward_oriented_pairs', 'integer'),
        ('samtools_pairs_with_other_orientation', 'integer'),
        ('het_snp_positions', 'integer'),
        ('het_snp_total_snps', 'integer'),
        ('het_snp_het_calls', 'integer'),
    ],

    'Read_counts': [
        ('seqrep_id', 'integer'),
        ('original_total', 'integer'),
        ('contamination', 'integer'),
        ('not_contamination', 'integer'),
        ('unmapped', 'integer'),
        ('total_after_remove_contam', 'integer'),
    ],

    'Reference': [
        ('reference_id', 'integer'),
        ('name', 'text'),
    ],

    'Sample': [
        ('sample_id', 'integer'),
        ('subject_id', 'text'),
        ('site_id', 'text'),
        ('sample_id_from_lab', 'text'),
        ('dataset_name', 'text'),
        ('ena_center_name', 'text'),
        ('ena_sample_accession', 'text'),
        ('ena_study_accession', 'text'),
    ],

    'Seqrep': [
        ('seqrep_id', 'integer'),
        ('isolate_id', 'integer'),
        ('sequence_replicate_number', 'integer'),
        ('original_reads_file_1_md5', 'text'),
        ('original_reads_file_2_md5', 'text'),
        ('remove_contam_reads_file_1_md5', 'text'),
        ('remove_contam_reads_file_2_md5', 'text'),
        ('withdrawn', 'integer'),
        ('import_status', 'integer'),
        ('submission_date', 'date'),
        ('submit_to_ena', 'integer'),
        ('instrument_model', 'text'),
        ('ena_run_accession', 'text'),
        ('ena_on_hold', 'integer'),
    ],

    'Version': [
        ('version', 'integer'),
    ]
}


primary_keys = {
    'Isolate': 'isolate_id',
    'Reference': 'reference_id',
    'Sample': 'sample_id',
    'Seqrep': 'seqrep_id',
}

