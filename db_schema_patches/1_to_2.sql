# Updates schema from version 1 to version 2
ALTER TABLE QC MODIFY fastqc1_max_sequence_length INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc1_min_sequence_length INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc1_sequences_flagged_as_poor_quality INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc1_total_sequences INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc2_max_sequence_length INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc2_min_sequence_length INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc2_sequences_flagged_as_poor_quality INT UNSIGNED;
ALTER TABLE QC MODIFY fastqc2_total_sequences INT UNSIGNED;
ALTER TABLE QC MODIFY samtools_raw_total_sequences INT UNSIGNED;
ALTER TABLE QC MODIFY samtools_reads_mapped INT UNSIGNED;
ALTER TABLE QC MODIFY samtools_reads_duplicated INT UNSIGNED;
ALTER TABLE QC MODIFY samtools_bases_mapped_cigar BIGINT UNSIGNED;
ALTER TABLE QC MODIFY samtools_bases_trimmed BIGINT UNSIGNED;
ALTER TABLE QC MODIFY samtools_inward_oriented_pairs INT UNSIGNED;
ALTER TABLE QC MODIFY samtools_outward_oriented_pairs INT UNSIGNED;
ALTER TABLE QC MODIFY samtools_pairs_with_other_orientation INT UNSIGNED;
ALTER TABLE QC MODIFY het_snp_positions INT UNSIGNED;
ALTER TABLE QC MODIFY het_snp_total_snps INT UNSIGNED;
ALTER TABLE QC MODIFY het_snp_het_calls INT UNSIGNED;

ALTER TABLE Read_counts MODIFY original_total INT UNSIGNED;
ALTER TABLE Read_counts MODIFY contamination INT UNSIGNED;
ALTER TABLE Read_counts MODIFY not_contamination INT UNSIGNED;
ALTER TABLE Read_counts MODIFY unmapped INT UNSIGNED;
ALTER TABLE Read_counts MODIFY total_after_remove_contam INT UNSIGNED;

UPDATE Version SET version=2 WHERE version=1;

