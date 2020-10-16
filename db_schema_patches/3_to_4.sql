# Updates schema from version 3 to version 4
ALTER TABLE QC ADD samtools_positions_with_depth_over_0 INT UNSIGNED, samtools_positions_with_depth_over_10 INT UNSIGNED, samtools_positions_with_depth_over_100 INT UNSIGNED;

UPDATE Version SET version=4 WHERE version=3;

