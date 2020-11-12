# Updates schema from version 3 to version 4
ALTER TABLE QC ADD samtools_positions_with_depth_of_0 INT UNSIGNED,
	samtools_positions_with_depth_atleast_2 INT UNSIGNED,
	samtools_positions_with_depth_atleast_5 INT UNSIGNED,
	samtools_positions_with_depth_atleast_10 INT UNSIGNED,
	samtools_positions_with_depth_atleast_20 INT UNSIGNED,
	samtools_positions_with_depth_atleast_100 INT UNSIGNED;

UPDATE Version SET version=4 WHERE version=3;

