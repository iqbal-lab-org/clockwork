# Updates schema from version 2 to version 3
ALTER TABLE Seqrep MODIFY sequence_replicate_number BIGINT UNSIGNED;
UPDATE Version SET version=3 WHERE version=2;

