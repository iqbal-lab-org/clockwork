#!/usr/bin/env bash

# This script was used to generate the files in this directory,
# which are used to test the nextflow QC pipeline.
# Expected to be run from here, ie: ./make_files.sh
set -vex


tmp_ref_dir=$PWD/Reference.tmp
pipeline_ref_dir=$PWD/Pipeline_refs
ref_dir_remove_contam=$PWD/Reference_remove_contam
ref_dir_qc=$PWD/Reference
dropbox=$PWD/Dropbox
import_tsv=$dropbox/import.tsv
pipeline_root=$PWD/Pipeline_root
db_pipeline_config=$PWD/../db.ini
mysql_cnf=db.cnf
spreadsheet_archive=$PWD/Spreadsheet_archive
reads_dir=$PWD/Reads
db_backup_dir=$PWD/db_backup
ref_metadata_tsv=$tmp_ref_dir/remove_contam.ref.metadata.tsv

sed 's/db_login/client/' $db_pipeline_config | egrep -v '^db' > $mysql_cnf
db_name=$(awk '$1~/^db/ {print $NF}' $db_pipeline_config)
mysql --defaults-file=$mysql_cnf -e "DROP DATABASE IF EXISTS $db_name"

rm -rf $pipeline_ref_dir $ref_dir_remove_contam $ref_dir_qc $tmp_ref_dir $dropbox $pipeline_root $db $spreadsheet_archive $reads_dir $db_backup_dir work .nextflow*
mkdir $pipeline_ref_dir $tmp_ref_dir $dropbox $pipeline_root $spreadsheet_archive $reads_dir $db_backup_dir

fastaq make_random_contigs --seed 42 --prefix ref. 1 1000 $tmp_ref_dir/ref.fa
fastaq make_random_contigs --seed 43 --prefix contam. 1 1000 $tmp_ref_dir/contam.fa
fastaq make_random_contigs --seed 44 --prefix other. 1 1000 $tmp_ref_dir/other.fa

fastaq to_perfect_reads --seed 45 $tmp_ref_dir/ref.fa - 200 1 10 75 | fastaq deinterleave - tmp.reads.1.ref.1.fq tmp.reads.1.ref.2.fq
fastaq to_perfect_reads --seed 46 $tmp_ref_dir/contam.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.1.contam.1.fq tmp.reads.1.contam.2.fq
fastaq to_perfect_reads --seed 47 $tmp_ref_dir/other.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.1.other.1.fq tmp.reads.1.other.2.fq

fastaq to_perfect_reads --seed 48 $tmp_ref_dir/ref.fa - 200 1 10 75 | fastaq deinterleave - tmp.reads.2.ref.1.fq tmp.reads.2.ref.2.fq
fastaq to_perfect_reads --seed 49 $tmp_ref_dir/contam.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.2.contam.1.fq tmp.reads.2.contam.2.fq
fastaq to_perfect_reads --seed 50 $tmp_ref_dir/other.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.2.other.1.fq tmp.reads.2.other.2.fq

fastaq to_perfect_reads --seed 51 $tmp_ref_dir/ref.fa - 200 1 10 75 | fastaq deinterleave - tmp.reads.3.ref.1.fq tmp.reads.3.ref.2.fq
fastaq to_perfect_reads --seed 52 $tmp_ref_dir/contam.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.3.contam.1.fq tmp.reads.3.contam.2.fq
fastaq to_perfect_reads --seed 53 $tmp_ref_dir/other.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.3.other.1.fq tmp.reads.3.other.2.fq

fastaq to_perfect_reads --seed 54 $tmp_ref_dir/ref.fa - 200 1 10 75 | fastaq deinterleave - tmp.reads.4.ref.1.fq tmp.reads.4.ref.2.fq
fastaq to_perfect_reads --seed 55 $tmp_ref_dir/contam.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.4.contam.1.fq tmp.reads.4.contam.2.fq
fastaq to_perfect_reads --seed 56 $tmp_ref_dir/other.fa - 200 1 1 75 | fastaq deinterleave - tmp.reads.4.other.1.fq tmp.reads.4.other.2.fq

for i in 1 2 3 4; do for j in 1 2; do cat tmp.reads.$i.*.$j.fq | gzip -c > $dropbox/reads.$i.$j.fq.gz; done; done
rm tmp.reads.*
cp $dropbox/reads.1.1.fq.gz $dropbox/reads.1.2.fq.gz $reads_dir/


echo "reference	0	ref.1
contaminants	1	contam.1" > $ref_metadata_tsv


for i in 1 2 3 4; do for j in 1 2; do md5sum $dropbox/reads.$i.$j.fq.gz > $dropbox/reads.$i.$j.fq.gz.md5; done; done

echo "subject_id,site_id,lab_id,isolate_number,sequence_replicate_number,submission_date,reads_file_1,reads_file_1_md5,reads_file_2,reads_file_2_md5,dataset_name,instrument_model,ena_center_name,submit_to_ena,ena_on_hold,ena_run_accession,ena_sample_accession
p1,s1,l1,42,43,20171225,$dropbox/reads.1.1.fq.gz,0,$dropbox/reads.1.2.fq.gz,0,g1,Illumina HiSeq 2000,Center 1,0,0,0,0
p2,s1,l1,44,45,20171225,$dropbox/reads.2.1.fq.gz,0,$dropbox/reads.2.2.fq.gz,0,g1,Illumina HiSeq 2000,Center 1,0,0,0,0
p3,s2,l1,46,47,20171225,$dropbox/reads.3.1.fq.gz,0,$dropbox/reads.3.2.fq.gz,0,g1,Illumina HiSeq 2000,Center 1,0,0,0,0
p4,s2,l1,48,49,20171225,$dropbox/reads.4.1.fq.gz,0,$dropbox/reads.4.2.fq.gz,0,g2,Illumina HiSeq 2000,Center 1,0,0,0,0" | sed 's/,/\t/g' > $import_tsv


cat $tmp_ref_dir/ref.fa $tmp_ref_dir/contam.fa > $tmp_ref_dir/remove_contam.ref.fa

clockwork make_empty_db $db_pipeline_config
clockwork reference_prepare --cortex_mem_height 17 --db_config_file $db_pipeline_config --pipeline_references_root $pipeline_ref_dir --contam_tsv $ref_metadata_tsv --name remove_contam_ref $tmp_ref_dir/remove_contam.ref.fa
clockwork reference_prepare --cortex_mem_height 17 --db_config_file $db_pipeline_config --pipeline_references_root $pipeline_ref_dir --name qc_ref $tmp_ref_dir/ref.fa
clockwork reference_prepare --cortex_mem_height 17 --outdir $ref_dir_qc $tmp_ref_dir/ref.fa

nextflow run --db_backup_dir $db_backup_dir --dropbox_dir $dropbox --pipeline_root $pipeline_root --db_config_file $db_pipeline_config --xlsx_archive_dir $spreadsheet_archive ../../../../../nextflow/import.nf

nextflow run --testing --ref_id 1 --references_root $pipeline_ref_dir --pipeline_root $pipeline_root --db_config_file $db_pipeline_config ../../../../../nextflow/remove_contam.nf

rm -rf work .nextflow* $dropbox $spreadsheet_archive $db_backup_dir $tmp_ref_dir

# Make a bad pair of fastq files, so pipeline will not work on them.
# Should get a pipeline status of -1 in the database
rm $pipeline_root/00/00/00/03/3/Reads/reads.remove_contam.47.1.fq.gz
rm $pipeline_root/00/00/00/03/3/Reads/reads.remove_contam.47.2.fq.gz
echo "whoops" | gzip -c > $pipeline_root/00/00/00/03/3/Reads/reads.remove_contam.47.1.fq.gz
echo "covfefe" | gzip -c > $pipeline_root/00/00/00/03/3/Reads/reads.remove_contam.47.2.fq.gz

mysqldump --defaults-file=$mysql_cnf $db_name > mysql.dump
