#!/usr/bin/env bash

# This script was used to generate the files in this directory,
# which are used to test the nextflow remove_contam pipeline.
# Expected to be run from here, ie: ./make_files.sh
set -vex


tmp_ref_dir=$PWD/Reference.tmp
pipeline_ref_dir=$PWD/Pipeline_refs
ref_dir=$PWD/Reference
dropbox=$PWD/Dropbox
import_tsv=$dropbox/import.tsv
pipeline_root=$PWD/Pipeline_root
db_pipeline_config=$PWD/../db.ini
mysql_cnf=db.cnf
spreadsheet_archive=$PWD/Spreadsheet_archive
reads_dir=$PWD/Reads
db_backup_dir=$PWD/db_backup

sed 's/db_login/client/' $db_pipeline_config | egrep -v '^db' > $mysql_cnf
db_name=$(awk '$1~/^db/ {print $NF}' $db_pipeline_config)
mysql --defaults-file=$mysql_cnf -e "DROP DATABASE IF EXISTS $db_name"


rm -rf $pipeline_ref_dir $ref_dir $tmp_ref_dir $dropbox $pipeline_root $spreadsheet_archive $reads_dir $db_backup_dir work .nextflow*
mkdir $pipeline_ref_dir $tmp_ref_dir $dropbox $pipeline_root $spreadsheet_archive $reads_dir $db_backup_dir

fastaq make_random_contigs --seed 42 --prefix ref. 1 1000 $tmp_ref_dir/ref.fa


for i in 1 2 3 4
do
    fastaq to_perfect_reads --seed $i $tmp_ref_dir/ref.fa - 200 1 $i 75 | fastaq deinterleave - $dropbox/reads.$i.1.fq.gz $dropbox/reads.$i.2.fq.gz
     md5sum $dropbox/reads.$i.1.fq.gz > $dropbox/reads.$i.1.fq.gz.md5
     md5sum $dropbox/reads.$i.2.fq.gz > $dropbox/reads.$i.2.fq.gz.md5
done


echo "subject_id,site_id,lab_id,isolate_number,sequence_replicate_number,submission_date,reads_file_1,reads_file_1_md5,reads_file_2,reads_file_2_md5,dataset_name,instrument_model,ena_center_name,submit_to_ena,ena_on_hold,ena_run_accession,ena_sample_accession
p1,s1,l1,1,1,20171225,$dropbox/reads.1.1.fq.gz,0,$dropbox/reads.1.2.fq.gz,0,g1,Illumina HiSeq 2000,Center 1,0,0,0,0
p2,s1,l1,2,1,20171225,$dropbox/reads.2.1.fq.gz,0,$dropbox/reads.2.2.fq.gz,0,g1,Illumina HiSeq 2000,Center 1,0,0,0,0
p3,s1,l1,3,1,20171225,$dropbox/reads.3.1.fq.gz,0,$dropbox/reads.3.2.fq.gz,0,g1,Illumina HiSeq 2000,Center 1,0,0,0,0
p4,s1,l1,4,1,20171225,$dropbox/reads.4.1.fq.gz,0,$dropbox/reads.4.2.fq.gz,0,g2,Illumina HiSeq 2000,Center 1,0,0,0,0" | sed 's/,/\t/g' > $import_tsv


clockwork make_empty_db $db_pipeline_config
nextflow run ../../../../../nextflow/import.nf --db_backup_dir $db_backup_dir --dropbox_dir $dropbox --pipeline_root $pipeline_root --db_config_file $db_pipeline_config --xlsx_archive_dir $spreadsheet_archive

rm -rf work .nextflow* $dropbox $spreadsheet_archive $db_backup_dir


# Make a bad pair of fastq files, so pipeline will not work on them.
# Should get a pipeline status of -1 in the database
echo "covfefe" | gzip -c > $pipeline_root/00/00/00/03/3/Reads/reads.original.1.1.fq.gz

mysqldump --defaults-file=$mysql_cnf $db_name > mysql.dump

rm -r $tmp_ref_dir

