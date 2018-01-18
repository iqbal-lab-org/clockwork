#!/usr/bin/env bash

./run.make_ref_fastas.py
samtools faidx run.ref.fa
bwa index run.ref.fa

fastaq to_perfect_reads run.ref_to_make_reads.fa - 180 1 15 75 | fastaq deinterleave - run.reads_1.fq.gz run.reads_2.fq.gz

bwa mem run.ref.fa run.reads_1.fq.gz run.reads_2.fq.gz  > tmp.$$.sam
samtools sort -o run.bam tmp.$$.sam
rm tmp.$$.sam


