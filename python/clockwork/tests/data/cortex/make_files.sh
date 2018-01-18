#!/usr/bin/env bash
set -e

refdir=Reference
rm -fr mutated.fa reads.fq $refdir


./make_ref_fastas.py
fastaq to_perfect_reads --seed 42 mutated.fa reads.fq 200 1 30 75
samtools faidx ref.fa
$CLOCKWORK_STAMPY_SCRIPT -G ref.stampy ref.fa
$CLOCKWORK_STAMPY_SCRIPT -g ref.stampy -H ref.stampy

echo $PWD/ref.fa > ref.cortex_se_list.fofn

$CLOCKWORK_CORTEX_DIR/bin/cortex_var_31_c1 \
    --kmer_size 31 \
    --mem_height 17 \
    --mem_width 100 \
    --se_list ref.cortex_se_list.fofn \
    --max_read_len 10000 \
    --dump_binary ref.k31.ctx \
    --sample_id REF

rm ref.cortex_se_list.fofn
mkdir $refdir
mv ref.* $refdir
