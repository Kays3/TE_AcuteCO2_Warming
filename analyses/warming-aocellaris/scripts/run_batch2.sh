#!/usr/bin/env bash
# Batch 2 (after the 5.1/E3 nulls and E1 Stages 1-2): single-threaded jobs, 9 at most at once, all within the
# 20 A4.6 slots, each under /usr/bin/time -v and taskset -c 0-20.
cd ~/jim_te_warm; mkdir -p logs/b2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R="/home/kaisar/KD/workspace/envs/te-ocean-acidification/bin/R --no-echo --no-restore"
P=~/jim_te_warm/mamba/envs/jim_te_e1/bin/python
RM=~/KD/workspace/TE_OA_XEON_TRANSFER/mobilome-comparative-2026-09-01/results/GCF_022539595.1/repeatmasker
T() { n=$1; shift; /usr/bin/time -v -o logs/b2/$n.time taskset -c 0-20 "$@" > logs/b2/$n.log 2>&1; echo "$n rc=$? $(date -u +%T)" >> logs/b2/status.log; }
T e1_stage3 $P e1_stage3_calls.py --stage2 e1/stage2 --gene-counts e1/genecounts --out e1/stage3 &
T e2 bash e2_match.sh &
T index $P index_build.py &
T young $P young_kimura.py $RM/genome.fa.align analysis/young_kimura.tsv &
for t in brain liver muscle digestive_tract; do T uniq_$t $R --file=deseq_fit.R --args $t uniq 0 0 analysis & done
T tissue_model $R --file=tissue_model.R --args analysis &
wait
echo "B2_DONE $(date -u +%FT%TZ)" >> logs/b2/status.log
