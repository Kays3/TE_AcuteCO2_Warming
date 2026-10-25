#!/usr/bin/env bash
# E2 (A1.6): 80-80-80 reciprocal best hits between the A. ocellaris library (16f27731) and the OA
# A. polyacanthus library (152a2a7e), rmblastn from dfam/tetools:2.00 (image a6179a355ba7), single-threaded.
set -euo pipefail
cd ~/jim_te_warm; W=~/jim_te_warm/e2; mkdir -p $W
AO=~/jim_te_warm/ref/aocel-families.fa
AP=~/KD/workspace/TE_OA_XEON_TRANSFER/mobilome-dfam-overlap-2026-09-01/apoly-families.fa
[ "$(sha256sum $AO | cut -c1-8)" = 16f27731 ] && [ "$(sha256sum $AP | cut -c1-8)" = 152a2a7e ] || { echo SHA MISMATCH; exit 2; }
cp $AO $W/aocel.fa; cp $AP $W/apoly.fa
D="docker run --rm --cpuset-cpus=0-20 -u $(id -u):$(id -g) -v $W:/w -w /w dfam/tetools:2.00"
F="6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
/usr/bin/time -v -o $W/e2.time bash -c "
$D makeblastdb -in aocel.fa -dbtype nucl -out aocel > /dev/null
$D makeblastdb -in apoly.fa -dbtype nucl -out apoly > /dev/null
$D rmblastn -query aocel.fa -db apoly -outfmt '$F' -reward 1 -penalty -2 -gapopen 5 -gapextend 2 -evalue 1e-10 -num_threads 1 -out ao_vs_ap.tsv
$D rmblastn -query apoly.fa -db aocel -outfmt '$F' -reward 1 -penalty -2 -gapopen 5 -gapextend 2 -evalue 1e-10 -num_threads 1 -out ap_vs_ao.tsv"
~/jim_te_warm/mamba/envs/jim_te_e1/bin/python ~/jim_te_warm/e2_rbh.py $W/ao_vs_ap.tsv $W/ap_vs_ao.tsv $W/rbh_808080.tsv
