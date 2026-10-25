#!/usr/bin/env bash
# Strandedness test (pre-registration section 4): 2M-pair subset of one brain library, seqtk seed 20261003.
set -euo pipefail
cd ~/jim_te_warm
acc=SRR24489122; F=~/KD/workspace/envs/te-ocean-acidification/bin; J=mamba/envs/jim_te_warm/bin; d=work/$acc
mkdir -p $d strand
T="/usr/bin/time -v -o"
[ -s $d/${acc}_2.fastq ] || $T logs/$acc.dump.time taskset -c 0-20 $J/fasterq-dump --split-files --threads 8 --temp $d --outdir $d sra/$acc/$acc.sra > logs/$acc.dump.log 2>&1
for r in 1 2; do $T logs/strand.seqtk$r.time taskset -c 0-20 $J/seqtk sample -s20261003 $d/${acc}_$r.fastq 2000000 > strand/sub_$r.fq; done
$T logs/strand.star.time taskset -c 0-20 $F/STAR --runThreadN 16 --genomeDir ref/star_149 --readFilesIn strand/sub_1.fq strand/sub_2.fq \
  --quantMode GeneCounts --outSAMtype None --outFileNamePrefix strand/ --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 > logs/strand.star.log 2>&1
awk "NR>4{c2+=\$2;c3+=\$3;c4+=\$4} END{printf \"unstranded=%d col3=%d col4=%d ratio_col3=%.4f\n\",c2,c3,c4,c3/(c3+c4)}" strand/ReadsPerGene.out.tab | tee strand/ratio.txt
rm -f strand/sub_1.fq strand/sub_2.fq
