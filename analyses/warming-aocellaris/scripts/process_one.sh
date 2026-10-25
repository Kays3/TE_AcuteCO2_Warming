#!/usr/bin/env bash
# One library through the frozen OA pipeline (pre-registration section 4, Amendments 1 and 3):
# fasterq-dump -> fastp 1.3.6 -> STAR 2.7.11b -> TEcount 2.2.4 multi and uniq -> delete FASTQ.
# Usage: process_one.sh <SRR> <stranded: no|forward|reverse>   (from ~/jim_te_warm, under taskset -c 0-20)
set -euo pipefail
acc=$1; STRANDED=$2
cd ~/jim_te_warm
F=~/KD/workspace/envs/te-ocean-acidification/bin; J=~/jim_te_warm/mamba/envs/jim_te_warm/bin
GENE_GTF=ref/GCF_022539595.1_ASM2253959v1_genomic.gtf; TE_GTF=ref/aocel.de-novo-te.gtf
check() { echo "$1  $2" | sha256sum -c --quiet || { echo "SHA MISMATCH $2"; exit 2; }; }
check c63a989abf5eff2a3ee27d99d88e86958bab96a63db10b0153c25d33529f450f $GENE_GTF
check d8b3625ebe12af7cb10a115c3bbbd573bdd1cb6dabdd1daa0ae050b8272badb5 $TE_GTF
d=work/$acc; o=out/$acc; L=logs/$acc; mkdir -p $d $o
T="/usr/bin/time -v -o"
[ -s $o/$acc.multi.cntTable ] && [ -s $o/$acc.uniq.cntTable ] && { echo "SKIP $acc"; exit 0; }
[ -s $d/${acc}_2.fastq ] || $T $L.dump.time $J/fasterq-dump --split-files --threads 8 --temp $d --outdir $d sra/$acc/$acc.sra > $L.dump.log 2>&1
$T $L.fastp.time $F/fastp -i $d/${acc}_1.fastq -I $d/${acc}_2.fastq -o $d/t_1.fastq.gz -O $d/t_2.fastq.gz \
  --detect_adapter_for_pe -w 8 -j $o/$acc.fastp.json -h /dev/null > $L.fastp.log 2>&1
rm -f $d/${acc}_1.fastq $d/${acc}_2.fastq
$T $L.star.time $F/STAR --runThreadN 16 --limitBAMsortRAM 30000000000 --genomeDir ref/star_149 \
  --readFilesIn $d/t_1.fastq.gz $d/t_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix $o/$acc. \
  --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD > $L.star.log 2>&1
rm -f $d/t_1.fastq.gz $d/t_2.fastq.gz
$T $L.index.time $F/samtools index -@ 8 $o/$acc.Aligned.sortedByCoord.out.bam
for m in multi uniq; do
  $T $L.tecount_$m.time $F/TEcount -b $o/$acc.Aligned.sortedByCoord.out.bam --GTF $GENE_GTF --TE $TE_GTF \
    --mode $m --stranded $STRANDED --project $o/$acc.$m --sortByPos > $L.tecount_$m.log 2>&1
done
rm -rf $d
# The .sra is regenerable (checksums in sra_checksums.tsv); drop it once both count tables exist.
[ -s $o/$acc.multi.cntTable ] && [ -s $o/$acc.uniq.cntTable ] && rm -f sra/$acc/$acc.sra
echo "OK $acc $(date -u +%FT%TZ)"
