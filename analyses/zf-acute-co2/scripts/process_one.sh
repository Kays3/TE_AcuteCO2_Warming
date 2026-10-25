#!/usr/bin/env bash
# Paths are placeholders: WORK (working folder), REF (reference folder), PIPE_ENV, TOOLS_ENV, R_ENV (software environments; see RUN_RECORD).
# One PRJNA1102381 library through the zebrafish archive pipeline (acute PREREGISTRATION.md section 3 at 77270d9;
# command lines from the archive scripts as pasted by Kevin, kev-zfcmd-1905):
# fasterq-dump -> fastp 1.3.6 (defaults, --detect_adapter_for_pe) -> STAR 2.7.11b (multimap 100) -> TEcount 2.2.4 multi, unstranded.
# BAMs are kept for the featureCounts cluster and genic/intergenic passes. Genome loaded per process (NoSharedMemory)
# instead of the archive LoadAndKeep; loading mode does not change alignments.
# Temporary FASTQ go to a fixed slot folder (work/slot<K>) that the next library in that slot overwrites; nothing is deleted.
# Usage: process_one.sh <SRR> <slot>   (from $WORK, under taskset -c 0-20)
set -euo pipefail
acc=$1; slot=$2; cd $WORK; N=${THREADS:-7}
F=$PIPE_ENV/bin; J=$TOOLS_ENV/bin; R=$REF
check() { echo "$1  $2" | sha256sum -c --quiet || { echo "SHA MISMATCH $2"; exit 2; }; }
check f2034791a53bd277a64fa6202b8403e7bd23fa44994d5cc1d7b125218dcfc937 $R/genes.gtf
check 4bf43918324146bb1a5391225f0e804b0aa3c0dcd7d79fd764f3a8ec20a44cbf $R/GRCz11_rmsk_TE.gtf
d=work/slot$slot; o=out/$acc; L=logs/$acc; mkdir -p $d $o
T="/usr/bin/time -v -o"
[ -s $o/$acc.multi.cntTable ] && { echo "SKIP $acc"; exit 0; }
$T $L.dump.time $J/fasterq-dump --force --split-files --threads $N --temp $d --outdir $d -o r.fastq sra/$acc/$acc.sra > $L.dump.log 2>&1
$T $L.fastp.time $F/fastp --in1 $d/r_1.fastq --in2 $d/r_2.fastq --out1 $d/t_1.fastq.gz --out2 $d/t_2.fastq.gz \
  --detect_adapter_for_pe --thread $N --json $o/$acc.fastp.json --html $o/$acc.fastp.html > $L.fastp.log 2>&1
$T $L.star.time $F/STAR --limitBAMsortRAM 20000000000 --runThreadN $N --genomeDir $R/star_index \
  --readFilesIn $d/t_1.fastq.gz $d/t_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix $o/$acc. \
  --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD > $L.star.log 2>&1
$T $L.index.time $F/samtools index -@ $N $o/$acc.Aligned.sortedByCoord.out.bam
$T $L.tecount_multi.time $F/TEcount --sortByPos --format BAM --mode multi --stranded no \
  --BAM $o/$acc.Aligned.sortedByCoord.out.bam --GTF $R/genes.gtf --TE $R/GRCz11_rmsk_TE.gtf \
  --project $acc.multi --outdir $o > $L.tecount_multi.log 2>&1
echo "OK $acc $(date -u +%FT%TZ)"
