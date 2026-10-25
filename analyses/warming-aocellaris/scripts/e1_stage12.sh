#!/usr/bin/env bash
# E1 (A1.5) Stages 1-2 on pSS with the committed mobilisation scripts, unchanged except for inputs.
# Sample ids are <tissue>_<run> with digestive_tract written digestivetract (stage 3 splits on "_").
# Slots (A4.6): Stage 1 and the region build are single-threaded; featureCounts -T 6, FCJ jobs at once (default 1).
set -uo pipefail
cd ~/jim_te_warm; W=~/jim_te_warm/e1; mkdir -p $W/{stage1,stage2/fc,stage2/fc_extra,bam,genecounts,logs/fc}
E=~/jim_te_warm/mamba/envs/jim_te_e1/bin
S=~/jim_te_warm/repos/TE_Ocean_Acidification/paper/results/te-mobilisation-transcription-20261002/scripts
RM=~/KD/workspace/TE_OA_XEON_TRANSFER/mobilome-comparative-2026-09-01/results/GCF_022539595.1/repeatmasker
REF=~/jim_te_warm/ref; G=$REF/GCF_022539595.1_ASM2253959v1_genomic.fna
# file identity (A1.11): exit on any mismatch
check() { [ "$(sha256sum "$1" | cut -c1-8)" = "$2" ] || { echo "SHA MISMATCH $1"; exit 2; }; }
check $REF/aocel-families.fa 16f27731; check $RM/genome.fa.out 9f51c6e7; check $RM/genome.fa.align 9d0105e5
check $REF/aocel.de-novo-te.gtf d8b3625e; check $G e139d17c
check $S/stage1_candidates.py 3acc6acd; check $S/stage2_regions.py ba99035a
[ -s $G.fai ] || ~/KD/workspace/envs/te-ocean-acidification/bin/samtools faidx $G
# inputs named as stage 3 expects
tail -n +2 design.tsv | while IFS=$'\t' read run name bs iso tissue temp clutch lab; do
  s="${tissue/_/}_$run"
  ln -sf ~/jim_te_warm/out/$run/$run.Aligned.sortedByCoord.out.bam $W/bam/$s.bam
  ln -sf ~/jim_te_warm/out/$run/$run.multi.cntTable $W/genecounts/$s.cntTable
  echo "$s"
done | sort > $W/stage2/samples.txt
T() { /usr/bin/time -v -o "$1" taskset -c 0-20 "${@:2}"; }
[ -s $W/stage1/copies_tiered.tsv.gz ] || T $W/logs/stage1.time $E/python $S/stage1_candidates.py $REF/aocel-families.fa $RM/genome.fa.out $RM/genome.fa.align $G $W/stage1 > $W/logs/stage1.log 2>&1 || { echo "STAGE1 FAIL"; exit 1; }
[ -s $W/stage2/regions_shuffle_extra.saf ] || T $W/logs/stage2_regions.time $E/python $S/stage2_regions.py --copies $W/stage1/copies_tiered.tsv.gz \
  --genes-gtf $REF/GCF_022539595.1_ASM2253959v1_genomic.gtf --te-gtf $REF/aocel.de-novo-te.gtf --fai $G.fai --genome $G \
  --out $W/stage2 --extra-shuffles 900 > $W/logs/stage2_regions.log 2>&1 || { echo "REGIONS FAIL"; exit 1; }
one() {  # one BAM: unique pass (all regions), 1/NH pass (loci), unique pass (900 extra shuffles), as s2one.sh / s2extra_one.sh
  s=$1
  for p in uniq multi extra; do
    case $p in uniq) saf=regions_unique.saf; o=fc/$s.uniq.txt; x="";; multi) saf=regions_loci.saf; o=fc/$s.multi.txt; x="-M --fraction";;
               extra) saf=regions_shuffle_extra.saf; o=fc_extra/$s.extra.txt; x="";; esac
    [ -s $W/stage2/$o.gz ] && continue
    /usr/bin/time -v -o $W/logs/fc/$s.$p.time taskset -c 0-20 $E/featureCounts -T 6 -F SAF -a $W/stage2/$saf -p --countReadPairs -s 0 -O --minOverlap 1 $x \
      -o $W/stage2/$o $W/bam/$s.bam > $W/logs/fc/$s.$p.log 2>&1 && gzip -f $W/stage2/$o
  done
}
export -f one; export W E
xargs -P "${FCJ:-1}" -I{} bash -c 'one {}' < $W/stage2/samples.txt
echo "E1_S12_DONE $(date -u +%FT%TZ) uniq=$(ls $W/stage2/fc/*.uniq.txt.gz | wc -l) multi=$(ls $W/stage2/fc/*.multi.txt.gz | wc -l) extra=$(ls $W/stage2/fc_extra/*.extra.txt.gz | wc -l)"
