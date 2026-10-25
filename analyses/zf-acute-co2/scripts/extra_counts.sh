#!/usr/bin/env bash
# Paths are placeholders: WORK (working folder), REF (reference folder), PIPE_ENV, TOOLS_ENV, R_ENV (software environments; see RUN_RECORD).
# Acute PREREGISTRATION.md (77270d9) section 3 and 5 count passes on the 26 BAMs, cores 0-20:
# (1) featureCounts 2.0.8 on the Tier A R1 cluster GTF (r1_shared573.gtf sha256 e4831c84), 1/NH fractional, -O, unstranded, paired fragments;
# (2) TEcount 2.2.4 --mode uniq (section 5 sensitivity), same GTFs as the multi pass. 6 BAMs at a time. Deletes nothing.
set -u
cd $WORK
F=$PIPE_ENV/bin; FC=$TOOLS_ENV/bin/featureCounts; R=$REF; GTF=$R/r1_shared573.gtf
echo "e4831c849554f8ee04a4da1c3ec9aa5dc610b368d7e34a4c7282493edcb5475a  $GTF" | sha256sum -c --quiet || { echo "GTF SHA MISMATCH"; exit 2; }
mkdir -p fc fc_tmp logs/extra
one() {
  r=$1; b=out/$r/$r.Aligned.sortedByCoord.out.bam
  [ -s fc/$r.frac.txt ] || /usr/bin/time -v -o logs/extra/$r.fc.time taskset -c 0-20 $FC -T 2 -s 0 -O --minOverlap 1 -p --countReadPairs -M --fraction \
      -t exon -g gene_id -a $GTF --tmpDir fc_tmp -o fc/$r.frac.txt $b > logs/extra/$r.fc.log 2>&1 || echo "FAIL $r fc" >> logs/extra/FAILED
  [ -s out/$r/$r.uniq.cntTable ] || /usr/bin/time -v -o logs/extra/$r.teu.time taskset -c 0-20 $F/TEcount --sortByPos --format BAM --mode uniq --stranded no \
      --BAM $b --GTF $R/genes.gtf --TE $R/GRCz11_rmsk_TE.gtf --project $r.uniq --outdir out/$r > logs/extra/$r.teu.log 2>&1 || echo "FAIL $r teu" >> logs/extra/FAILED
  echo "$r done $(date -u +%FT%TZ)" >> logs/extra/progress
}
export -f one; export F FC R GTF
tail -n +2 runs.tsv | xargs -P 6 -n 1 bash -c 'one "$0"'
echo "EXTRA ALL_DONE $(date -u +%FT%TZ)" >> logs/extra/progress
