#!/usr/bin/env bash
# Amendment 5 runner. STAGE=observed: the 4 M1 observed fits (family level, multi), plus the M1 unique-mode observed fits.
# STAGE=null: one 924-labeling null per tissue for each JOB in "$JOBS" (default "family:m1 family:m2 superfamily:m1"),
# in that order, chunks of 84, single-threaded R, at most MAXJ at once, inside taskset -c 0-20, each under /usr/bin/time -v.
# Resumable: skips chunks whose output exists. A file STOP stops new launches. Guards as A4.5 (MemAvailable >= 20 %).
cd "$WORK"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R="$PIPE_ENV/bin/R --no-echo --no-restore"; MAXJ=${MAXJ:-18}; O=analysis_mtcorr
mkdir -p logs/mtcorr $O
run() { # tissue mode first last level model
  tag=$(printf '%s_%s_%s_%s_%04d' "$1" "$2" "$5" "$6" "$3")
  /usr/bin/time -v -o logs/mtcorr/$tag.time taskset -c 0-20 $R --file=deseq_fit_mtcorr.R --args "$1" "$2" "$3" "$4" $O "$5" "$6" > logs/mtcorr/$tag.log 2>&1
}
if [ "$STAGE" = observed ]; then
  for t in brain liver muscle digestive_tract; do run $t multi 0 0 family m1 & done; wait
  for t in brain liver muscle digestive_tract; do run $t uniq 0 0 family m1 & done; wait
  echo "OBSERVED_DONE $(date -u +%FT%TZ)" >> logs/mtcorr/status.log; exit 0
fi
for job in ${JOBS:-family:m1 family:m2 superfamily:m1}; do
  lv=${job%%:*}; md=${job##*:}; sfx=$([ "$lv" = superfamily ] && echo multi_sf || echo multi)
  for t in brain liver muscle digestive_tract; do run $t multi 0 0 $lv $md; done   # builds each prep file first
  for t in brain liver muscle digestive_tract; do for s in $(seq 1 84 924); do echo "$t $s $((s + 83))"; done; done | {
    while read t a b; do
      [ -e STOP ] && break
      [ -s "$(printf '%s/perm_%s_%s_%s_%04d_%04d.tsv' $O "$t" "$sfx" "$md" "$a" "$b")" ] && continue
      while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do sleep 15; done
      avail=$(awk '/MemAvailable/{a=$2}/MemTotal/{t=$2}END{print int(100*a/t)}' /proc/meminfo)
      [ "$avail" -ge 20 ] || { echo "GUARD STOP mem=${avail}% $(date -u +%FT%TZ)" >> logs/mtcorr/status.log; break; }
      run "$t" multi "$a" "$b" "$lv" "$md" & sleep 1
    done; wait; }
  echo "NULL_DONE $job $(date -u +%FT%TZ)" >> logs/mtcorr/status.log
done
