#!/usr/bin/env bash
# 5.1 exact labeling null and E3 (observed + null): 924 labelings x 4 tissues per level, multi mode.
# Amendment 4 A4.2/A4.6: single-threaded R processes, at most MAXJ at once (default 14, leaving 6 of the
# 20 thread slots for one E1 featureCounts -T 6 job), all inside taskset -c 0-20, each under /usr/bin/time -v.
# Chunks of 84 labelings; E3 observed fits (index 0) first. Resumable: skips chunks whose output exists.
cd ~/jim_te_warm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R="/home/kaisar/KD/workspace/envs/te-ocean-acidification/bin/R --no-echo --no-restore"
MAXJ=${MAXJ:-14}
mkdir -p logs/analysis/perm
{
  for t in brain liver muscle digestive_tract; do echo "superfamily $t 0 0"; done
  for lv in family superfamily; do
    for t in brain liver muscle digestive_tract; do
      for s in $(seq 1 84 924); do echo "$lv $t $s $((s + 83))"; done
    done
  done
} | { while read lv t a b; do
  sfx=$([ "$lv" = superfamily ] && echo multi_sf || echo multi)
  f=$(printf 'analysis/perm_%s_%s_%04d_%04d.tsv' "$t" "$sfx" "$a" "$b")
  [ -s "$f" ] && continue
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do sleep 15; done
  avail=$(awk '/MemAvailable/{a=$2}/MemTotal/{t=$2}END{print int(100*a/t)}' /proc/meminfo)
  floor=$([ "$(date -u +%s)" -lt "$(date -u -d 2026-10-05T05:31:00Z +%s)" ] && echo 10 || echo 20)  # A4.5
  [ "$avail" -ge "$floor" ] || { echo "GUARD STOP mem=${avail}% $(date -u +%FT%TZ)" >> logs/analysis/status.log; exit 3; }
  tag=$(printf '%s_%s_%04d' "$t" "$sfx" "$a")
  /usr/bin/time -v -o logs/analysis/perm/$tag.time taskset -c 0-20 $R --file=deseq_fit.R --args "$t" multi "$a" "$b" analysis "$lv" > logs/analysis/perm/$tag.log 2>&1 &
  # the first superfamily fit of a tissue builds its prep file; let it finish before that tissue's chunks start
  [ "$a" = 0 ] && wait $!
  sleep 1
done
wait  # inside the pipeline subshell that owns the jobs (fix 2026-10-04 09:13Z: the old outer wait returned at once)
}
echo "PERM_DONE $(date -u +%FT%TZ)" >> logs/analysis/status.log
