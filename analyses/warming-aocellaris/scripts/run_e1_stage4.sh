#!/usr/bin/env bash
# E1 Stage 4 (A1.5): observed fits (index 0, with the descriptive controls) and the 924-labeling null per tissue,
# both measures; single-threaded R, at most MAXJ (default 18) at once inside taskset -c 0-20, /usr/bin/time -v each.
cd ~/jim_te_warm; mkdir -p logs/e1s4 e1/stage4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
R="/home/kaisar/KD/workspace/envs/te-ocean-acidification/bin/R --no-echo --no-restore"
MAXJ=${MAXJ:-18}
{ for t in brain liver muscle digestive_tract; do echo "$t 0 0"; for s in $(seq 1 84 924); do echo "$t $s $((s + 83))"; done; done; } | {
while read t a b; do
  [ -s "$(printf 'e1/stage4/e1_perm_%s_%04d_%04d.tsv' $t $a $b)" ] && continue
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do sleep 10; done
  tag=$(printf '%s_%04d' $t $a)
  /usr/bin/time -v -o logs/e1s4/$tag.time taskset -c 0-20 $R --file=e1_stage4_warming.R --args e1/stage3 e1/genecounts $t $a $b e1/stage4 > logs/e1s4/$tag.log 2>&1 &
  sleep 1
done
wait  # inside the pipeline subshell that owns the jobs (fix 2026-10-04 09:13Z: the old outer wait returned at once)
}
echo "E1S4_DONE $(date -u +%FT%TZ)" >> logs/e1s4/status.log
