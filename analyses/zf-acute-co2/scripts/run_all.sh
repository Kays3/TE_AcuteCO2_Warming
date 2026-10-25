#!/usr/bin/env bash
# Paths are placeholders: WORK (working folder), REF (reference folder), PIPE_ENV, TOOLS_ENV, R_ENV (software environments; see RUN_RECORD).
# All 26 libraries, 3 at a time in slots 1-3. After the first 3 finish, waits for the file go_rest (rate check, 35 CPU-h cap).
# Guards: stop launching if / has < 15 % free or MemAvailable < 20 %. A file STOP stops new launches.
cd $WORK
guard() {
  local used=$(df --output=pcent / | tail -1 | tr -dc 0-9)
  local mem_pct=$(awk "/MemTotal/{t=\$2}/MemAvailable/{a=\$2}END{printf \"%d\", 100*a/t}" /proc/meminfo)
  [ $((100-used)) -ge 15 ] && [ "$mem_pct" -ge 20 ] || { echo "GUARD STOP free=$((100-used))% mem=${mem_pct}% $(date -u +%FT%TZ)"; return 1; }
}
take_slot() { while :; do for s in 1 2 3; do mkdir work/slot$s.lock 2>/dev/null && { echo $s; return; }; done; sleep 30; done; }
n=0
for a in $(tail -n +2 runs.tsv); do
  [ -s out/$a/$a.multi.cntTable ] && continue
  if [ $n -eq 3 ]; then wait; echo "FIRST_DONE $(date -u +%FT%TZ)"; while [ ! -e go_rest ]; do [ -e STOP ] && exit 4; sleep 60; done; fi
  [ -e STOP ] && { echo "STOP file $(date -u +%FT%TZ)"; break; }
  s=$(take_slot); guard || exit 3
  ( /usr/bin/time -v -o logs/run_$a.time taskset -c 0-20 ./process_one.sh $a $s || echo "FAIL $a"; rmdir work/slot$s.lock ) &
  n=$((n+1)); sleep 5
done
wait
echo "ALL_DONE $(date -u +%FT%TZ)"
