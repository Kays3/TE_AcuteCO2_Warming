#!/usr/bin/env bash
# All 48 libraries, MAXJ at a time (default 2), as their downloads complete. Resource guards (Michael 57e9ab):
# stop launching if / has < 300 GB free or MemAvailable < 20 % of RAM; FASTQ and .sra deleted per library.
cd ~/jim_te_warm
guard() {
  local free_gb=$(df -BG --output=avail / | tail -1 | tr -dc 0-9)
  local mem_pct=$(awk '/MemTotal/{t=$2}/MemAvailable/{a=$2}END{printf "%d", 100*a/t}' /proc/meminfo)
  [ "$free_gb" -ge 300 ] && [ "$mem_pct" -ge 20 ] || { echo "GUARD STOP free=${free_gb}G mem=${mem_pct}% $(date -u +%FT%TZ)"; return 1; }
}
pending() { cut -f1 design.tsv | tail -n +2 | while read a; do [ -s out/$a/$a.uniq.cntTable ] || echo $a; done; }
while [ -n "$(pending)" ]; do
  launched=0
  for a in $(pending); do
    [ -s logs/download/$a.done ] || continue
    [ -e work/$a.running ] && continue
    while [ "$(ls work/*.running 2>/dev/null | wc -l)" -ge "${MAXJ:-2}" ]; do sleep 30; done
    guard || exit 3
    touch work/$a.running
    ( /usr/bin/time -v -o logs/run_$a.time taskset -c 0-20 ./process_one.sh $a reverse >> logs/run_all.log 2>&1 || echo "FAIL $a" >> logs/run_all.log; rm -f work/$a.running ) &
    launched=1; sleep 5
  done
  [ $launched -eq 0 ] && sleep 60
done
wait
echo "ALL_DONE $(date -u +%FT%TZ)" >> logs/run_all.log
