#!/usr/bin/env bash
cd ~/jim_te_warm
for a in SRR24489122 SRR24489123; do /usr/bin/time -v -o logs/pilot_$a.time taskset -c 0-20 ./process_one.sh $a reverse; done
echo PILOT_DONE $(date -u +%FT%TZ)
