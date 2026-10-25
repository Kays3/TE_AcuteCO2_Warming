#!/usr/bin/env bash
cd ~/jim_te_warm
P=4 exec /usr/bin/time -v -o logs/download_total.time taskset -c 0-20 ./download_sra.sh design.tsv PRJNA970523_runinfo.csv > logs/download_total.log 2>&1 < /dev/null
