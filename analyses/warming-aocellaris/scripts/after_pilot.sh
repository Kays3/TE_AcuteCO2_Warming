#!/usr/bin/env bash
# Start the bulk run once the registered 2-library pilot has finished.
cd ~/jim_te_warm
until grep -q PILOT_DONE logs/pilot.log 2>/dev/null; do sleep 60; done
MAXJ=6 exec ./run_all.sh > logs/run_all.driver.log 2>&1 < /dev/null
