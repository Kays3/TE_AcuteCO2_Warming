#!/usr/bin/env bash
# Paths are placeholders: WORK (working folder), REF (reference folder), PIPE_ENV, TOOLS_ENV, R_ENV (software environments; see RUN_RECORD).
# Acute downstream fits in the registered order: TE permutations (brain 600, gill 210), brain gene control (200), then sensitivities. Cores 0-20.
cd $WORK; R=$R_ENV/bin/R
for j in "brain perm_te" "gill perm_te" "brain gene_null" "brain sens" "gill sens"; do
  set -- $j
  /usr/bin/time -v -o logs/an_$1_$2.time taskset -c 0-20 $R --no-echo --no-restore --file=acute_deseq.R --args analysis samples.tsv analysis/assign $1 $2 8 > logs/an_$1_$2.log 2>&1 || echo "FAIL $1 $2" >> logs/an_FAILED
  echo "$1 $2 done $(date -u +%FT%TZ)" >> logs/an_progress
done
echo "AN ALL_DONE $(date -u +%FT%TZ)" >> logs/an_progress
