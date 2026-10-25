"""Acute check: enumerate every within-cohort treatment assignment per tissue (PREREGISTRATION s4) and draw the
brain gene-control null (RUN_RECORD choice: 200 of the 599 non-observed brain assignments, default_rng(20261005),
without replacement). Order: cohorts 1..3, libraries sorted by run accession within a cohort, CO2 positions from
itertools.combinations, cohorts combined with itertools.product (cohort 1 varies slowest). Row 0 of each output is
always the observed assignment. Usage: acute_assign.py SAMPLES_TSV OUT_DIR"""
import sys, os, itertools
import numpy as np, pandas as pd

S, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
s = pd.read_csv(S, sep="\t")
for t, n_exp in (("brain", 600), ("gill", 210)):
    d = s[s.tissue == t].sort_values(["cohort", "run_accession"])
    runs = list(d.run_accession); obs = tuple(d.co2_treatment == "CO2")
    per = []
    for c, g in d.groupby("cohort", sort=True):
        idx = [runs.index(r) for r in g.run_accession]; k = int((g.co2_treatment == "CO2").sum())
        per.append([set(idx[i] for i in comb) for comb in itertools.combinations(range(len(idx)), k)])
    alls = [tuple(i in set().union(*p) for i in range(len(runs))) for p in itertools.product(*per)]
    assert len(alls) == n_exp and len(set(alls)) == n_exp and obs in alls, t
    alls = [obs] + [a for a in alls if a != obs]
    m = pd.DataFrame([[ "CO2" if x else "Control" for x in a] for a in alls], columns=runs)
    m.insert(0, "assignment", range(len(m))); m.to_csv(f"{O}/{t}_assignments.tsv", sep="\t", index=False)
    if t == "brain":
        draw = np.random.default_rng(20261005).choice(len(alls) - 1, size=200, replace=False) + 1
        pd.DataFrame({"assignment": np.sort(draw)}).to_csv(f"{O}/brain_gene_null_draws.tsv", sep="\t", index=False)
    print(t, len(alls), "observed CO2:", sum(obs))
