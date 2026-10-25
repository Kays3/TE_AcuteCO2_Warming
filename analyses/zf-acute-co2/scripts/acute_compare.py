"""Acute check s5, exploratory comparison (RUN_RECORD choice 7): Spearman rho between zebrafish CO2 apeglm log2FC and damselfish
OA cluster log2FC, per side; 95 % percentile bootstrap over clusters (5,000, seed 20261005); same-sign share; the same for the MLE log2FC. No p, no reading.
Usage: acute_compare.py ANALYSIS_DIR"""
import sys
import numpy as np, pandas as pd
from scipy import stats
O = sys.argv[1]; d = pd.read_csv(f"{O}/acute_damsel_clusters.tsv", sep="\t"); rows = []
for side, tis in (("gill_developmental", "gill"), ("gill_intergenerational", "gill"), ("brain_intergenerational", "brain")):
    z = pd.read_csv(f"{O}/{tis}_clusters_zf.tsv", sep="\t"); j = d[d.side == side].merge(z, on="cluster")
    j.to_csv(f"{O}/compare_{side}.tsv", sep="\t", index=False)
    x, y = j.d_lfc.values, j.lfc_apeglm.values; rho = stats.spearmanr(x, y)[0]
    rng = np.random.default_rng(20261005); bs = []
    for _ in range(5000):
        i = rng.integers(0, len(j), len(j)); bs.append(stats.spearmanr(x[i], y[i])[0])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5]); ym = j.lfc_mle.values
    rng = np.random.default_rng(20261005); bm = []  # MLE sensitivity: apeglm shrank every effect to ~1e-5 in these data
    for _ in range(5000):
        i = rng.integers(0, len(j), len(j)); bm.append(stats.spearmanr(x[i], ym[i])[0])
    mlo, mhi = np.nanpercentile(bm, [2.5, 97.5])
    rows.append(dict(side=side, n_matched=len(j), rho=rho, ci_lo=lo, ci_hi=hi, same_sign=float(np.mean(np.sign(x) == np.sign(y))),
                     rho_mle=stats.spearmanr(x, ym)[0], mle_ci_lo=mlo, mle_ci_hi=mhi, mle_same_sign=float(np.mean(np.sign(x) == np.sign(ym)))))
r = pd.DataFrame(rows); r.to_csv(f"{O}/compare_results.tsv", sep="\t", index=False); print(r.round(3).to_string())
