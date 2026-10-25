"""Section 5.1 outcomes, 5.2.1 sign test, 5.2.3 and 5.3.1 against the exact labeling null, and E3, per tissue.
Input: the analysis/ folder from pSS (observed_*.tsv, perm_*.tsv). Output: tables/ (small TSVs).
p = (labelings in 1..924 with count >= observed) / 924; the observed labeling is one of the 924 (obs_idx),
and its index-0 refit must give the same count (checked). Outcome threshold 0.05 / 4 = 0.0125.
Usage: summarize_nulls.py ANALYSIS_DIR TABLES_DIR"""
import sys, glob, os
import numpy as np, pandas as pd
from scipy.stats import binomtest

A, T = sys.argv[1], sys.argv[2]
os.makedirs(T, exist_ok=True)
TISSUES = ["brain", "digestive_tract", "liver", "muscle"]
rows, signs, nulls = [], [], []
for level, sfx in (("family", "multi"), ("superfamily", "multi_sf")):
    for t in TISSUES:
        fs = sorted(glob.glob(f"{A}/perm_{t}_{sfx}_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].tsv"))
        P = pd.concat([pd.read_csv(f, sep="\t") for f in fs])
        obs0 = P[P.labeling == 0]
        perm = P[P.labeling > 0].drop_duplicates("labeling").sort_values("labeling")
        assert len(perm) == 924 and set(perm.labeling) == set(range(1, 925)), (t, sfx, len(perm))
        assert len(obs0) == 1
        o = obs0.iloc[0]
        oi = perm[perm.observed]
        assert len(oi) == 1 and int(oi.n_te_primary_fdr05.iloc[0]) == int(o.n_te_primary_fdr05), (t, sfx)
        n_obs = int(o.n_te_primary_fdr05)
        p = (perm.n_te_primary_fdr05 >= n_obs).sum() / 924
        outcome = "TANK_TEMPERATURE_SHIFT_DETECTED" if (n_obs >= 1 and p < 0.0125) else "NO_DETECTABLE_TANK_SHIFT"
        g_obs = int(o.n_gene_fdr05); g95 = np.percentile(perm.n_gene_fdr05, 95)
        r_obs = float(o.te_gene_matched_ratio)
        rows.append(dict(level=level, tissue=t, n_primary_fdr05=n_obs, perm_p=round(p, 4),
                         null_median=float(perm.n_te_primary_fdr05.median()),
                         null_p95=float(np.percentile(perm.n_te_primary_fdr05, 95)),
                         outcome=outcome, n_unknown_fdr05=int(o.n_te_unknown_fdr05),
                         n_gene_fdr05=g_obs, gene_null_p95=g95,
                         gene_perm_p=round((perm.n_gene_fdr05 >= g_obs).sum() / 924, 4),
                         gene_exceeds_p95=bool(g_obs > g95),
                         te_gene_ratio=round(r_obs, 3),
                         te_gene_ratio_perm_p_upper=round((perm.te_gene_matched_ratio >= r_obs).sum() / 924, 4),
                         te_gene_ratio_null_median=round(float(perm.te_gene_matched_ratio.median()), 3)))
        nulls.append(perm.assign(level=level, tissue=t))
        R = pd.read_csv(f"{A}/observed_{t}_{sfx}.tsv", sep="\t")
        hit = R[(R.set == "te_primary") & (R.padj_set < 0.05)]
        up, dn = int((hit.log2FoldChange > 0).sum()), int((hit.log2FoldChange < 0).sum())
        bp = binomtest(up, up + dn, 0.5).pvalue if up + dn else float("nan")
        signs.append(dict(level=level, tissue=t, n_pass=up + dn, up=up, down=dn, sign_test_p=round(bp, 4),
                          underpowered=(up + dn) < 10))
pd.DataFrame(rows).to_csv(f"{T}/outcomes_5_1_E3.tsv", sep="\t", index=False)
pd.DataFrame(signs).to_csv(f"{T}/sign_test_5_2_1.tsv", sep="\t", index=False)
pd.concat(nulls).to_csv(f"{T}/null_distributions.tsv.gz", sep="\t", index=False)
print(pd.DataFrame(rows).to_string()); print(pd.DataFrame(signs).to_string())
