"""5.2.2 index, 5.3.2 unique-only sensitivity, 5.3.3 young families, E2 (A1.6, A1.13), per tissue.
Usage: summarize_secondary.py LOCAL_DIR(analysis/, rbh_808080.tsv) DESIGN.tsv OA_RESPONSE_PATTERNS.tsv.gz TABLES_DIR"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats

L, DESIGN, OA, T = sys.argv[1:5]
A = f"{L}/analysis"
TIS = ["brain", "digestive_tract", "liver", "muscle"]
d = pd.read_csv(DESIGN, sep="\t")

# ---- 5.2.2 global TE index (strict TE % of all counts), 31 minus 28 per tissue ----
ix = pd.read_csv(f"{A}/te_index_per_sample.tsv", sep="\t").merge(d[["run", "tissue", "temperature_c"]],
                                                                 left_on="sample_id", right_on="run")
rows = []
for t in TIS:
    a = ix[(ix.tissue == t) & (ix.temperature_c == 31)].te_index_pct.values
    b = ix[(ix.tissue == t) & (ix.temperature_c == 28)].te_index_pct.values
    diff = a.mean() - b.mean()
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    dfw = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    half = stats.t.ppf(0.975, dfw) * np.sqrt(va + vb)
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    g = diff / sp * (1 - 3 / (4 * (len(a) + len(b)) - 9))
    rows.append(dict(tissue=t, mean_31=a.mean(), mean_28=b.mean(), diff_pp=diff, ci_lo=diff - half, ci_hi=diff + half,
                     hedges_g=g, welch_p=stats.ttest_ind(a, b, equal_var=False).pvalue,
                     mwu_p=stats.mannwhitneyu(a, b, alternative="two-sided").pvalue, n31=len(a), n28=len(b)))
pd.DataFrame(rows).to_csv(f"{T}/index_5_2_2.tsv", sep="\t", index=False, float_format="%.4g")
ix.to_csv(f"{T}/te_index_per_sample.tsv", sep="\t", index=False)

# ---- 5.3.2 unique-only sensitivity and 5.3.3 young families, on the primary multi-mode hits ----
yk = pd.read_csv(f"{A}/young_kimura.tsv", sep="\t").set_index("family")
fam = lambda s: s.str.strip('"').str.split(":").str[0]
rows, hits_all = [], []
for t in TIS:
    M = pd.read_csv(f"{A}/observed_{t}_multi.tsv", sep="\t"); U = pd.read_csv(f"{A}/observed_{t}_uniq.tsv", sep="\t")
    M["family"] = fam(M.feature_id); U["family"] = fam(U.feature_id)
    hm = M[(M.set == "te_primary") & (M.padj_set < 0.05)]
    hu = U[(U.set == "te_primary") & (U.padj_set < 0.05)]
    young = M[M.set == "te_primary"].family.map(yk.median_kimura_cpg) < 5
    hm = hm.assign(young=hm.family.map(yk.median_kimura_cpg) < 5, kimura=hm.family.map(yk.median_kimura_cpg),
                   passes_uniq=hm.family.isin(hu.family), tissue=t)
    hits_all.append(hm)
    rows.append(dict(tissue=t, n_multi=len(hm), n_uniq=len(hu), multi_only=int((~hm.passes_uniq).sum()),
                     uniq_only=int((~hu.family.isin(hm.family)).sum()),
                     n_young_tested=int(young.sum()), n_young_hits=int(hm.young.sum())))
pd.DataFrame(rows).to_csv(f"{T}/uniq_young_5_3_2_3.tsv", sep="\t", index=False)
H = pd.concat(hits_all)
H[["tissue", "family", "class", "baseMean", "log2FoldChange", "ci95_low", "ci95_high", "padj_set", "passes_uniq",
   "kimura", "young"]].round(4).to_csv(f"{T}/primary_hits_5_1.tsv", sep="\t", index=False)

# ---- E2: matched families, warming vs CO2 (brain, liver) ----
rbh = pd.read_csv(f"{L}/rbh_808080.tsv", sep="\t")
rbh["ao"] = rbh.aocel_family.str.split("#").str[0]; rbh["ap"] = rbh.apoly_family.str.split("#").str[0]
oa = pd.read_csv(OA, sep="\t")
oa = oa[oa.analysis_set == "full"]
oa["fam"] = oa.feature_id.str.split(":").str[0]
rows = []
for t in ["brain", "liver"]:
    W = pd.read_csv(f"{A}/observed_{t}_multi.tsv", sep="\t"); W["fam"] = fam(W.feature_id)
    W = W[W.set.isin(["te_primary", "te_unknown"])]
    o = oa[oa.tissue == t].groupby("fam")[["developmental_vs_control_lfc", "intergenerational_vs_control_lfc"]].mean()
    m = rbh.merge(W[["fam", "log2FoldChange"]], left_on="ao", right_on="fam").merge(o, left_on="ap", right_index=True)
    for c in ["developmental_vs_control_lfc", "intergenerational_vs_control_lfc"]:
        x, y = m.log2FoldChange.values, m[c].values
        rho, p = stats.spearmanr(x, y); n = len(x)
        z, se = np.arctanh(rho), 1 / np.sqrt(n - 3)
        rows.append(dict(tissue=t, oa_contrast=c.replace("_lfc", ""), n_matched_tested=n,
                         share_of_warming_tested=n / len(W), rho=rho, rho_ci_lo=np.tanh(z - 1.96 * se),
                         rho_ci_hi=np.tanh(z + 1.96 * se), p_descriptive=p,
                         sign_concordance=float(np.mean(np.sign(x) == np.sign(y)))))
e2 = pd.DataFrame(rows).round(4)
e2.to_csv(f"{T}/e2_warming_vs_co2.tsv", sep="\t", index=False)
print(pd.read_csv(f"{T}/index_5_2_2.tsv", sep="\t").to_string()); print(pd.read_csv(f"{T}/uniq_young_5_3_2_3.tsv", sep="\t").to_string())
print(e2.to_string()); print("rbh pairs", len(rbh))
