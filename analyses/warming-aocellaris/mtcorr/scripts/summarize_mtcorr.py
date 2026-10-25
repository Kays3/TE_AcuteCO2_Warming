"""Amendment 5 (POST-HOC) summary: shares (A5.2), M1 readings with the relaxed per-tissue control (A5.5, A5.6), the relaxed
control on the registered fits (A5.5 a, descriptive), M2 and E3 under M1 against their own 924-labeling nulls, unique-mode
M1 observed counts, the family overlap with the registered result, and the global TE index with rRNA removed (A5.7).
p = (labelings 1..924 with count >= observed) / 924. Control percentile = value at rank 832 of the 924 ordered counts.
Usage: summarize_mtcorr.py WORK_DIR OUT_DIR   (WORK_DIR holds design.tsv, out/, analysis/, analysis_mtcorr/, rrna_genes.txt)"""
import sys, os, glob, math, re
import numpy as np, pandas as pd
from scipy import stats
W, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
T = ["brain", "digestive_tract", "liver", "muscle"]; R832 = math.ceil(0.9 * 924); R878 = math.ceil(0.95 * 924)
PRIM = {"DNA", "LINE", "SINE", "LTR", "RC", "PLE", "Retroposon"}
des = pd.read_csv(f"{W}/design.tsv", sep="\t")

def null(d, tis, tag):
    fs = [f for f in glob.glob(f"{W}/{d}/perm_{tis}_{tag}_*.tsv") if re.search(rf"perm_{tis}_{tag}_\d+_\d+\.tsv$", f) and not f.endswith("_0000_0000.tsv")]
    x = pd.concat(pd.read_csv(f, sep="\t") for f in fs).sort_values("labeling")
    assert list(x.labeling) == list(range(1, 925)), (d, tis, tag, len(x))
    o = pd.read_csv(f"{W}/{d}/perm_{tis}_{tag}_0000_0000.tsv", sep="\t").iloc[0]
    return o, x

def pctl(v, r): return float(np.sort(v)[r - 1])

shares = pd.concat(pd.read_csv(f"{W}/analysis_mtcorr/rrna_share_{t}_multi_m1.tsv", sep="\t").assign(tissue=t) for t in T)
shares = shares.merge(des[["run", "temperature_c"]], on="run"); shares.to_csv(f"{O}/a52_shares.tsv", sep="\t", index=False)
shares.groupby(["tissue", "temperature_c"])[["mt16S", "mt_rRNA", "all_rRNA_genes", "mtrRNA_centred"]].mean().reset_index().to_csv(f"{O}/a52_shares_by_group.tsv", sep="\t", index=False)

rows = []
for t in T:
    r = dict(tissue=t)
    o, x = null("analysis", t, "multi")  # registered fits, stored nulls
    r.update(reg_te=int(o.n_te_primary_fdr05), reg_p=float((x.n_te_primary_fdr05 >= o.n_te_primary_fdr05).mean()),
             reg_gene=int(o.n_gene_fdr05), reg_gene_p95=pctl(x.n_gene_fdr05, R878), reg_gene_p90=pctl(x.n_gene_fdr05, R832))
    r["reg_relaxed_control_pass"] = r["reg_gene"] > r["reg_gene_p90"]
    for lab, tag in (("m1", "multi_m1"), ("m2", "multi_m2"), ("e3m1", "multi_sf_m1")):
        try: o, x = null("analysis_mtcorr", t, tag)
        except (AssertionError, ValueError): continue
        r[f"{lab}_te"] = int(o.n_te_primary_fdr05); r[f"{lab}_p"] = float((x.n_te_primary_fdr05 >= o.n_te_primary_fdr05).mean())
        r[f"{lab}_gene"] = int(o.n_gene_fdr05); r[f"{lab}_gene_p90"] = pctl(x.n_gene_fdr05, R832); r[f"{lab}_gene_p95"] = pctl(x.n_gene_fdr05, R878)
    u = pd.read_csv(f"{W}/analysis_mtcorr/perm_{t}_uniq_m1_0000_0000.tsv", sep="\t").iloc[0]; r["m1_uniq_te"] = int(u.n_te_primary_fdr05)
    if "m1_te" in r:
        ctrl = r["m1_gene"] > r["m1_gene_p90"]; r["m1_relaxed_control_pass"] = ctrl
        det = r["m1_te"] >= 1 and r["m1_p"] < 0.0125
        r["posthoc_reading"] = "POSTHOC_TANK_TEMPERATURE_SHIFT_DETECTED" if det else ("POSTHOC_NO_DETECTABLE_SHIFT" if ctrl else "POSTHOC_UNINFORMATIVE")
    reg = pd.read_csv(f"{W}/analysis/observed_{t}_multi.tsv", sep="\t"); m1 = pd.read_csv(f"{W}/analysis_mtcorr/observed_{t}_multi_m1.tsv", sep="\t")
    a = set(reg[(reg.set == "te_primary") & (reg.padj_set < .05)].feature_id); b = m1[(m1.set == "te_primary") & (m1.padj_set < .05)]
    r.update(m1_up=int((b.log2FoldChange > 0).sum()), m1_down=int((b.log2FoldChange < 0).sum()), shared_with_registered=len(a & set(b.feature_id)),
             only_registered=len(a - set(b.feature_id)), only_m1=len(set(b.feature_id) - a))
    j = reg[reg.set == "te_primary"].merge(m1[m1.set == "te_primary"], on="feature_id", suffixes=("_reg", "_m1"))
    r["lfc_spearman_reg_vs_m1"] = float(stats.spearmanr(j.log2FoldChange_reg, j.log2FoldChange_m1, nan_policy="omit")[0])
    b.assign(tissue=t).to_csv(f"{O}/m1_families_{t}.tsv", sep="\t", index=False)
    rows.append(r)
res = pd.DataFrame(rows); res.to_csv(f"{O}/a5_results.tsv", sep="\t", index=False)

# A5.7 global index: strict-TE (A1.9 classes) counts as % of all counts, with rRNA gene rows and RM rRNA rows removed
rg = set(open(f"{W}/rrna_genes.txt").read().split()); ix = []
for run, tis, temp in des[["run", "tissue", "temperature_c"]].itertuples(index=False):
    c = pd.read_csv(f"{W}/out/{run}/{run}.multi.cntTable", sep="\t", index_col=0).iloc[:, 0]; c.index = c.index.str.strip('"')
    p = c.index.str.split(":"); te = p.str.len() == 3; cl = p.str[-1].where(te)
    keep = ~(c.index.isin(rg) | (te & (cl == "rRNA"))); st = te & cl.isin(PRIM)
    ix.append(dict(run=run, tissue=tis, temperature_c=temp, index_pct=100 * c[st & keep].sum() / c[keep].sum()))
ix = pd.DataFrame(ix); ix.to_csv(f"{O}/a57_index_per_library.tsv", sep="\t", index=False); out = []
for t, d in ix.groupby("tissue"):
    a, b = d[d.temperature_c == 31].index_pct.values, d[d.temperature_c == 28].index_pct.values
    diff = a.mean() - b.mean(); se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)); w = stats.ttest_ind(a, b, equal_var=False)
    dfw = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / ((a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1) + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)); g = diff / sp * (1 - 3 / (4 * (len(a) + len(b)) - 9))
    out.append(dict(tissue=t, mean_31=a.mean(), mean_28=b.mean(), diff_pp=diff, ci_lo=diff - stats.t.ppf(.975, dfw) * se, ci_hi=diff + stats.t.ppf(.975, dfw) * se,
                    hedges_g=g, welch_p=w.pvalue, mwu_p=stats.mannwhitneyu(a, b, alternative="two-sided").pvalue))
pd.DataFrame(out).to_csv(f"{O}/a57_index.tsv", sep="\t", index=False)
print(res.T.to_string()); print(pd.DataFrame(out).round(3).to_string())
