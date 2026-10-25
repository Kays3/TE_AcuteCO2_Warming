"""Acute check s4: permutation p, Holm over two tissues, brain gene control and readings from the fit outputs.
p = share of the enumerated assignments (observed included) with a statistic at least as large. Gene control: observed brain
gene count at FDR 0.05 must be strictly greater than the 190th of the 200 ordered null counts (RUN_RECORD choice).
Also writes the sensitivity table against the primary null. Usage: acute_primary.py ANALYSIS_DIR"""
import sys, os
import pandas as pd
O = sys.argv[1]; rows = []
for t in ("brain", "gill"):
    p = pd.read_csv(f"{O}/{t}_perm_te.tsv", sep="\t"); obs = int(p.n_fdr05.iloc[0]); assert p.assignment.iloc[0] == 0
    te = pd.read_csv(f"{O}/{t}_te_observed.tsv", sep="\t"); sig = te[te.padj < 0.05]
    rows.append(dict(tissue=t, te_rows=len(te), n_fdr05=obs, up=int((sig.lfc_mle > 0).sum()), down=int((sig.lfc_mle < 0).sum()),
                     n_assign=len(p), perm_p=float((p.n_fdr05 >= obs).mean()), null_median=float(p.n_fdr05.median()),
                     null_95=float(p.n_fdr05.quantile(0.95)), null_max=int(p.n_fdr05.max())))
r = pd.DataFrame(rows).sort_values("perm_p").reset_index(drop=True)
r["holm_p"] = [min(1.0, max((2 - i) * r.perm_p.iloc[i] for i in range(k + 1))) for k in range(len(r))]
g = pd.read_csv(f"{O}/brain_gene_null.tsv", sep="\t"); gobs = int(g.n_fdr05.iloc[0]); gnull = sorted(g.n_fdr05.iloc[1:])
assert g.assignment.iloc[0] == 0 and len(gnull) == 200
thr = gnull[189]; gpass = gobs > thr
def reading(x):
    if x.tissue == "brain" and not gpass: return "UNINFORMATIVE"
    return "CO2_TE_SHIFT_DETECTED" if x.holm_p < 0.05 else "NO_DETECTABLE_TE_SHIFT"
r["reading"] = r.apply(reading, axis=1)
r.to_csv(f"{O}/primary_results.tsv", sep="\t", index=False)
pd.DataFrame([dict(gene_rows=len(pd.read_csv(f"{O}/brain_gene_observed.tsv", sep="\t")), observed_fdr05=gobs, null_190th_of_200=thr,
                   null_median=float(pd.Series(gnull).median()), null_max=max(gnull), passes=gpass)]).to_csv(f"{O}/brain_gene_control.tsv", sep="\t", index=False)
print(r.to_string()); print("gene control", gobs, ">", thr, gpass)
sens = []
for t in ("brain", "gill"):
    f = f"{O}/{t}_sensitivity.tsv"
    if not os.path.exists(f): continue
    p = pd.read_csv(f"{O}/{t}_perm_te.tsv", sep="\t")
    for x in pd.read_csv(f, sep="\t").itertuples():
        sens.append(dict(tissue=t, sensitivity=x.sensitivity, n_fdr05=x.n_fdr05, share_primary_null_ge=float((p.n_fdr05 >= x.n_fdr05).mean())))
if sens: pd.DataFrame(sens).to_csv(f"{O}/sensitivity_results.tsv", sep="\t", index=False); print(pd.DataFrame(sens).to_string())
