"""Report figures for the warming study. Usage: make_figures.py TABLES_DIR LOCAL_DIR OUT_PNG_DIR"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T, L, O = sys.argv[1:4]
os.makedirs(O, exist_ok=True)
A = f"{L}/analysis"
TIS = ["brain", "digestive_tract", "liver", "muscle"]
LAB = {"brain": "Brain", "digestive_tract": "Digestive tract", "liver": "Liver", "muscle": "Muscle"}
C28, C31 = "#3b6ea8", "#c0504d"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
out = pd.read_csv(f"{T}/outcomes_5_1_E3.tsv", sep="\t")
nul = pd.read_csv(f"{T}/null_distributions.tsv.gz", sep="\t")

# Fig 1: observed count of classified TE families at FDR < 0.05 against the 924-labeling null (family, superfamily)
fig, ax = plt.subplots(2, 4, figsize=(10, 4.6), sharey="row")
for i, lv in enumerate(["family", "superfamily"]):
    for j, t in enumerate(TIS):
        a = ax[i, j]; n = nul[(nul.level == lv) & (nul.tissue == t)].n_te_primary_fdr05
        o = out[(out.level == lv) & (out.tissue == t)].iloc[0]
        mx = max(n.max(), o.n_primary_fdr05) + 1
        a.hist(n, bins=np.arange(-0.5, mx + 1.5, 1 if mx < 30 else max(1, mx // 25)), color="#9aa5b1")
        a.axvline(o.n_primary_fdr05, color=C31, lw=2)
        a.set_title(f"{LAB[t]} ({lv})\nobserved {o.n_primary_fdr05}, p = {o.perm_p:.4f}", fontsize=8)
        if j == 0: a.set_ylabel("labelings (of 924)")
        if i == 1: a.set_xlabel("families at FDR < 0.05")
fig.tight_layout(); fig.savefig(f"{O}/fig1_null.png"); plt.close(fig)

# Fig 2: passing classified families, log2 fold change (31 vs 28) with 95 % CI, by tissue
H = pd.read_csv(f"{T}/primary_hits_5_1.tsv", sep="\t")
H = H[H.tissue != "muscle"].sort_values(["tissue", "log2FoldChange"])
cls_col = {"DNA": "#4c72b0", "LINE": "#dd8452", "LTR": "#55a868", "SINE": "#c44e52", "RC": "#8172b3",
           "PLE": "#937860", "Retroposon": "#da8bc3"}
fig, ax = plt.subplots(1, 3, figsize=(10, 5.2), sharex=True)
for a, t in zip(ax, ["brain", "digestive_tract", "liver"]):
    h = H[H.tissue == t].reset_index(drop=True)
    for k, r in h.iterrows():
        a.plot([r.ci95_low, r.ci95_high], [k, k], color="#888", lw=1)
        a.plot(r.log2FoldChange, k, "o", color=cls_col.get(r["class"], "k"), ms=5,
               mfc="white" if not r.passes_uniq else cls_col.get(r["class"], "k"))
    a.set_yticks(range(len(h))); a.set_yticklabels(h.family, fontsize=6)
    a.axvline(0, color="k", lw=0.6); a.set_title(f"{LAB[t]} (n = {len(h)})"); a.set_xlabel("log2 fold change, 31 vs 28 °C tank")
hs = [plt.Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in cls_col.items() if k in set(H["class"])]
hs.append(plt.Line2D([], [], marker="o", ls="", mfc="white", color="#555", label="fails unique-only"))
fig.legend(handles=hs, loc="lower center", ncol=len(hs), frameon=False, fontsize=7)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(f"{O}/fig2_families.png"); plt.close(fig)

# Fig 3: positive control - gene count vs its null (A) and heat-shock genes (B)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [1, 1.4]})
f = out[out.level == "family"]
x = np.arange(4)
for j, t in enumerate(TIS):
    n = nul[(nul.level == "family") & (nul.tissue == t)].n_gene_fdr05
    ax[0].boxplot(n + 1, positions=[j], widths=0.5, showfliers=False)
    ax[0].plot([j - 0.3, j + 0.3], [np.percentile(n, 95) + 1] * 2, color="k", lw=2)
    ax[0].plot(j, f[f.tissue == t].n_gene_fdr05.iloc[0] + 1, "D", color=C31)
ax[0].set_yscale("log")
ax[0].set_xticks(x); ax[0].set_xticklabels([LAB[t] for t in TIS], fontsize=8)
ax[0].set_ylabel("genes at FDR < 0.05 (+1, log)"); ax[0].set_title("A  Genes: observed (red), null box and 95th pct (bar)")
hsp = []
for t in TIS:
    R = pd.read_csv(f"{A}/observed_{t}_multi.tsv", sep="\t"); g = R[R.set == "gene"].copy()
    g["fid"] = g.feature_id.str.strip('"')
    s = g[g.fid.str.match(r"^(hsp|dnaj)", case=False) & (g.padj_set < 0.05)]
    hsp.append(s.assign(tissue=t))
hsp = pd.concat(hsp)
top = hsp.groupby("fid").size().sort_values(ascending=False).index[:12]
mat = pd.DataFrame({t: hsp[hsp.tissue == t].set_index("fid").log2FoldChange for t in TIS}).reindex(top)
im = ax[1].imshow(mat.values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto")
ax[1].set_xticks(range(4)); ax[1].set_xticklabels([LAB[t] for t in TIS], fontsize=8)
ax[1].set_yticks(range(len(top))); ax[1].set_yticklabels(top, fontsize=7)
ax[1].set_title("B  Heat-shock genes at FDR < 0.05 (log2 FC; blank = not passing)")
fig.colorbar(im, ax=ax[1], shrink=0.8); fig.tight_layout(); fig.savefig(f"{O}/fig3_control.png"); plt.close(fig)

# Fig 4: global TE index by tissue and tank
ix = pd.read_csv(f"{T}/te_index_per_sample.tsv", sep="\t")
fig, ax = plt.subplots(1, 4, figsize=(10, 2.8))
for a, t in zip(ax, TIS):
    for k, (temp, c) in enumerate(((28, C28), (31, C31))):
        v = ix[(ix.tissue == t) & (ix.temperature_c == temp)].te_index_pct
        a.plot(np.full(len(v), k) + np.linspace(-0.08, 0.08, len(v)), v, "o", color=c, ms=4)
        a.plot([k - 0.2, k + 0.2], [v.mean()] * 2, color="k", lw=1)
    a.set_xticks([0, 1]); a.set_xticklabels(["28 °C", "31 °C"]); a.set_title(LAB[t]); a.set_xlim(-0.6, 1.6)
ax[0].set_ylabel("strict TE % of counts")
fig.tight_layout(); fig.savefig(f"{O}/fig4_index.png"); plt.close(fig)

# Fig 5: E2 warming vs CO2 on matched families (developmental contrast shown; both in the table)
e2 = pd.read_csv(f"{T}/e2_warming_vs_co2.tsv", sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
rbh = pd.read_csv(f"{L}/rbh_808080.tsv", sep="\t")
rbh["ao"] = rbh.aocel_family.str.split("#").str[0]; rbh["ap"] = rbh.apoly_family.str.split("#").str[0]
oa = pd.read_csv(sys.argv[4] if len(sys.argv) > 4 else "", sep="\t") if len(sys.argv) > 4 else None
for a, t in zip(ax, ["brain", "liver"]):
    W = pd.read_csv(f"{A}/observed_{t}_multi.tsv", sep="\t"); W["fam"] = W.feature_id.str.split(":").str[0]
    if oa is not None:
        o = oa[(oa.analysis_set == "full") & (oa.tissue == t)].copy(); o["fam"] = o.feature_id.str.split(":").str[0]
        o = o.groupby("fam").developmental_vs_control_lfc.mean()
        m = rbh.merge(W[["fam", "log2FoldChange"]], left_on="ao", right_on="fam").merge(o, left_on="ap", right_index=True)
        a.plot(m.log2FoldChange, m.developmental_vs_control_lfc, "o", ms=3, color="#555", alpha=0.7)
    r = e2[(e2.tissue == t) & (e2.oa_contrast == "developmental_vs_control")].iloc[0]
    a.axhline(0, color="#bbb", lw=0.6); a.axvline(0, color="#bbb", lw=0.6)
    a.set_title(f"{LAB[t]}: rho = {r.rho:.2f} [{r.rho_ci_lo:.2f}, {r.rho_ci_hi:.2f}], n = {r.n_matched_tested}", fontsize=8)
    a.set_xlabel("warming log2 FC (A. ocellaris)"); a.set_ylabel("CO2 log2 FC, developmental (A. polyacanthus)")
fig.tight_layout(); fig.savefig(f"{O}/fig5_e2.png"); plt.close(fig)

# Fig 6: tissue structure of TE expression (5.2.4, descriptive PCA over the 500 most variable classified TE families)
pc = pd.read_csv(f"{A}/tissue_model_te_pca.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(4.6, 3.6))
mk = {"brain": "o", "digestive_tract": "s", "liver": "^", "muscle": "D"}
for t in TIS:
    for temp, c in ((28, C28), (31, C31)):
        s = pc[(pc.tissue == t) & (pc.temperature == temp)]
        ax.plot(s.PC1, s.PC2, mk[t], color=c, ms=5, ls="", label=f"{LAB[t]} {temp} °C")
ax.set_xlabel(f"PC1 ({pc.var1.iloc[0]*100:.0f} %)"); ax.set_ylabel(f"PC2 ({pc.var2.iloc[0]*100:.0f} %)")
ax.legend(fontsize=6, frameon=False, ncol=2); fig.tight_layout(); fig.savefig(f"{O}/fig6_tissue_pca.png"); plt.close(fig)

# Fig 7: E1 per tissue, observed count of classified intact-tier families at FDR < 0.05 vs its null (if present)
if os.path.exists(f"{T}/e1_outcomes.tsv"):
    e1 = pd.read_csv(f"{T}/e1_outcomes.tsv", sep="\t"); en = pd.read_csv(f"{T}/e1_null.tsv.gz", sep="\t")
    fig, ax = plt.subplots(2, 4, figsize=(10, 4.4), sharey="row")
    for i, m in enumerate(["unique", "multi"]):
        for j, t in enumerate(TIS):
            a = ax[i, j]; n = en[(en.measure == m) & (en.tissue == t)].n_classified_fdr05
            o = e1[(e1.measure == m) & (e1.tissue == t)].iloc[0]
            a.hist(n, bins=np.arange(-0.5, max(n.max(), o.n_obs) + 1.5, 1), color="#9aa5b1")
            a.axvline(o.n_obs, color=C31, lw=2)
            a.set_title(f"{LAB[t]} ({m})\nobserved {o.n_obs}, p = {o.perm_p:.4f}", fontsize=8)
            if j == 0: a.set_ylabel("labelings (of 924)")
            if i == 1: a.set_xlabel("families at FDR < 0.05")
    fig.tight_layout(); fig.savefig(f"{O}/fig7_e1.png"); plt.close(fig)
print("figures written to", O)
