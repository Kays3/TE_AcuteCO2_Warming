"""Acute check report figures from the committed tables. Usage: acute_figs.py TABLES_DIR OUT_DIR"""
import sys, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
T, O = sys.argv[1:3]; os.makedirs(O, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
INK, MUTE, CTRL, CO2 = "#1f2a33", "#8a949c", "#5f7f8c", "#c0652f"
CLS = {"LTR": "#2f6f8f", "DNA": "#b5653a", "LINE": "#6a8f3a", "SINE": "#8a5fa0", "RC": "#c9a227"}

# Fig 1: permutation nulls with the observed value
fig, axs = plt.subplots(1, 3, figsize=(7.8, 2.6))
for ax, f, lab, thr in zip(axs, ["brain_perm_te", "gill_perm_te", "brain_gene_null"],
                           ["Brain TE subfamilies\n600 assignments", "Gill TE subfamilies\n210 assignments", "Brain genes (positive control)\n200 drawn assignments"],
                           [None, None, "brain_gene_control"]):
    d = pd.read_csv(f"{T}/primary/{f}.tsv", sep="\t"); obs = d.n_fdr05.iloc[0]; null = d.n_fdr05.iloc[1:]
    bins = np.r_[-0.5, 0.5, np.logspace(0, np.log10(max(null.max(), 10) + 1), 18)]
    ax.hist(null + 0.0, bins=bins, color=MUTE, edgecolor="white", lw=.4); ax.set_xscale("symlog", linthresh=1)
    ax.axvline(obs, color=CO2, lw=1.6); ax.text(obs, ax.get_ylim()[1] * .92, f"  observed {obs}", color=CO2, fontsize=8, va="top")
    if thr:
        g = pd.read_csv(f"{T}/primary/{thr}.tsv", sep="\t").iloc[0]; ax.axvline(g.null_190th_of_200, color=INK, ls="--", lw=.9)
        ax.text(g.null_190th_of_200, ax.get_ylim()[1] * .75, f" 95th pct {int(g.null_190th_of_200)}", fontsize=8)
    ax.set_title(lab, fontsize=8.5); ax.set_xlabel("features at FDR 0.05")
axs[0].set_ylabel("assignments")
fig.tight_layout(); fig.savefig(f"{O}/fig1_nulls.png", dpi=200); plt.close(fig)

# Fig 2: TE subfamily effects (MLE; apeglm shrank every effect to ~1e-5) against mean count
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
for ax, t in zip(axs, ["brain", "gill"]):
    d = pd.read_csv(f"{T}/primary/{t}_te_observed.tsv", sep="\t"); d["c"] = d["class"].str.rstrip("?")
    for c, col in CLS.items():
        x = d[d.c == c]; ax.scatter(x.baseMean, x.lfc_mle, s=7, color=col, alpha=.7, lw=0, label=f"{c} ({len(x)})")
    ax.set_xscale("log"); ax.axhline(0, color=MUTE, lw=.6); ax.set_xlabel("mean normalised count")
    ax.set_title(f"{t.capitalize()}: {len(d)} subfamilies, 0 at FDR 0.05", fontsize=8.5); ax.legend(fontsize=6.5, frameon=False, markerscale=1.5, loc="lower right")
axs[0].set_ylabel("CO2 vs control log2FC (MLE)")
fig.tight_layout(); fig.savefig(f"{O}/fig2_te_effects.png", dpi=200); plt.close(fig)

# Fig 3: library composition. mt 16S rRNA share and gene PCA
p = pd.read_csv(f"{T}/descriptive/gene_pca.tsv", sep="\t")
fig, axs = plt.subplots(2, 2, figsize=(7.6, 5.2), gridspec_kw={"height_ratios": [1, 1.15]})
for j, t in enumerate(["brain", "gill"]):
    d = p[p.tissue == t].reset_index(drop=True); ax = axs[0, j]
    ax.bar(range(len(d)), 100 * d.top_gene_share, color=[CO2 if x == "CO2" else CTRL for x in d.treatment])
    ax.set_xticks(range(len(d))); ax.set_xticklabels([f"{a}\nc{c}" for a, c in zip(d.alias, d.cohort)], fontsize=5.5)
    ax.set_ylim(0, 100); ax.set_title(f"{t.capitalize()}: mt 16S rRNA share of gene counts", fontsize=8.5)
    ax = axs[1, j]
    for tr, col in (("Control", CTRL), ("CO2", CO2)):
        for c, mk in zip([1, 2, 3], "os^"):
            x = d[(d.treatment == tr) & (d.cohort == c)]; ax.scatter(x.PC1, x.PC2, color=col, marker=mk, s=28, label=f"{tr}, cohort {c}")
    for r in d.itertuples(): ax.annotate(r.alias, (r.PC1, r.PC2), fontsize=6, xytext=(3, 2), textcoords="offset points", color=MUTE)
    ax.set_xlabel(f"PC1 ({d.PC1_pct[0]:.0f} %)"); ax.set_ylabel(f"PC2 ({d.PC2_pct[0]:.0f} %)"); ax.set_title(f"{t.capitalize()}: gene PCA (2,000 most variable)", fontsize=8.5)
axs[0, 0].set_ylabel("% of gene counts"); axs[1, 1].legend(fontsize=6, frameon=False, ncol=2)
fig.tight_layout(); fig.savefig(f"{O}/fig3_composition.png", dpi=200); plt.close(fig)

# Fig 4: TE share per library
fig, axs = plt.subplots(1, 2, figsize=(6.2, 2.7)); m = pd.read_csv(f"{T}/descriptive/te_share_model.tsv", sep="\t")
for ax, t in zip(axs, ["brain", "gill"]):
    d = pd.read_csv(f"{T}/descriptive/{t}_te_share.tsv", sep="\t"); r = m[m.tissue == t].iloc[0]
    for i, (tr, col) in enumerate((("Control", CTRL), ("CO2", CO2))):
        for c, mk in zip([1, 2, 3], "os^"):
            x = d[(d.treatment == tr) & (d.cohort == c)]; ax.scatter(i + (c - 2) * .12 + np.zeros(len(x)), x.share, color=col, marker=mk, s=22)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "CO2"]); ax.set_xlim(-.5, 1.5)
    ax.set_title(f"{t.capitalize()}: {r.effect_pp:+.2f} pp [{r.ci_lo:.2f}, {r.ci_hi:.2f}]", fontsize=8.5)
axs[0].set_ylabel("TE share of assigned counts (%)")
fig.tight_layout(); fig.savefig(f"{O}/fig4_te_share.png", dpi=200); plt.close(fig)

# Fig 5: cluster comparison with damselfish OA
c = pd.read_csv(f"{T}/compare/compare_results.tsv", sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(7.8, 2.7), sharey=False)
for ax, side, lab in zip(axs, ["gill_developmental", "gill_intergenerational", "brain_intergenerational"],
                         ["Gill, developmental", "Gill, intergenerational", "Brain, intergenerational"]):
    j = pd.read_csv(f"{T}/compare/compare_{side}.tsv", sep="\t"); r = c[c.side == side].iloc[0]
    ax.scatter(j.d_lfc, j.lfc_mle, s=10, color=CO2 if side.startswith("brain") else "#2f6f8f", alpha=.6, lw=0)
    ax.axhline(0, color=MUTE, lw=.6); ax.axvline(0, color=MUTE, lw=.6)
    ax.set_title(f"{lab}, n = {int(r.n_matched)}\napeglm ρ {r.rho:.2f} [{r.ci_lo:.2f}, {r.ci_hi:.2f}]\nMLE ρ {r.rho_mle:.2f} [{r.mle_ci_lo:.2f}, {r.mle_ci_hi:.2f}]", fontsize=8.5); ax.set_xlabel("damselfish OA cluster log2FC")
axs[0].set_ylabel("zebrafish CO2 cluster log2FC\n(MLE)")
fig.tight_layout(); fig.savefig(f"{O}/fig5_compare.png", dpi=200); plt.close(fig)
print("ok")
