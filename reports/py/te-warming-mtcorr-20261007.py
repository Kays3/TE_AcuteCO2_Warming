"""Figures for te-warming-mtcorr-20261007 (Amendment 5, post-hoc mt-rRNA correction of the A. ocellaris warming study).

Usage: python te-warming-mtcorr-20261007.py [TE_OA repo, default ~/workspace/github/TE_Ocean_Acidification]
Run with hive/tools/.venv-report/bin/python. Reads committed tables only; writes
  hive/reports/tsv/te-warming-mtcorr-20261007/   tables copied unchanged from TE_OA, plus derived TSVs
  hive/reports/png/te-warming-mtcorr-20261007/   figures 1 to 5
"""
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import report_figs as rf  # noqa: E402

REPO = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/workspace/github/TE_Ocean_Acidification")
REF = "analysis/te-warming-mtcorr-20261006"
P = "paper/results/te-warming-aocellaris-20261003"
NAME = "te-warming-mtcorr-20261007"
TSV = os.path.join(HERE, "..", "tsv", NAME)
PNG = os.path.join(HERE, "..", "png", NAME)
SCRIPT = os.path.abspath(__file__)
T = ["brain", "digestive_tract", "liver", "muscle"]
TL = {"brain": "Brain", "digestive_tract": "Digestive tract", "liver": "Liver", "muscle": "Muscle"}
R832 = math.ceil(0.9 * 924)
rf.style()
C28, C31 = rf.TREAT["control"], rf.TREAT["treated"]


def table(rel, gz=False):
    dest = os.path.join(TSV, os.path.basename(rel))
    src, sha = rf.export_table(REPO, REF, f"{P}/{rel}", dest)
    return pd.read_csv(dest, sep="\t", compression="gzip" if gz else None), f"TE_OA:{src}", sha


res, res_src, res_sha = table("mtcorr/tables/a5_results.tsv")
sh, sh_src, sh_sha = table("mtcorr/tables/a52_shares.tsv")
nm, nm_src, nm_sha = table("mtcorr/tables/null_mtcorr.tsv.gz", gz=True)
nr, nr_src, nr_sha = table("tables/null_distributions.tsv.gz", gz=True)
ix, ix_src, ix_sha = table("mtcorr/tables/a57_index.tsv")
ixr, ixr_src, ixr_sha = table("tables/index_5_2_2.tsv")
res = res.set_index("tissue")


def nullset(df, tissue, fit=None):
    d = df[(df.tissue == tissue) & (df.labeling >= 1)]
    if fit is not None:
        d = d[(d.fit == fit) & (d["mode"] == "multi")]
    else:
        d = d[d.level == "family"]
    d = d.sort_values("labeling")
    assert len(d) == 924, (tissue, fit, len(d))
    return d


# Figure 1: per-library mt rRNA share of gene counts (a) and the M2 covariate (b), by tissue and temperature.
fig, axs = rf.figure(1, 2, width="full", height=2.7)
rows = []
for k, t in enumerate(T):
    d = sh[sh.tissue == t]
    for temp, col, off in [(28, C28, -0.17), (31, C31, 0.17)]:
        g = d[d.temperature_c == temp]
        x = k + off + (np.arange(len(g)) - (len(g) - 1) / 2) * 0.035
        axs[0].scatter(x, 100 * g.mt_rRNA, s=14, color=col, edgecolor="white", lw=0.3, zorder=3)
        axs[0].hlines(100 * g.mt_rRNA.mean(), k + off - 0.12, k + off + 0.12, color=rf.INK, lw=1)
        axs[1].scatter(x, g.mtrRNA_centred, s=14, color=col, edgecolor="white", lw=0.3, zorder=3)
        axs[1].hlines(g.mtrRNA_centred.mean(), k + off - 0.12, k + off + 0.12, color=rf.INK, lw=1)
    r = np.corrcoef(d.mtrRNA_centred, (d.temperature_c == 31).astype(float))[0, 1]
    rows.append(dict(tissue=t, r_pointbiserial=r, mt_rRNA_pct_min=100 * d.mt_rRNA.min(), mt_rRNA_pct_max=100 * d.mt_rRNA.max()))
    axs[1].text(k, 0.62, f"r = {r:.2f}", ha="center", fontsize=6.5, color=rf.MUTE)
pd.DataFrame(rows).to_csv(os.path.join(TSV, "fig1_covariate_vs_temperature.tsv"), sep="\t", index=False)
for ax in axs:
    ax.set_xticks(range(4), [TL[t] for t in T])
axs[0].set_ylabel("mt 12S + 16S rRNA, % of gene counts")
axs[1].set_ylabel("M2 covariate (centred log ratio)")
axs[1].set_ylim(-0.75, 0.7)
axs[0].scatter([], [], color=C28, label="28 °C tank")
axs[0].scatter([], [], color=C31, label="31 °C tank")
axs[0].legend(loc="upper left", frameon=False, fontsize=6.5)
rf.panel_labels(axs)
rf.save(fig, os.path.join(PNG, "fig1_shares.png"), sh_src, sh_sha, n=48, script=SCRIPT)

# Figure 2: M1 TE-family null strips, one panel per tissue.
fig, axs = rf.figure(2, 2, width="full", height=3.6, layout="constrained")
axs = axs.ravel()
for ax, t in zip(axs, T):
    d = nullset(nm, t, "m1_family")
    obs = int(res.loc[t, "m1_te"])
    rf.null_strip(ax, d.n_te_primary_fdr05, obs)
    ax.set_title(f"{TL[t]}: p = {res.loc[t, 'm1_p']:.4f}", loc="left")
    ax.set_xlim(-0.3, 100)
for ax in axs[2:]:
    ax.set_xlabel("classified TE families at FDR 0.05 (symlog)")
rf.panel_labels(axs)
rf.save(fig, os.path.join(PNG, "fig2_te_null_m1.png"), nm_src, nm_sha, n="924 labelings per tissue", script=SCRIPT)

# Figure 3: relaxed gene control. Row a: registered stored nulls (descriptive); row b: M1 nulls.
fig, axs = rf.figure(2, 4, width="full", height=3.9, layout="constrained")
ctl = []
for j, t in enumerate(T):
    for i, (lab, d, obs) in enumerate([("registered", nullset(nr, t), int(res.loc[t, "reg_gene"])),
                                       ("M1", nullset(nm, t, "m1_family"), int(res.loc[t, "m1_gene"]))]):
        ax = axs[i, j]
        v = np.sort(d.n_gene_fdr05.values)
        thr = v[R832 - 1]
        rf.null_strip(ax, v, obs)
        for tx in list(ax.texts):
            tx.remove()
        ge = int((v >= obs).sum())
        ax.axvline(thr, color=rf.THRESHOLD, ls="--", lw=0.9)
        ax.text(0.01, 0.98, f"observed {obs}\n{ge} of 924 \u2265 obs.\nrank 832: {thr:g}", transform=ax.transAxes, fontsize=6.5, va="top", color=rf.INK)
        ax.set_xlim(-0.3, 20000)
        ax.set_ylim(-0.5, 1.15)
        ax.set_title(f"{TL[t]}, {lab}", fontsize=7.5, loc="left")
        ctl.append(dict(tissue=t, fit=lab, observed=obs, rank832=thr, margin=obs - thr, n_ge_obs=int((v >= obs).sum())))
for ax in axs[-1]:
    ax.set_xlabel("genes at FDR 0.05 (symlog)")
pd.DataFrame(ctl).to_csv(os.path.join(TSV, "fig3_gene_control.tsv"), sep="\t", index=False)
for c in ctl:
    col = "reg_gene_p90" if c["fit"] == "registered" else "m1_gene_p90"
    assert c["rank832"] == res.loc[c["tissue"], col], c
rf.panel_labels(axs.ravel())
rf.save(fig, os.path.join(PNG, "fig3_gene_control.png"), f"{nr_src} ; {nm_src}", f"{nr_sha} ; {nm_sha}",
        n="924 labelings per panel", script=SCRIPT)

# Figure 4: global strict-TE index, 31 minus 28 °C with Welch 95 % CI; registered (filled) and rRNA removed (open).
fig, ax = rf.figure(1, 1, width="full", height=2.9)
rows = []
for t in T:
    a = ixr.set_index("tissue").loc[t]
    b = ix.set_index("tissue").loc[t]
    rows.append(dict(group=TL[t], label="registered", est=a.diff_pp, lo=a.ci_lo, hi=a.ci_hi, color=rf.INK, n=12))
    rows.append(dict(group=TL[t], label="rRNA removed", est=b.diff_pp, lo=b.ci_lo, hi=b.ci_hi, color=rf.INK, filled=False, n=12))
rf.forest(ax, rows, xlabel="31 minus 28 °C, percentage points")
rf.save(fig, os.path.join(PNG, "fig4_index.png"), f"{ixr_src} ; {ix_src}", f"{ixr_sha} ; {ix_sha}",
        n="12 libraries per tissue", script=SCRIPT)

# Figure 5: compute against the 250 CPU-h ceiling (numbers from RUN_RECORD.md, written to a TSV first).
rr_sha = rf.subprocess.check_output(["git", "-C", REPO, "log", "-1", "--format=%h", REF, "--", f"{P}/RUN_RECORD.md"], text=True).strip()
comp = pd.DataFrame([("Registered study (to 4 October)", 220.04), ("A5 M1 family", 5.911), ("A5 M2 family", 4.534),
                     ("A5 M1 superfamily (E3)", 5.526), ("A5 M1 uniq observed", 0.012)], columns=["item", "cpu_h"])
comp.to_csv(os.path.join(TSV, "fig5_compute.tsv"), sep="\t", index=False)
fig, ax = rf.figure(1, 1, width="full", height=1.5)
left = 0
cols = [rf.MUTE, rf.OI["blue"], rf.OI["orange"], rf.OI["green"], rf.INK]
for (item, v), c in zip(comp.itertuples(index=False), cols):
    ax.barh(0, v, left=left, color=c, height=0.5, label=f"{item}: {v:g}")
    left += v
ax.axvline(250, color=rf.THRESHOLD, ls="--", lw=0.9)
ax.text(250, 0.33, " ceiling 250", fontsize=6.5)
ax.axvspan(237, 242, color=rf.GRID, zorder=0)
ax.text(239.5, 0.33, "A5 projection", fontsize=6, ha="center")
ax.set_xlim(200, 255)
ax.set_yticks([])
ax.set_xlabel(f"cumulative CPU-h (total {left:.2f})")
ax.legend(loc="upper left", bbox_to_anchor=(0, -0.65), ncol=2, frameon=False, fontsize=6.3)
rf.save(fig, os.path.join(PNG, "fig5_compute.png"), f"TE_OA:{P}/RUN_RECORD.md", rr_sha, script=SCRIPT)
print("ok", res_sha, nm_sha, nr_sha, ix_sha, ixr_sha, rr_sha)
print(pd.DataFrame(ctl).to_string())
print(pd.read_csv(os.path.join(TSV, "fig1_covariate_vs_temperature.tsv"), sep="\t").round(3).to_string())
