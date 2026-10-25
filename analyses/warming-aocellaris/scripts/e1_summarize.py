"""E1 (A1.5) outcomes per tissue and measure against the exact labeling null, plus Stage 3 family outcomes and the
shuffle p (descriptive). TANK_TEMPERATURE_SHIFT_INTACT if >= 1 classified family at FDR < 0.05 and p < 0.0125 on
that measure; otherwise NO_DETECTABLE_INTACT_SHIFT. Measures are reported separately and never combined.
Usage: e1_summarize.py STAGE4_DIR STAGE3_DIR TABLES_DIR"""
import sys, glob
import pandas as pd
S4, S3, T = sys.argv[1:4]
TIS = ["brain", "digestive_tract", "liver", "muscle"]
rows, nulls, hits = [], [], []
for t in TIS:
    P = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{S4}/e1_perm_{t}_*.tsv")])
    for m in ["unique", "multi"]:
        q = P[P.measure == m]
        o = q[q.labeling == 0]; perm = q[q.labeling > 0].drop_duplicates("labeling")
        assert len(o) == 1 and len(perm) == 924, (t, m, len(o), len(perm))
        n = int(o.n_classified_fdr05.iloc[0]); p = (perm.n_classified_fdr05 >= n).sum() / 924
        rows.append(dict(tissue=t, measure=m, n_obs=n, perm_p=round(p, 4), n_unknown_obs=int(o.n_unknown_fdr05.iloc[0]),
                         class_not_found=int(o.n_class_not_found.iloc[0]),
                         outcome="TANK_TEMPERATURE_SHIFT_INTACT" if (n >= 1 and p < 0.0125) else "NO_DETECTABLE_INTACT_SHIFT"))
        nulls.append(perm.assign(tissue=t))
        R = pd.read_csv(f"{S4}/e1_observed_{t}_{m}.tsv", sep="\t")
        hits.append(R[(R.set == "classified") & (R.padj_set < 0.05)].assign(tissue=t, measure=m))
pd.DataFrame(rows).to_csv(f"{T}/e1_outcomes.tsv", sep="\t", index=False)
pd.concat(nulls).to_csv(f"{T}/e1_null.tsv.gz", sep="\t", index=False)
H = pd.concat(hits)
fo = pd.read_csv(f"{S3}/family_outcomes.tsv", sep="\t")
H = H.merge(fo[fo.intact_tier == "T1"][["family", "outcome", "n_intact"]], left_on="feature", right_on="family", how="left")
H[["tissue", "measure", "feature", "class", "baseMean", "log2FoldChange", "ci_lo", "ci_hi", "padj_set", "outcome", "n_intact"]] \
    .round(4).to_csv(f"{T}/e1_hits.tsv", sep="\t", index=False)
fo.groupby(["intact_tier", "outcome"]).size().rename("n").reset_index().to_csv(f"{T}/e1_family_outcomes_counts.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string()); print(H[["tissue", "measure", "feature", "class", "log2FoldChange", "padj_set", "outcome"]].to_string())
print(fo.groupby(["intact_tier", "outcome"]).size().to_string())
