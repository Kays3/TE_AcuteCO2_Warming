"""A5.5 descriptive: the registered heat-shock and DnaJ gene rows (tables/heatshock_dnaj_genes.tsv) read again in the M1 observed fits.

Usage: heatshock_m1.py WORK_DIR REG_TABLE OUT_DIR   (WORK_DIR holds analysis_mtcorr/observed_<tissue>_multi_m1.tsv)"""
import sys
import pandas as pd

W, REG, O = sys.argv[1:]
reg = pd.read_csv(REG, sep="\t"); rows, summ = [], []
for t, r in reg.groupby("tissue"):
    m = pd.read_csv(f"{W}/analysis_mtcorr/observed_{t}_multi_m1.tsv", sep="\t")
    m = m[m.set == "gene"].set_index("feature_id")
    x = r.set_index("fid").join(m[["log2FoldChange", "padj_set"]], rsuffix="_m1", how="left")
    x["tissue"] = t; rows.append(x.reset_index().rename(columns={"index": "fid"}))
    summ.append(dict(tissue=t, n_genes=len(x), n_missing_m1=int(x.log2FoldChange_m1.isna().sum()),
                     n_fdr05_reg=int((x.padj_set < 0.05).sum()), n_fdr05_m1=int((x.padj_set_m1 < 0.05).sum()),
                     n_up_fdr05_m1=int(((x.padj_set_m1 < 0.05) & (x.log2FoldChange_m1 > 0)).sum()),
                     max_abs_lfc_change=float((x.log2FoldChange_m1 - x.log2FoldChange).abs().max())))
pd.concat(rows).to_csv(f"{O}/heatshock_m1_genes.tsv", sep="\t", index=False)
s = pd.DataFrame(summ); s.to_csv(f"{O}/heatshock_m1_summary.tsv", sep="\t", index=False); print(s.to_string())
