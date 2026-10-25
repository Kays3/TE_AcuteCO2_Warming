"""Acute check s5, damselfish side of the exploratory cluster comparison (RUN_RECORD choice 7; Tier A Q3 rule):
full-cohort OA results, analysis_set full, strict_te TRUE, test_status tested; family value = mean log2FoldChange over the
two parental phenotypes; cluster value = unweighted mean over its tested families (A. polyacanthus cluster members).
Sides: gill developmental, gill intergenerational, brain intergenerational.
Usage: acute_damsel.py SHARED_TSV CLUSTERS_TSV OA_ALL_TE_GZ OUT_TSV"""
import sys
import pandas as pd
SH, CL, OA, O = sys.argv[1:5]
APOLY = "GCF_021347895.1"
sh = pd.read_csv(SH, sep="\t")[["cluster", "class"]]
cl = pd.read_csv(CL, sep="\t"); cl = cl[cl.cluster.isin(sh.cluster) & (cl.accession == APOLY)][["family", "cluster"]]
oa = pd.read_csv(OA, sep="\t", usecols=["family_id", "strict_te", "analysis_set", "tissue", "parental_phenotype",
                                         "contrast", "baseMean", "log2FoldChange", "padj", "test_status"])
oa = oa[(oa.analysis_set == "full") & (oa.strict_te.astype(str).str.upper() == "TRUE") & (oa.test_status == "tested")].rename(columns={"family_id": "family"})
out = []
for tis, con in (("gill", "developmental_vs_control"), ("gill", "intergenerational_vs_control"), ("brain", "intergenerational_vs_control")):
    r = oa[(oa.tissue == tis) & (oa.contrast == con)]
    assert r.parental_phenotype.nunique() == 2, (tis, con, r.parental_phenotype.unique())
    f = r.groupby("family").agg(lfc=("log2FoldChange", "mean"), bm=("baseMean", "mean")).reset_index()
    m = cl.merge(f, on="family").groupby("cluster").agg(d_lfc=("lfc", "mean"), d_bm=("bm", "sum"), n_fam=("family", "size")).reset_index()
    m.insert(0, "side", f"{tis}_{con.split('_')[0]}"); out.append(m.merge(sh, on="cluster"))
    print(tis, con, len(r), "rows", len(m), "clusters")
pd.concat(out).to_csv(O, sep="\t", index=False)
