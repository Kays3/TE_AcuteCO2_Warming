## Acute check, descriptive diagnostic: PCA of variance-stabilised gene counts (filtered gene rows, primary size factors, blind vst)
## per tissue, plus each library's share of gene counts in its single largest gene. Usage: R --file=acute_pca.R --args ANALYSIS_DIR SAMPLES
a <- commandArgs(trailingOnly = TRUE); suppressPackageStartupMessages(library(DESeq2)); O <- a[1]; s0 <- read.delim(a[2]); out <- NULL
for (t in c("brain", "gill")) {
  s <- s0[s0$tissue == t, ]; s <- s[order(s$cohort, s$run_accession), ]; runs <- s$run_accession
  l <- lapply(runs, function(r) { x <- read.delim(sprintf("out/%s/%s.multi.cntTable", r, r), check.names = FALSE); setNames(x[[2]], gsub('"', "", x[[1]])) })
  m <- round(do.call(cbind, l)); colnames(m) <- runs; g <- m[!grepl(":", rownames(m)), ]
  rows <- read.delim(file.path(O, sprintf("%s_gene_observed.tsv", t)))$feature; sf <- read.delim(file.path(O, sprintf("%s_size_factors.tsv", t)))
  dds <- DESeqDataSetFromMatrix(g[rows, ], data.frame(row.names = runs, x = factor(rep(1, length(runs)))), ~ 1); sizeFactors(dds) <- sf$size_factor
  v <- assay(vst(dds, blind = TRUE)); v <- v[order(-apply(v, 1, var))[1:2000], ]; pc <- prcomp(t(v)); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  top <- apply(g, 2, max) / colSums(g); topid <- rownames(g)[apply(g, 2, which.max)]
  out <- rbind(out, data.frame(tissue = t, run = runs, alias = s$sample_alias, cohort = s$cohort, treatment = s$co2_treatment, size_factor = sf$size_factor,
                               PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC1_pct = ve[1], PC2_pct = ve[2], top_gene = topid, top_gene_share = top))
}
write.table(out, file.path(O, "gene_pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(out[, c(1, 3, 4, 5, 6, 7, 8, 11, 12)], digits = 3)
