## Acute check s5, zebrafish side of the exploratory cluster comparison (RUN_RECORD choice 7): featureCounts 1/NH cluster counts
## (R1 GTF), clusters with >= 10 counts in >= 4 libraries of the tissue, DESeq2 ~ cohort + treatment with the tissue's primary size
## factors, apeglm log2FC. Usage: R --file=acute_clusters.R --args ANALYSIS_DIR SAMPLES
a <- commandArgs(trailingOnly = TRUE); suppressPackageStartupMessages(library(DESeq2)); O <- a[1]; s0 <- read.delim(a[2])
for (t in c("brain", "gill")) {
  s <- s0[s0$tissue == t, ]; s <- s[order(s$cohort, s$run_accession), ]; runs <- s$run_accession
  l <- lapply(runs, function(r) { x <- read.delim(sprintf("fc/%s.frac.txt", r), comment.char = "#", check.names = FALSE); setNames(x[[7]], x$Geneid) })
  m <- round(do.call(cbind, l)); colnames(m) <- runs; m <- m[rowSums(m >= 10) >= 4, ]
  sf <- read.delim(file.path(O, sprintf("%s_size_factors.tsv", t))); stopifnot(identical(sf$run, runs))
  cd <- data.frame(row.names = runs, cohort = factor(s$cohort), treatment = factor(s$co2_treatment, levels = c("Control", "CO2")))
  dds <- DESeqDataSetFromMatrix(m, cd, ~ cohort + treatment); sizeFactors(dds) <- sf$size_factor; dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = "treatment_CO2_vs_Control", independentFiltering = FALSE)
  sh <- lfcShrink(dds, coef = "treatment_CO2_vs_Control", type = "apeglm", quiet = TRUE)
  out <- data.frame(cluster = rownames(m), baseMean = r$baseMean, lfc_mle = r$log2FoldChange, lfc_apeglm = sh$log2FoldChange, padj = r$padj)
  write.table(out, file.path(O, sprintf("%s_clusters_zf.tsv", t)), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(t, "covered clusters", nrow(m), "padj<.05", sum(out$padj < .05, na.rm = TRUE), "\n")
}
