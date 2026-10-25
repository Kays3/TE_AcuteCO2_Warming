## Acute check (PREREGISTRATION s4-5 at 77270d9; RUN_RECORD implementation choices): DESeq2 1.46.0 per tissue, ~ cohort + treatment.
## One procedure for observed and permuted fits: rows and size factors fixed once per tissue (08 rule: >= 10 counts in at least as
## many libraries as the smaller group, on genes + TE-set rows), dispersions re-estimated per fit, DESeq2 defaults otherwise,
## independentFiltering = FALSE so BH runs over all tested rows of the set. Statistic = number of rows at padj < 0.05.
## Usage: R --no-echo --file=acute_deseq.R --args OUT_DIR SAMPLES ASSIGN_DIR TISSUE JOB [CORES]
##   JOB: observed | perm_te | gene_null | sens   (sens = unique counts, label exclusions, no cohort term; observed labels only)
a <- commandArgs(trailingOnly = TRUE); suppressPackageStartupMessages({library(DESeq2); library(parallel)})
O <- a[1]; s <- read.delim(a[2]); A <- a[3]; tis <- a[4]; job <- a[5]; nc <- if (length(a) > 5) as.integer(a[6]) else 1
TECL <- c("DNA", "DNA?", "LINE", "LTR", "SINE", "RC")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
s <- s[s$tissue == tis, ]; s <- s[order(s$cohort, s$run_accession), ]; runs <- s$run_accession
rd <- function(r, mode) { x <- read.delim(sprintf("out/%s/%s.%s.cntTable", r, r, mode), check.names = FALSE); setNames(x[[2]], gsub('"', "", x[[1]])) }
mat <- function(mode, rr = runs) { l <- lapply(rr, rd, mode = mode); ids <- names(l[[1]]); stopifnot(all(sapply(l, function(v) identical(names(v), ids))))
  m <- do.call(cbind, l); dimnames(m) <- list(ids, rr); round(m) }
cls <- function(ids) sapply(strsplit(ids, ":"), function(v) v[length(v)])
prep <- function(mode, rr = runs) {
  m <- mat(mode, rr); te <- grepl(":", rownames(m)); tes <- te & cls(rownames(m)) %in% TECL; gene <- !te
  lab <- s$co2_treatment[match(rr, runs)]; k <- min(table(lab)); keep <- rowSums(m >= 10) >= k
  sf <- estimateSizeFactors(DESeqDataSetFromMatrix(m[keep & (gene | tes), ], data.frame(row.names = rr, x = factor(rep(1, length(rr)))), ~ 1))
  list(m = m, sf = sizeFactors(sf), te_rows = rownames(m)[keep & tes], gene_rows = rownames(m)[keep & gene], k = k, all_te = rownames(m)[tes])
}
fit <- function(p, rows, labels, rr = colnames(p$m), cohort = TRUE, shrink = FALSE) {
  cd <- data.frame(row.names = rr, cohort = factor(s$cohort[match(rr, runs)]), treatment = factor(labels, levels = c("Control", "CO2")))
  dds <- DESeqDataSetFromMatrix(p$m[rows, rr], cd, if (cohort) ~ cohort + treatment else ~ treatment)
  sizeFactors(dds) <- p$sf[rr]; dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, name = "treatment_CO2_vs_Control", independentFiltering = FALSE)
  out <- data.frame(feature = rows, baseMean = r$baseMean, lfc_mle = r$log2FoldChange, lfcSE = r$lfcSE, pvalue = r$pvalue, padj = r$padj)
  if (shrink) out$lfc_apeglm <- lfcShrink(dds, coef = "treatment_CO2_vs_Control", type = "apeglm", quiet = TRUE)$log2FoldChange
  out
}
nsig <- function(r) sum(r$padj < 0.05, na.rm = TRUE)
p <- prep("multi")
asg <- read.delim(file.path(A, sprintf("%s_assignments.tsv", tis)), check.names = FALSE)
lab_of <- function(i) as.character(unlist(asg[asg$assignment == i, runs]))
stopifnot(identical(lab_of(0), s$co2_treatment))
t0 <- proc.time()
if (job == "observed") {
  write.table(data.frame(run = runs, size_factor = p$sf, te_rows = length(p$te_rows), gene_rows = length(p$gene_rows), filter_k = p$k),
              file.path(O, sprintf("%s_size_factors.tsv", tis)), sep = "\t", quote = FALSE, row.names = FALSE)
  te <- fit(p, p$te_rows, lab_of(0), shrink = TRUE); te$class <- cls(te$feature)
  write.table(te, file.path(O, sprintf("%s_te_observed.tsv", tis)), sep = "\t", quote = FALSE, row.names = FALSE)
  g <- fit(p, p$gene_rows, lab_of(0))
  write.table(g, file.path(O, sprintf("%s_gene_observed.tsv", tis)), sep = "\t", quote = FALSE, row.names = FALSE)
  sh <- data.frame(run = runs, cohort = s$cohort, treatment = s$co2_treatment, te_set = colSums(p$m[p$all_te, ]), assigned = colSums(p$m))
  write.table(sh, file.path(O, sprintf("%s_te_share_input.tsv", tis)), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(tis, "TE rows", length(p$te_rows), "TE sig", nsig(te), "gene rows", length(p$gene_rows), "gene sig", nsig(g), "\n")
} else if (job %in% c("perm_te", "gene_null")) {
  ids <- if (job == "perm_te") asg$assignment else c(0, read.delim(file.path(A, "brain_gene_null_draws.tsv"))$assignment)
  rows <- if (job == "perm_te") p$te_rows else p$gene_rows
  st <- unlist(mclapply(ids, function(i) nsig(fit(p, rows, lab_of(i))), mc.cores = nc, mc.preschedule = FALSE))
  write.table(data.frame(assignment = ids, n_fdr05 = st), file.path(O, sprintf("%s_%s.tsv", tis, job)), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(tis, job, "observed", st[1], "n", length(st), "\n")
} else if (job == "sens") {
  res <- list()
  pu <- prep("uniq"); res$unique_counts <- nsig(fit(pu, pu$te_rows, lab_of(0)))
  ex <- if (tis == "brain") s$run_accession[s$sample_alias %in% c("EB19", "EB21")] else s$run_accession[s$sample_alias %in% c("EG19", "EG22")]
  rr <- setdiff(runs, ex); lab <- s$co2_treatment[match(rr, runs)]; k <- min(table(lab))
  rows <- rownames(p$m)[rowSums(p$m[, rr] >= 10) >= k & rownames(p$m) %in% p$all_te]
  res[[sprintf("exclude_%s", paste(s$sample_alias[match(ex, runs)], collapse = "_"))]] <- nsig(fit(p, rows, lab, rr = rr))
  res$no_cohort_term <- nsig(fit(p, p$te_rows, lab_of(0), cohort = FALSE))
  write.table(data.frame(sensitivity = names(res), n_fdr05 = unlist(res)), file.path(O, sprintf("%s_sensitivity.tsv", tis)), sep = "\t", quote = FALSE, row.names = FALSE)
  print(unlist(res))
}
tt <- proc.time() - t0; cat(sprintf("%s %s CPU-s self %.0f children %.0f elapsed %.0f\n", tis, job, tt[1] + tt[2], tt[4] + tt[5], tt[3]))
