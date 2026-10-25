## E1 Stage 4 (A1.5): warming test per tissue on T1/T1u/SINE-T/T2 family sums, unique and rounded 1/NH,
## size factors from gene counts (poscounts), BH within tissue, with the 5.1 exact labeling null.
## Same structure as the mobilisation stage4_deseq.R fit(); design ~ temperature replaces the CO2 designs.
## BH sets (A1.9): classified families (outcome set) and Unknown-class families (separate, never in the outcome).
## Usage: R --file=e1_stage4_warming.R --args <stage3_dir> <genecounts_dir> <tissue> <first> <last> <outdir>
##   labelings 1..924 = combn(12, 6) over the tissue's runs sorted by run id (31 C set); 0 = observed.
## Single-threaded (A4.2). Run from ~/jim_te_warm.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
s3 <- args[[1]]; gene_dir <- args[[2]]; tissue <- args[[3]]
first <- as.integer(args[[4]]); last <- as.integer(args[[5]]); out <- args[[6]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(DESeq2))
BiocParallel::register(BiocParallel::SerialParam())

design <- read.delim("design.tsv", stringsAsFactors = FALSE)
design <- design[design$tissue == tissue, ]
design <- design[order(design$run), ]
stopifnot(nrow(design) == 12L)
design$sample <- paste0(gsub("_", "", design$tissue), "_", design$run)

gl <- lapply(design$sample, function(s) {
  g <- read.delim(file.path(gene_dir, paste0(s, ".cntTable")), header = TRUE, check.names = FALSE)
  setNames(g[[2]], gsub('"', "", g[[1]]))
})
ids <- Reduce(union, lapply(gl, names))
A <- sapply(gl, function(x) { v <- x[ids]; v[is.na(v)] <- 0; v }); rownames(A) <- ids
parts <- strsplit(ids, ":", fixed = TRUE)
is_te <- lengths(parts) == 3L
G <- A[!grepl(":", ids), ]; storage.mode(G) <- "integer"           # genes, as in stage4_deseq.R
fam_class <- setNames(vapply(parts[is_te], `[`, "", 3L), vapply(parts[is_te], `[`, "", 1L))
sf <- sizeFactors(estimateSizeFactors(DESeqDataSetFromMatrix(
  G, data.frame(row.names = design$sample, x = factor(rep(1:2, 6))), ~ 1), type = "poscounts"))

rd <- function(f) { m <- as.matrix(read.delim(file.path(s3, f), row.names = 1, check.names = FALSE)); m[, design$sample] }
mats <- list(unique = rd("family_sums_unique.tsv.gz"), multi = round(rd("family_sums_multi.tsv.gz")))
mats <- lapply(mats, function(m) { storage.mode(m) <- "integer"; m[rowSums(m) > 0, , drop = FALSE] })

combos <- combn(12L, 6L)
obs31 <- which(design$temperature_c == 31)

fit_one <- function(M, idx) {
  hot <- if (idx == 0L) obs31 else combos[, idx]
  cd <- data.frame(row.names = design$sample,
                   temperature = factor(ifelse(seq_len(12) %in% hot, "31", "28"), levels = c("28", "31")))
  dds <- DESeqDataSetFromMatrix(M, cd, ~ temperature)
  sizeFactors(dds) <- sf
  dds <- DESeq(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, name = "temperature_31_vs_28"))
  r$class <- unname(fam_class[rownames(r)])
  r$set <- ifelse(is.na(r$class), "class_not_found", ifelse(r$class == "Unknown", "unknown", "classified"))
  r$padj_set <- NA_real_
  for (s in c("classified", "unknown")) {
    w <- r$set == s & !is.na(r$pvalue)
    r$padj_set[w] <- p.adjust(r$pvalue[w], "BH")
  }
  r
}

rows <- list()
for (m in names(mats)) for (idx in seq(first, last)) {
  r <- fit_one(mats[[m]], idx)
  rows[[length(rows) + 1L]] <- data.frame(measure = m, labeling = idx,
    n_classified_fdr05 = sum(r$padj_set[r$set == "classified"] < 0.05, na.rm = TRUE),
    n_unknown_fdr05 = sum(r$padj_set[r$set == "unknown"] < 0.05, na.rm = TRUE),
    n_class_not_found = sum(r$set == "class_not_found"))
  if (idx == 0L) {
    r$feature <- rownames(r); r$ci_lo <- r$log2FoldChange - 1.96 * r$lfcSE; r$ci_hi <- r$log2FoldChange + 1.96 * r$lfcSE
    write.table(r, file.path(out, sprintf("e1_observed_%s_%s.tsv", tissue, m)), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
write.table(do.call(rbind, rows), file.path(out, sprintf("e1_perm_%s_%04d_%04d.tsv", tissue, first, last)),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Observed labeling only, descriptive, as in the mobilisation Stage 4: C1 family sums (1/NH) fold-change
## distribution, detected T1 loci (unique), host genes of detected intronic/exonic T1 loci.
if (first == 0L) {
  desc <- function(M, label) {
    M <- M[rowSums(M) > 0, , drop = FALSE]
    if (!nrow(M)) return(NULL)
    r <- fit_one(M, 0L); r$feature <- rownames(r); r$label <- label; r
  }
  C1 <- round(rd("family_sums_multi_C1.tsv.gz")); storage.mode(C1) <- "integer"
  L <- rd("locus_unique_noC1.tsv.gz"); storage.mode(L) <- "integer"
  calls <- read.delim(file.path(s3, "locus_calls.tsv.gz"), stringsAsFactors = FALSE, colClasses = c(rm_id = "character"))
  tcol <- paste0(gsub("_", "", tissue), "_detected")
  det <- calls[calls$tier == "T1" & calls[[tcol]] %in% c(TRUE, "True"), ]
  hosts <- intersect(unique(na.omit(det$host_gene)), rownames(G))
  d <- list(desc(C1, "C1_multi"),
            if (nrow(det)) desc(L[intersect(det$rm_id, rownames(L)), , drop = FALSE], "T1loci_unique"),
            if (length(hosts)) desc(G[hosts, , drop = FALSE], "host_gene"))
  d <- d[!vapply(d, is.null, TRUE)]
  if (length(d)) {
    cols <- Reduce(intersect, lapply(d, names))
    write.table(do.call(rbind, lapply(d, `[`, cols)), file.path(out, sprintf("e1_descriptive_%s.tsv", tissue)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
}
