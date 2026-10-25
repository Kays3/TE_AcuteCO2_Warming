## Amendment 5 (POST-HOC) variant of deseq_fit.R: identical procedure, except that the 95 gene_biotype rRNA gene rows
## (A5.1, listed in rrna_genes.txt) are dropped before the filter, size factors and every fit (M1), and model m2 adds the
## per-library covariate mtrRNA = log((KEG47_r01 + KEG47_r02 + 0.5) / (non-rRNA gene counts + 0.5)), centred within tissue (A5.4).
## Also writes per-library shares (A5.2). Usage: R --file=deseq_fit_mtcorr.R --args <tissue> <mode multi|uniq> <first> <last>
##   <outdir> <level family|superfamily> <model m1|m2>
##   level superfamily (E3, A1.7): TE family rows are summed to class:superfamily before the filter;
##   gene and other rows are unchanged; the same filter, poscounts and BH sets apply to that table.
##   labelings are indexed 1..924 over combn(12, 6) of the tissue's runs sorted by run id
##   (the set chosen is the 31 C group); index 0 means the observed labeling.
## Single-threaded (Amendment 4, A4.2). Run from ~/jim_te_warm.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7L)
tissue <- args[[1]]; mode <- args[[2]]
first <- as.integer(args[[3]]); last <- as.integer(args[[4]]); outdir <- args[[5]]
level <- args[[6]]; model <- args[[7]]
stopifnot(level %in% c("family", "superfamily"), model %in% c("m1", "m2"))
tag <- paste0(if (level == "family") mode else paste0(mode, "_sf"), "_", model)
rrna_genes <- readLines("rrna_genes.txt"); stopifnot(length(rrna_genes) == 95L)
MT <- c("KEG47_r01", "KEG47_r02")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(DESeq2))
BiocParallel::register(BiocParallel::SerialParam())

primary_classes <- c("DNA", "LINE", "SINE", "LTR", "RC", "PLE", "Retroposon")  # A1.9
design <- read.delim("design.tsv", stringsAsFactors = FALSE)
design <- design[design$tissue == tissue, ]
design <- design[order(design$run), ]
stopifnot(nrow(design) == 12L)

## Prepared once per tissue and mode: joint table, rRNA rows dropped (A3 addendum),
## filter >= 10 counts in >= 6 libraries, poscounts size factors; all label-free.
prep_path <- file.path(outdir, sprintf("prep_%s_%s.rds", tissue, tag))
if (file.exists(prep_path)) {
  prep <- readRDS(prep_path)
} else {
  read_one <- function(run) {
    tab <- read.delim(sprintf("out/%s/%s.%s.cntTable", run, run, mode),
                      header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    setNames(as.integer(round(tab[[2]])), tab[[1]])
  }
  cols <- lapply(design$run, read_one)
  ids <- unique(unlist(lapply(cols, names)))
  m <- matrix(0L, length(ids), 12L, dimnames = list(ids, design$run))
  for (i in seq_along(cols)) m[names(cols[[i]]), i] <- cols[[i]]
  parts <- strsplit(rownames(m), ":", fixed = TRUE)
  is_te <- lengths(parts) == 3L
  cls <- ifelse(is_te, vapply(parts, function(p) p[length(p)], ""), NA_character_)
  rrna <- is_te & cls == "rRNA"
  rrna_share <- colSums(m[rrna, , drop = FALSE]) / colSums(m)
  stopifnot(all(rrna_genes %in% rownames(m)))
  gene_tot <- colSums(m[!is_te, , drop = FALSE]); rg <- rownames(m) %in% rrna_genes
  shares <- data.frame(run = colnames(m), mt16S = m["KEG47_r01", ] / gene_tot, mt_rRNA = colSums(m[MT, ]) / gene_tot,
                       all_rRNA_genes = colSums(m[rg, ]) / gene_tot, gene_counts = gene_tot)
  mtcov <- log((colSums(m[MT, ]) + 0.5) / (gene_tot - colSums(m[rg, ]) + 0.5)); mtcov <- mtcov - mean(mtcov)
  drop <- rrna | rg
  m <- m[!drop, ]; is_te <- is_te[!drop]; cls <- cls[!drop]
  if (level == "superfamily") {
    sf_id <- ifelse(is_te, vapply(strsplit(rownames(m), ":", fixed = TRUE),
                                  function(p) paste0("SF:", p[2], ":", p[3]), ""), rownames(m))
    m <- rowsum(m, sf_id, reorder = FALSE)
    is_te <- startsWith(rownames(m), "SF:")
    cls <- ifelse(is_te, vapply(strsplit(rownames(m), ":", fixed = TRUE), function(p) p[3], ""), NA_character_)
  }
  keep <- rowSums(m >= 10L) >= 6L
  m <- m[keep, ]; is_te <- is_te[keep]; cls <- cls[keep]
  sf <- estimateSizeFactors(DESeqDataSetFromMatrix(
    m, data.frame(row.names = colnames(m), x = factor(rep(1:2, 6))), ~ 1),
    type = "poscounts")
  feature_set <- ifelse(!is_te, "gene",
                 ifelse(cls %in% primary_classes, "te_primary",
                 ifelse(cls == "Unknown", "te_unknown", "repeat_other")))
  prep <- list(m = m, sf = sizeFactors(sf), cls = cls, set = feature_set,
               rrna_share = rrna_share, runs = design$run, shares = shares, mtcov = mtcov)
  saveRDS(prep, prep_path)
}

combos <- combn(12L, 6L)
obs31 <- which(design$temperature_c == 31)
obs_idx <- which(apply(combos, 2, function(cc) setequal(cc, obs31)))
stopifnot(length(obs_idx) == 1L)

fit_one <- function(idx) {
  hot <- if (idx == 0L) obs31 else combos[, idx]
  temp <- factor(ifelse(seq_len(12) %in% hot, "31", "28"), levels = c("28", "31"))
  cd <- data.frame(row.names = colnames(prep$m), temperature = temp, mtrRNA = prep$mtcov[colnames(prep$m)])
  dds <- DESeqDataSetFromMatrix(prep$m, cd, if (model == "m2") ~ mtrRNA + temperature else ~ temperature)
  sizeFactors(dds) <- prep$sf
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "temperature_31_vs_28", independentFiltering = FALSE)
  res <- as.data.frame(res)
  res$set <- prep$set; res$class <- prep$cls
  res$padj_set <- NA_real_
  for (s in c("te_primary", "te_unknown", "gene")) {
    w <- res$set == s & !is.na(res$pvalue)
    res$padj_set[w] <- p.adjust(res$pvalue[w], "BH")
  }
  ## 5.2.3: median |LFC| of primary TE families over the median |LFC| of genes in the
  ## same gene-baseMean decile, matched decile by decile (weighted by TE count per decile).
  g <- res$set == "gene" & !is.na(res$log2FoldChange)
  t <- res$set == "te_primary" & !is.na(res$log2FoldChange)
  br <- unique(quantile(res$baseMean[g], probs = seq(0, 1, 0.1)))
  dg <- cut(res$baseMean[g], br, include.lowest = TRUE)
  dt <- cut(pmin(pmax(res$baseMean[t], min(br)), max(br)), br, include.lowest = TRUE)
  gmed <- tapply(abs(res$log2FoldChange[g]), dg, median)
  matched_gene <- median(gmed[as.character(dt)], na.rm = TRUE)
  stat <- median(abs(res$log2FoldChange[t])) / matched_gene
  list(res = res, row = data.frame(
    labeling = idx, observed = (idx == 0L || idx == obs_idx),
    n_te_primary_fdr05 = sum(res$padj_set[res$set == "te_primary"] < 0.05, na.rm = TRUE),
    n_te_unknown_fdr05 = sum(res$padj_set[res$set == "te_unknown"] < 0.05, na.rm = TRUE),
    n_gene_fdr05 = sum(res$padj_set[res$set == "gene"] < 0.05, na.rm = TRUE),
    te_gene_matched_ratio = stat))
}

rows <- list()
for (idx in seq(first, last)) {
  t0 <- proc.time()
  f <- fit_one(idx)
  el <- (proc.time() - t0)[["elapsed"]]
  f$row$fit_seconds <- el
  rows[[length(rows) + 1L]] <- f$row
  if (idx == 0L) {
    out <- f$res
    out$ci95_low <- out$log2FoldChange - qnorm(0.975) * out$lfcSE
    out$ci95_high <- out$log2FoldChange + qnorm(0.975) * out$lfcSE
    out$feature_id <- rownames(out)
    write.table(out, file.path(outdir, sprintf("observed_%s_%s.tsv", tissue, tag)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(cbind(prep$shares, mtrRNA_centred = prep$mtcov, rm_rRNA_class_share = prep$rrna_share),
                file.path(outdir, sprintf("rrna_share_%s_%s.tsv", tissue, tag)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
write.table(do.call(rbind, rows),
            file.path(outdir, sprintf("perm_%s_%s_%04d_%04d.tsv", tissue, tag, first, last)),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("observed labeling index in combn order: %d\n", obs_idx))
