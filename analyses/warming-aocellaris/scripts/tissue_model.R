## 5.2.4 tissue contrasts (descriptive): ~ tissue + temperature over all 48 libraries, multi mode.
## The fish term cannot be fitted (libraries cannot be mapped to fish), so the 48 libraries are treated as
## independent and p-values are too small: only effect sizes are written. Same row rules as 5.1
## (rRNA rows dropped, >= 10 counts in >= 6 libraries, poscounts). Reference tissue: brain. Single-threaded.
## Usage: R --file=tissue_model.R --args <outdir>   (run from ~/jim_te_warm)
args <- commandArgs(trailingOnly = TRUE); out <- args[[1]]
suppressPackageStartupMessages(library(DESeq2))
BiocParallel::register(BiocParallel::SerialParam())
d <- read.delim("design.tsv", stringsAsFactors = FALSE); d <- d[order(d$tissue, d$run), ]
cols <- lapply(d$run, function(r) {
  tab <- read.delim(sprintf("out/%s/%s.multi.cntTable", r, r), header = TRUE, check.names = FALSE)
  setNames(as.integer(round(tab[[2]])), tab[[1]])
})
ids <- unique(unlist(lapply(cols, names)))
m <- matrix(0L, length(ids), nrow(d), dimnames = list(ids, d$run))
for (i in seq_along(cols)) m[names(cols[[i]]), i] <- cols[[i]]
parts <- strsplit(ids, ":", fixed = TRUE); is_te <- lengths(parts) == 3L
cls <- ifelse(is_te, vapply(parts, function(p) p[length(p)], ""), NA_character_)
keep <- !(is_te & cls == "rRNA") & rowSums(m >= 10L) >= 6L
m <- m[keep, ]; cls <- cls[keep]
cd <- data.frame(row.names = d$run, tissue = relevel(factor(d$tissue), "brain"),
                 temperature = factor(d$temperature_c, levels = c(28, 31)))
dds <- DESeqDataSetFromMatrix(m, cd, ~ tissue + temperature)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, quiet = TRUE)
res <- do.call(rbind, lapply(grep("^tissue_", resultsNames(dds), value = TRUE), function(w) {
  r <- results(dds, name = w, independentFiltering = FALSE)
  data.frame(contrast = w, feature = rownames(r), class = cls, baseMean = r$baseMean,
             log2FC = r$log2FoldChange, lfcSE = r$lfcSE)
}))
res <- res[!is.na(res$class) & res$class %in% c("DNA", "LINE", "SINE", "LTR", "RC", "PLE", "Retroposon", "Unknown"), ]
write.table(res, file.path(out, "tissue_model_te_effects.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
vs <- assay(vst(dds, blind = TRUE))
te <- !is.na(cls) & cls %in% c("DNA", "LINE", "SINE", "LTR", "RC", "PLE", "Retroposon")
pc <- prcomp(t(vs[te, ][order(-apply(vs[te, ], 1, var))[1:500], ]))
write.table(data.frame(run = rownames(pc$x), tissue = d$tissue, temperature = d$temperature_c, pc$x[, 1:2],
                       var1 = summary(pc)$importance[2, 1], var2 = summary(pc)$importance[2, 2]),
            file.path(out, "tissue_model_te_pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
