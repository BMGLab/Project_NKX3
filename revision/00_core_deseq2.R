#!/usr/bin/env Rscript
# Revision core: one DESeq2 model on all 12 samples for two quantifications.
#   fc  = STAR genome alignments counted with featureCounts (unique reads, reverse-stranded, exon/gene_id)
#   sal = original manuscript input (nf-core STAR->Salmon, tximport lengthScaledTPM, protein-coding tx2gene; v15)
# Design ~ NKX + CM + NKX:CM; per-context contrasts are taken as model coefficients after re-levelling,
# so every contrast can be shrunk with apeglm. Additive ~ NKX + CM gives the main effects.
# Sign convention everywhere: treatment minus reference (CM - ctrl, NKX3.1 - vector).
# Run with: /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript scripts/revision/00_core_deseq2.R

suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(BiocParallel)
  library(readr); library(dplyr); library(AnnotationDbi); library(org.Hs.eg.db)
})

proj <- "/mnt/ssd1/Projects/NKX3_wd"
out  <- file.path(proj, "results/revision_2026/core")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
register(MulticoreParam(16))

min_count <- 10   # pre-filter: >= min_count reads ...
min_samps <- 3    # ... in >= min_samps samples (smallest group size)

# ---------------------- samples ----------------------
st <- read_csv(file.path(proj, "results/sample_table.csv"), show_col_types = FALSE)
cd <- data.frame(sample_id = st$sample_id, condition = st$condition, stringsAsFactors = FALSE)
cd$NKX   <- factor(ifelse(grepl("NKX", cd$condition), "NKX3.1", "vector"), levels = c("vector", "NKX3.1"))
cd$CM    <- factor(ifelse(grepl("CM$", cd$condition), "CM", "ctrl"), levels = c("ctrl", "CM"))
cd$group <- factor(paste(cd$NKX, cd$CM, sep = "_"),
                   levels = c("vector_ctrl", "vector_CM", "NKX3.1_ctrl", "NKX3.1_CM"))
cd$label <- c(vector_ctrl = "Control-basal", vector_CM = "Control-inflammation",
              NKX3.1_ctrl = "NKX3.1-basal", NKX3.1_CM = "NKX3.1-inflammation")[as.character(cd$group)]
cd$num   <- as.integer(sub("^24-?004-(\\d+)_.*$", "\\1", cd$sample_id))
cd <- cd[order(cd$group, cd$num), ]
cd$rep   <- ave(seq_len(nrow(cd)), cd$group, FUN = seq_along)
rownames(cd) <- cd$sample_id

# ---------------------- counts ----------------------
fc_tab <- read.delim(file.path(proj, "results/revision_checks_2026-10/featureCounts_unique_gene_counts.tsv"),
                     check.names = FALSE)
ann <- fc_tab[, c("Geneid", "gene_name", "gene_biotype")]
colnames(ann) <- c("gene_id", "symbol", "biotype")
ann$entrez <- suppressMessages(mapIds(org.Hs.eg.db, keys = ann$gene_id, column = "ENTREZID",
                                      keytype = "ENSEMBL", multiVals = "first"))
rownames(ann) <- ann$gene_id
counts_fc <- as.matrix(fc_tab[, cd$sample_id]); rownames(counts_fc) <- fc_tab$Geneid

dds_v15 <- readRDS(file.path(proj, "results/rna_downstreams_out_v15/dds_prepped_for_deseq2.rds"))
counts_sal <- counts(dds_v15)[, cd$sample_id]
vsd_v15 <- readRDS(file.path(proj, "results/rna_downstreams_out_v15/vsd.rds"))
vst_sal_v15 <- SummarizedExperiment::assay(vsd_v15)[, cd$sample_id]   # the matrix behind Fig 1B-D, 3, GSEA

# ---------------------- model fitting ----------------------
contrast_specs <- list(
  CM_in_vector = list(levels = list(NKX = c("vector", "NKX3.1"), CM = c("ctrl", "CM")), coef = "CM_CM_vs_ctrl"),
  CM_in_NKX    = list(levels = list(NKX = c("NKX3.1", "vector"), CM = c("ctrl", "CM")), coef = "CM_CM_vs_ctrl"),
  NKX_in_ctrl  = list(levels = list(NKX = c("vector", "NKX3.1"), CM = c("ctrl", "CM")), coef = "NKX_NKX3.1_vs_vector"),
  NKX_in_CM    = list(levels = list(NKX = c("vector", "NKX3.1"), CM = c("CM", "ctrl")), coef = "NKX_NKX3.1_vs_vector"),
  interaction  = list(levels = list(NKX = c("vector", "NKX3.1"), CM = c("ctrl", "CM")), coef = "NKXNKX3.1.CMCM")
)

tidy_res <- function(dds, coef, contrast_name) {
  r  <- results(dds, name = coef, alpha = 0.05)
  rT <- results(dds, name = coef, alpha = 0.05, lfcThreshold = log2(1.5), altHypothesis = "greaterAbs")
  s  <- lfcShrink(dds, coef = coef, type = "apeglm", res = r, quiet = TRUE, parallel = TRUE)
  data.frame(gene_id = rownames(r), contrast = contrast_name,
             baseMean = r$baseMean, log2FC_MLE = r$log2FoldChange, lfcSE_MLE = r$lfcSE,
             stat = r$stat, pvalue = r$pvalue, padj = r$padj,
             log2FC_shrunk = s$log2FoldChange, lfcSE_shrunk = s$lfcSE,
             padj_LFC1.5 = rT$padj, stringsAsFactors = FALSE)
}

fit_all <- function(counts, label) {
  keep <- rowSums(counts >= min_count) >= min_samps
  message(sprintf("[%s] pre-filter keeps %d / %d genes (>= %d reads in >= %d samples)",
                  label, sum(keep), nrow(counts), min_count, min_samps))
  m <- round(counts[keep, cd$sample_id])
  res <- list(); dds_full <- NULL
  fits <- list()
  for (cn in names(contrast_specs)) {
    sp <- contrast_specs[[cn]]
    key <- paste(sp$levels$NKX[1], sp$levels$CM[1])
    if (is.null(fits[[key]])) {
      cdx <- cd
      cdx$NKX <- factor(as.character(cdx$NKX), levels = sp$levels$NKX)
      cdx$CM  <- factor(as.character(cdx$CM),  levels = sp$levels$CM)
      d <- DESeqDataSetFromMatrix(m, cdx, ~ NKX + CM + NKX:CM)
      fits[[key]] <- DESeq(d, quiet = TRUE, parallel = TRUE)
    }
    res[[cn]] <- tidy_res(fits[[key]], sp$coef, cn)
  }
  dds_full <- fits[["vector ctrl"]]
  d_add <- DESeq(DESeqDataSetFromMatrix(m, cd, ~ NKX + CM), quiet = TRUE, parallel = TRUE)
  res$NKX_main <- tidy_res(d_add, "NKX_NKX3.1_vs_vector", "NKX_main")
  res$CM_main  <- tidy_res(d_add, "CM_CM_vs_ctrl", "CM_main")
  list(dds_full = dds_full, dds_add = d_add, res = res, keep = keep,
       vst = assay(vst(dds_full, blind = FALSE)),
       vst_blind = assay(vst(dds_full, blind = TRUE)),
       norm = counts(dds_full, normalized = TRUE))
}

fc  <- fit_all(counts_fc,  "featureCounts")
sal <- fit_all(counts_sal, "Salmon v15")

# ---------------------- per-group summaries ----------------------
group_summary <- function(norm) {
  g <- cd[colnames(norm), "group"]
  mu <- sapply(levels(cd$group), function(l) rowMeans(norm[, g == l, drop = FALSE]))
  sdv <- sapply(levels(cd$group), function(l) apply(norm[, g == l, drop = FALSE], 1, sd))
  colnames(mu) <- paste0("mean_", colnames(mu)); colnames(sdv) <- paste0("sd_", colnames(sdv))
  data.frame(gene_id = rownames(norm), mu, sdv, check.names = FALSE)
}

write_tables <- function(obj, tag) {
  gs <- group_summary(obj$norm)
  for (cn in names(obj$res)) {
    t <- obj$res[[cn]] %>% left_join(ann[, c("gene_id", "symbol", "biotype", "entrez")], by = "gene_id") %>%
      left_join(gs, by = "gene_id") %>% arrange(pvalue) %>% relocate(symbol, biotype, .after = gene_id)
    write_tsv(t, file.path(out, sprintf("deseq2_%s_%s.tsv", tag, cn)))
  }
  write_tsv(data.frame(gene_id = rownames(obj$vst), symbol = ann[rownames(obj$vst), "symbol"], obj$vst, check.names = FALSE),
            file.path(out, sprintf("vst_%s.tsv", tag)))
  write_tsv(data.frame(gene_id = rownames(obj$norm), symbol = ann[rownames(obj$norm), "symbol"], obj$norm, check.names = FALSE),
            file.path(out, sprintf("normalized_counts_%s.tsv", tag)))
}
write_tables(fc, "fc")
write_tables(sal, "sal")

sink(file.path(out, "summary_counts.txt"))
for (tag in c("fc", "sal")) {
  obj <- get(tag)
  cat(sprintf("== %s: %d genes tested\n", tag, sum(obj$keep)))
  for (cn in names(obj$res)) {
    r <- obj$res[[cn]]
    cat(sprintf("  %-13s padj<0.05: %5d   padj<0.1: %5d   padj(|LFC|>log2 1.5)<0.05: %5d\n", cn,
                sum(r$padj < 0.05, na.rm = TRUE), sum(r$padj < 0.1, na.rm = TRUE),
                sum(r$padj_LFC1.5 < 0.05, na.rm = TRUE)))
  }
}
sink()

saveRDS(list(cd = cd, ann = ann, counts_fc = counts_fc, counts_sal = counts_sal, vst_sal_v15 = vst_sal_v15,
             fc = fc, sal = sal, min_count = min_count, min_samps = min_samps),
        file.path(out, "core.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat(readLines(file.path(out, "summary_counts.txt")), sep = "\n")
