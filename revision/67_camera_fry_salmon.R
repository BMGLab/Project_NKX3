#!/usr/bin/env Rscript
# 67_camera_fry_salmon.R - correlation-aware Hallmark tests (limma fry; camera with the inter-gene correlation estimated
# from the data) on the Salmon VST (core$sal$vst, plasmid-flagged genes removed); same design and contrasts as
# section E of 50_verify_fig3_02_recompute.R. Output: results/revision_2026/salmon_v3/camera_fry_hallmark.tsv (+ _summary)
suppressPackageStartupMessages({library(DESeq2); library(limma); library(msigdbr); library(data.table)})
proj <- "/mnt/ssd1/Projects/NKX3_wd"; rv <- file.path(proj, "results/revision_2026")
core <- readRDS(file.path(rv, "core/core.rds")); cd <- core$cd; ann <- core$ann; S <- cd$sample_id
FLAG <- fread(file.path(rv, "salmon_v3/plasmid_flagged_genes.tsv"))$gene_id
V <- core$sal$vst[, S]; V <- V[!rownames(V) %in% FLAG, ]
H <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, .(term = gs_name, entrez = as.character(ncbi_gene))]
dY <- data.table(id = rownames(V), entrez = ann[rownames(V), "entrez"], m = rowMeans(V))[!is.na(entrez)][order(-m)][!duplicated(entrez)]
Y <- V[dY$id, S]
idx <- lapply(split(H$entrez, H$term), function(s) which(dY$entrez %in% s)); idx <- idx[lengths(idx) >= 15]
gfac <- factor(as.character(cd$group), levels = c("vector_ctrl", "vector_CM", "NKX3.1_ctrl", "NKX3.1_CM"))
DES <- model.matrix(~ 0 + gfac); colnames(DES) <- levels(gfac)
CON <- list(CM_in_vector = c(-1, 1, 0, 0), NKX_in_ctrl = c(-1, 0, 1, 0), NKX_in_CM = c(0, -1, 0, 1),
            CM_in_NKX = c(0, 0, -1, 1), CM_main = c(-.5, .5, -.5, .5), NKX_main = c(-.5, -.5, .5, .5))
all <- list(); summ <- list()
for (cn in names(CON)) {
  fr <- fry(Y, idx, DES, CON[[cn]]); ce <- camera(Y, idx, DES, CON[[cn]], inter.gene.cor = NA)
  all[[cn]] <- data.table(contrast = cn, set = rownames(fr), n_genes = fr$NGenes, fry_direction = fr$Direction, fry_p = fr$PValue,
                          fry_FDR = fr$FDR)[data.table(set = rownames(ce), camera_direction = ce$Direction, camera_p = ce$PValue,
                                                       camera_FDR = ce$FDR, camera_inter_gene_cor = ce$Correlation), on = "set"]
  summ[[cn]] <- data.table(contrast = cn, n_sets = length(idx), fry_FDR_lt_0.05 = sum(fr$FDR < 0.05), camera_FDR_lt_0.05 = sum(ce$FDR < 0.05),
                           min_fry_FDR = signif(min(fr$FDR), 2), min_camera_FDR = signif(min(ce$FDR), 2))
}
fwrite(rbindlist(all), file.path(rv, "salmon_v3/camera_fry_hallmark.tsv"), sep = "\t")
s <- rbindlist(summ); fwrite(s, file.path(rv, "salmon_v3/camera_fry_hallmark_summary.tsv"), sep = "\t"); print(s)
