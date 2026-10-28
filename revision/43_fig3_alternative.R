#!/usr/bin/env Rscript
# =============================================================================
# 43_fig3_alternative.R -- task "fig3_alternative"
#
# Question from a co-author: can the original Figure 3 heatmaps be kept, with
# their "modularity" supported by WGCNA etc.?  This script
#
# Part 1  audits the ORIGINAL Fig 3 (nf-core Salmon counts; rows = Welch-z
#         'significant' genes; boxed NKX3-1 "modules")
#   1a  reproduces the 200 rows of each original heatmap (A: CM, B: NKX3.1 basal,
#       C: NKX3.1 under CM) with the original procedure (rna_downstreams_salmon.R,
#       as re-implemented in 14_modules.R) and checks them against the published PDFs
#   1b  status of every row in the corrected analysis (unique reads, DESeq2)
#   1c  applies the original row filter to the nine "nonsense" groupings of the
#       same six libraries (one library swapped between the two groups) and draws
#       the original display for the real and for one nonsense grouping
# Part 2  new angle: WGCNA modules versus the pathways that GSEA found
#   2a  round-1 network (corrected VST, 5,000 most variable genes, signed bicor;
#       14_modules.R): module x Hallmark over-representation (one-sided Fisher,
#       BH over modules x 50 sets) and eigengene x design association (exact
#       permutation over all 924 balanced labelings of the 12 libraries, BH over
#       modules x 3 terms)
#   2b  sensitivity network, pre-specified because the 5,000-gene network holds only
#       ~25-35% of the genes of each focal Hallmark set: identical settings on all
#       expressed genes (mean normalised count >= 10)
#   2c  a module is a supportable "pathway module" only if it is over-represented
#       for a focal set (BH < 0.05) AND its eigengene is associated with CM or NKX3.1
#       (BH-adjusted exact permutation P < 0.05); reported whatever the outcome
# Part 3  candidate replacement Fig 3: leading-edge heatmaps of the corrected GSEA
#   A  CM vs control (vector cells): HALLMARK_MYC_TARGETS_V1 + _OXIDATIVE_PHOSPHORYLATION
#   B  NKX3.1 vs control (basal):    HALLMARK_E2F_TARGETS + _G2M_CHECKPOINT
#   C  NKX3.1 vs control (CM):       HALLMARK_APICAL_JUNCTION
#   rows = leading-edge (core_enrichment) genes, top ~30 by |Wald statistic|
#   (2-set panels: top k of each set, k raised until the union reaches 30);
#   corrected VST (core$fc$vst), row z-scores, row annotation = log2FC (MLE).
#   Variants: all 12 libraries; unscaled (VST difference from the control mean).
#   Replicate consistency, with a selection-matched null: the same top-k rule
#   applied to data with no group effect (random rotations of the residuals of the
#   12-library group-means model; gene-gene correlation preserved).
# Part 4  correlation-aware gene-set tests of the same Hallmark sets (limma camera with
#   fixed or estimated inter-gene correlation; fry/roast rotation tests; voom input;
#   exploratory replicate-factor block; leave-one-library-out), because preranked GSEA
#   permutes genes and so treats the genes of a set as independent.
#
# Sign convention everywhere: treatment minus reference (CM - control; NKX3.1 - vector).
# Inputs are read-only (core.rds, round-1 modules/gsea/markers_qc outputs, v15 PDFs).
# Run from a clean shell:
#   /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/43_fig3_alternative.R
# =============================================================================

# plotmath minus signs / Unicode need a UTF-8 character locale (a clean env -i shell starts in C)
if (!isTRUE(l10n_info()[["UTF-8"]])) for (loc in c("en_US.UTF-8", "C.UTF-8"))
  if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break

suppressPackageStartupMessages({
  library(DESeq2)            # classes stored in core.rds
  library(data.table)
  library(matrixStats)
  library(msigdbr)
  library(WGCNA)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
options(stringsAsFactors = FALSE, warn = 1)
invisible(allowWGCNAThreads(16))
ht_opt$message <- FALSE
t_start <- Sys.time()

proj <- "/mnt/ssd1/Projects/NKX3_wd"
revd <- file.path(proj, "results/revision_2026")
out  <- file.path(revd, "fig3_alternative")
work <- "/mnt/ssd1/Projects/NKX3_wd/work/fig3_alternative"
v15  <- file.path(proj, "results/rna_downstreams_out_v15")
gsea_tab <- file.path(revd, "gsea/tables")
for (d in c(out, work)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEED   <- 42L
N_ROT  <- 2000L                     # residual rotations for the selection-matched null
N_ROWS <- 30L                       # rows per leading-edge heatmap
NKX_ID <- "ENSG00000167034"         # NKX3-1
FONT   <- "Liberation Sans"         # Arial-metric; Arial is not installed
GRP_COL <- c("Control-basal" = "#909FC9", "Control-inflammation" = "#E28DBB",
             "NKX3.1-basal" = "#57C5D1", "NKX3.1-inflammation" = "#F5906A")
PAL <- c("#2166AC", "white", "#B2182B")       # round-1 blue-white-red
set.seed(SEED)

log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
wtsv <- function(x, f) fwrite(as.data.table(x), file.path(out, f), sep = "\t", quote = FALSE, na = "NA")
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, formatC(p, format = "e", digits = 1),
                                                    formatC(signif(p, 2), format = "fg", digits = 2)))
save_pdf_png <- function(draw_fun, name, w, h) {
  cairo_pdf(file.path(out, paste0(name, ".pdf")), width = w, height = h, family = FONT)
  draw_fun(); invisible(dev.off())
  png(file.path(out, paste0(name, ".png")), width = w, height = h, units = "in", res = 300, type = "cairo", family = FONT)
  draw_fun(); invisible(dev.off())
}

# =============================================================================
# 0. Data
# =============================================================================
log_msg("reading core.rds and round-1 outputs")
core <- readRDS(file.path(revd, "core/core.rds"))
r1   <- readRDS(file.path(revd, "modules/modules_results.rds"))
cd   <- core$cd; ann <- core$ann
S    <- cd$sample_id
cd$short <- paste0("#", cd$num)
grp_lab <- setNames(cd$label, S)
V    <- core$fc$vst[, S]                       # corrected VST (blind = FALSE), 18,704 genes
NORM <- core$fc$norm[, S]
RES  <- core$fc$res
sym_of <- function(id) { s <- ann[id, "symbol"]; ifelse(is.na(s) | s == "", id, s) }

qc <- fread(file.path(revd, "markers_qc/qc_sample_metrics.tsv"))[match(S, sample_id)]
stopifnot(identical(qc$sample_id, S))
batchB <- setNames(qc$hidden_factor_split == "B", S)           # replicate-level factor (markers_qc T4)
fcsum  <- as.data.table(r1$fc_sum)[match(S, sample_id)]
pct_out <- setNames(fcsum$no_feature_pct, S)                    # % reads outside annotated exons (round 1)
log_msg("replicate-level factor (markers_qc): %s", paste(cd$short[batchB], collapse = ", "))

# ranked lists exactly as 12_gsea.R: Wald statistic, Ensembl -> Entrez, duplicated Entrez -> max |stat|
rank_table <- function(cn) {
  r <- as.data.table(RES[[cn]])[, .(gene_id, stat, log2FC_MLE, lfcSE_MLE, log2FC_shrunk, pvalue, padj, baseMean)]
  r[, entrez := ann[gene_id, "entrez"]]
  r <- r[!is.na(entrez) & is.finite(stat)][order(-abs(stat))][!duplicated(entrez)]
  r[, rank_in_list := frank(-stat, ties.method = "first")]
  r
}

# Hallmark gene sets (same source as 12_gsea.R)
H_tab <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
MSIG_VERSION <- paste(unique(H_tab$db_version), collapse = ",")
H <- unique(H_tab[!is.na(ncbi_gene) & ncbi_gene != "", .(term = gs_name, entrez = as.character(ncbi_gene))])
log_msg("msigdbr %s (MSigDB %s): %d Hallmark sets", as.character(packageVersion("msigdbr")), MSIG_VERSION, uniqueN(H$term))
FOCAL <- c(MYC = "HALLMARK_MYC_TARGETS_V1", OXPHOS = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           E2F = "HALLMARK_E2F_TARGETS", G2M = "HALLMARK_G2M_CHECKPOINT", AJ = "HALLMARK_APICAL_JUNCTION")
SET_LAB <- c(HALLMARK_MYC_TARGETS_V1 = "MYC targets V1", HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OXPHOS",
             HALLMARK_E2F_TARGETS = "E2F targets", HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
             HALLMARK_APICAL_JUNCTION = "Apical junction")
read_gsea <- function(q, cn) fread(file.path(gsea_tab, sprintf("gsea_%s_hallmark_%s.tsv", q, cn)),
                                   colClasses = list(character = c("core_enrichment")))

# =============================================================================
# Part 1. Original Fig 3 audit
# =============================================================================
log_msg("Part 1: original Fig 3 audit")
st   <- read.csv(file.path(proj, "results/sample_table.csv"))
V0   <- core$vst_sal_v15[, st$sample_id]                         # published blind VST (Salmon), original column order
ids0 <- rownames(V0); s0 <- ann[ids0, "symbol"]
rn0  <- ifelse(!is.na(s0) & s0 != "", s0, ids0)
cond0 <- st$condition
welch <- function(M, i1, i2) {                                   # = run_comparison() Welch block
  mu1 <- rowMeans(M[, i1, drop = FALSE]); mu2 <- rowMeans(M[, i2, drop = FALSE])
  v1 <- rowVars(M[, i1, drop = FALSE]); v2 <- rowVars(M[, i2, drop = FALSE])
  lfc <- mu1 - mu2
  z <- lfc / sqrt(pmax(v1, 1e-6) / length(i1) + pmax(v2, 1e-6) / length(i2))
  list(lfc = lfc, z = z, sig = pmax(mu1, mu2) >= 2.0 & abs(lfc) >= 0.25 & abs(z) >= 2.0)
}
orig_rows <- function(M, rn, i1, i2) {                          # = PLOT 3 of run_comparison(): top 200 by |z|
  w <- welch(M, i1, i2)
  o <- order(-abs(w$z), -abs(w$lfc)); so <- o[w$sig[o]]
  keep <- unique(rn[so]); keep <- keep[seq_len(min(200, length(keep)))]
  rows <- which(rn %in% keep)
  mat <- t(scale(t(M[rows, c(i1, i2)]))); rownames(mat) <- make.unique(rn[rows])
  hc <- hclust(dist(mat), method = "complete")                  # pheatmap default
  list(rows = rows, mat = mat, hc = hc, order = rownames(mat)[hc$order], n_sig = sum(w$sig), w = w)
}
pdf_row_labels <- function(pdf) {                                # published row labels, top to bottom
  tmp <- file.path(work, paste0(basename(pdf), ".bbox.html"))
  system2("pdftotext", c("-bbox", shQuote(pdf), shQuote(tmp)))
  x <- grep("<word ", readLines(tmp, warn = FALSE, encoding = "UTF-8"), value = TRUE)
  xmin <- sub('.*xMin="([0-9.]+)".*', "\\1", x)
  ymin <- as.numeric(sub('.*yMin="([0-9.]+)".*', "\\1", x))
  w <- gsub("−", "-", sub('.*">(.*)</word>.*', "\\1", x))
  col <- names(which.max(table(xmin)))
  sel <- xmin == col & w != "condition"
  unlink(tmp)
  w[sel][order(ymin[sel])]
}
ORIG <- list(
  A = list(g1 = "HM-vector-C",  g2 = "HM-vector-CM", sign = -1, contrast = "CM_in_vector", module = NA_character_,
           label = "Fig 3A: CM vs control (vector cells)", trt = "vector_CM", ref = "vector_ctrl"),
  B = list(g1 = "HM-NKX3.1-C",  g2 = "HM-vector-C",  sign = 1, contrast = "NKX_in_ctrl", module = "basal",
           label = "Fig 3B/D: NKX3.1 vs control (basal)", trt = "NKX3.1_ctrl", ref = "vector_ctrl"),
  C = list(g1 = "HM-NKX3.1-CM", g2 = "HM-vector-CM", sign = 1, contrast = "NKX_in_CM", module = "inflammation",
           label = "Fig 3C/E: NKX3.1 vs control (CM)", trt = "NKX3.1_CM", ref = "vector_CM"))
ORIG_HM <- list(); row_status <- list(); orig_summary <- list()
for (k in names(ORIG)) {
  o  <- ORIG[[k]]
  i1 <- which(cond0 == o$g1); i2 <- which(cond0 == o$g2)
  h  <- orig_rows(V0, rn0, i1, i2)
  pdf <- file.path(v15, sprintf("cmp_%s_vs_%s_heatmap_sig_genes.pdf", o$g1, o$g2))
  pl <- pdf_row_labels(pdf)
  exact <- identical(h$order, pl)
  lab2id <- setNames(ids0[h$rows], rownames(h$mat))
  ids <- unname(lab2id[h$order])
  rr <- as.data.frame(RES[[o$contrast]]); rownames(rr) <- rr$gene_id
  tested <- ids %in% rownames(rr)
  boxed <- if (is.na(o$module)) character(0) else unname(r1$modules[[o$module]])
  stopifnot(all(ids %in% rownames(core$counts_fc)))
  cnt <- core$counts_fc[ids, S, drop = FALSE]
  d <- data.table(panel = k, figure = o$label, heatmap_row = seq_along(ids), symbol = h$order, gene_id = ids,
                  in_boxed_module = ids %in% boxed,
                  orig_welch_diff_vst_treat_minus_ref = round(o$sign * h$w$lfc[match(ids, ids0)], 3),
                  orig_welch_z_treat_minus_ref = round(o$sign * h$w$z[match(ids, ids0)], 2),
                  unique_reads_total = rowSums(cnt), unique_reads_max_per_library = rowMaxs(cnt),
                  corrected_tested = tested,
                  corrected_baseMean = round(rr[ids, "baseMean"], 1),
                  corrected_log2FC_MLE = round(rr[ids, "log2FC_MLE"], 3),
                  corrected_log2FC_shrunk = signif(rr[ids, "log2FC_shrunk"], 3),
                  corrected_stat = round(rr[ids, "stat"], 2),
                  corrected_padj = signif(rr[ids, "padj"], 3))
  d[, same_direction := ifelse(corrected_tested, sign(corrected_log2FC_MLE) == sign(orig_welch_diff_vst_treat_minus_ref), NA)]
  row_status[[k]] <- d
  ORIG_HM[[k]] <- c(h, list(i1 = i1, i2 = i2, ids = ids, exact = exact))
  sigc <- d[corrected_padj < 0.05]
  orig_summary[[k]] <- data.table(
    panel = k, figure = o$label, corrected_contrast = o$contrast, rows = nrow(d),
    reproduced_row_order_equals_published_pdf = exact, welch_filter_genes_passing = h$n_sig,
    rows_in_boxed_module = sum(d$in_boxed_module),
    rows_not_testable_from_unique_reads = sum(!d$corrected_tested),
    rows_max_5_unique_reads_per_library = sum(d$unique_reads_max_per_library <= 5),
    rows_corrected_padj_lt_0.05 = nrow(sigc),
    which_rows_corrected_padj_lt_0.05 = paste(sigc$symbol, collapse = ","),
    rows_tested_same_direction = sum(d$same_direction, na.rm = TRUE),
    rows_tested_opposite_direction = sum(!d$same_direction, na.rm = TRUE),
    median_abs_orig_diff = round(median(abs(d$orig_welch_diff_vst_treat_minus_ref)), 3),
    median_abs_corrected_log2FC_MLE = round(median(abs(d$corrected_log2FC_MLE), na.rm = TRUE), 3),
    median_abs_corrected_log2FC_MLE_excl_NKX3_1 = round(median(abs(d[gene_id != NKX_ID]$corrected_log2FC_MLE), na.rm = TRUE), 3),
    boxed_rows_not_testable = sum(d$in_boxed_module & !d$corrected_tested),
    boxed_rows_padj_lt_0.05 = sum(d$in_boxed_module & d$corrected_padj < 0.05, na.rm = TRUE))
  log_msg("  %s: %d rows (published order reproduced: %s); not testable %d; corrected padj<0.05: %s",
          k, nrow(d), exact, sum(!d$corrected_tested), paste(sigc$symbol, collapse = ",") )
}
row_status <- rbindlist(row_status); orig_summary <- rbindlist(orig_summary)

# ---- 1c. the original filter on the nine nonsense groupings of the same six libraries
split_stats <- function(M, idx6, label) {
  sp <- combn(6, 3); sp <- sp[, sp[1, ] == 1, drop = FALSE]      # 10 distinct 3+3 splits; first = real
  rbindlist(lapply(seq_len(ncol(sp)), function(j) {
    a <- idx6[sp[, j]]; b <- setdiff(idx6, a)
    w <- welch(M, a, b)
    zs <- sort(abs(w$z[w$sig]), decreasing = TRUE)
    data.table(quantification = label, split = j, real_grouping = j == 1,
               group1 = paste(a, collapse = ","), group2 = paste(b, collapse = ","),
               libraries_swapped = as.integer(min(sum(a %in% idx6[1:3]), 3 - sum(a %in% idx6[1:3]))),
               genes_passing_filter = sum(w$sig), genes_abs_z_ge_2 = sum(abs(w$z) >= 2),
               z_of_200th_row = if (length(zs) >= 200) round(zs[200], 2) else NA_real_,
               rows_with_abs_z_ge_5_among_top200 = sum(head(zs, 200) >= 5))
  }))
}
nonsense <- list()
for (k in names(ORIG)) {
  o <- ORIG[[k]]
  idx_o <- c(which(cond0 == o$g1), which(cond0 == o$g2))                      # positions in V0 (st order)
  idx_c <- match(st$sample_id[idx_o], colnames(V))                            # same libraries in V (cd order)
  s1 <- split_stats(V0, idx_o, "original (Salmon VST, published)")
  s2 <- split_stats(V, idx_c, "corrected (featureCounts VST)")
  lib_lab <- function(x) vapply(strsplit(x, ","), function(v) paste(cd$short[match(st$sample_id[as.integer(v)], S)], collapse = " "), "")
  s1[, `:=`(group1 = lib_lab(group1), group2 = lib_lab(group2))]
  lib_lab2 <- function(x) vapply(strsplit(x, ","), function(v) paste(cd$short[as.integer(v)], collapse = " "), "")
  s2[, `:=`(group1 = lib_lab2(group1), group2 = lib_lab2(group2))]
  nonsense[[k]] <- cbind(panel = k, figure = o$label, rbind(s1, s2))
}
nonsense <- rbindlist(nonsense)
nonsense_summary <- nonsense[, .(real_genes_passing = genes_passing_filter[real_grouping],
                                 nonsense_median_genes_passing = median(genes_passing_filter[!real_grouping]),
                                 nonsense_min_genes_passing = min(genes_passing_filter[!real_grouping]),
                                 nonsense_max_genes_passing = max(genes_passing_filter[!real_grouping]),
                                 real_z_200th = z_of_200th_row[real_grouping],
                                 nonsense_median_z_200th = median(z_of_200th_row[!real_grouping], na.rm = TRUE),
                                 nonsense_groupings_with_200_rows = sum(!real_grouping & genes_passing_filter >= 200)),
                             by = .(panel, figure, quantification)]
print(nonsense_summary)

# =============================================================================
# Part 2. WGCNA modules vs the GSEA pathways
# =============================================================================
log_msg("Part 2: WGCNA module x Hallmark over-representation and eigengene x design")
traits <- data.frame(NKX = ifelse(cd$NKX == "NKX3.1", 1, -1), CM = ifelse(cd$CM == "CM", 1, -1))
traits$NKXxCM <- traits$NKX * traits$CM
LAB <- combn(12, 6)                                               # all 924 balanced labelings
TALL <- apply(LAB, 2, function(L) { t <- rep(-1, 12); t[L] <- 1; t })   # 12 x 924

module_trait <- function(MEs, label, sizes) {
  MEs <- as.matrix(MEs[S, , drop = FALSE])
  rn <- cor(MEs, TALL)                                            # modules x 924
  rows <- lapply(colnames(traits), function(tr) {
    r_obs <- cor(MEs, traits[[tr]])[, 1]
    data.table(network = label, module = sub("^ME", "", colnames(MEs)), trait = tr, r = r_obs,
               p_exact = rowMeans(abs(rn) >= abs(r_obs) - 1e-12))
  })
  mt <- rbindlist(rows)
  mt[, padj_exact := p.adjust(p_exact, "BH")]
  gm <- t(apply(MEs, 2, function(x) tapply(x, cd$group, mean)))
  extra <- data.table(module = sub("^ME", "", colnames(MEs)), size = as.integer(sizes[sub("^ME", "", colnames(MEs))]),
                      r_pct_reads_outside_exons = cor(MEs, pct_out)[, 1],
                      r_replicate_factor_3_6_7_10 = cor(MEs, as.numeric(batchB))[, 1],
                      group_R2 = apply(MEs, 2, function(x) summary(lm(x ~ cd$group))$r.squared),
                      ME_mean_Control_basal = gm[, "vector_ctrl"], ME_mean_Control_inflammation = gm[, "vector_CM"],
                      ME_mean_NKX3.1_basal = gm[, "NKX3.1_ctrl"], ME_mean_NKX3.1_inflammation = gm[, "NKX3.1_CM"])
  merge(mt, extra, by = "module", sort = FALSE)
}
module_enrichment <- function(colors, label) {
  U <- data.table(gene_id = names(colors), module = unname(colors))
  U[, entrez := ann[gene_id, "entrez"]]
  U <- U[!is.na(entrez)][!duplicated(entrez)]
  mods <- setdiff(sort(unique(U$module)), "grey")
  sets <- sort(unique(H$term))
  inS_all <- sapply(sets, function(s) U$entrez %in% H[term == s, entrez])
  res <- rbindlist(lapply(mods, function(m) {
    inM <- U$module == m
    rbindlist(lapply(sets, function(s) {
      inS <- inS_all[, s]
      a <- sum(inM & inS); b <- sum(inM & !inS); c <- sum(!inM & inS); d <- sum(!inM & !inS)
      ft <- fisher.test(matrix(c(a, c, b, d), 2), alternative = "greater")
      data.table(network = label, module = m, set = s, overlap = a, module_genes_with_entrez = a + b,
                 set_genes_in_network = a + c, universe = nrow(U), expected = round((a + b) * (a + c) / nrow(U), 2),
                 odds_ratio = round(unname(ft$estimate), 2), p = ft$p.value,
                 overlap_genes = paste(sort(sym_of(U$gene_id[inM & inS])), collapse = ","))
    }))
  }))
  res[, padj := p.adjust(p, "BH")]
  res[]
}
set_coverage <- function(colors, label) {
  U <- data.table(gene_id = names(colors), module = unname(colors))
  U[, entrez := ann[gene_id, "entrez"]]
  expr_entrez <- unique(na.omit(ann[rownames(V), "entrez"]))
  rbindlist(lapply(FOCAL, function(s) {
    se <- H[term == s, entrez]; inn <- U[entrez %in% se]
    tb <- sort(table(inn$module), decreasing = TRUE)
    data.table(network = label, set = s, set_size = length(se), expressed = sum(se %in% expr_entrez),
               in_network = nrow(inn), pct_of_expressed_in_network = round(100 * nrow(inn) / sum(se %in% expr_entrez), 1),
               module_distribution = paste(names(tb), tb, sep = ":", collapse = " "))
  }))
}

# ---- 2a. round-1 network
w1 <- r1$wg$corrected
net1 <- list(label = "round-1: 5,000 most variable genes", short = "round1_top5000",
             colors = setNames(w1$assign$module, w1$assign$gene_id), MEs = w1$MEs, power = w1$beta,
             rule = w1$rule, sizes = table(w1$assign$module), n_genes = nrow(w1$assign))

# ---- 2b. sensitivity network on all expressed genes
expressed <- rownames(NORM)[rowMeans(NORM) >= 10]
log_msg("  sensitivity WGCNA on %d expressed genes (mean normalised count >= 10)", length(expressed))
run_wgcna <- function(genes, label, short) {
  datExpr <- t(V[genes, S])
  powers <- c(1:10, seq(12, 30, by = 2))
  invisible(capture.output(sft <- suppressWarnings(pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                           corFnc = bicor, corOptions = list(maxPOutliers = 0.05), blockSize = 6000, verbose = 0))))
  fi <- sft$fitIndices; fi$signed_R2 <- -sign(fi$slope) * fi$SFT.R.sq
  beta <- fi$Power[which(fi$signed_R2 >= 0.80)[1]]
  rule <- "lowest power with signed scale-free R^2 >= 0.80 (as round 1)"
  if (is.na(beta)) { beta <- 18; rule <- "no power reached signed R^2 >= 0.80 -> WGCNA FAQ value for signed networks with <20 samples (18)" }
  assign("cor", WGCNA::cor, envir = globalenv())                  # as 14_modules.R: avoid stats::cor masking
  on.exit(rm("cor", envir = globalenv()))
  net <- suppressWarnings(blockwiseModules(datExpr, power = beta, networkType = "signed", TOMType = "signed",
                          corType = "bicor", maxPOutliers = 0.05, minModuleSize = 30, mergeCutHeight = 0.25,
                          deepSplit = 2, maxBlockSize = ncol(datExpr) + 1, numericLabels = FALSE, saveTOMs = FALSE,
                          randomSeed = SEED, nThreads = 16, verbose = 0))
  colors <- setNames(net$colors, colnames(datExpr))
  MEs <- orderMEs(net$MEs); MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
  kme <- suppressWarnings(bicor(datExpr, as.matrix(MEs), maxPOutliers = 0.05))
  own <- paste0("ME", colors)
  kme_own <- ifelse(own %in% colnames(kme), kme[cbind(seq_along(colors), match(own, colnames(kme)))], NA)
  list(label = label, short = short, colors = colors, MEs = MEs, power = beta, rule = rule, sizes = table(colors),
       n_genes = length(genes), sft = data.table(network = label, fi, chosen = fi$Power == beta),
       assign = data.table(network = label, gene_id = names(colors), symbol = sym_of(names(colors)),
                           entrez = ann[names(colors), "entrez"], module = colors, kME_own_module = round(kme_own, 3)))
}
# FIG3ALT_DEV_CACHE=1 reuses a cached sensitivity network during development only; the default
# (and the recorded run) recomputes it from scratch.
cache_f <- file.path(work, "dev_cache_net2.rds")
if (Sys.getenv("FIG3ALT_DEV_CACHE") == "1" && file.exists(cache_f)) {
  net2 <- readRDS(cache_f); log_msg("  DEV: sensitivity network read from cache")
} else {
  net2 <- run_wgcna(expressed, sprintf("sensitivity: all %s expressed genes", format(length(expressed), big.mark = ",")),
                    "sensitivity_all_expressed")
  if (Sys.getenv("FIG3ALT_DEV_CACHE") == "1") saveRDS(net2, cache_f)
}
log_msg("  sensitivity network: power %d (%s); %d modules (+ grey %d genes)", net2$power, net2$rule,
        length(setdiff(names(net2$sizes), "grey")), as.integer(net2$sizes["grey"]))

NETS <- list(net1, net2)
enr  <- rbindlist(lapply(NETS, function(n) module_enrichment(n$colors, n$label)))
mtr  <- rbindlist(lapply(NETS, function(n) module_trait(n$MEs, n$label, n$sizes)))
cov  <- rbindlist(lapply(NETS, function(n) set_coverage(n$colors, n$label)))
# round-1 r values must be reproduced exactly from the saved eigengenes
chk <- merge(mtr[network == net1$label, .(module, trait, r)],
             as.data.table(w1$module_trait)[, .(module, trait, r_round1 = r)], by = c("module", "trait"))
stopifnot(max(abs(chk$r - chk$r_round1)) < 1e-3)

# verdict per module: focal over-representation AND design association
focal_best <- enr[set %in% FOCAL, .SD[which.min(padj)], by = .(network, module)][
  , .(network, module, best_focal_set = set, best_focal_overlap = overlap, best_focal_set_in_network = set_genes_in_network,
      best_focal_OR = odds_ratio, best_focal_padj = padj)]
focal_sig <- enr[set %in% FOCAL & padj < 0.05, .(focal_sets_padj_lt_0.05 = paste(sprintf("%s (%d/%d, OR %.1f, padj %s)",
  SET_LAB[set], overlap, set_genes_in_network, odds_ratio, fmt_p(padj)), collapse = "; ")), by = .(network, module)]
wide_tr <- dcast(mtr, network + module + size + r_pct_reads_outside_exons + r_replicate_factor_3_6_7_10 + group_R2 ~ trait,
                 value.var = c("r", "p_exact", "padj_exact"))
verdict <- merge(merge(wide_tr, focal_best, by = c("network", "module"), all.x = TRUE),
                 focal_sig, by = c("network", "module"), all.x = TRUE)
verdict[is.na(focal_sets_padj_lt_0.05), focal_sets_padj_lt_0.05 := ""]
verdict[, focal_enriched := best_focal_padj < 0.05]
verdict[, design_associated := pmin(padj_exact_CM, padj_exact_NKX) < 0.05]
verdict[, qualifies_as_pathway_module := focal_enriched & design_associated]
verdict[, nominal_design_p_lt_0.05 := pmin(p_exact_CM, p_exact_NKX) < 0.05]
# descriptive: does a focal-enriched module's eigengene follow CM/NKX3.1 once the library-level covariates are in the model?
adj_fit <- function(netlab, mod) {
  n <- NETS[[which(vapply(NETS, `[[`, "", "label") == netlab)]]
  me <- n$MEs[S, paste0("ME", mod)]
  f <- summary(lm(me ~ traits$NKX + traits$CM + pct_out + as.numeric(batchB)))$coefficients
  data.table(network = netlab, module = mod, adj_t_NKX = round(f[2, 3], 2), adj_p_NKX = signif(f[2, 4], 2),
             adj_t_CM = round(f[3, 3], 2), adj_p_CM = signif(f[3, 4], 2),
             adj_t_pct_outside = round(f[4, 3], 2), adj_t_replicate_factor = round(f[5, 3], 2),
             adj_p_replicate_factor = signif(f[5, 4], 2))
}
fe <- verdict[focal_enriched == TRUE]
if (nrow(fe)) {
  adj <- rbindlist(lapply(seq_len(nrow(fe)), function(i) adj_fit(fe$network[i], fe$module[i])))
  verdict <- merge(verdict, adj, by = c("network", "module"), all.x = TRUE)
}
setorder(verdict, network, best_focal_padj)
print(verdict[focal_enriched == TRUE, .(network, module, size, focal_sets_padj_lt_0.05, r_CM, padj_exact_CM, r_NKX, padj_exact_NKX,
                                        r_replicate_factor_3_6_7_10, r_pct_reads_outside_exons, qualifies_as_pathway_module)])
log_msg("  modules qualifying as pathway modules (focal BH<0.05 AND CM/NKX3.1 BH<0.05): %d",
        sum(verdict$qualifies_as_pathway_module, na.rm = TRUE))

# =============================================================================
# Part 3. Leading-edge heatmaps (candidate replacement Fig 3)
# =============================================================================
log_msg("Part 3: leading-edge heatmaps")
PANELS <- list(
  A = list(contrast = "CM_in_vector", sets = FOCAL[c("MYC", "OXPHOS")], ref = "vector_ctrl", trt = "vector_CM",
           other = "CM_in_NKX", other_label = "CM vs control in NKX3.1 cells", title = "CM vs control (vector cells)"),
  B = list(contrast = "NKX_in_ctrl", sets = FOCAL[c("E2F", "G2M")], ref = "vector_ctrl", trt = "NKX3.1_ctrl",
           other = "NKX_in_CM", other_label = "NKX3.1 vs control under CM", title = "NKX3.1 vs control (basal)"),
  C = list(contrast = "NKX_in_CM", sets = FOCAL[c("AJ")], ref = "vector_CM", trt = "NKX3.1_CM",
           other = "NKX_in_ctrl", other_label = "NKX3.1 vs control, basal", title = "NKX3.1 vs control (CM)"))
G_FC  <- lapply(setNames(nm = unique(c(sapply(PANELS, `[[`, "contrast"), sapply(PANELS, `[[`, "other")))), function(cn) read_gsea("fc", cn))
G_ADJ <- lapply(setNames(nm = sapply(PANELS, `[[`, "contrast")), function(cn) read_gsea("fcadj", cn))
G_SAL <- lapply(setNames(nm = sapply(PANELS, `[[`, "contrast")), function(cn) read_gsea("sal", cn))
v15pub <- function(cn, set) {
  f <- file.path(gsea_tab, sprintf("gsea_v15published_hallmark_%s.tsv", cn))
  if (!file.exists(f)) return(c(NA, NA))
  p <- fread(f); p <- p[ID == set]
  if (!nrow(p)) return(c(NA, NA)) else c(p$NES, p$p.adjust)
}
le_of <- function(g, set) { x <- g[ID == set]; if (!nrow(x)) character(0) else strsplit(x$core_enrichment, "/", fixed = TRUE)[[1]] }
select_top <- function(le_ids_by_set, n = N_ROWS) {              # le_ids_by_set: gene ids ordered by |stat| desc
  if (length(le_ids_by_set) == 1) return(list(ids = head(le_ids_by_set[[1]], n), k = n))
  k <- ceiling(n / length(le_ids_by_set))
  repeat {
    sel <- unique(unlist(lapply(le_ids_by_set, head, k)))
    if (length(sel) >= n || k >= max(lengths(le_ids_by_set))) break
    k <- k + 1
  }
  list(ids = sel, k = k)
}
w1_mod <- setNames(net1$colors, names(net1$colors)); w2_mod <- setNames(net2$colors, names(net2$colors))
PAN <- list(); le_tables <- list(); set_tables <- list()
for (pn in names(PANELS)) {
  p  <- PANELS[[pn]]
  rk <- rank_table(p$contrast)
  ro <- as.data.table(RES[[p$other]])[, .(gene_id, other_log2FC_MLE = log2FC_MLE, other_stat = stat, other_padj = padj)]
  le_by_set <- list(); set_rows <- list()
  for (s in p$sets) {
    g <- G_FC[[p$contrast]][ID == s]
    le <- le_of(G_FC[[p$contrast]], s)
    m  <- rk[entrez %in% le][order(-abs(stat))]
    if (nrow(m) != length(le)) stop("leading-edge Entrez IDs missing from the ranking for ", s)
    le_by_set[[s]] <- m$gene_id
    pub <- v15pub(p$contrast, s)
    go <- G_FC[[p$other]][ID == s]
    le_adj <- le_of(G_ADJ[[p$contrast]], s); le_sal <- le_of(G_SAL[[p$contrast]], s)
    set_rows[[s]] <- data.table(panel = pn, contrast = p$contrast, set = s, set_label = SET_LAB[[s]], set_size_tested = g$setSize,
      NES = round(g$NES, 2), FDR = signif(g$p.adjust, 2), leading_edge_size = length(le),
      NES_trend_adjusted = round(G_ADJ[[p$contrast]][ID == s]$NES, 2), FDR_trend_adjusted = signif(G_ADJ[[p$contrast]][ID == s]$p.adjust, 2),
      NES_salmon_refit = round(G_SAL[[p$contrast]][ID == s]$NES, 2), FDR_salmon_refit = signif(G_SAL[[p$contrast]][ID == s]$p.adjust, 2),
      NES_published_sign_corrected = round(pub[1], 2), FDR_published = signif(pub[2], 2),
      other_context = p$other_label, NES_other_context = round(go$NES, 2), FDR_other_context = signif(go$p.adjust, 2),
      le_size_trend_adjusted = length(le_adj), le_size_salmon = length(le_sal))
  }
  sel <- select_top(le_by_set)
  all_le <- unique(unlist(le_by_set))
  tab <- rk[gene_id %in% all_le]
  tab[, symbol := sym_of(gene_id)]
  for (s in p$sets) {
    tab[, (paste0("LE_", names(p$sets)[p$sets == s])) := gene_id %in% le_by_set[[s]]]
    tab[, (paste0("LEtrendadj_", names(p$sets)[p$sets == s])) := entrez %in% le_of(G_ADJ[[p$contrast]], s)]
    tab[, (paste0("LEsalmon_", names(p$sets)[p$sets == s])) := entrez %in% le_of(G_SAL[[p$contrast]], s)]
  }
  tab[, displayed := gene_id %in% sel$ids]
  tab <- merge(tab, ro, by = "gene_id", all.x = TRUE)
  gm <- t(apply(NORM[tab$gene_id, , drop = FALSE], 1, function(x) tapply(x, cd$group, mean)))
  tab[, `:=`(mean_norm_Control_basal = round(gm[, "vector_ctrl"], 1), mean_norm_Control_inflammation = round(gm[, "vector_CM"], 1),
             mean_norm_NKX3.1_basal = round(gm[, "NKX3.1_ctrl"], 1), mean_norm_NKX3.1_inflammation = round(gm[, "NKX3.1_CM"], 1))]
  tab[, `:=`(wgcna_module_round1 = ifelse(gene_id %in% names(w1_mod), w1_mod[gene_id], "not in network"),
             wgcna_module_sensitivity = ifelse(gene_id %in% names(w2_mod), w2_mod[gene_id], "not in network"))]
  # membership category and display order: set1-only, both, set2-only; within category by Wald statistic (most negative first)
  if (length(p$sets) == 2) {
    n1 <- names(p$sets)[1]; n2 <- names(p$sets)[2]
    tab[, membership := fifelse(get(paste0("LE_", n1)) & get(paste0("LE_", n2)), "both",
                                fifelse(get(paste0("LE_", n1)), paste(n1, "only"), paste(n2, "only")))]
    tab[, mem_order := match(membership, c(paste(n1, "only"), "both", paste(n2, "only")))]
  } else { tab[, membership := names(p$sets)[1]]; tab[, mem_order := 1L] }
  setorder(tab, -displayed, mem_order, stat)
  tab[, display_row := ifelse(displayed, seq_len(.N), NA_integer_)]
  le_tables[[pn]] <- cbind(panel = pn, contrast = p$contrast, k_per_set = sel$k, tab)
  set_tables[[pn]] <- rbindlist(set_rows)
  PAN[[pn]] <- c(p, list(rank = rk, le_by_set = le_by_set, k = sel$k, genes = tab[displayed == TRUE]$gene_id,
                         tab = tab[displayed == TRUE], sets_tab = rbindlist(set_rows),
                         ref_s = S[cd$group == p$ref], trt_s = S[cd$group == p$trt]))
  log_msg("  panel %s: %d displayed genes (k = %d per set); median log2FC %.3f; max gene-level padj %s, min %s",
          pn, length(PAN[[pn]]$genes), sel$k, median(PAN[[pn]]$tab$log2FC_MLE),
          fmt_p(max(PAN[[pn]]$tab$padj, na.rm = TRUE)), fmt_p(min(PAN[[pn]]$tab$padj, na.rm = TRUE)))
}
le_all <- rbindlist(le_tables, fill = TRUE)
sets_all <- rbindlist(set_tables)

# ---- replicate consistency + selection-matched null (residual rotation)
consistency <- function(Xref, Xtrt) {
  mr <- rowMeans(Xref); mt <- rowMeans(Xtrt); d <- sign(mt - mr)
  trt_ok <- rowSums(sign(Xtrt - mr) == d) == ncol(Xtrt)
  ref_ok <- rowSums(sign(Xref - mt) == -d) == ncol(Xref)
  sep <- ifelse(d < 0, rowMaxs(Xtrt) < rowMins(Xref), rowMins(Xtrt) > rowMaxs(Xref))
  c(pct_treated_all3_agree = 100 * mean(trt_ok), pct_reference_all3_agree = 100 * mean(ref_ok),
    pct_both_agree = 100 * mean(trt_ok & ref_ok), pct_complete_separation = 100 * mean(sep),
    median_abs_diff_vst = median(abs(mt - mr)))
}
Dm <- model.matrix(~ 0 + cd$group)
Pm <- Dm %*% solve(crossprod(Dm)) %*% t(Dm)
eB <- eigen(diag(12) - Pm, symmetric = TRUE)
B8 <- eB$vectors[, eB$values > 0.5]                                # 12 x 8 orthonormal residual basis
stopifnot(ncol(B8) == 8)
haar <- function(n) { q <- qr(matrix(rnorm(n * n), n)); Q <- qr.Q(q); Q %*% diag(sign(diag(qr.R(q)))) }
tstat_sel <- function(X, ref_i, trt_i, s2, set_ids) {     # top-k most negative t per set (= leading edge top by |stat|)
  tt <- (rowMeans(X[, trt_i, drop = FALSE]) - rowMeans(X[, ref_i, drop = FALSE])) / sqrt(s2 * (2 / 3))
  by_set <- lapply(set_ids, function(ids) ids[order(tt[ids])])
  select_top(by_set)$ids
}
cons_rows <- list(); per_sample <- list(); null_store <- list()
for (pn in names(PAN)) {
  p <- PAN[[pn]]
  X <- V[p$genes, S]
  obs <- consistency(X[, p$ref_s], X[, p$trt_s])
  # per library: share of displayed genes on its own group's side of the 6-library mean
  m6 <- rowMeans(X[, c(p$ref_s, p$trt_s)]); dsg <- sign(rowMeans(X[, p$trt_s]) - rowMeans(X[, p$ref_s]))
  med_d <- function(tr, rf) median(rowMeans(X[, tr, drop = FALSE]) - rowMeans(X[, rf, drop = FALSE]))
  ps <- rbindlist(lapply(c(p$ref_s, p$trt_s), function(s) {
    expect <- if (s %in% p$trt_s) dsg else -dsg
    data.table(panel = pn, sample_id = s, library = cd$short[S == s], group = grp_lab[[s]],
               replicate_factor_B = batchB[[s]], pct_reads_outside_exons = pct_out[[s]],
               pct_genes_on_own_group_side = round(100 * mean(sign(X[, s] - m6) == expect), 1),
               median_diff_from_own_group_mean = round(median(X[, s] - rowMeans(X[, if (s %in% p$trt_s) p$trt_s else p$ref_s])), 3),
               median_group_diff_all_6 = round(med_d(p$trt_s, p$ref_s), 3),
               median_group_diff_without_this_library = round(med_d(setdiff(p$trt_s, s), setdiff(p$ref_s, s)), 3))
  }))
  per_sample[[pn]] <- ps
  # null: universe = all genes of the panel's sets present in the ranked list
  set_ids <- lapply(p$sets, function(s) p$rank[entrez %in% H[term == s, entrez]]$gene_id)
  U <- unique(unlist(set_ids))
  XU <- V[U, S]; mu <- rowMeans(XU); Zr <- XU %*% B8
  q0 <- rep(0, 12); q0[S %in% p$trt_s] <- 1 / 3; q0[S %in% p$ref_s] <- -1 / 3; q0 <- q0 / sqrt(sum(q0^2))
  ref_i <- match(p$ref_s, S); trt_i <- match(p$trt_s, S)
  s2_obs <- rowSums(Zr^2) / 8
  sel_obs_t <- tstat_sel(XU, ref_i, trt_i, s2_obs, set_ids)
  obs_t <- consistency(XU[sel_obs_t, ref_i], XU[sel_obs_t, trt_i])
  set.seed(SEED + match(pn, names(PAN)))
  nul <- t(vapply(seq_len(N_ROT), function(b) {
    Zs <- Zr %*% haar(8)
    Xs <- mu + outer(Zs[, 1], q0) + Zs[, 2:8] %*% t(B8[, 2:8])
    s2 <- rowSums(Zs[, 2:8]^2) / 7
    sel <- tstat_sel(Xs, ref_i, trt_i, s2, set_ids)
    consistency(Xs[sel, ref_i, drop = FALSE], Xs[sel, trt_i, drop = FALSE])
  }, numeric(5)))
  null_store[[pn]] <- nul
  qs <- function(v) sprintf("%.0f (%.0f-%.0f)", mean(v), quantile(v, 0.025), quantile(v, 0.975))
  tb <- p$tab
  cons_rows[[pn]] <- data.table(
    panel = pn, contrast = p$contrast, gene_sets = paste(SET_LAB[p$sets], collapse = " + "),
    libraries = paste(c(cd$short[S %in% p$ref_s], "vs", cd$short[S %in% p$trt_s]), collapse = " "),
    n_genes_displayed = length(p$genes), k_per_set = p$k,
    median_log2FC_MLE = round(median(tb$log2FC_MLE), 3), IQR_log2FC_MLE = paste(round(quantile(tb$log2FC_MLE, c(.25, .75)), 3), collapse = " to "),
    range_log2FC_MLE = paste(round(range(tb$log2FC_MLE), 3), collapse = " to "),
    median_fold_change_pct = round(100 * (2^median(tb$log2FC_MLE) - 1), 1),
    median_log2FC_shrunk_apeglm = signif(median(tb$log2FC_shrunk), 2),
    Wald_stat_range = paste(round(range(tb$stat), 2), collapse = " to "),
    genes_padj_lt_0.05 = sum(tb$padj < 0.05, na.rm = TRUE), min_gene_padj = signif(min(tb$padj, na.rm = TRUE), 2),
    pct_treated_all3_agree = round(obs[1], 1), pct_reference_all3_agree = round(obs[2], 1),
    pct_both_groups_all_replicates_agree = round(obs[3], 1), pct_complete_separation = round(obs[4], 1),
    null_universe_genes = length(U),
    tsel_observed_pct_both_agree = round(obs_t[3], 1), tsel_observed_pct_complete_separation = round(obs_t[4], 1),
    tsel_observed_median_abs_diff = round(obs_t[5], 3),
    null_pct_both_agree_mean_95 = qs(nul[, 3]), null_pct_complete_separation_mean_95 = qs(nul[, 4]),
    null_median_abs_diff_mean_95 = sprintf("%.3f (%.3f-%.3f)", mean(nul[, 5]), quantile(nul[, 5], .025), quantile(nul[, 5], .975)),
    null_frac_ge_observed_complete_sep = round(mean(nul[, 4] >= obs_t[4]), 3),
    null_frac_ge_observed_median_abs_diff = round(mean(nul[, 5] >= obs_t[5]), 4),
    other_context = p$other_label,
    pct_same_sign_other_context = round(100 * mean(sign(tb$other_log2FC_MLE) == sign(tb$log2FC_MLE), na.rm = TRUE), 1),
    r_log2FC_with_other_context = round(cor(tb$log2FC_MLE, tb$other_log2FC_MLE, use = "complete.obs"), 2),
    pct_in_trend_adjusted_leading_edge = round(100 * mean(Reduce(`|`, tb[, .SD, .SDcols = patterns("^LEtrendadj_")])), 1),
    pct_in_salmon_leading_edge = round(100 * mean(Reduce(`|`, tb[, .SD, .SDcols = patterns("^LEsalmon_")])), 1),
    wgcna_modules_round1 = paste(names(sort(table(tb$wgcna_module_round1), decreasing = TRUE)), sort(table(tb$wgcna_module_round1), decreasing = TRUE), sep = ":", collapse = " "),
    wgcna_modules_sensitivity = paste(names(sort(table(tb$wgcna_module_sensitivity), decreasing = TRUE)), sort(table(tb$wgcna_module_sensitivity), decreasing = TRUE), sep = ":", collapse = " "))
}
cons <- rbindlist(cons_rows); per_sample <- rbindlist(per_sample)
print(t(cons))
print(per_sample)

# =============================================================================
# Part 4. Correlation-aware gene-set tests
#   Preranked GSEA (Fig 2; fgsea gene permutation) treats the genes of a set as
#   independent. The same Hallmark sets are tested here with methods that keep the
#   gene-gene correlation (limma): camera, competitive, inter-gene correlation fixed at
#   0.01 (= cameraPR default used in verify_gsea) or estimated from the data; fry and
#   roast, self-contained rotation tests of the residual space. Corrected VST, design
#   ~ 0 + group (cell means; same contrasts as DESeq2); fry also on voom(counts);
#   exploratory: fry with the replicate-level factor (#3, #6, #7, #10) as a block.
#   Leave-one-library-out fry for the three panel contrasts.
# =============================================================================
log_msg("Part 4: correlation-aware gene-set tests (limma camera / fry / roast)")
suppressPackageStartupMessages({ library(limma); library(edgeR) })
dY <- data.table(id = rownames(V), entrez = ann[rownames(V), "entrez"], m = rowMeans(V))[!is.na(entrez)][order(-m)][!duplicated(entrez)]
Y  <- V[dY$id, S]                                                  # one gene per Entrez ID (highest mean)
idxY <- lapply(split(H$entrez, H$term), function(s) which(dY$entrez %in% s)); idxY <- idxY[lengths(idxY) >= 15]
gfac <- factor(as.character(cd$group), levels = c("vector_ctrl", "vector_CM", "NKX3.1_ctrl", "NKX3.1_CM"))
DES  <- model.matrix(~ 0 + gfac); colnames(DES) <- levels(gfac)
DESB <- cbind(DES, replicate_factor = as.numeric(batchB))
CON <- list(CM_in_vector = c(-1, 1, 0, 0), NKX_in_ctrl = c(-1, 0, 1, 0), NKX_in_CM = c(0, -1, 0, 1),
            CM_in_NKX = c(0, 0, -1, 1), CM_main = c(-.5, .5, -.5, .5), NKX_main = c(-.5, -.5, .5, .5))
CON_LAB <- c(CM_in_vector = "CM vs control, vector cells", NKX_in_ctrl = "NKX3.1 vs control, basal",
             NKX_in_CM = "NKX3.1 vs control, CM", CM_in_NKX = "CM vs control, NKX3.1 cells",
             CM_main = "CM, both cell types (12 libraries)", NKX_main = "NKX3.1, both conditions (12 libraries)")
VO <- voom(calcNormFactors(DGEList(core$counts_fc[dY$id, S])), DES)
set_tests <- list()
for (cn in names(CON)) {
  cc <- CON[[cn]]; sets <- names(idxY)
  f_v  <- fry(Y, idxY, DES, cc); f_vo <- fry(VO, idxY, DES, cc); f_b <- fry(Y, idxY, DESB, c(cc, 0))
  c01  <- camera(Y, idxY, DES, cc, inter.gene.cor = 0.01)
  cest <- camera(Y, idxY, DES, cc, inter.gene.cor = NA)
  set.seed(SEED); ro <- mroast(Y, idxY[unname(FOCAL)], DES, cc, nrot = 9999, set.statistic = "mean")
  gs <- read_gsea("fc", cn)
  set_tests[[cn]] <- data.table(
    contrast = cn, contrast_label = CON_LAB[[cn]], set = sets, n_genes = lengths(idxY)[sets],
    GSEA_NES = round(gs$NES[match(sets, gs$ID)], 2), GSEA_FDR = gs$p.adjust[match(sets, gs$ID)],
    camera_r0.01_direction = c01[sets, "Direction"], camera_r0.01_p = c01[sets, "PValue"], camera_r0.01_FDR = c01[sets, "FDR"],
    camera_estimated_r = round(cest[sets, "Correlation"], 3), camera_est_direction = cest[sets, "Direction"],
    camera_est_p = cest[sets, "PValue"], camera_est_FDR = cest[sets, "FDR"],
    fry_direction = f_v[sets, "Direction"], fry_p = f_v[sets, "PValue"], fry_FDR = f_v[sets, "FDR"],
    fry_voom_p = f_vo[sets, "PValue"], fry_voom_FDR = f_vo[sets, "FDR"],
    fry_block_direction = f_b[sets, "Direction"], fry_block_p = f_b[sets, "PValue"], fry_block_FDR = f_b[sets, "FDR"],
    roast_mean_9999rot_p = ro[match(sets, rownames(ro)), "PValue"])
}
set_tests <- rbindlist(set_tests)
# within-set residual correlation (group-means model) and the variance inflation it implies for a set mean
fitY <- lmFit(Y, DES); RY <- residuals(fitY, Y)
fitB <- lmFit(Y, DESB); RB <- residuals(fitB, Y)
mean_r <- function(R) { cc <- cor(t(R)); mean(cc[upper.tri(cc)]) }
set.seed(SEED); bg <- sample(nrow(RY), 2000)
rho_tab <- data.table(set = unname(FOCAL), set_label = unname(SET_LAB[FOCAL]), n_genes = lengths(idxY[unname(FOCAL)]),
                      mean_residual_r = sapply(unname(FOCAL), function(s) round(mean_r(RY[idxY[[s]], ]), 3)),
                      mean_residual_r_replicate_factor_blocked = sapply(unname(FOCAL), function(s) round(mean_r(RB[idxY[[s]], ]), 3)))
rho_tab[, variance_inflation_of_set_mean := round(1 + (n_genes - 1) * mean_residual_r, 1)]
rho_bg <- round(mean_r(RY[bg, ]), 3); rho_bg_b <- round(mean_r(RB[bg, ]), 3)
rho_tab <- rbind(rho_tab, data.table(set = "2,000 random expressed genes", set_label = "background", n_genes = 2000L,
                                     mean_residual_r = rho_bg, mean_residual_r_replicate_factor_blocked = rho_bg_b,
                                     variance_inflation_of_set_mean = NA_real_))
print(rho_tab)
n_sig_methods <- set_tests[, .(GSEA = sum(GSEA_FDR < 0.05, na.rm = TRUE), camera_r0.01 = sum(camera_r0.01_FDR < 0.05),
                               camera_estimated = sum(camera_est_FDR < 0.05), fry = sum(fry_FDR < 0.05),
                               fry_voom = sum(fry_voom_FDR < 0.05), fry_block = sum(fry_block_FDR < 0.05),
                               min_fry_FDR = signif(min(fry_FDR), 2), min_camera_est_FDR = signif(min(camera_est_FDR), 2)),
                           by = .(contrast, contrast_label)]
print(n_sig_methods)
# leave one library out (panel contrasts, panel sets)
loo <- rbindlist(lapply(names(PANELS), function(pn) {
  p <- PANELS[[pn]]; cc <- CON[[p$contrast]]
  rbindlist(lapply(which(cd$group %in% c(p$ref, p$trt)), function(l) {
    keep <- seq_len(12) != l
    f <- fry(Y[, keep], idxY[unname(p$sets)], DES[keep, ], cc)
    data.table(panel = pn, contrast = p$contrast, left_out = cd$short[l], left_out_group = cd$label[l],
               left_out_in_replicate_factor_B = batchB[[l]], set = rownames(f), direction = f$Direction, fry_p = signif(f$PValue, 3))
  }))
}))
focal_tests <- set_tests[set %in% FOCAL]
focal_tests <- merge(focal_tests, loo[, .(loo_fry_p_min = min(fry_p), loo_fry_p_max = max(fry_p),
                                          loo_library_giving_max_p = left_out[which.max(fry_p)]), by = .(contrast, set)],
                     by = c("contrast", "set"), all.x = TRUE)
print(focal_tests[contrast %in% c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM", "CM_main"),
                  .(contrast, set = SET_LAB[set], GSEA_FDR = signif(GSEA_FDR, 2), cam01 = signif(camera_r0.01_FDR, 2),
                    cam_est_p = signif(camera_est_p, 2), fry_p = signif(fry_p, 2), fry_FDR = signif(fry_FDR, 2),
                    fry_block_p = signif(fry_block_p, 2), roast = signif(roast_mean_9999rot_p, 2), loo_fry_p_max)])

# =============================================================================
# Figures
# =============================================================================
log_msg("figures")
LFC_LO <- -1                                                       # shared log2FC axis (-1 to 0); bars below are truncated

panel_ht <- function(p, samples, mode = c("z", "unscaled"), cell_w = 3.8, cell_h = 3.2) {
  mode <- match.arg(mode)
  X <- V[p$genes, samples]
  M <- if (mode == "z") t(scale(t(X))) else X - rowMeans(V[p$genes, p$ref_s])
  col_fun <- if (mode == "z") colorRamp2(c(-2, 0, 2), PAL) else colorRamp2(c(-0.5, 0, 0.5), PAL)
  glab <- factor(grp_lab[samples], levels = names(GRP_COL)[names(GRP_COL) %in% grp_lab[samples]])
  top <- HeatmapAnnotation(Group = glab, col = list(Group = GRP_COL[levels(glab)]), show_legend = FALSE,
                           show_annotation_name = FALSE, simple_anno_size = unit(2.2, "mm"), gp = gpar(col = NA))
  tb <- p$tab
  lfc_bar <- pmin(pmax(tb$log2FC_MLE, LFC_LO), 0)
  gene_lab <- ifelse(tb$log2FC_MLE < LFC_LO,
                     sprintf("%s (%s)", tb$symbol, sub("-", "−", sprintf("%.2f", tb$log2FC_MLE), fixed = TRUE)), tb$symbol)
  right <- rowAnnotation(
    log2FC = anno_barplot(lfc_bar, baseline = 0, ylim = c(LFC_LO, 0), bar_width = 0.72, border = FALSE,
                          gp = gpar(fill = "grey45", col = NA), width = unit(10, "mm"),
                          axis_param = list(side = "bottom", at = c(-1, -0.5, 0), labels = c("−1", "−0.5", "0"),
                                            labels_rot = 0, gp = gpar(fontsize = 5.5))),
    gene = anno_text(gene_lab, gp = gpar(fontsize = 6, fontface = "italic"), just = "left", location = unit(0.6, "mm")),
    show_annotation_name = c(log2FC = TRUE, gene = FALSE), annotation_name_side = "top", annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 6), annotation_name_offset = unit(0.6, "mm"), gap = unit(0.6, "mm"))
  rs <- NULL; rt <- NULL
  if (length(p$sets) == 2) {                                        # row blocks by leading-edge membership
    nm <- names(p$sets)
    rs <- droplevels(factor(tb$membership, levels = c(paste(nm[1], "only"), "both", paste(nm[2], "only"))))
    disp <- c(setNames(unname(SET_LAB[p$sets]), paste(nm, "only")), both = "Both")
    rt <- unname(disp[levels(rs)])
  }
  Heatmap(M, name = paste0("ht", mode), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
          row_split = rs, row_title = rt, row_title_rot = 90, row_title_gp = gpar(fontsize = 6),
          row_gap = unit(1.2, "mm"), cluster_row_slices = FALSE,
          column_split = glab, column_gap = unit(0.8, "mm"), column_title = NULL,
          show_row_names = FALSE, show_column_names = TRUE, column_labels = cd$short[match(samples, S)],
          column_names_rot = if (length(samples) > 6) 90 else 0, column_names_centered = TRUE,
          column_names_gp = gpar(fontsize = 5), top_annotation = top, right_annotation = right,
          width = unit(cell_w * length(samples), "mm"), height = unit(cell_h * nrow(M), "mm"),
          rect_gp = gpar(col = "white", lwd = 0.25), border = FALSE, show_heatmap_legend = FALSE)
}
ht_dims <- function(ht) {                                          # drawn size in inches
  pdf(NULL); on.exit(dev.off())
  d <- draw(ht, padding = unit(c(0, 0, 0, 0), "mm"))
  c(w = convertWidth(ComplexHeatmap:::width(d), "inch", valueOnly = TRUE),
    h = convertHeight(ComplexHeatmap:::height(d), "inch", valueOnly = TRUE))
}
header_lines <- function(pn) {                                    # set name + NES (statistics live in Fig 2 / legend)
  p <- PAN[[pn]]
  sprintf("%s (NES %s)", p$sets_tab$set_label, sub("-", "−", sprintf("%.2f", p$sets_tab$NES), fixed = TRUE))
}
draw_figure <- function(pns, samples_of, mode, legend_title, file, fig_w, cell_w, cell_h = 3.2, gap_in = 0.16,
                        letters = pns) {
  hts <- lapply(pns, function(pn) panel_ht(PAN[[pn]], samples_of(PAN[[pn]]), mode, cell_w, cell_h))
  dims <- t(sapply(hts, ht_dims))
  head_h <- 0.50; leg_h <- if (grepl("\n", legend_title)) 0.52 else 0.42; mar <- 0.08
  fig_h <- mar + head_h + max(dims[, "h"]) + leg_h + mar
  need_w <- 2 * mar + sum(dims[, "w"]) + gap_in * (length(pns) - 1)
  if (need_w > fig_w + 1e-6) log_msg("  note: %s needs %.2f in (> %.2f in requested); width enlarged", file, need_w, fig_w)
  fig_w <- max(fig_w, need_w)
  extra <- (fig_w - need_w) / length(pns)
  col_fun <- if (mode == "z") colorRamp2(c(-2, 0, 2), PAL) else colorRamp2(c(-0.5, 0, 0.5), PAL)
  grp_present <- names(GRP_COL)[names(GRP_COL) %in% unique(unlist(lapply(pns, function(pn) grp_lab[samples_of(PAN[[pn]])])))]
  draw_fun <- function() {
    grid.newpage()
    pushViewport(viewport(width = unit(fig_w, "in"), height = unit(fig_h, "in")))
    x <- mar
    for (i in seq_along(pns)) {
      pn <- pns[i]; w <- dims[i, "w"]; h <- dims[i, "h"]
      x0 <- x + extra / 2
      ytop <- fig_h - mar
      grid.text(letters[i], x = unit(x0, "in"), y = unit(ytop, "in"), just = c("left", "top"),
                gp = gpar(fontsize = 9, fontface = "bold"))
      grid.text(PAN[[pn]]$title, x = unit(x0 + 0.14, "in"), y = unit(ytop - 0.005, "in"), just = c("left", "top"),
                gp = gpar(fontsize = 7.5, fontface = "bold"))
      hl <- header_lines(pn)
      for (j in seq_along(hl)) grid.text(hl[[j]], x = unit(x0 + 0.14, "in"), y = unit(ytop - 0.165 - 0.105 * (j - 1), "in"),
                                         just = c("left", "top"), gp = gpar(fontsize = 6.5))
      pushViewport(viewport(x = unit(x0, "in"), y = unit(ytop - head_h, "in"), width = unit(w, "in"), height = unit(h, "in"),
                            just = c("left", "top")))
      draw(hts[[i]], newpage = FALSE, padding = unit(c(0, 0, 0, 0), "mm"))
      popViewport()
      x <- x + w + gap_in + extra
    }
    lg1 <- Legend(col_fun = col_fun, title = legend_title, direction = "horizontal",
                  at = if (mode == "z") -2:2 else c(-0.5, -0.25, 0, 0.25, 0.5),
                  labels = sub("-", "−", if (mode == "z") as.character(-2:2) else c("-0.5", "-0.25", "0", "0.25", "0.5"), fixed = TRUE),
                  legend_width = unit(26, "mm"), grid_height = unit(2.2, "mm"), title_position = "topcenter",
                  title_gp = gpar(fontsize = 6.5, fontface = "plain"), labels_gp = gpar(fontsize = 6))
    lg2 <- Legend(labels = grp_present, legend_gp = gpar(fill = GRP_COL[grp_present], col = NA), title = "Group",
                  title_position = "topleft", title_gp = gpar(fontsize = 6.5, fontface = "plain"),
                  nrow = 1, labels_gp = gpar(fontsize = 6.5), grid_width = unit(2.6, "mm"), grid_height = unit(2.6, "mm"),
                  column_gap = unit(2.5, "mm"))
    pk <- packLegend(lg1, lg2, direction = "horizontal", column_gap = unit(7, "mm"))
    draw(pk, x = unit(fig_w / 2, "in"), y = unit(mar + 0.02, "in"), just = "bottom")   # c("centre","bottom") is not drawn
    popViewport()
  }
  save_pdf_png(draw_fun, file, fig_w, fig_h)
  log_msg("  wrote %s (%.2f x %.2f in)", file, fig_w, fig_h)
  c(width = fig_w, height = fig_h)
}
six  <- function(p) c(p$ref_s, p$trt_s)
all12 <- function(p) S
fig_sizes <- list()
fig_sizes$main6 <- draw_figure(c("A", "B", "C"), six, "z", "Row z-score (VST)", "fig3alt_candidate_leading_edge_heatmaps_6samples", 7.0, 3.8)
fig_sizes$all12 <- draw_figure(c("A", "B", "C"), all12, "z", "Row z-score (VST, 12 libraries)", "fig3alt_candidate_leading_edge_heatmaps_12samples", 7.0, 2.25)
fig_sizes$unscaled <- draw_figure(c("A", "B", "C"), six, "unscaled", "log2 expression relative\nto the control-group mean",
                                  "fig3alt_candidate_leading_edge_unscaled_log2_vs_control", 7.0, 3.8)
fig_sizes$AC <- draw_figure(c("A", "C"), six, "z", "Row z-score (VST)", "fig3alt_candidate_leading_edge_heatmaps_AC_only", 4.8, 3.8,
                            letters = c("A", "B"))

# ---- gene-set evidence by method (focal sets of the three panels + the CM main effect)
okabe <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")
meth <- c(GSEA_FDR = "GSEA, genes permuted (as Fig 2)", camera_r0.01_FDR = "camera, inter-gene r fixed at 0.01",
          camera_est_FDR = "camera, inter-gene r estimated", fry_FDR = "fry rotation test",
          fry_block_FDR = "fry + replicate-factor block (exploratory)")
show <- rbind(data.table(contrast = "CM_in_vector", set = FOCAL[c("MYC", "OXPHOS")], panel = "A"),
              data.table(contrast = "NKX_in_ctrl", set = FOCAL[c("E2F", "G2M")], panel = "B"),
              data.table(contrast = "NKX_in_CM", set = FOCAL["AJ"], panel = "C"),
              data.table(contrast = "CM_main", set = FOCAL[c("MYC", "OXPHOS")], panel = "CM main"))
fd <- merge(show, set_tests, by = c("contrast", "set"), sort = FALSE)
fd[, row_lab := sprintf("%s: %s, %s", panel, SET_LAB[set], CON_LAB[contrast])]
fl <- melt(fd, id.vars = c("row_lab"), measure.vars = names(meth), variable.name = "method", value.name = "FDR")
fl[, method := factor(meth[as.character(method)], levels = meth)]
fl[, row_lab := factor(row_lab, levels = rev(unique(fd$row_lab)))]
fl[, x := pmax(FDR, 1e-16)]
pS <- ggplot(fl, aes(x, row_lab, colour = method, shape = method)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
  geom_point(size = 1.9, stroke = 0.6, position = position_dodge(width = 0.55)) +
  scale_x_continuous(trans = scales::compose_trans("log10", "reverse"), breaks = c(1, 0.05, 1e-4, 1e-8, 1e-12, 1e-16),
                     labels = c("1", "0.05", expression(10^-4), expression(10^-8), expression(10^-12), expression(""<=10^-16)),
                     limits = c(1, 1e-16)) +
  scale_colour_manual(values = setNames(okabe, meth), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 17, 15, 18, 4), meth), name = NULL) +
  labs(x = "FDR across the 50 Hallmark sets (log scale)", y = NULL) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.text = element_text(colour = "black", size = 6.5), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.height = unit(3, "mm"), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        plot.background = element_rect(fill = "white", colour = NA))
save_pdf_png(function() print(pS), "fig3alt_set_tests_by_method", 7.0, 3.1)

# ---- WGCNA figure: focal over-representation and eigengene correlations, both networks
minus_lab <- function(x) sub("-", "−", formatC(x, format = "f", digits = 1), fixed = TRUE)
wg_parts <- function(n, title) {
  vv <- verdict[network == n$label]
  mods <- sub("^ME", "", colnames(n$MEs))                          # orderMEs order
  vv <- vv[match(mods, module)]
  ylab <- sprintf("%s (%d)", vv$module, vv$size)
  e <- enr[network == n$label & set %in% FOCAL]
  e[, set_lab := factor(SET_LAB[set], levels = SET_LAB[FOCAL])]
  e[, mod_lab := factor(ylab[match(module, mods)], levels = rev(ylab))]
  e[, nl := pmin(-log10(padj), 12)]
  e[, lab := ifelse(padj < 0.05, sprintf("%d/%d", overlap, set_genes_in_network), "")]
  th <- theme_minimal(base_size = 6.5, base_family = FONT) +
    theme(panel.grid = element_blank(), axis.text.x.top = element_text(angle = 40, hjust = 0, vjust = 0, size = 6, colour = "black"),
          legend.position = "bottom", legend.key.height = unit(2.2, "mm"), legend.key.width = unit(6, "mm"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
          plot.subtitle = element_text(size = 6.5, face = "bold"), plot.title.position = "plot",
          plot.title = element_text(size = 7.5, face = "bold", margin = margin(b = 4)))
  pE <- ggplot(e, aes(set_lab, mod_lab, fill = nl)) +
    geom_tile(colour = "grey85", linewidth = 0.25) +
    geom_text(aes(label = lab), size = 1.9, colour = ifelse(e$nl > 6, "white", "black")) +
    scale_fill_gradient(low = "white", high = "#3F007D", limits = c(0, 12), breaks = c(0, 4, 8, 12),
                        labels = c("0", "4", "8", "≥12"), name = expression(-log[10]~FDR)) +
    scale_x_discrete(position = "top") + labs(x = NULL, y = NULL, title = title, subtitle = "Hallmark over-representation") +
    th + theme(axis.text.y = element_text(size = 6, colour = "black"))
  tr <- melt(vv[, .(module, r_NKX, r_CM, r_NKXxCM, r_pct_reads_outside_exons, r_replicate_factor_3_6_7_10,
                    padj_exact_NKX, padj_exact_CM, padj_exact_NKXxCM)],
             id.vars = "module", measure.vars = list(r = c("r_NKX", "r_CM", "r_NKXxCM", "r_pct_reads_outside_exons", "r_replicate_factor_3_6_7_10"),
                                                     padj = c("padj_exact_NKX", "padj_exact_CM", "padj_exact_NKXxCM", "padj_exact_NKX", "padj_exact_NKX")))
  tr_lab <- c("NKX3.1", "CM", "NKX3.1×CM", "% reads outside exons", "Replicate factor\n(#3, #6, #7, #10)")
  tr[, trait := factor(tr_lab[as.integer(variable)], levels = tr_lab)]
  tr[as.integer(variable) > 3, padj := NA]
  tr[, mod_lab := factor(ylab[match(module, mods)], levels = rev(ylab))]
  tr[, lab := sprintf("%s%s", sub("-", "−", sprintf("%.2f", r), fixed = TRUE), ifelse(!is.na(padj) & padj < 0.05, "*", ""))]
  pT <- ggplot(tr, aes(trait, mod_lab, fill = r)) +
    geom_tile(colour = "grey85", linewidth = 0.25) +
    geom_text(aes(label = lab), size = 1.9) +
    geom_vline(xintercept = 3.5, linewidth = 0.4, colour = "grey30") +
    scale_fill_gradient2(low = PAL[1], mid = PAL[2], high = PAL[3], midpoint = 0, limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1), labels = minus_lab, name = "Eigengene r") +
    scale_x_discrete(position = "top") + labs(x = NULL, y = NULL, title = " ", subtitle = "Eigengene correlation (12 libraries)") +
    th + theme(axis.text.y = element_blank(), plot.margin = margin(4, 18, 4, 4))
  list(pE = pE, pT = pT)
}
nm1 <- ncol(net1$MEs); nm2 <- ncol(net2$MEs)
w1p <- wg_parts(net1, sprintf("A  WGCNA network of the %s most variable genes (soft power %d)", format(net1$n_genes, big.mark = ","), net1$power))
w2p <- wg_parts(net2, sprintf("B  WGCNA network of all %s expressed genes (soft power %d)", format(net2$n_genes, big.mark = ","), net2$power))
pW <- (w1p$pE | w1p$pT) / (w2p$pE | w2p$pT) +
  plot_layout(heights = c(nm1, nm2), guides = "collect") & theme(legend.position = "bottom")
save_pdf_png(function() print(pW), "fig3alt_wgcna_modules_vs_hallmark", 7.0, 2.6 + 0.135 * (nm1 + nm2))

# ---- the original display on the real and on one nonsense grouping (panel B, Salmon VST, as published)
nb <- nonsense[panel == "B" & quantification == "original (Salmon VST, published)"]
pick <- nb[real_grouping == FALSE][order(abs(genes_passing_filter - median(genes_passing_filter)), split)][1]
oB <- ORIG$B; idx6 <- c(which(cond0 == oB$g1), which(cond0 == oB$g2))
sp <- combn(6, 3); sp <- sp[, sp[1, ] == 1, drop = FALSE]
orig_display <- function(j) {
  a <- idx6[sp[, j]]; b <- setdiff(idx6, a)
  h <- orig_rows(V0, rn0, a, b)
  samp <- st$sample_id[c(a, b)]
  glab <- factor(grp_lab[samp], levels = c("Control-basal", "NKX3.1-basal"))
  sel_lab <- factor(rep(c("grouping 1", "grouping 2"), each = 3), levels = c("grouping 1", "grouping 2"))
  top <- HeatmapAnnotation(`True group` = glab, `Grouping used` = sel_lab,
                           col = list(`True group` = GRP_COL[levels(glab)], `Grouping used` = c(`grouping 1` = "grey20", `grouping 2` = "grey75")),
                           show_legend = FALSE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 6),
                           simple_anno_size = unit(2.2, "mm"), gap = unit(0.5, "mm"))
  Heatmap(h$mat, col = colorRamp2(c(-1.5, 0, 1.5), PAL), cluster_rows = h$hc, cluster_columns = FALSE,
          show_row_dend = TRUE, row_dend_width = unit(5, "mm"), row_dend_gp = gpar(lwd = 0.3),
          column_split = sel_lab, column_gap = unit(0.8, "mm"), column_title = NULL,
          show_row_names = FALSE, column_labels = cd$short[match(samp, S)], column_names_rot = 0,
          column_names_centered = TRUE, column_names_gp = gpar(fontsize = 5.5), top_annotation = top,
          width = unit(4 * 6, "mm"), height = unit(100, "mm"), show_heatmap_legend = FALSE, name = paste0("o", j))
}
htR <- orig_display(1); htN <- orig_display(pick$split)
draw_orig <- function() {
  grid.newpage()
  lay <- grid.layout(2, 2, heights = unit(c(0.42, 1), c("in", "null")))
  pushViewport(viewport(layout = lay))
  txt <- list(c("Real grouping", sprintf("%d genes pass the original filter; 200 shown", nb[real_grouping == TRUE]$genes_passing_filter)),
              c("Nonsense grouping (one library swapped)", sprintf("%d genes pass the original filter; 200 shown", pick$genes_passing_filter)))
  for (j in 1:2) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j))
    grid.text(c("A", "B")[j], x = unit(0.06, "in"), y = unit(1, "npc") - unit(0.05, "in"), just = c("left", "top"), gp = gpar(fontsize = 9, fontface = "bold"))
    grid.text(txt[[j]][1], x = unit(0.22, "in"), y = unit(1, "npc") - unit(0.06, "in"), just = c("left", "top"), gp = gpar(fontsize = 7, fontface = "bold"))
    grid.text(txt[[j]][2], x = unit(0.22, "in"), y = unit(1, "npc") - unit(0.20, "in"), just = c("left", "top"), gp = gpar(fontsize = 6))
    popViewport()
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j))
    draw(list(htR, htN)[[j]], newpage = FALSE, padding = unit(c(10, 2, 1, 2), "mm"))
    popViewport()
  }
  popViewport()
  lg1 <- Legend(col_fun = colorRamp2(c(-1.5, 0, 1.5), PAL), title = "Row z-score (original Salmon VST)", direction = "horizontal",
                at = c(-1.5, 0, 1.5), labels = c("−1.5", "0", "1.5"), legend_width = unit(24, "mm"), grid_height = unit(2.2, "mm"),
                title_position = "topcenter", title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5))
  lg2 <- Legend(labels = c("Control-basal", "NKX3.1-basal"), legend_gp = gpar(fill = GRP_COL[c("Control-basal", "NKX3.1-basal")], col = NA),
                nrow = 1, labels_gp = gpar(fontsize = 6), title = NULL, grid_width = unit(2.6, "mm"), grid_height = unit(2.6, "mm"))
  draw(packLegend(lg1, lg2, direction = "horizontal", column_gap = unit(6, "mm")), x = unit(0.5, "npc"), y = unit(0.04, "in"), just = "bottom")
}
save_pdf_png(draw_orig, "fig3alt_original_procedure_real_vs_nonsense_grouping", 4.6, 5.25)

# =============================================================================
# Tables, objects, session
# =============================================================================
log_msg("writing tables")
wtsv(orig_summary, "fig3alt_original_fig3_summary.tsv")
wtsv(row_status, "fig3alt_original_fig3_rows_corrected_status.tsv")
wtsv(nonsense, "fig3alt_original_filter_real_vs_nonsense_groupings.tsv")
wtsv(nonsense_summary, "fig3alt_original_filter_real_vs_nonsense_summary.tsv")
wtsv(enr[order(network, padj)], "fig3alt_wgcna_module_hallmark_enrichment.tsv")
wtsv(mtr, "fig3alt_wgcna_module_trait_exact_permutation.tsv")
wtsv(verdict, "fig3alt_wgcna_module_verdict.tsv")
wtsv(cov, "fig3alt_wgcna_focal_set_coverage.tsv")
wtsv(net2$sft, "fig3alt_wgcna_sensitivity_soft_threshold.tsv")
fwrite(net2$assign, file.path(out, "fig3alt_wgcna_sensitivity_gene_assignment.tsv.gz"), sep = "\t", quote = FALSE)
le_out <- copy(le_all)
setcolorder(le_out, c("panel", "contrast", "displayed", "display_row", "symbol", "gene_id", "entrez", "membership"))
wtsv(le_out, "fig3alt_leading_edge_genes.tsv")
wtsv(sets_all, "fig3alt_panel_gene_sets.tsv")
wtsv(cons, "fig3alt_replicate_consistency.tsv")
wtsv(per_sample, "fig3alt_replicate_consistency_per_library.tsv")
wtsv(set_tests[order(contrast, fry_p)], "fig3alt_set_tests_correlation_aware_all_hallmark.tsv")
wtsv(focal_tests[order(contrast, set)], "fig3alt_set_tests_correlation_aware_focal.tsv")
wtsv(n_sig_methods, "fig3alt_set_tests_n_significant_by_method.tsv")
wtsv(rho_tab, "fig3alt_set_tests_within_set_residual_correlation.tsv")
wtsv(loo, "fig3alt_set_tests_leave_one_library_out.tsv")
saveRDS(list(orig_summary = orig_summary, row_status = row_status, nonsense = nonsense, nonsense_summary = nonsense_summary,
             enrichment = enr, module_trait = mtr, verdict = verdict, coverage = cov,
             net2 = net2[c("label", "power", "rule", "sizes", "n_genes", "MEs", "sft")], net2_colors = net2$colors,
             panels = lapply(PAN, function(p) p[c("contrast", "sets", "genes", "k", "tab", "sets_tab", "ref_s", "trt_s")]),
             consistency = cons, per_library = per_sample, null = null_store, fig_sizes = fig_sizes,
             set_tests = set_tests, focal_tests = focal_tests, n_sig_methods = n_sig_methods, rho = rho_tab, loo = loo,
             msigdb_version = MSIG_VERSION, n_rot = N_ROT, seed = SEED),
        file.path(out, "fig3alt_results.rds"))
writeLines(c(sprintf("43_fig3_alternative.R finished %s (%.1f min)", format(Sys.time()), as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
             sprintf("msigdbr %s, MSigDB %s", as.character(packageVersion("msigdbr")), MSIG_VERSION), "",
             capture.output(sessionInfo())), file.path(out, "sessionInfo_43_fig3_alternative.txt"))
log_msg("done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
