#!/usr/bin/env Rscript
# =============================================================================
# 50_verify_fig3_02_recompute.R -- task verify_round2/fig3 (adversarial check of fig3_alternative)
#
# Recomputes, independently of 43_fig3_alternative.R where possible:
#   A. original Fig 3 row filter on the nine mislabelled groupings (which groupings pass most genes;
#      do they line up with the replicate-level factor #3/#6/#7/#10?)
#   B. WGCNA: eigengenes recomputed from the saved module assignments (both networks); the sensitivity
#      network re-run from scratch with identical settings; module x design association by
#      (i) unrestricted exact permutation over all 924 balanced labelings (as fig3alt),
#      (ii) stratified exact permutation (CM labels permuted within NKX3.1 strata and vice versa, 400 each),
#      (iii) the parametric 2x2 linear model; resolution floor of the exact test and the smallest BH value
#      a single module can reach; module x Hallmark Fisher tests for the focal sets; technical covariates
#      of the MYC/OXPHOS module (within-group correlations).
#   C. GC-content confound of the candidate panels: GC of displayed genes, per-library panel score vs
#      library GC-bias slope, and a GC-corrected DESeq2 + Hallmark GSEA (per-library loess GC offsets as
#      normalisation factors), plus leave-one-library-out GSEA (#12 for panel C, #1 for panels A/B).
#      The pipeline is first validated by reproducing the published corrected GSEA tables.
#   D. selection-matched null for the replicate-consistency claims with a corrected construction
#      (contrast noise + the four within-group residual dimensions of the six panel libraries drawn
#      from one random rotation of the eight residual dimensions).
#   E. correlation-aware set tests (limma camera / fry) recomputed for all 50 Hallmark sets.
#   (the diagnostic figure and the full leave-one-library-out GSEA are in 50_verify_fig3_03_loo_figure.R)
# Inputs read-only. Outputs: results/revision_2026/verify_round2/fig3/v02_*
# Run: /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/50_verify_fig3_02_recompute.R
# =============================================================================
if (!isTRUE(l10n_info()[["UTF-8"]])) for (loc in c("en_US.UTF-8", "C.UTF-8"))
  if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break

suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(matrixStats); library(WGCNA)
  library(limma); library(edgeR); library(msigdbr); library(clusterProfiler); library(BiocParallel)
})
options(stringsAsFactors = FALSE, warn = 1)
invisible(allowWGCNAThreads(16))
register(SerialParam())                         # GSEA (fgsea inside clusterProfiler) serial, as 12_gsea.R
BPP <- MulticoreParam(16)                       # DESeq2 only
t0 <- Sys.time()

proj <- "/mnt/ssd1/Projects/NKX3_wd"
revd <- file.path(proj, "results/revision_2026")
out  <- file.path(revd, "verify_round2/fig3")
work <- "/mnt/ssd1/Projects/NKX3_wd/work/verify_round2_fig3"
for (d in c(out, work)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
SEED <- 42L
FONT <- "Liberation Sans"
GRP_COL <- c("Control-basal" = "#909FC9", "Control-inflammation" = "#E28DBB",
             "NKX3.1-basal" = "#57C5D1", "NKX3.1-inflammation" = "#F5906A")
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
wtsv <- function(x, f) fwrite(as.data.table(x), file.path(out, f), sep = "\t", quote = FALSE, na = "NA")
scor <- function(...) stats::cor(...)

# -----------------------------------------------------------------------------
# data
# -----------------------------------------------------------------------------
core <- readRDS(file.path(revd, "core/core.rds"))
r1   <- readRDS(file.path(revd, "modules/modules_results.rds"))
f3   <- readRDS(file.path(revd, "fig3_alternative/fig3alt_results.rds"))
cd <- core$cd; ann <- core$ann; S <- cd$sample_id
cd$short <- paste0("#", cd$num)
V    <- core$fc$vst[, S]
NORM <- core$fc$norm[, S]
qc   <- fread(file.path(revd, "markers_qc/qc_sample_metrics.tsv"))[match(S, sample_id)]
stopifnot(identical(qc$sample_id, S))
batchB  <- setNames(qc$hidden_factor_split == "B", S)
fcsum   <- as.data.table(r1$fc_sum)[match(S, sample_id)]
pct_out <- setNames(fcsum$no_feature_pct, S)                       # as used by fig3alt
gcslope <- setNames(qc$gc_bias_slope_per10pct, S)
gtab <- fread(file.path(revd, "markers_qc/gene_gc_length.tsv"))
GC   <- setNames(gtab$gc, gtab$gene_id)
H_tab <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
H <- unique(H_tab[!is.na(ncbi_gene) & ncbi_gene != "", .(term = gs_name, entrez = as.character(ncbi_gene))])
T2G <- as.data.frame(H[, .(term, gene = entrez)])
FOCAL <- c(MYC = "HALLMARK_MYC_TARGETS_V1", OXPHOS = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           E2F = "HALLMARK_E2F_TARGETS", G2M = "HALLMARK_G2M_CHECKPOINT", AJ = "HALLMARK_APICAL_JUNCTION")
PANELS <- list(A = list(contrast = "CM_in_vector", sets = FOCAL[c("MYC", "OXPHOS")], ref = "vector_ctrl", trt = "vector_CM"),
               B = list(contrast = "NKX_in_ctrl",  sets = FOCAL[c("E2F", "G2M")],   ref = "vector_ctrl", trt = "NKX3.1_ctrl"),
               C = list(contrast = "NKX_in_CM",    sets = FOCAL["AJ"],              ref = "vector_CM",   trt = "NKX3.1_CM"))
le_tab <- fread(file.path(revd, "fig3_alternative/fig3alt_leading_edge_genes.tsv"))
disp <- lapply(setNames(nm = names(PANELS)), function(pn) le_tab[panel == pn & displayed == TRUE][order(display_row)]$gene_id)
grp_of <- setNames(as.character(cd$group), S)
lab_of <- setNames(cd$label, S)
gcen <- function(x) x - ave(x, cd$group)                           # group-centred (within-group) values

# =============================================================================
# A. original filter on mislabelled groupings: which groupings pass most genes?
# =============================================================================
log_msg("A: original Welch filter on the mislabelled groupings")
st <- read.csv(file.path(proj, "results/sample_table.csv"))
V0 <- core$vst_sal_v15[, st$sample_id]; cond0 <- st$condition
short0 <- cd$short[match(st$sample_id, S)]
B0 <- batchB[st$sample_id]
welch_n <- function(M, i1, i2) {
  mu1 <- rowMeans(M[, i1]); mu2 <- rowMeans(M[, i2])
  z <- (mu1 - mu2) / sqrt(pmax(rowVars(M[, i1]), 1e-6) / 3 + pmax(rowVars(M[, i2]), 1e-6) / 3)
  sum(pmax(mu1, mu2) >= 2 & abs(mu1 - mu2) >= 0.25 & abs(z) >= 2)
}
OG <- list(A = c("HM-vector-C", "HM-vector-CM"), B = c("HM-NKX3.1-C", "HM-vector-C"), C = c("HM-NKX3.1-CM", "HM-vector-CM"))
mis <- rbindlist(lapply(names(OG), function(k) {
  idx <- c(which(cond0 == OG[[k]][1]), which(cond0 == OG[[k]][2]))
  sp <- combn(6, 3); sp <- sp[, sp[1, ] == 1]
  rbindlist(lapply(seq_len(ncol(sp)), function(j) {
    a <- idx[sp[, j]]; b <- setdiff(idx, a)
    data.table(panel = k, split = j, real = j == 1, group1 = paste(short0[a], collapse = " "), group2 = paste(short0[b], collapse = " "),
               replicate_factor_libraries_together = (sum(B0[a]) == 2 || sum(B0[b]) == 2),
               genes_passing = welch_n(V0, a, b))
  }))
}))
mis_sum <- mis[, .(real = genes_passing[real], mislabelled_median = median(genes_passing[!real]),
                   mislabelled_ge200 = sum(genes_passing[!real] >= 200),
                   median_when_factor_libraries_together = median(genes_passing[!real & replicate_factor_libraries_together]),
                   median_when_apart = median(genes_passing[!real & !replicate_factor_libraries_together]),
                   max_grouping = paste(group1[which.max(genes_passing)], "vs", group2[which.max(genes_passing)]),
                   max_genes = max(genes_passing)), by = panel]
print(mis_sum)
wtsv(mis, "v02_A_mislabelled_groupings.tsv"); wtsv(mis_sum, "v02_A_mislabelled_groupings_summary.tsv")

# =============================================================================
# B. WGCNA
# =============================================================================
log_msg("B: WGCNA eigengenes, re-run, permutations")
w1 <- r1$wg$corrected
col1 <- setNames(w1$assign$module, w1$assign$gene_id)
col2 <- f3$net2_colors
recompute_ME <- function(colors) {
  me <- moduleEigengenes(t(V[names(colors), S]), colors, excludeGrey = TRUE)$eigengenes
  rownames(me) <- S; me
}
ME1 <- recompute_ME(col1); ME2 <- recompute_ME(col2)
cmp_me <- function(new, old, lab) {
  old <- as.data.frame(old)[S, , drop = FALSE]
  m <- intersect(colnames(new), colnames(old))
  data.table(network = lab, module = sub("^ME", "", m),
             r_recomputed_vs_saved = vapply(m, function(x) scor(new[[x]], old[[x]]), 0))
}
me_check <- rbind(cmp_me(ME1, w1$MEs, "round-1 (5,000 genes)"), cmp_me(ME2, f3$net2$MEs, "sensitivity (17,405 genes)"))
log_msg("  eigengene recomputation: min r with saved = %.6f", min(me_check$r_recomputed_vs_saved))

# re-run the sensitivity network with the settings of 43_fig3_alternative.R
expressed <- rownames(NORM)[rowMeans(NORM) >= 10]
stopifnot(setequal(expressed, names(col2)))
datExpr <- t(V[expressed, S])
powers <- c(1:10, seq(12, 30, by = 2))
invisible(capture.output(sft <- suppressWarnings(pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
  corFnc = bicor, corOptions = list(maxPOutliers = 0.05), blockSize = 6000, verbose = 0))))
fi <- sft$fitIndices; fi$signed_R2 <- -sign(fi$slope) * fi$SFT.R.sq
beta2 <- fi$Power[which(fi$signed_R2 >= 0.80)[1]]; if (is.na(beta2)) beta2 <- 18
assign("cor", WGCNA::cor, envir = globalenv())
net_re <- suppressWarnings(blockwiseModules(datExpr, power = beta2, networkType = "signed", TOMType = "signed",
  corType = "bicor", maxPOutliers = 0.05, minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
  maxBlockSize = ncol(datExpr) + 1, numericLabels = FALSE, saveTOMs = FALSE, randomSeed = SEED, nThreads = 16, verbose = 0))
rm("cor", envir = globalenv())
col_re <- setNames(net_re$colors, colnames(datExpr))
xt <- table(saved = col2[expressed], rerun = col_re[expressed])
rerun_check <- data.table(power_rerun = beta2, max_signed_R2 = round(max(fi$signed_R2), 3),
                          identical_partition = all(col2[expressed] == col_re[expressed]),
                          n_modules_saved = length(setdiff(unique(col2), "grey")), n_modules_rerun = length(setdiff(unique(col_re), "grey")),
                          genes_same_label = sum(col2[expressed] == col_re[expressed]))
print(rerun_check)

# module x design
traits <- data.frame(NKX = ifelse(cd$NKX == "NKX3.1", 1, -1), CM = ifelse(cd$CM == "CM", 1, -1)); traits$NKXxCM <- traits$NKX * traits$CM
LAB <- combn(12, 6); TALL <- apply(LAB, 2, function(L) { t <- rep(-1, 12); t[L] <- 1; t })
strat_labels <- function(stratum) {                             # balanced relabelings within the two strata of 'stratum'
  s1 <- which(stratum == 1); s2 <- which(stratum == -1)
  c1 <- combn(s1, 3); c2 <- combn(s2, 3)
  do.call(cbind, lapply(seq_len(ncol(c1)), function(i) vapply(seq_len(ncol(c2)), function(j) {
    t <- rep(-1, 12); t[c(c1[, i], c2[, j])] <- 1; t }, numeric(12))))
}
T_CM_in_NKX <- strat_labels(traits$NKX)                          # CM labels permuted within NKX3.1 strata (400)
T_NKX_in_CM <- strat_labels(traits$CM)                           # NKX3.1 labels permuted within CM strata (400)
mod_tests <- function(MEs, net) {
  MEs <- as.matrix(MEs[S, , drop = FALSE])
  rA <- scor(MEs, TALL); rC <- scor(MEs, T_CM_in_NKX); rN <- scor(MEs, T_NKX_in_CM)
  res <- rbindlist(lapply(colnames(MEs), function(m) {
    x <- MEs[, m]
    f <- summary(lm(x ~ traits$NKX + traits$CM + traits$NKXxCM))$coefficients
    xo <- data.table(network = net, module = sub("^ME", "", m), size = NA_integer_,
      r_CM = scor(x, traits$CM), r_NKX = scor(x, traits$NKX), r_NKXxCM = scor(x, traits$NKXxCM),
      p_exact924_CM = mean(abs(rA[m, ]) >= abs(scor(x, traits$CM)) - 1e-12),
      p_exact924_NKX = mean(abs(rA[m, ]) >= abs(scor(x, traits$NKX)) - 1e-12),
      p_exact924_NKXxCM = mean(abs(rA[m, ]) >= abs(scor(x, traits$NKXxCM)) - 1e-12),
      p_strat400_CM = mean(abs(rC[m, ]) >= abs(scor(x, traits$CM)) - 1e-12),
      p_strat400_NKX = mean(abs(rN[m, ]) >= abs(scor(x, traits$NKX)) - 1e-12),
      t_lm_CM = f[3, 3], p_lm_CM = f[3, 4], t_lm_NKX = f[2, 3], p_lm_NKX = f[2, 4], t_lm_NKXxCM = f[4, 3], p_lm_NKXxCM = f[4, 4],
      r_pct_reads_outside_exons = scor(x, pct_out), r_gc_slope = scor(x, gcslope),
      r_within_group_pct_outside = scor(gcen(x), gcen(pct_out)), r_within_group_gc_slope = scor(gcen(x), gcen(gcslope)),
      r_replicate_factor = scor(x, as.numeric(batchB)),
      cm_separation_complete = (max(x[traits$CM == 1]) < min(x[traits$CM == -1])) || (min(x[traits$CM == 1]) > max(x[traits$CM == -1])))
    xo
  }))
  m3 <- nrow(res) * 3
  res[, `:=`(padj_exact924 = NA_real_)]
  pe <- c(res$p_exact924_CM, res$p_exact924_NKX, res$p_exact924_NKXxCM); pa <- p.adjust(pe, "BH")
  res[, `:=`(padj_exact924_CM = pa[seq_len(.N)], padj_exact924_NKX = pa[.N + seq_len(.N)], padj_exact924_NKXxCM = pa[2 * .N + seq_len(.N)])]
  ps <- c(res$p_strat400_CM, res$p_strat400_NKX); pas <- p.adjust(ps, "BH", n = m3)
  res[, `:=`(padj_strat400_CM = pas[seq_len(.N)], padj_strat400_NKX = pas[.N + seq_len(.N)])]
  pl <- c(res$p_lm_CM, res$p_lm_NKX, res$p_lm_NKXxCM); pal <- p.adjust(pl, "BH")
  res[, `:=`(padj_lm_CM = pal[seq_len(.N)], padj_lm_NKX = pal[.N + seq_len(.N)], padj_lm_NKXxCM = pal[2 * .N + seq_len(.N)])]
  res[, padj_exact924 := NULL]
  res[, `:=`(n_tests_in_BH = m3, exact924_floor_p = 2 / 924, best_BH_single_test_at_floor = min(1, (2 / 924) * m3),
             strat400_floor_p = 2 / 400, best_BH_single_strat_test_at_floor = min(1, (2 / 400) * m3))]
  res[]
}
mt <- rbind(mod_tests(ME1, "round-1 (5,000 genes)"), mod_tests(ME2, "sensitivity (17,405 genes)"))
mt[network == "round-1 (5,000 genes)", size := as.integer(table(col1)[module])]
mt[network == "sensitivity (17,405 genes)", size := as.integer(table(col2)[module])]
# compare the unrestricted exact p with fig3alt
fx <- as.data.table(f3$module_trait)[, .(network = ifelse(grepl("round-1", network), "round-1 (5,000 genes)", "sensitivity (17,405 genes)"),
                                          module, trait, p_exact_fig3alt = p_exact, padj_fig3alt = padj_exact)]
chk <- merge(melt(mt[, .(network, module, CM = p_exact924_CM, NKX = p_exact924_NKX, NKXxCM = p_exact924_NKXxCM)],
                  id.vars = c("network", "module"), variable.name = "trait", value.name = "p_mine"), fx, by = c("network", "module", "trait"))
log_msg("  exact-permutation p: max |mine - fig3alt| = %.2e over %d tests", max(abs(chk$p_mine - chk$p_exact_fig3alt)), nrow(chk))
print(mt[order(network, pmin(p_exact924_CM, p_exact924_NKX))][, .(network, module, size, r_CM = round(r_CM, 2), p924_CM = signif(p_exact924_CM, 2),
      padj924_CM = signif(padj_exact924_CM, 2), pstrat_CM = signif(p_strat400_CM, 2), p_lm_CM = signif(p_lm_CM, 2), padj_lm_CM = signif(padj_lm_CM, 2),
      r_NKX = round(r_NKX, 2), p924_NKX = signif(p_exact924_NKX, 2), p_lm_NKX = signif(p_lm_NKX, 2), padj_lm_NKX = signif(padj_lm_NKX, 2),
      r_out = round(r_pct_reads_outside_exons, 2), r_out_within = round(r_within_group_pct_outside, 2), cm_sep = cm_separation_complete)][1:14])

# module x focal Hallmark (Fisher, one-sided; BH over modules x 50 sets per network), recomputed
fisher_tab <- function(colors, net) {
  U <- data.table(gene_id = names(colors), module = unname(colors))
  U[, entrez := ann[gene_id, "entrez"]]; U <- U[!is.na(entrez)][!duplicated(entrez)]
  mods <- setdiff(sort(unique(U$module)), "grey"); sets <- sort(unique(H$term))
  res <- rbindlist(lapply(mods, function(m) rbindlist(lapply(sets, function(s) {
    inM <- U$module == m; inS <- U$entrez %in% H[term == s, entrez]
    a <- sum(inM & inS); b <- sum(inM & !inS); c <- sum(!inM & inS); d <- sum(!inM & !inS)
    data.table(network = net, module = m, set = s, overlap = a, set_in_network = a + c,
               p = fisher.test(matrix(c(a, c, b, d), 2), alternative = "greater")$p.value)
  }))))
  res[, padj := p.adjust(p, "BH")][]
}
enr <- rbind(fisher_tab(col1, "round-1 (5,000 genes)"), fisher_tab(col2, "sensitivity (17,405 genes)"))
enr_f <- enr[set %in% FOCAL & padj < 0.05][order(network, padj)]
fe <- as.data.table(f3$enrichment)[set %in% FOCAL & padj < 0.05]
log_msg("  focal enrichment (BH<0.05): mine %d module-set pairs, fig3alt %d; overlaps identical: %s", nrow(enr_f), nrow(fe),
        isTRUE(all.equal(sort(paste(enr_f$module, enr_f$set, enr_f$overlap)), sort(paste(fe$module, fe$set, fe$overlap)))))
print(enr_f[, .(network, module, set = sub("HALLMARK_", "", set), overlap, set_in_network, padj = signif(padj, 2))])

# MYC/OXPHOS (pink) eigengene per library
pinkME <- ME2[, "MEpink"]
pink_lib <- data.table(library = cd$short, group = cd$label, CM = cd$CM, pink_eigengene = round(pinkME, 3),
                       pct_reads_outside_exons = pct_out, gc_slope = round(gcslope, 3), replicate_factor_B = batchB)
print(pink_lib)

# =============================================================================
# C. GC confound of the candidate panels; GC-corrected and leave-one-out GSEA
# =============================================================================
log_msg("C: GC content of displayed genes and GC-corrected re-analysis")
expr_ids <- rownames(V)
gc_rows <- rbindlist(lapply(names(PANELS), function(pn) {
  ids <- disp[[pn]]; le_all <- le_tab[panel == pn]$gene_id
  sc <- colMeans(V[ids, S] - rowMeans(V[ids, S]))
  data.table(panel = pn, n_displayed = length(ids),
             median_GC_displayed = round(median(GC[ids], na.rm = TRUE), 3), median_GC_all_leading_edge = round(median(GC[le_all], na.rm = TRUE), 3),
             median_GC_expressed = round(median(GC[expr_ids], na.rm = TRUE), 3),
             median_GC_percentile_displayed = round(median(ecdf(GC[expr_ids])(GC[ids])), 2),
             wilcox_p_displayed_vs_expressed = signif(wilcox.test(GC[ids], GC[setdiff(expr_ids, ids)])$p.value, 2),
             r_score_vs_gc_slope_12lib = round(scor(sc, gcslope), 2), r_score_vs_gc_slope_within_group = round(scor(gcen(sc), gcen(gcslope)), 2),
             r_score_vs_pct_out_within_group = round(scor(gcen(sc), gcen(pct_out)), 2),
             score_by_library = paste(sprintf("%s %.2f", cd$short, sc), collapse = "; "))
}))
print(gc_rows[, -"score_by_library"])
panel_scores <- as.data.table(sapply(names(PANELS), function(pn) colMeans(V[disp[[pn]], S] - rowMeans(V[disp[[pn]], S]))))
panel_scores[, `:=`(sample_id = S, library = cd$short, group = cd$label, gc_slope = gcslope, pct_out = pct_out)]

# --- DESeq2 refits: plain (validation), GC-corrected (two variants), leave-one-library-out
dds0 <- core$fc$dds_full
cnt  <- counts(dds0)[, S]; sf <- sizeFactors(dds0)[S]
stopifnot(identical(colnames(cnt), S))
coldat <- data.frame(row.names = S, NKX = factor(as.character(cd$NKX), levels = c("vector", "NKX3.1")),
                     CM = factor(as.character(cd$CM), levels = c("ctrl", "CM")))
gc_nf <- function(reference = c("overall", "group")) {
  reference <- match.arg(reference)
  g <- GC[rownames(cnt)]; g[is.na(g)] <- median(g, na.rm = TRUE)
  nm <- sweep(cnt, 2, sf, "/"); y <- log2(nm + 1)
  dev <- if (reference == "overall") y - rowMeans(y) else y - t(apply(y, 1, function(v) ave(v, cd$group)))
  train <- rowMeans(nm) >= 50 & !is.na(GC[rownames(cnt)])
  rng <- range(g[train])
  fit <- sapply(S, function(s) {
    lo <- loess(dev[train, s] ~ g[train], span = 0.3, degree = 1, family = "symmetric")
    predict(lo, newdata = pmin(pmax(g, rng[1]), rng[2]))
  })
  nf <- sweep(2^fit, 2, sf, "*")
  nf <- nf / exp(rowMeans(log(nf)))
  list(nf = nf, fit = fit)
}
fit_contrasts <- function(samples, nf = NULL) {
  d <- DESeqDataSetFromMatrix(cnt[, samples], coldat[samples, ], ~ NKX + CM + NKX:CM)
  if (!is.null(nf)) normalizationFactors(d) <- nf[, samples]
  d <- DESeq(d, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  list(CM_in_vector = results(d, name = "CM_CM_vs_ctrl"),
       NKX_in_ctrl  = results(d, name = "NKX_NKX3.1_vs_vector"),
       NKX_in_CM    = results(d, contrast = list(c("NKX_NKX3.1_vs_vector", "NKXNKX3.1.CMCM"))))
}
ranks_of <- function(r) {
  d <- data.table(gene_id = rownames(r), stat = r$stat)[, entrez := ann[gene_id, "entrez"]][!is.na(entrez) & is.finite(stat)]
  d <- d[order(-abs(stat))][!duplicated(entrez)]
  sort(setNames(d$stat, d$entrez), decreasing = TRUE)
}
run_gsea <- function(r) {
  set.seed(42)
  g <- GSEA(geneList = ranks_of(r), TERM2GENE = T2G, exponent = 1, minGSSize = 15, maxGSSize = 500, eps = 0,
            pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, seed = FALSE, by = "fgsea")
  as.data.table(as.data.frame(g))[, .(ID, NES, p.adjust, core_enrichment)]
}
NFo <- gc_nf("overall"); NFg <- gc_nf("group")
log_msg("  DESeq2 refits (plain, GC overall, GC within-group, leave-one-out)")
FITS <- list(plain = fit_contrasts(S), gc_overall = fit_contrasts(S, NFo$nf), gc_within_group = fit_contrasts(S, NFg$nf),
             drop_12 = fit_contrasts(setdiff(S, "24004-12_S15")), drop_1 = fit_contrasts(setdiff(S, "24-004-1_S56")))
# validation: plain refit stats vs core results
val <- rbindlist(lapply(names(PANELS), function(pn) {
  cn <- PANELS[[pn]]$contrast; a <- FITS$plain[[cn]]; b <- as.data.frame(core$fc$res[[cn]]); rownames(b) <- b$gene_id
  data.table(contrast = cn, max_abs_stat_diff = max(abs(a$stat - b[rownames(a), "stat"]), na.rm = TRUE))
}))
print(val)
gsea_rows <- list()
for (fn in names(FITS)) for (pn in names(PANELS)) {
  cn <- PANELS[[pn]]$contrast
  if (fn == "drop_12" && pn != "C") next
  if (fn == "drop_1" && pn == "C") next
  g <- run_gsea(FITS[[fn]][[cn]])
  r <- FITS[[fn]][[cn]]
  for (s in PANELS[[pn]]$sets) {
    x <- g[ID == s]
    gsea_rows[[length(gsea_rows) + 1]] <- data.table(analysis = fn, panel = pn, contrast = cn, set = s, NES = round(x$NES, 3),
      FDR = signif(x$p.adjust, 3), n_sets_FDR_lt_0.05 = sum(g$p.adjust < 0.05),
      median_log2FC_displayed = round(median(r[intersect(disp[[pn]], rownames(r)), "log2FoldChange"], na.rm = TRUE), 3))
  }
}
gsea_cmp <- rbindlist(gsea_rows)
pub <- rbindlist(lapply(names(PANELS), function(pn) {
  cn <- PANELS[[pn]]$contrast
  g <- fread(file.path(revd, sprintf("gsea/tables/gsea_fc_hallmark_%s.tsv", cn)))
  g[ID %in% PANELS[[pn]]$sets, .(panel = pn, set = ID, NES_published_corrected = round(NES, 3), FDR_published_corrected = signif(p.adjust, 3))]
}))
gsea_cmp <- merge(gsea_cmp, pub, by = c("panel", "set"), all.x = TRUE)
print(gsea_cmp[, .(analysis, panel, set = sub("HALLMARK_", "", set), NES, FDR, NES_published_corrected, median_log2FC_displayed, n_sets_FDR_lt_0.05)])

# --- GC-corrected log-expression for rotation tests: VST minus the per-library GC trend (overall reference)
gc_trend_vst <- function(M) {
  g <- GC[rownames(M)]; g[is.na(g)] <- median(g, na.rm = TRUE)
  dev <- M - rowMeans(M); train <- rowMeans(NORM[rownames(M), S]) >= 50 & !is.na(GC[rownames(M)])
  rng <- range(g[train])
  fit <- sapply(colnames(M), function(s) predict(loess(dev[train, s] ~ g[train], span = 0.3, degree = 1, family = "symmetric"),
                                                 newdata = pmin(pmax(g, rng[1]), rng[2])))
  M - fit
}
Vgc <- gc_trend_vst(V)

# =============================================================================
# D. selection-matched null with the corrected construction
# =============================================================================
log_msg("D: selection-matched null (corrected construction)")
N_ROT <- 2000L
consistency <- function(Xref, Xtrt) {
  mr <- rowMeans(Xref); mt_ <- rowMeans(Xtrt); d <- sign(mt_ - mr)
  trt_ok <- rowSums(sign(Xtrt - mr) == d) == ncol(Xtrt)
  ref_ok <- rowSums(sign(Xref - mt_) == -d) == ncol(Xref)
  sep <- ifelse(d < 0, rowMaxs(Xtrt) < rowMins(Xref), rowMins(Xtrt) > rowMaxs(Xref))
  c(both_agree = 100 * mean(trt_ok & ref_ok), complete_sep = 100 * mean(sep), median_abs_diff = median(abs(mt_ - mr)))
}
select_top <- function(by_set, n = 30L) {
  if (length(by_set) == 1) return(head(by_set[[1]], n))
  k <- ceiling(n / length(by_set))
  repeat { sel <- unique(unlist(lapply(by_set, head, k))); if (length(sel) >= n || k >= max(lengths(by_set))) break; k <- k + 1 }
  sel
}
Dm <- model.matrix(~ 0 + cd$group); Pm <- Dm %*% solve(crossprod(Dm)) %*% t(Dm)
eB <- eigen(diag(12) - Pm, symmetric = TRUE); B8 <- eB$vectors[, eB$values > 0.5]; stopifnot(ncol(B8) == 8)
haar <- function(n) { q <- qr(matrix(rnorm(n * n), n)); Q <- qr.Q(q); Q %*% diag(sign(diag(qr.R(q)))) }
helm3 <- cbind(c(1, -1, 0) / sqrt(2), c(1, 1, -2) / sqrt(6))
null_rows <- list()
for (pn in names(PANELS)) {
  p <- PANELS[[pn]]
  rk <- as.data.table(core$fc$res[[p$contrast]])[, .(gene_id, stat)][, entrez := ann[gene_id, "entrez"]][!is.na(entrez) & is.finite(stat)]
  rk <- rk[order(-abs(stat))][!duplicated(entrez)]
  set_ids <- lapply(p$sets, function(s) rk[entrez %in% H[term == s, entrez]]$gene_id)
  U <- unique(unlist(set_ids))
  ref_s <- S[cd$group == p$ref]; trt_s <- S[cd$group == p$trt]; six <- c(ref_s, trt_s)
  XU <- V[U, S]; Zr <- XU %*% B8
  W4 <- rbind(cbind(helm3, matrix(0, 3, 2)), cbind(matrix(0, 3, 2), helm3))       # 6 x 4 within-group basis (ref rows, trt rows)
  q6 <- c(rep(-1, 3), rep(1, 3)) / sqrt(6)
  mu6 <- rowMeans(XU[, six])
  sel_t <- function(X6, s2) {
    tt <- (rowMeans(X6[, 4:6]) - rowMeans(X6[, 1:3])) / sqrt(s2 * 2 / 3)
    select_top(lapply(set_ids, function(ids) ids[order(tt[ids])]))
  }
  X6o <- XU[, six]; s2o <- rowSums(Zr^2) / 8
  so <- sel_t(X6o, s2o); obs <- consistency(X6o[so, 1:3], X6o[so, 4:6])
  # observed for the displayed genes (DESeq2 |Wald| selection)
  od <- consistency(V[disp[[pn]], ref_s], V[disp[[pn]], trt_s])
  set.seed(SEED + match(pn, names(PANELS)))
  nul <- t(vapply(seq_len(N_ROT), function(b) {
    Zs <- Zr %*% haar(8)
    X6 <- mu6 + outer(Zs[, 1], q6) + Zs[, 2:5] %*% t(W4)
    s2 <- rowSums(Zs[, 2:8]^2) / 7
    sl <- sel_t(X6, s2)
    consistency(X6[sl, 1:3, drop = FALSE], X6[sl, 4:6, drop = FALSE])
  }, numeric(3)))
  fnull <- as.data.table(f3$consistency)[panel == pn]
  q95 <- function(v) sprintf("%.0f (%.0f-%.0f)", mean(v), quantile(v, .025), quantile(v, .975))
  null_rows[[pn]] <- data.table(panel = pn, universe = length(U),
    observed_displayed_both_agree = round(od[1], 1), observed_tsel_both_agree = round(obs[1], 1),
    null_both_agree_mean_95_corrected = q95(nul[, 1]), null_both_agree_mean_95_fig3alt = fnull$null_pct_both_agree_mean_95,
    frac_null_ge_observed_tsel_both_agree = round(mean(nul[, 1] >= obs[1]), 3),
    observed_tsel_complete_sep = round(obs[2], 1), null_complete_sep_mean_95_corrected = q95(nul[, 2]),
    null_complete_sep_mean_95_fig3alt = fnull$null_pct_complete_separation_mean_95,
    frac_null_ge_observed_complete_sep = round(mean(nul[, 2] >= obs[2]), 3),
    observed_tsel_median_abs_diff = round(obs[3], 3),
    null_median_abs_diff_mean_95_corrected = sprintf("%.3f (%.3f-%.3f)", mean(nul[, 3]), quantile(nul[, 3], .025), quantile(nul[, 3], .975)),
    null_median_abs_diff_mean_95_fig3alt = fnull$null_median_abs_diff_mean_95,
    frac_null_ge_observed_median_abs_diff_corrected = round(mean(nul[, 3] >= obs[3]), 3),
    frac_null_ge_observed_median_abs_diff_fig3alt = fnull$null_frac_ge_observed_median_abs_diff)
}
null_cmp <- rbindlist(null_rows)
print(t(null_cmp))

# =============================================================================
# E. correlation-aware set tests (all 50 Hallmark sets), recomputed; fry on GC-corrected VST
# =============================================================================
log_msg("E: camera / fry recomputed")
dY <- data.table(id = rownames(V), entrez = ann[rownames(V), "entrez"], m = rowMeans(V))[!is.na(entrez)][order(-m)][!duplicated(entrez)]
Y <- V[dY$id, S]; Ygc <- Vgc[dY$id, S]
idx <- lapply(split(H$entrez, H$term), function(s) which(dY$entrez %in% s)); idx <- idx[lengths(idx) >= 15]
gfac <- factor(as.character(cd$group), levels = c("vector_ctrl", "vector_CM", "NKX3.1_ctrl", "NKX3.1_CM"))
DES <- model.matrix(~ 0 + gfac); colnames(DES) <- levels(gfac)
CON <- list(CM_in_vector = c(-1, 1, 0, 0), NKX_in_ctrl = c(-1, 0, 1, 0), NKX_in_CM = c(0, -1, 0, 1),
            CM_in_NKX = c(0, 0, -1, 1), CM_main = c(-.5, .5, -.5, .5), NKX_main = c(-.5, -.5, .5, .5))
st_rows <- list(); foc_rows <- list()
for (cn in names(CON)) {
  fr <- fry(Y, idx, DES, CON[[cn]]); ce <- camera(Y, idx, DES, CON[[cn]], inter.gene.cor = NA); c1 <- camera(Y, idx, DES, CON[[cn]], inter.gene.cor = 0.01)
  fg <- fry(Ygc, idx, DES, CON[[cn]])
  st_rows[[cn]] <- data.table(contrast = cn, n_sets = length(idx), fry_FDR_lt_0.05 = sum(fr$FDR < 0.05), camera_est_FDR_lt_0.05 = sum(ce$FDR < 0.05),
                              camera_0.01_FDR_lt_0.05 = sum(c1$FDR < 0.05), fry_gc_corrected_FDR_lt_0.05 = sum(fg$FDR < 0.05),
                              min_fry_FDR = signif(min(fr$FDR), 2), min_camera_est_FDR = signif(min(ce$FDR), 2))
  foc_rows[[cn]] <- data.table(contrast = cn, set = unname(FOCAL), fry_p = signif(fr[unname(FOCAL), "PValue"], 3), fry_FDR = signif(fr[unname(FOCAL), "FDR"], 3),
                               fry_dir = fr[unname(FOCAL), "Direction"], camera_est_p = signif(ce[unname(FOCAL), "PValue"], 3),
                               camera_est_r = round(ce[unname(FOCAL), "Correlation"], 3),
                               fry_gc_corrected_p = signif(fg[unname(FOCAL), "PValue"], 3), fry_gc_corrected_dir = fg[unname(FOCAL), "Direction"])
}
st_cmp <- merge(rbindlist(st_rows), as.data.table(f3$n_sig_methods)[, .(contrast, fig3alt_fry = fry, fig3alt_camera_est = camera_estimated, fig3alt_camera_0.01 = camera_r0.01)], by = "contrast")
print(st_cmp)
foc <- merge(rbindlist(foc_rows), as.data.table(f3$set_tests)[, .(contrast, set, fig3alt_fry_p = signif(fry_p, 3), fig3alt_roast_p = signif(roast_mean_9999rot_p, 3))],
             by = c("contrast", "set"), all.x = TRUE)
print(foc[contrast %in% c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM", "CM_main")][order(contrast, set)])

# =============================================================================
# tables, session
# =============================================================================
wtsv(me_check, "v02_B_eigengene_recomputation_check.tsv")
wtsv(rerun_check, "v02_B_sensitivity_network_rerun_check.tsv")
wtsv(mt, "v02_B_module_design_tests.tsv")
wtsv(chk, "v02_B_exact_permutation_vs_fig3alt.tsv")
wtsv(enr_f, "v02_B_module_focal_enrichment.tsv")
wtsv(pink_lib, "v02_B_pink_module_by_library.tsv")
wtsv(gc_rows, "v02_C_panel_gc_confound.tsv")
wtsv(panel_scores, "v02_C_panel_scores_by_library.tsv")
wtsv(val, "v02_C_refit_validation.tsv")
wtsv(gsea_cmp, "v02_C_gsea_gc_corrected_and_leave_one_out.tsv")
wtsv(null_cmp, "v02_D_selection_null_corrected_vs_fig3alt.tsv")
wtsv(st_cmp, "v02_E_set_tests_n_significant.tsv")
wtsv(foc, "v02_E_set_tests_focal.tsv")
writeLines(c(sprintf("50_verify_fig3_02_recompute.R finished %s (%.1f min)", format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "mins"))),
             sprintf("msigdbr %s, MSigDB %s", as.character(packageVersion("msigdbr")), paste(unique(H_tab$db_version), collapse = ",")), "",
             capture.output(sessionInfo())), file.path(out, "sessionInfo_50_verify_fig3_02.txt"))
log_msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
