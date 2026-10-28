#!/usr/bin/env Rscript
# =============================================================================
# 14_modules.R -- Reviewer 1, major point 5: formal tests of the
# "NKX3.1 co-expression module" claims (manuscript Fig. 3D/E and the Results
# subsection "NKX3.1 co-expression modules exhibit context-dependent rewiring").
#
# Claims tested
#   (i)   basal: NKX3.1 is co-expressed with a module incl. SLX1A, BCL11B, ZNF331
#   (ii)  inflammation: the module shifts to a distinct one incl. TPGS1 and SRMS
#   (iii) BCL11B stays associated with NKX3.1 across conditions
#   (iv)  LRIG1 is induced (by NKX3.1) under inflammation
#
# Every test is run on two quantifications of the same 12 libraries:
#   corrected = STAR genome alignments counted with featureCounts (unique reads,
#               reverse-stranded, Ensembl 99 exons); VST = core$fc$vst
#   original  = nf-core STAR->Salmon, tximport lengthScaledTPM, blind VST used
#               for manuscript Fig 1B-D/3 = core$vst_sal_v15
#
# Steps
#   0. Reproduce the original Welch-z heatmaps (Fig 3B/C; script
#      rna_downstreams_salmon.R) and define each claimed module formally as the
#      dendrogram clade containing NKX3-1 whose leaves are exactly the boxed rows.
#   1. Fig-3 view: correlation with NKX3-1 across the 6 samples of one comparison;
#      exact permutation p over the 20 balanced 3-vs-3 relabelings.
#   2. Co-expression beyond group membership
#      (a) dose-response among the 6 NKX3.1-transfected samples:
#          lm(gene ~ NKX3-1 + CM) on VST (3 residual df) + exact within-CM permutation
#      (b) partial correlation with NKX3-1 across all 12 samples after removing
#          the 4 group means; NKX3-1 residuals permuted within groups
#          (10,000 Monte-Carlo draws, plus exact enumeration of all 6^4 = 1,296).
#   3. Module coherence vs 10,000 random gene sets matched on mean-expression decile.
#   4. WGCNA (signed, bicor) on the 5,000 most variable genes.
#   5. Read support, NKX3.1 x CM interaction ("rewiring") and a claim-by-claim verdict.
#
# Sign convention everywhere: treatment minus reference (NKX3.1 - vector; CM - control).
# Run: /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/14_modules.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(WGCNA)
})
options(stringsAsFactors = FALSE, warn = 1)
invisible(allowWGCNAThreads(16))

proj   <- "/mnt/ssd1/Projects/NKX3_wd"
out    <- file.path(proj, "results/revision_2026/modules")
work   <- "/mnt/ssd1/Projects/NKX3_wd/work/modules"
v15    <- file.path(proj, "results/rna_downstreams_out_v15")
salmon_merged <- "/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon/salmon.merged.gene_counts.tsv"
dir.create(out,  recursive = TRUE, showWarnings = FALSE)
dir.create(work, recursive = TRUE, showWarnings = FALSE)

SEED   <- 42L
N_PERM <- 10000L                 # Monte-Carlo permutations / null gene sets
NKX_ID <- "ENSG00000167034"      # NKX3-1
set.seed(SEED)

# ----------------------------------------------------------------------------- helpers
wtsv <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                   row.names = FALSE, na = "NA")
bh <- function(p) p.adjust(p, method = "BH")
rowcor <- function(X, y) {                      # Pearson r of every row of X with y
  Xc <- X - rowMeans(X); yc <- y - mean(y)
  as.vector(Xc %*% yc) / sqrt(rowSums(Xc^2) * sum(yc^2))
}
rowranks <- function(X) t(apply(X, 1, rank))
center_by <- function(X, g) {                   # subtract group means row-wise
  g <- droplevels(factor(g)); G <- model.matrix(~ 0 + g)
  X - (X %*% G %*% diag(1 / colSums(G), ncol(G))) %*% t(G)
}
perms <- function(v) if (length(v) <= 1) matrix(v, 1) else
  do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], perms(v[-i]))))
within_group_perms <- function(g) {             # every within-group shuffle; column 1 = identity
  idx <- split(seq_along(g), droplevels(factor(g)))
  P <- lapply(idx, perms)
  grid <- expand.grid(lapply(P, function(m) seq_len(nrow(m))))
  sapply(seq_len(nrow(grid)), function(k) {
    o <- seq_along(g)
    for (j in seq_along(idx)) o[idx[[j]]] <- P[[j]][grid[k, j], ]
    o
  })
}
fmt_p <- function(p) trimws(ifelse(is.na(p), "NA", ifelse(p < 1e-3, formatC(p, format = "e", digits = 1),
                                                          formatC(signif(p, 2), format = "fg", digits = 2))))

# ----------------------------------------------------------------------------- data
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
cd   <- core$cd; ann <- core$ann
S    <- cd$sample_id
grp  <- droplevels(cd$group)
is_nkx <- cd$NKX == "NKX3.1"; is_cm <- cd$CM == "CM"
symbol_of <- function(id) ann[id, "symbol"]

# =============================================================================
# 0. Reproduce the original heatmaps and define the claimed modules
# =============================================================================
st  <- read.csv(file.path(proj, "results/sample_table.csv"))
V0  <- core$vst_sal_v15[, st$sample_id]                   # original column order
ids0 <- rownames(V0); s0 <- ann[ids0, "symbol"]
rn0 <- ifelse(!is.na(s0) & s0 != "", s0, ids0)             # original pipeline: symbol, else Ensembl ID
cond0 <- st$condition

orig_heatmap <- function(g1, g2) {                         # = run_comparison() PLOT 3 in rna_downstreams_salmon.R
  i1 <- which(cond0 == g1); i2 <- which(cond0 == g2)
  mu1 <- rowMeans(V0[, i1]); mu2 <- rowMeans(V0[, i2])
  v1 <- rowVars(V0[, i1]); v2 <- rowVars(V0[, i2])
  lfc <- mu1 - mu2
  z   <- lfc / sqrt(pmax(v1, 1e-6) / length(i1) + pmax(v2, 1e-6) / length(i2))
  sig <- pmax(mu1, mu2) >= 2.0 & abs(lfc) >= 0.25 & abs(z) >= 2.0
  o   <- order(-abs(z), -abs(lfc)); so <- o[sig[o]]
  keep <- unique(rn0[so]); keep <- keep[seq_len(min(200, length(keep)))]
  rows <- which(rn0 %in% keep)
  mat <- t(scale(t(V0[rows, c(i1, i2)])))
  rownames(mat) <- make.unique(rn0[rows])
  hc <- hclust(dist(mat), method = "complete")              # pheatmap default: euclidean, complete
  list(hc = hc, ids = setNames(ids0[rows], rownames(mat)), n_sig = sum(sig),
       welch = data.frame(gene_id = ids0, lfc_vst = lfc, z = z))
}
pdf_row_labels <- function(pdf) {                          # row labels of the published heatmap PDF, top to bottom
  tmp <- file.path(work, paste0(basename(pdf), ".bbox.html"))
  system2("pdftotext", c("-bbox", shQuote(pdf), shQuote(tmp)))
  x <- grep("<word ", readLines(tmp, warn = FALSE, encoding = "UTF-8"), value = TRUE)
  xmin <- sub('.*xMin="([0-9.]+)".*', "\\1", x)
  ymin <- as.numeric(sub('.*yMin="([0-9.]+)".*', "\\1", x))
  w <- gsub("−", "-", sub('.*">(.*)</word>.*', "\\1", x))
  col <- names(which.max(table(xmin)))                     # the 200 row labels share one x position
  sel <- xmin == col & w != "condition"
  unlink(tmp)
  w[sel][order(ymin[sel])]
}
clade_walk <- function(hc, leaf) {                         # ancestors of a leaf, bottom-up
  m <- hc$merge; lab <- hc$labels
  leaves <- function(node) if (node < 0) lab[-node] else c(leaves(m[node, 1]), leaves(m[node, 2]))
  cur <- -which(lab == leaf); res <- list()
  repeat {
    node <- which(m[, 1] == cur | m[, 2] == cur)
    if (!length(node)) break
    res[[length(res) + 1]] <- list(node = node, height = hc$height[node], leaves = leaves(node))
    cur <- node
  }
  res
}

boxed <- list(   # rows inside the dashed boxes of Fig 3D and 3E, top to bottom
  basal = c("EHD1", "SLX1A", "NKX3-1", "BCL11B", "SSC4D", "ZNF331", "NPHP4", "KIAA0895", "AGAP6",
            "TRPM4", "D2HGDH", "ZC3H12A", "KLHL3", "ULK4", "TPGS1", "FANCA", "TRPV3", "DIPK1B",
            "MCTS1", "IGLL5", "SRMS", "ARHGEF19"),
  inflammation = c("DIPK1B", "SSC4D", "FRMD4A", "IGLL5", "CNIH4", "LRIG1", "SRMS", "CCN1", "NKX3-1",
                   "BCL11B", "TPGS1", "PPIP5K1", "HSPA6", "KLRC3"))
hm_spec <- list(basal = c("HM-NKX3.1-C", "HM-vector-C"), inflammation = c("HM-NKX3.1-CM", "HM-vector-CM"))
ctx_label <- c(basal = "basal (Control-basal vs NKX3.1-basal; Fig 3D)",
               inflammation = "inflammation (Control-inflammation vs NKX3.1-inflammation; Fig 3E)")

modules <- list(); module_def <- list(); hm_store <- list()
for (m in names(hm_spec)) {
  h  <- orig_heatmap(hm_spec[[m]][1], hm_spec[[m]][2])
  pl <- pdf_row_labels(file.path(v15, sprintf("cmp_%s_vs_%s_heatmap_sig_genes.pdf", hm_spec[[m]][1], hm_spec[[m]][2])))
  ord <- h$hc$labels[h$hc$order]
  if (!identical(ord, pl)) stop("Reproduced row order differs from the published heatmap for ", m)
  walk <- clade_walk(h$hc, "NKX3-1")
  hit  <- which(vapply(walk, function(w) setequal(w$leaves, boxed[[m]]), logical(1)))
  if (length(hit) != 1) stop("No dendrogram clade equals the boxed rows for ", m)
  cl <- walk[[hit]]
  h_next <- walk[[hit + 1]]$height
  k_mid  <- length(unique(cutree(h$hc, h = (cl$height + h_next) / 2)))
  leaves_ord <- ord[ord %in% cl$leaves]
  ids <- unname(h$ids[leaves_ord])
  modules[[m]] <- setNames(ids, leaves_ord)
  wl <- h$welch[match(ids, h$welch$gene_id), ]
  module_def[[m]] <- data.frame(
    module = m, heatmap = sprintf("cmp_%s_vs_%s_heatmap_sig_genes.pdf", hm_spec[[m]][1], hm_spec[[m]][2]),
    heatmap_row = match(leaves_ord, ord), symbol = leaves_ord, gene_id = ids,
    orig_welch_lfc_vst_NKX_minus_ctrl = round(wl$lfc_vst, 3), orig_welch_z = round(wl$z, 2),
    clade_node = cl$node, clade_merge_height = round(cl$height, 3),
    parent_merge_height = round(h_next, 3), n_clusters_if_cut_inside_clade_range = k_mid,
    n_welch_sig_genes = h$n_sig, heatmap_rows = length(ord),
    walk_up_from_NKX3_1 = paste(vapply(walk[seq_len(hit + 1)], function(w)
      sprintf("%d genes@h=%.2f", length(w$leaves), w$height), ""), collapse = " > "))
  hm_store[[m]] <- list(order = ord, ids = h$ids)
  message(sprintf("[module %s] published heatmap reproduced exactly (%d rows); boxed module = clade node %d, %d genes, merge height %.3f (parent %.3f)",
                  m, length(ord), cl$node, length(ids), cl$height, h_next))
}
module_def <- do.call(rbind, module_def)
union_ids <- unique(c(modules$basal, modules$inflammation))
names(union_ids) <- symbol_of(union_ids)
shared_ids <- intersect(modules$basal, modules$inflammation)
gene_class <- setNames(ifelse(union_ids %in% shared_ids, "both",
                              ifelse(union_ids %in% modules$basal, "Fig3D_only", "Fig3E_only")), union_ids)
module_def$fc_unique_reads_total <- rowSums(core$counts_fc[module_def$gene_id, ])
module_def$fc_quantifiable <- module_def$gene_id %in% rownames(core$fc$vst)
module_def$in_both_modules <- module_def$gene_id %in% shared_ids
wtsv(module_def, "modules_definition.tsv")

# =============================================================================
# 1. Expression matrices
# =============================================================================
vst_all_fc <- function() {   # VST of ALL featureCounts genes with exactly the transform behind core$fc$vst
  d  <- core$fc$dds_full
  bm <- rowMeans(counts(d, normalized = TRUE)); sub <- d[bm > 5, ]; bm <- bm[bm > 5]
  o  <- order(bm); sub <- sub[o[round(seq(from = 1, to = length(o), length = 1000))], ]   # = DESeq2::vst() subsetting
  sub <- estimateDispersionsGeneEst(sub, quiet = TRUE)
  sub <- estimateDispersionsFit(sub, fitType = "parametric", quiet = TRUE)
  dall <- DESeqDataSetFromMatrix(core$counts_fc[, colnames(d)], colData(d), design(d))
  sizeFactors(dall) <- sizeFactors(d)
  suppressMessages(dispersionFunction(dall) <- dispersionFunction(sub))
  va <- assay(varianceStabilizingTransformation(dall, blind = FALSE))
  dev <- max(abs(va[rownames(core$fc$vst), ] - core$fc$vst))
  if (dev > 1e-8) stop("all-gene VST does not reproduce core$fc$vst (max diff ", dev, ")")
  va[, S]
}
Q <- list(
  corrected = list(V = core$fc$vst[, S], name = "Corrected (featureCounts, unique reads)"),
  original  = list(V = core$vst_sal_v15[, S], name = "Original (Salmon; manuscript)"))
vst_fc_all <- vst_all_fc()

# library-level check: reads outside annotated exons per sample (featureCounts summaries of the recount)
fc_sum <- do.call(rbind, lapply(S, function(s) {
  x <- read.delim(file.path(proj, "results/revision_checks_2026-10", sprintf("fc_%s.txt.summary", s)), check.names = FALSE)
  v <- setNames(x[[2]], x[[1]]); tot <- sum(v)
  data.frame(sample_id = s, group = cd$label[cd$sample_id == s], total = tot,
             assigned_pct = round(100 * v[["Assigned"]] / tot, 1),
             no_feature_pct = round(100 * v[["Unassigned_NoFeatures"]] / tot, 1),
             multimapping_pct = round(100 * v[["Unassigned_MultiMapping"]] / tot, 1),
             ambiguity_pct = round(100 * v[["Unassigned_Ambiguity"]] / tot, 1))
}))
wtsv(fc_sum, "library_featureCounts_assignment.tsv")
nofeat_rng <- function(lab) sprintf("%.1f-%.1f%%", min(fc_sum$no_feature_pct[fc_sum$group == lab]), max(fc_sum$no_feature_pct[fc_sum$group == lab]))
nofeat_other <- sprintf("%.1f-%.1f%%", min(fc_sum$no_feature_pct[fc_sum$group != "NKX3.1-inflammation"]), max(fc_sum$no_feature_pct[fc_sum$group != "NKX3.1-inflammation"]))

message(sprintf("[data] corrected VST: %d genes; original VST: %d genes; NKX3-1 VST range among NKX3.1-transfected: corrected %.2f-%.2f, original %.2f-%.2f",
                nrow(Q$corrected$V), nrow(Q$original$V),
                min(Q$corrected$V[NKX_ID, is_nkx]), max(Q$corrected$V[NKX_ID, is_nkx]),
                min(Q$original$V[NKX_ID, is_nkx]), max(Q$original$V[NKX_ID, is_nkx])))

# =============================================================================
# 2. Test 1 -- Fig 3 view: correlation with NKX3-1 across the 6 samples of one comparison
# =============================================================================
splits <- combn(6, 3)                                       # the 20 balanced 3-vs-3 relabelings
test1_genomewide <- function(V, ctx) {
  s <- which(cd$CM == ifelse(ctx == "basal", "ctrl", "CM"))
  X <- V[, s]; X <- X[rowVars(X) > 0, , drop = FALSE]
  nk <- X[NKX_ID, ]; ind <- as.numeric(is_nkx[s])
  r_p <- rowcor(X, nk); r_s <- rowcor(rowranks(X), rank(nk)); r_g <- rowcor(X, ind)
  Rperm <- vapply(seq_len(ncol(splits)), function(j) { v <- numeric(6); v[splits[, j]] <- 1; rowcor(X, v) },
                  numeric(nrow(X)))
  p20 <- rowMeans(abs(Rperm) >= abs(r_g) - 1e-9)
  data.frame(gene_id = rownames(X), context = ctx, r_pearson = r_p, r_spearman = r_s,
             r_group_indicator = r_g, p_perm20 = p20)
}
t1_all <- list(); t1_summary <- list()
for (q in names(Q)) for (ctx in c("basal", "inflammation")) {
  t <- test1_genomewide(Q[[q]]$V, ctx); t$quant <- q
  t1_all[[paste(q, ctx)]] <- t
  tt <- t[t$gene_id != NKX_ID, ]
  r_crit <- 0.9; tcrit <- r_crit * sqrt(4 / (1 - r_crit^2))
  t1_summary[[paste(q, ctx)]] <- data.frame(
    quant = q, context = ctx, n_genes = nrow(tt),
    cor_r_NKX3_1_vs_r_group = round(cor(tt$r_pearson, tt$r_group_indicator, use = "complete.obs"), 4),
    median_abs_diff_r = round(median(abs(tt$r_pearson - tt$r_group_indicator), na.rm = TRUE), 4),
    min_p_perm20 = min(tt$p_perm20, na.rm = TRUE),
    n_genes_p_perm20_eq_0.1 = sum(tt$p_perm20 <= 0.1 + 1e-9, na.rm = TRUE),
    frac_abs_r_ge_0.9 = round(mean(abs(tt$r_pearson) >= r_crit, na.rm = TRUE), 4),
    n_abs_r_ge_0.9 = sum(abs(tt$r_pearson) >= r_crit, na.rm = TRUE),
    expected_frac_abs_r_ge_0.9_if_independent = round(2 * pt(-tcrit, df = 4), 4))
}
t1_summary <- do.call(rbind, t1_summary)
t1_mod <- do.call(rbind, lapply(names(t1_all), function(k) {
  t <- t1_all[[k]]; q <- t$quant[1]; ctx <- t$context[1]
  data.frame(quant = q, context = ctx, symbol = names(union_ids), gene_id = unname(union_ids),
             gene_class = gene_class[union_ids],
             in_module_of_this_context = unname(union_ids) %in% modules[[ctx]],
             t[match(union_ids, t$gene_id), c("r_pearson", "r_spearman", "r_group_indicator", "p_perm20")])
}))
wtsv(t1_mod, "test1_fig3view_correlation_module_genes.tsv")
wtsv(t1_summary, "test1_fig3view_genomewide_summary.tsv")

# =============================================================================
# 3. Test 2 -- co-expression with NKX3-1 beyond group membership
# =============================================================================
test2a_genomewide <- function(V) {                          # dose-response among NKX3.1-transfected samples
  s <- which(is_nkx); X <- V[, s]; X <- X[rowVars(X) > 0, , drop = FALSE]
  nk <- X[NKX_ID, ]; cmv <- as.numeric(is_cm[s])
  fit <- function(nkv) {
    D <- cbind(1, nkv, cmv); XtXi <- solve(crossprod(D))
    B <- X %*% D %*% XtXi; E <- X - B %*% t(D)
    list(b = B[, 2], se = sqrt(rowSums(E^2) / 3 * XtXi[2, 2]))
  }
  f  <- fit(nk); tt <- f$b / f$se
  PI <- within_group_perms(cmv)                              # 36 within-CM shuffles of NKX3-1
  Tp <- vapply(seq_len(ncol(PI)), function(k) { fk <- fit(nk[PI[, k]]); fk$b / fk$se }, numeric(nrow(X)))
  res <- data.frame(gene_id = rownames(X), dose_slope = f$b, dose_se = f$se, dose_t = tt,
                    dose_partial_r = tt / sqrt(tt^2 + 3), dose_p = 2 * pt(-abs(tt), df = 3),
                    dose_p_perm36 = rowMeans(abs(Tp) >= abs(tt) - 1e-9))
  res[res$gene_id == NKX_ID, -1] <- NA
  res$dose_padj_genomewide <- bh(res$dose_p)
  res
}
PI12 <- within_group_perms(grp)                             # 1,296 within-group shuffles (all 4 groups)
set.seed(SEED)
MC12 <- replicate(N_PERM, {                                  # 10,000 random within-group shuffles
  o <- seq_len(12)
  for (ix in split(seq_len(12), grp)) o[ix] <- ix[sample.int(length(ix))]
  o
})
test2b_genomewide <- function(V) {                          # partial correlation, 12 samples, group-centred
  R <- center_by(V, grp); R <- R[rowSums(R^2) > 1e-12, , drop = FALSE]
  Rz <- R / sqrt(rowSums(R^2)); nz <- R[NKX_ID, ] / sqrt(sum(R[NKX_ID, ]^2))
  r  <- as.vector(Rz %*% nz)
  Rp <- Rz %*% matrix(nz[PI12], nrow = 12)                  # genes x 1296
  p_exact <- rowMeans(abs(Rp) >= abs(r) - 1e-9)
  df <- 12 - 4 - 1
  tt <- r * sqrt(df / pmax(1 - r^2, 1e-12))
  res <- data.frame(gene_id = rownames(R), partial_r = r, partial_p_exact1296 = p_exact,
                    partial_p_parametric_df7 = 2 * pt(-abs(tt), df))
  res[res$gene_id == NKX_ID, -1] <- NA
  res$partial_padj_genomewide <- bh(res$partial_p_exact1296)
  res$partial_rank_abs_r <- rank(-abs(res$partial_r), na.last = "keep", ties.method = "min")
  res$partial_percentile <- 100 * (1 - (res$partial_rank_abs_r - 1) / sum(!is.na(res$partial_r)))
  attr(res, "Rz") <- Rz; attr(res, "nz") <- nz
  res
}
heldout_partial <- function(V, ids, m) {   # partial r with NKX3-1 in the 6 samples NOT used to define module m
  s <- which(cd$CM == ifelse(m == "basal", "CM", "ctrl"))
  R <- center_by(V[unique(c(NKX_ID, intersect(ids, rownames(V)))), s, drop = FALSE], grp[s])
  R <- R[rowSums(R^2) > 1e-12, , drop = FALSE]
  Rz <- R / sqrt(rowSums(R^2)); nz <- Rz[NKX_ID, ]
  r  <- as.vector(Rz %*% nz)
  PI <- within_group_perms(droplevels(grp[s]))              # 36 within-group shuffles
  Rp <- Rz %*% matrix(nz[PI], nrow = length(s))
  data.frame(gene_id = rownames(R), heldout_partial_r = r, heldout_p_exact36 = rowMeans(abs(Rp) >= abs(r) - 1e-9))
}
t2 <- list(); t2_gw <- list(); resid_store <- list()
for (q in names(Q)) {
  a <- test2a_genomewide(Q[[q]]$V)
  b <- test2b_genomewide(Q[[q]]$V)
  Rz <- attr(b, "Rz"); nz <- attr(b, "nz")
  # Monte-Carlo permutation p (10,000 within-group shuffles), as specified, for the claimed genes
  ids <- intersect(union_ids, rownames(Rz)); ids <- setdiff(ids, NKX_ID)
  rmc <- Rz[ids, , drop = FALSE] %*% matrix(nz[MC12], nrow = 12)
  robs <- b$partial_r[match(ids, b$gene_id)]
  p_mc <- (1 + rowSums(abs(rmc) >= abs(robs) - 1e-9)) / (N_PERM + 1)
  gw <- merge(a, b, by = "gene_id", all = TRUE)
  for (ctx in c("basal", "inflammation")) {
    t1x <- t1_all[[paste(q, ctx)]]
    gw[[paste0("r6_", ctx)]] <- t1x$r_pearson[match(gw$gene_id, t1x$gene_id)]
  }
  gw$symbol <- symbol_of(gw$gene_id); gw$quant <- q
  gw <- gw[order(gw$partial_p_exact1296), c("quant", "gene_id", "symbol", setdiff(names(gw), c("quant", "gene_id", "symbol")))]
  t2_gw[[q]] <- gw
  resid_store[[q]] <- list(R = center_by(Q[[q]]$V[unique(c(NKX_ID, intersect(union_ids, rownames(Q[[q]]$V)))), ], grp))
  for (m in names(modules)) {
    mid <- setdiff(modules[[m]], NKX_ID)
    x <- data.frame(quant = q, module = m, symbol = symbol_of(mid), gene_id = mid,
                    gene_class = gene_class[mid], quantifiable = mid %in% rownames(Q[[q]]$V))
    x <- cbind(x, gw[match(mid, gw$gene_id), setdiff(names(gw), c("quant", "gene_id", "symbol"))])
    x$partial_p_mc10000 <- p_mc[match(mid, ids)]
    x$dose_padj_module <- NA; x$partial_padj_module <- NA
    ok <- x$quantifiable
    x$dose_padj_module[ok]    <- bh(x$dose_p[ok])
    x$partial_padj_module[ok] <- bh(x$partial_p_mc10000[ok])
    ho <- heldout_partial(Q[[q]]$V, mid, m)
    x$heldout_samples <- ifelse(m == "basal", "inflammation context (not used to define the basal module)",
                                "basal context (not used to define the inflammation module)")
    x$heldout_partial_r <- ho$heldout_partial_r[match(mid, ho$gene_id)]
    x$heldout_p_exact36 <- ho$heldout_p_exact36[match(mid, ho$gene_id)]
    x$heldout_padj_module <- NA
    x$heldout_padj_module[ok] <- bh(x$heldout_p_exact36[ok])
    t2[[paste(q, m)]] <- x
  }
}
t2 <- do.call(rbind, t2); rownames(t2) <- NULL
wtsv(t2, "test2_beyond_group_module_genes.tsv")
for (q in names(t2_gw)) {
  gz <- gzfile(file.path(out, sprintf("test1_2_genomewide_%s.tsv.gz", q)), "w")
  write.table(t2_gw[[q]], gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(gz)
}
t2_gw_summary <- do.call(rbind, lapply(names(t2_gw), function(q) {
  g <- t2_gw[[q]]
  data.frame(quant = q, n_genes = sum(!is.na(g$partial_r)),
             dose_n_padj_lt_0.1 = sum(g$dose_padj_genomewide < 0.1, na.rm = TRUE),
             dose_min_p = min(g$dose_p, na.rm = TRUE),
             partial_n_padj_lt_0.1 = sum(g$partial_padj_genomewide < 0.1, na.rm = TRUE),
             partial_min_p_exact = min(g$partial_p_exact1296, na.rm = TRUE),
             partial_n_p_exact_lt_0.05 = sum(g$partial_p_exact1296 < 0.05, na.rm = TRUE),
             partial_expected_n_p_lt_0.05 = round(0.05 * sum(!is.na(g$partial_r))),
             partial_top10 = paste(head(g$symbol[order(g$partial_p_exact1296, -abs(g$partial_r))], 10), collapse = ","),
             critical_abs_r_dose_df3 = round(qt(0.975, 3) / sqrt(qt(0.975, 3)^2 + 3), 3),
             critical_abs_r_partial_df7 = round(qt(0.975, 7) / sqrt(qt(0.975, 7)^2 + 7), 3),
             NKX3_1_vst_range_transfected = sprintf("%.2f-%.2f", min(Q[[q]]$V[NKX_ID, is_nkx]), max(Q[[q]]$V[NKX_ID, is_nkx])),
             NKX3_1_within_group_sd_transfected = round(sqrt(mean(c(var(Q[[q]]$V[NKX_ID, is_nkx & !is_cm]),
                                                                    var(Q[[q]]$V[NKX_ID, is_nkx & is_cm])))), 3))
}))
wtsv(t2_gw_summary, "test2_genomewide_summary.tsv")

# =============================================================================
# 4. Test 3 -- module coherence vs expression-matched random gene sets
# =============================================================================
sample_sets <- list(
  NKX3.1_transfected_CMcentred = list(idx = which(is_nkx), center = cd$CM[is_nkx]),
  NKX3.1_transfected_raw       = list(idx = which(is_nkx), center = NULL),
  all12_groupcentred           = list(idx = seq_len(12), center = grp))
coherence <- function(V, partners, set, seed) {
  X <- V[, set$idx]
  if (!is.null(set$center)) X <- center_by(X, set$center)
  X <- X - rowMeans(X); ss <- rowSums(X^2); X <- X[ss > 1e-10, , drop = FALSE]
  Z <- X / sqrt(rowSums(X^2))
  pid <- intersect(partners, rownames(Z))
  stat <- function(ids) {
    C <- tcrossprod(Z[ids, , drop = FALSE])
    c(coherence = mean(C[upper.tri(C)]), assoc_NKX3_1 = mean(Z[ids, , drop = FALSE] %*% Z[NKX_ID, ]))
  }
  obs <- stat(pid)
  mu  <- rowMeans(V[rownames(Z), ])                         # mean VST over all 12 samples
  dec <- setNames(cut(mu, unique(quantile(mu, 0:10 / 10)), include.lowest = TRUE, labels = FALSE), rownames(Z))
  pool <- setdiff(rownames(Z), c(union_ids, NKX_ID))
  pool_by <- split(pool, dec[pool])
  need <- table(dec[pid])
  set.seed(seed)
  null <- vapply(seq_len(N_PERM), function(i)
    stat(unlist(lapply(names(need), function(d) sample(pool_by[[d]], need[[d]])))), numeric(2))
  list(obs = obs, null = null, n = length(pid), deciles = paste(names(need), need, sep = ":", collapse = ","))
}
t3 <- list(); t3_null <- list(); k <- 0
for (q in names(Q)) for (m in names(modules)) {
 ho_s <- which(cd$CM == ifelse(m == "basal", "CM", "ctrl"))  # context NOT used to define the module
 sets <- c(sample_sets, list(heldout_context_groupcentred = list(idx = ho_s, center = grp[ho_s])))
 for (ss in names(sets)) {
  k <- k + 1
  cr <- coherence(Q[[q]]$V, setdiff(modules[[m]], NKX_ID), sets[[ss]], seed = SEED + k)
  t3_null[[paste(q, m, ss)]] <- cr
  for (st_ in c("coherence", "assoc_NKX3_1")) {
    nl <- cr$null[st_, ]
    t3[[paste(q, m, ss, st_)]] <- data.frame(
      quant = q, module = m, sample_set = ss,
      statistic = ifelse(st_ == "coherence", "mean pairwise r among module genes", "mean r of module genes with NKX3-1"),
      n_module_genes_used = cr$n, observed = round(cr$obs[[st_]], 4),
      null_mean = round(mean(nl), 4), null_q95 = round(quantile(nl, 0.95), 4),
      z_vs_null = round((cr$obs[[st_]] - mean(nl)) / sd(nl), 2),
      p_empirical_upper = (1 + sum(nl >= cr$obs[[st_]])) / (N_PERM + 1),
      decile_composition = cr$deciles)
  }
 }
}
t3 <- do.call(rbind, t3); rownames(t3) <- NULL
t3$selection_note <- ifelse(t3$sample_set == "heldout_context_groupcentred",
  "independent of module selection (other context)",
  "includes samples used to select the module (clade selection favours co-varying genes)")
wtsv(t3, "test3_module_coherence.tsv")

# =============================================================================
# 5. Test 4 -- WGCNA
# =============================================================================
traits <- data.frame(NKX = ifelse(is_nkx, 1, -1), CM = ifelse(is_cm, 1, -1))
traits$NKXxCM <- traits$NKX * traits$CM                     # orthogonal +-1 coding (balanced 2x2)
run_wgcna <- function(V, q) {
  top <- names(sort(rowVars(V, useNames = TRUE), decreasing = TRUE))[1:5000]
  datExpr <- t(V[top, ])
  powers <- c(1:10, seq(12, 30, by = 2))
  # bicor falls back to Pearson for genes with zero MAD (WGCNA default pearsonFallback = "individual");
  # the resulting warnings are expected for genes at the VST floor in >= half of the samples.
  invisible(capture.output(sft <- suppressWarnings(pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",
                           corFnc = bicor, corOptions = list(maxPOutliers = 0.05), blockSize = 6000, verbose = 0))))
  fi <- sft$fitIndices; fi$signed_R2 <- -sign(fi$slope) * fi$SFT.R.sq
  beta <- fi$Power[which(fi$signed_R2 >= 0.80)[1]]
  rule <- "lowest power with signed scale-free R^2 >= 0.80"
  if (is.na(beta)) { beta <- 18; rule <- "no power reached signed R^2 >= 0.80 -> WGCNA FAQ value for signed networks with <20 samples (18)" }
  bwm <- function(pw) {
    assign("cor", WGCNA::cor, envir = globalenv())            # avoid stats::cor masking inside WGCNA
    on.exit(rm("cor", envir = globalenv()))
    suppressWarnings(blockwiseModules(datExpr, power = pw, networkType = "signed", TOMType = "signed",
                     corType = "bicor", maxPOutliers = 0.05, minModuleSize = 30,
                     mergeCutHeight = 0.25, deepSplit = 2, maxBlockSize = 6000,
                     numericLabels = FALSE, saveTOMs = FALSE, randomSeed = SEED,
                     nThreads = 16, verbose = 0))
  }
  net <- bwm(beta)
  # sensitivity: the WGCNA-FAQ power for signed networks with < 20 samples (18)
  sens <- NULL
  if (beta != 18) {
    n18 <- bwm(18); c18 <- setNames(n18$colors, top); k18 <- c18[NKX_ID]
    me18 <- n18$MEs[[paste0("ME", k18)]]
    sens <- data.frame(quant = q, power = 18, n_modules_excl_grey = length(setdiff(unique(c18), "grey")),
                       NKX3_1_module = k18, NKX3_1_module_size = sum(c18 == k18),
                       NKX3_1_module_r_NKX = if (k18 != "grey") round(stats::cor(me18, traits$NKX), 3) else NA,
                       claimed_genes_in_NKX3_1_module = paste(symbol_of(setdiff(intersect(union_ids, top)[c18[intersect(union_ids, top)] == k18], NKX_ID)), collapse = ","),
                       overlap_with_main_NKX3_1_module = sum(c18 == k18 & net$colors == net$colors[match(NKX_ID, top)]))
  }
  colors <- setNames(net$colors, top)
  MEs <- orderMEs(net$MEs); MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
  r  <- stats::cor(MEs, traits); p <- corPvalueStudent(r, nSamples = 12)
  set.seed(SEED)
  rperm <- replicate(N_PERM, stats::cor(MEs, traits[sample.int(12), ]))
  pperm <- (1 + apply(abs(rperm) >= array(abs(r), dim(rperm)) - 1e-12, c(1, 2), sum)) / (N_PERM + 1)
  r2grp <- vapply(colnames(MEs), function(me) summary(lm(MEs[[me]] ~ grp))$r.squared, numeric(1))
  nofeat <- fc_sum$no_feature_pct[match(S, fc_sum$sample_id)]
  r_nofeat <- vapply(colnames(MEs), function(me) stats::cor(MEs[[me]], nofeat), numeric(1))
  me_gm <- t(vapply(colnames(MEs), function(me) tapply(MEs[[me]], grp, mean), numeric(nlevels(grp))))
  sizes <- table(colors)
  mt <- do.call(rbind, lapply(colnames(MEs), function(me) {
    col <- sub("^ME", "", me)
    data.frame(quant = q, module = col, size = as.integer(sizes[col]), trait = colnames(traits),
               r = round(r[me, ], 3), p_student = p[me, ], p_perm10000 = pperm[me, ],
               group_R2 = round(r2grp[[me]], 3),
               ME_mean_Control_basal = round(me_gm[me, "vector_ctrl"], 3), ME_mean_Control_inflammation = round(me_gm[me, "vector_CM"], 3),
               ME_mean_NKX3.1_basal = round(me_gm[me, "NKX3.1_ctrl"], 3), ME_mean_NKX3.1_inflammation = round(me_gm[me, "NKX3.1_CM"], 3),
               r_with_pct_reads_outside_exons = round(r_nofeat[[me]], 3))
  }))
  mt$padj_student <- bh(mt$p_student); mt$padj_perm <- bh(mt$p_perm10000)
  # kME of NKX3-1 and of every claimed gene against every module eigengene
  ids <- intersect(c(NKX_ID, union_ids), rownames(V))
  kme <- suppressWarnings(bicor(t(V[ids, , drop = FALSE]), as.matrix(MEs), maxPOutliers = 0.05))
  if (is.null(dim(kme))) kme <- matrix(kme, nrow = length(ids), dimnames = list(ids, colnames(MEs)))
  rownames(kme) <- ids
  nkx_col <- colors[NKX_ID]
  genes <- data.frame(quant = q, symbol = symbol_of(ids), gene_id = ids,
                      gene_class = ifelse(ids == NKX_ID, "anchor", gene_class[ids]),
                      in_top5000 = ids %in% top,
                      module = ifelse(ids %in% top, colors[ids], "not in network"),
                      same_module_as_NKX3_1 = ids %in% top & colors[ids] == nkx_col,
                      kME_NKX3_1_module = if (nkx_col != "grey") round(kme[, paste0("ME", nkx_col)], 3) else NA,
                      best_module_by_signed_kME = sub("^ME", "", colnames(kme)[max.col(kme, ties.method = "first")]),
                      best_signed_kME = round(apply(kme, 1, max), 3))
  assign_tab <- data.frame(quant = q, gene_id = top, symbol = symbol_of(top), module = colors[top])
  if (nkx_col != "grey") {
    kme_all <- suppressWarnings(bicor(datExpr, MEs[[paste0("ME", nkx_col)]], maxPOutliers = 0.05))[, 1]
    assign_tab$kME_NKX3_1_module <- round(kme_all[top], 3)
    mem <- top[colors[top] == nkx_col]
    hubs <- head(mem[order(-kme_all[mem])], 25)
  } else hubs <- character(0)
  cin <- setdiff(intersect(union_ids, top), NKX_ID)          # claimed genes present in the network
  Kn  <- sum(colors[setdiff(top, NKX_ID)] == nkx_col); kin <- sum(colors[cin] == nkx_col)
  hyp <- data.frame(claimed_in_network = length(cin), claimed_in_NKX3_1_module = kin,
                    expected = round(length(cin) * Kn / 4999, 2),
                    p_hypergeometric = phyper(kin - 1, Kn, 4999 - Kn, length(cin), lower.tail = FALSE))
  fi_out <- data.frame(quant = q, fi[, c("Power", "SFT.R.sq", "slope", "signed_R2", "truncated.R.sq", "mean.k.", "median.k.", "max.k.")],
                       chosen = fi$Power == beta, rule = rule)
  list(sft = fi_out, beta = beta, rule = rule, module_trait = mt, genes = genes, assign = assign_tab,
       nkx_module = nkx_col, nkx_module_size = if (nkx_col != "grey") as.integer(sizes[nkx_col]) else NA,
       nkx_hubs = symbol_of(hubs), MEs = MEs, n_modules = ncol(MEs),
       n_grey = as.integer(sizes["grey"]), sizes = sizes, sensitivity18 = sens, hyper = hyp,
       claimed_in_nkx_module = sum(genes$same_module_as_NKX3_1 & genes$gene_id != NKX_ID))
}
wg <- list()
for (q in names(Q)) {
  wg[[q]] <- run_wgcna(Q[[q]]$V, q)
  message(sprintf("[WGCNA %s] power %d (%s); %d modules (+grey %d genes); NKX3-1 in '%s' (%s genes); claimed genes in NKX3-1 module: %d",
                  q, wg[[q]]$beta, wg[[q]]$rule, wg[[q]]$n_modules, wg[[q]]$n_grey, wg[[q]]$nkx_module,
                  wg[[q]]$nkx_module_size, wg[[q]]$claimed_in_nkx_module))
}
wtsv(do.call(rbind, lapply(wg, `[[`, "sft")), "test4_wgcna_soft_threshold.tsv")
wtsv(do.call(rbind, lapply(wg, `[[`, "module_trait")), "test4_wgcna_module_trait.tsv")
wtsv(do.call(rbind, lapply(wg, `[[`, "genes")), "test4_wgcna_claimed_genes.tsv")
for (q in names(wg)) wtsv(wg[[q]]$assign, sprintf("test4_wgcna_gene_assignment_%s.tsv", q))
sens_tab <- do.call(rbind, lapply(wg, `[[`, "sensitivity18"))
if (!is.null(sens_tab)) wtsv(sens_tab, "test4_wgcna_sensitivity_power18.tsv")
wg_summary <- do.call(rbind, lapply(names(wg), function(q) {
  w <- wg[[q]]; mt <- w$module_trait
  data.frame(quant = q, genes_in_network = 5000, soft_power = w$beta, power_rule = w$rule,
             n_modules_excl_grey = w$n_modules, n_grey = w$n_grey,
             median_module_group_R2 = round(median(unique(mt[, c("module", "group_R2")])$group_R2), 3),
             frac_modules_group_R2_ge_0.8 = round(mean(unique(mt[, c("module", "group_R2")])$group_R2 >= 0.8), 3),
             group_R2_expected_without_group_effect = round(3 / 11, 3),
             n_modules_NKX_padj_perm_lt_0.05 = sum(mt$trait == "NKX" & mt$padj_perm < 0.05),
             n_modules_CM_padj_perm_lt_0.05 = sum(mt$trait == "CM" & mt$padj_perm < 0.05),
             n_modules_NKXxCM_padj_perm_lt_0.05 = sum(mt$trait == "NKXxCM" & mt$padj_perm < 0.05),
             n_modules_abs_r_pct_reads_outside_exons_ge_0.7 = sum(abs(unique(mt[, c("module", "r_with_pct_reads_outside_exons")])$r_with_pct_reads_outside_exons) >= 0.7),
             NKX3_1_module = w$nkx_module, NKX3_1_module_size = w$nkx_module_size,
             NKX3_1_module_r_pct_reads_outside_exons = if (w$nkx_module != "grey") mt$r_with_pct_reads_outside_exons[mt$module == w$nkx_module & mt$trait == "NKX"] else NA,
             NKX3_1_module_r_NKX = if (w$nkx_module != "grey") mt$r[mt$module == w$nkx_module & mt$trait == "NKX"] else NA,
             NKX3_1_module_padj_perm_NKX = if (w$nkx_module != "grey") mt$padj_perm[mt$module == w$nkx_module & mt$trait == "NKX"] else NA,
             NKX3_1_module_r_NKXxCM = if (w$nkx_module != "grey") mt$r[mt$module == w$nkx_module & mt$trait == "NKXxCM"] else NA,
             claimed_genes_in_NKX3_1_module = w$claimed_in_nkx_module,
             claimed_genes_in_network = w$hyper$claimed_in_network,
             claimed_expected_in_NKX3_1_module = w$hyper$expected,
             claimed_enrichment_p_hypergeometric = signif(w$hyper$p_hypergeometric, 3),
             NKX3_1_module_top_hubs = paste(w$nkx_hubs, collapse = ","))
}))
wtsv(wg_summary, "test4_wgcna_summary.tsv")

# =============================================================================
# 6. Read support, per-context DESeq2 evidence and NKX3.1 x CM interaction ("rewiring")
# =============================================================================
grp_means <- function(M) sapply(levels(grp), function(l) rowMeans(M[, grp == l, drop = FALSE]))
norm_fc_all <- t(t(core$counts_fc[, S]) / sizeFactors(core$fc$dds_full)[S])
norm_sal    <- core$sal$norm[, S]
sm <- read.delim(salmon_merged, check.names = FALSE)
smm <- as.matrix(sm[, -(1:2)]); rownames(smm) <- sm$gene_id
colnames(smm) <- gsub(".", "-", sub("^X", "", colnames(smm)), fixed = TRUE); smm <- smm[, S]
norm_smm <- t(t(smm) / estimateSizeFactorsForMatrix(round(smm[rowSums(smm >= 10) >= 3, ])))
slx <- c(SLX1A = "ENSG00000132207", SLX1B = "ENSG00000181625")
rs_ids <- c(NKX_ID, setdiff(union_ids, NKX_ID))
deseq_cols <- function(qq, ids) {
  res <- core[[qq]]$res; o <- list()
  for (cn in c("NKX_in_ctrl", "NKX_in_CM", "interaction", "CM_in_NKX")) {
    r <- res[[cn]][match(ids, res[[cn]]$gene_id), ]
    o[[paste0(cn, "_LFC_MLE")]]    <- round(r$log2FC_MLE, 3)
    o[[paste0(cn, "_LFC_shrunk")]] <- round(r$log2FC_shrunk, 3)
    o[[paste0(cn, "_padj")]]       <- signif(r$padj, 3)
  }
  o <- as.data.frame(o); names(o) <- paste0(ifelse(qq == "fc", "corrected_", "original_"), names(o)); o
}
mk_rs <- function(ids, labels) {
  gm_fc <- grp_means(norm_fc_all[ids, , drop = FALSE]); gm_sal <- grp_means(norm_sal[ids, , drop = FALSE])
  gm_smm <- grp_means(norm_smm[ids, , drop = FALSE])
  lr <- function(g, a, b) round(log2((g[, a] + 1) / (g[, b] + 1)), 2)
  data.frame(symbol = labels, gene_id = ids, gene_class = ifelse(ids == NKX_ID, "anchor", gene_class[ids]),
             corrected_quantifiable = ids %in% rownames(core$fc$vst),
             corrected_max_unique_reads_per_sample = apply(core$counts_fc[ids, S, drop = FALSE], 1, max),
             setNames(as.data.frame(round(gm_fc, 1)), paste0("corrected_mean_", colnames(gm_fc))),
             setNames(as.data.frame(round(gm_sal, 1)), paste0("original_mean_", colnames(gm_sal))),
             setNames(as.data.frame(round(gm_smm, 1)), paste0("salmonNumReads_mean_", colnames(gm_smm))),
             corrected_log2ratio_NKX_basal = lr(gm_fc, "NKX3.1_ctrl", "vector_ctrl"),
             corrected_log2ratio_NKX_inflam = lr(gm_fc, "NKX3.1_CM", "vector_CM"),
             original_log2ratio_NKX_basal = lr(gm_sal, "NKX3.1_ctrl", "vector_ctrl"),
             original_log2ratio_NKX_inflam = lr(gm_sal, "NKX3.1_CM", "vector_CM"),
             salmonNumReads_log2ratio_NKX_basal = lr(gm_smm, "NKX3.1_ctrl", "vector_ctrl"),
             salmonNumReads_log2ratio_NKX_inflam = lr(gm_smm, "NKX3.1_CM", "vector_CM"),
             deseq_cols("fc", ids), deseq_cols("sal", ids), check.names = FALSE)
}
rs <- mk_rs(rs_ids, symbol_of(rs_ids))
# Artifact rule (fixed before inspecting outcomes): the original NKX3.1 effect in a context is significant
# (original DESeq2 padj < 0.05) but is not reproduced by genome-aligned unique reads: the gene fails the
# corrected expression filter, or the corrected MLE log2FC is < 1/4 of the original one.
art_ctx <- function(ctx) {
  cn <- ifelse(ctx == "basal", "NKX_in_ctrl", "NKX_in_CM")
  o_sig <- !is.na(rs[[paste0("original_", cn, "_padj")]]) & rs[[paste0("original_", cn, "_padj")]] < 0.05
  o_l <- rs[[paste0("original_", cn, "_LFC_MLE")]]; c_l <- rs[[paste0("corrected_", cn, "_LFC_MLE")]]
  o_sig & (!rs$corrected_quantifiable | is.na(c_l) | (c_l < 0.25 * o_l))
}
rs$artifact_basal <- art_ctx("basal"); rs$artifact_inflammation <- art_ctx("inflammation")
rs$original_signif_basal <- !is.na(rs$original_NKX_in_ctrl_padj) & rs$original_NKX_in_ctrl_padj < 0.05
rs$original_signif_inflammation <- !is.na(rs$original_NKX_in_CM_padj) & rs$original_NKX_in_CM_padj < 0.05
rs$corrected_signif_basal <- !is.na(rs$corrected_NKX_in_ctrl_padj) & rs$corrected_NKX_in_ctrl_padj < 0.05
rs$corrected_signif_inflammation <- !is.na(rs$corrected_NKX_in_CM_padj) & rs$corrected_NKX_in_CM_padj < 0.05
wtsv(rs, "read_support_and_deseq2_module_genes.tsv")
slx_block <- function(M, label) {                          # paralogue pair: per-gene and summed group means
  gm <- t(vapply(slx, function(id) vapply(levels(grp), function(l)
    if (id %in% rownames(M)) mean(M[id, grp == l]) else NA_real_, numeric(1)), numeric(nlevels(grp))))
  data.frame(feature = c(names(slx), "SLX1A+SLX1B"), quant = label, round(rbind(gm, colSums(gm)), 1), row.names = NULL)
}
slx_all <- rbind(slx_block(norm_sal, "original (manuscript counts; tximport lengthScaledTPM)"),
                 slx_block(norm_smm, "Salmon NumReads (nf-core salmon.merged.gene_counts)"),
                 slx_block(norm_fc_all, "corrected (featureCounts unique reads)"))
wtsv(slx_all, "read_support_SLX1A_SLX1B_paralogues.tsv")

# Rewiring: does the NKX3.1 effect differ between contexts?
rw_gene <- data.frame(symbol = rs$symbol, gene_id = rs$gene_id, gene_class = rs$gene_class,
                      corrected_NKX_basal_LFC = rs$corrected_NKX_in_ctrl_LFC_MLE, corrected_NKX_inflam_LFC = rs$corrected_NKX_in_CM_LFC_MLE,
                      corrected_interaction_LFC = rs$corrected_interaction_LFC_MLE, corrected_interaction_padj = rs$corrected_interaction_padj,
                      corrected_interaction_p = core$fc$res$interaction$pvalue[match(rs$gene_id, core$fc$res$interaction$gene_id)],
                      original_NKX_basal_LFC = rs$original_NKX_in_ctrl_LFC_MLE, original_NKX_inflam_LFC = rs$original_NKX_in_CM_LFC_MLE,
                      original_interaction_LFC = rs$original_interaction_LFC_MLE, original_interaction_padj = rs$original_interaction_padj,
                      original_interaction_p = core$sal$res$interaction$pvalue[match(rs$gene_id, core$sal$res$interaction$gene_id)])
wtsv(rw_gene, "test5_rewiring_gene_interaction.tsv")
module_score_test <- function(V, ids) {
  ids <- intersect(ids, rownames(V))
  sc <- colMeans(t(scale(t(V[ids, , drop = FALSE]))))
  f <- summary(lm(sc ~ cd$NKX * cd$CM))$coefficients
  data.frame(n_genes = length(ids), NKX_effect_basal = round(f[2, 1], 3), NKX_effect_basal_p = f[2, 4],
             interaction = round(f[4, 1], 3), interaction_t = round(f[4, 3], 2), interaction_p = f[4, 4])
}
gsets <- list(Fig3D_only = names(gene_class)[gene_class == "Fig3D_only"],
              Fig3E_only = names(gene_class)[gene_class == "Fig3E_only"],
              both_excl_NKX3_1 = setdiff(names(gene_class)[gene_class == "both"], NKX_ID))
rw_mod <- do.call(rbind, lapply(names(Q), function(q) do.call(rbind, lapply(names(gsets), function(s)
  data.frame(quant = q, gene_set = s, prediction_if_rewired = c(Fig3D_only = "interaction < 0", Fig3E_only = "interaction > 0", both_excl_NKX3_1 = "none")[s],
             module_score_test(Q[[q]]$V, gsets[[s]]))))))
wtsv(rw_mod, "test5_rewiring_module_score_interaction.tsv")
ctx_cor <- sapply(c("corrected", "original"), function(q) {
  a <- rw_gene[[paste0(q, "_NKX_basal_LFC")]]; b <- rw_gene[[paste0(q, "_NKX_inflam_LFC")]]
  ok <- !is.na(a) & !is.na(b) & rw_gene$gene_id != NKX_ID
  c(r = cor(a[ok], b[ok]), n = sum(ok))
})

# --- selection-aware null for "rewiring" -------------------------------------
# The published modules were picked from the same data whose context-specificity is then claimed.
# Re-run the published selection with randomly re-assigned contexts (3 of the 6 NKX3.1-transfected and
# 3 of the 6 vector samples called "basal"; 20 x 20 = 400 relabelings, the true one included).
# Rule that reproduces both published boxes exactly: on the top-200 Welch genes of a context, the
# smallest dendrogram clade that contains NKX3-1 and has >= 14 genes.
derive_modules <- function(V, pseudo_basal) {
  res <- list()
  for (ctx in c("basal", "inflammation")) {
    s <- if (ctx == "basal") which(pseudo_basal) else which(!pseudo_basal)
    i1 <- s[is_nkx[s]]; i2 <- s[!is_nkx[s]]
    mu1 <- rowMeans(V[, i1]); mu2 <- rowMeans(V[, i2]); lfc <- mu1 - mu2
    z <- lfc / sqrt(pmax(rowVars(V[, i1]), 1e-6) / 3 + pmax(rowVars(V[, i2]), 1e-6) / 3)
    sig <- pmax(mu1, mu2) >= 2 & abs(lfc) >= 0.25 & abs(z) >= 2
    o <- order(-abs(z), -abs(lfc)); so <- o[sig[o]]
    keep <- sort(so[seq_len(min(200, length(so)))])
    mat <- t(scale(t(V[keep, c(i1, i2)])))
    hc <- hclust(dist(mat), method = "complete")
    if (!NKX_ID %in% hc$labels) { res[[ctx]] <- character(0); next }
    walk <- clade_walk(hc, NKX_ID)
    res[[ctx]] <- walk[[which(vapply(walk, function(w) length(w$leaves) >= 14, logical(1)))[1]]]$leaves
  }
  res
}
rewire_stats <- function(V, mods, pseudo_basal) {
  b_only <- setdiff(mods$basal, c(mods$inflammation, NKX_ID))
  i_only <- setdiff(mods$inflammation, c(mods$basal, NKX_ID))
  pcm <- !pseudo_basal
  tint <- function(ids) {
    if (length(ids) < 2) return(NA_real_)
    sc <- colMeans(t(scale(t(V[ids, , drop = FALSE]))))
    summary(lm(sc ~ is_nkx * pcm))$coefficients[4, 3]
  }
  tb <- tint(b_only); ti <- tint(i_only)
  c(t_basal_only = tb, t_inflammation_only = ti, T = ti - tb,
    n_basal_only = length(b_only), n_inflammation_only = length(i_only),
    jaccard = length(intersect(mods$basal, mods$inflammation)) / length(union(mods$basal, mods$inflammation)))
}
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
cmb <- combn(6, 3); relab <- expand.grid(a = seq_len(ncol(cmb)), b = seq_len(ncol(cmb)))
nk_idx <- which(is_nkx); vec_idx <- which(!is_nkx); true_pb <- cd$CM == "ctrl"
pb_of <- function(k) { pb <- rep(FALSE, 12); pb[nk_idx[cmb[, relab$a[k]]]] <- TRUE; pb[vec_idx[cmb[, relab$b[k]]]] <- TRUE; pb }
rew_null <- list(); rew_mods <- list()
for (q in names(Q)) {
  V <- Q[[q]]$V
  rew_mods[[q]] <- derive_modules(V, true_pb)
  stt <- t(vapply(seq_len(nrow(relab)), function(k) { pb <- pb_of(k); rewire_stats(V, derive_modules(V, pb), pb) }, numeric(6)))
  rew_null[[q]] <- data.frame(quant = q, relab_nkx = relab$a, relab_vector = relab$b, stt,
                              observed_labeling = vapply(seq_len(nrow(relab)), function(k) identical(pb_of(k), true_pb), logical(1)))
}
rule_reproduces <- setequal(rew_mods$original$basal, modules$basal) && setequal(rew_mods$original$inflammation, modules$inflammation)
message(sprintf("[rewiring] rule-based selection reproduces the published modules on the original counts: %s", rule_reproduces))
rew_summary <- do.call(rbind, lapply(names(rew_null), function(q) {
  x <- rew_null[[q]]; o <- x[x$observed_labeling, ]
  data.frame(quant = q, n_relabelings = nrow(x),
             modules_basal = paste(symbol_of(rew_mods[[q]]$basal), collapse = ","),
             modules_inflammation = paste(symbol_of(rew_mods[[q]]$inflammation), collapse = ","),
             jaccard_with_published_basal = round(jac(rew_mods[[q]]$basal, modules$basal), 3),
             jaccard_with_published_inflammation = round(jac(rew_mods[[q]]$inflammation, modules$inflammation), 3),
             rule_reproduces_published = if (q == "original") rule_reproduces else NA,
             observed_t_basal_only = round(o$t_basal_only, 2), observed_t_inflammation_only = round(o$t_inflammation_only, 2),
             observed_T = round(o$T, 2), null_median_T = round(median(x$T, na.rm = TRUE), 2),
             null_q95_T = round(quantile(x$T, 0.95, na.rm = TRUE), 2),
             p_selection_aware_T = mean(x$T >= o$T - 1e-9, na.rm = TRUE),
             observed_jaccard_basal_vs_inflammation = round(o$jaccard, 3),
             null_median_jaccard = round(median(x$jaccard), 3),
             p_jaccard_as_low_or_lower = mean(x$jaccard <= o$jaccard + 1e-9))
}))
# composition of the context-only gene sets of the rule-derived modules (true contexts)
comp <- function(q, ids) {
  ids <- setdiff(ids, NKX_ID); res <- if (q == "corrected") core$fc$res$NKX_in_CM else core$sal$res$NKX_in_CM
  c(n = length(ids), noncoding = sum(ann[ids, "biotype"] != "protein_coding"),
    median_baseMean = median(res$baseMean[match(ids, res$gene_id)]))
}
for (q in names(rew_mods)) {
  m <- rew_mods[[q]]
  cb <- comp(q, setdiff(m$basal, m$inflammation)); ci <- comp(q, setdiff(m$inflammation, m$basal))
  rew_summary[rew_summary$quant == q, "basal_only_n_noncoding_median_baseMean"] <- sprintf("%d genes, %d non-coding, median baseMean %.0f", cb[["n"]], cb[["noncoding"]], cb[["median_baseMean"]])
  rew_summary[rew_summary$quant == q, "inflammation_only_n_noncoding_median_baseMean"] <- sprintf("%d genes, %d non-coding, median baseMean %.0f", ci[["n"]], ci[["noncoding"]], ci[["median_baseMean"]])
  V <- Q[[q]]$V; io <- setdiff(m$inflammation, c(m$basal, NKX_ID))
  sc <- colMeans(t(scale(t(V[io, , drop = FALSE]))))
  rew_summary[rew_summary$quant == q, "inflammation_only_score_group_means"] <- paste(sprintf("%s %.2f", levels(grp), tapply(sc, grp, mean)), collapse = "; ")
}
wtsv(rew_summary, "test5_rewiring_selection_aware_permutation.tsv")

wtsv(do.call(rbind, rew_null), "test5_rewiring_selection_aware_null_distribution.tsv")

# =============================================================================
# 7. Claim-by-claim verdicts
# =============================================================================
rsv  <- function(sym_, col) rs[rs$symbol == sym_, col]
t2v  <- function(q, m, sym_, col) t2[t2$quant == q & t2$module == m & t2$symbol == sym_, col]
t3v  <- function(q, m, ss, st_, col) t3[t3$quant == q & t3$module == m & t3$sample_set == ss &
                                         t3$statistic == c(coh = "mean pairwise r among module genes", assoc = "mean r of module genes with NKX3-1")[[st_]], col]
gm_str <- function(sym_, pre) paste(round(unlist(rs[rs$symbol == sym_, paste0(pre, "_mean_", levels(grp))])), collapse = "/")
nsig_mod <- function(q, m) {
  x <- t2[t2$quant == q & t2$module == m & t2$quantifiable, ]
  c(n = nrow(x), dose = sum(x$dose_padj_module < 0.1, na.rm = TRUE), partial = sum(x$partial_padj_module < 0.1, na.rm = TRUE),
    heldout = sum(x$heldout_padj_module < 0.1, na.rm = TRUE))
}
mod_text <- function(q, m) {
  n <- nsig_mod(q, m)
  sprintf("%s module, %d %spartner genes: module-BH<0.1 for dose-response %d, partial r (12 samples) %d, held-out partial r %d; held-out coherence r=%.2f (matched null %.2f, p=%s), held-out mean r with NKX3-1 %.2f (null %.2f, p=%s); all-12 group-centred coherence r=%.2f (p=%s), mean r with NKX3-1 %.2f (p=%s)",
          m, n[["n"]], ifelse(q == "corrected", "quantifiable ", ""), n[["dose"]], n[["partial"]], n[["heldout"]],
          t3v(q, m, "heldout_context_groupcentred", "coh", "observed"), t3v(q, m, "heldout_context_groupcentred", "coh", "null_mean"),
          fmt_p(t3v(q, m, "heldout_context_groupcentred", "coh", "p_empirical_upper")),
          t3v(q, m, "heldout_context_groupcentred", "assoc", "observed"), t3v(q, m, "heldout_context_groupcentred", "assoc", "null_mean"),
          fmt_p(t3v(q, m, "heldout_context_groupcentred", "assoc", "p_empirical_upper")),
          t3v(q, m, "all12_groupcentred", "coh", "observed"), fmt_p(t3v(q, m, "all12_groupcentred", "coh", "p_empirical_upper")),
          t3v(q, m, "all12_groupcentred", "assoc", "observed"), fmt_p(t3v(q, m, "all12_groupcentred", "assoc", "p_empirical_upper")))
}
wg_text <- function(q) {
  w <- wg[[q]]; mt <- w$module_trait[w$module_trait$module == w$nkx_module & w$module_trait$trait == "NKX", ]
  sprintf("WGCNA: NKX3-1 in module '%s' (%d genes; NKX3-1 kME %.2f); eigengene r with NKX3.1 status %.2f (permutation padj %s), with CM %.2f, with %% reads outside exons %.2f; %d/%d claimed genes in the network share NKX3-1's module (expected %.1f, hypergeometric p=%s); no module eigengene associated with NKX3.1, CM or NKX3.1xCM at padj<0.05",
          w$nkx_module, w$nkx_module_size, w$genes$kME_NKX3_1_module[w$genes$gene_id == NKX_ID], mt$r, fmt_p(mt$padj_perm),
          w$module_trait$r[w$module_trait$module == w$nkx_module & w$module_trait$trait == "CM"], mt$r_with_pct_reads_outside_exons,
          w$hyper$claimed_in_NKX3_1_module, w$hyper$claimed_in_network, w$hyper$expected, fmt_p(w$hyper$p_hypergeometric))
}
gene_text <- function(sym_, m, q) {
  cn <- ifelse(m == "basal", "NKX_in_ctrl", "NKX_in_CM"); pre <- ifelse(q == "corrected", "corrected", "original")
  if (q == "corrected" && !rsv(sym_, "corrected_quantifiable"))
    return(sprintf("%s: at most %d unique read%s per sample (group means %s); below expression filter, not testable",
                   sym_, rsv(sym_, "corrected_max_unique_reads_per_sample"),
                   ifelse(rsv(sym_, "corrected_max_unique_reads_per_sample") == 1, "", "s"), gm_str(sym_, "corrected")))
  sprintf("%s: group means %s; NKX3.1 effect (%s) log2FC %.2f, padj %s; dose slope %.2f (p %s); partial r %.2f (perm p %s); held-out partial r %.2f (p %s)",
          sym_, gm_str(sym_, pre), m, rsv(sym_, paste0(pre, "_", cn, "_LFC_MLE")), fmt_p(rsv(sym_, paste0(pre, "_", cn, "_padj"))),
          t2v(q, m, sym_, "dose_slope"), fmt_p(t2v(q, m, sym_, "dose_p")), t2v(q, m, sym_, "partial_r"),
          fmt_p(t2v(q, m, sym_, "partial_p_mc10000")), t2v(q, m, sym_, "heldout_partial_r"), fmt_p(t2v(q, m, sym_, "heldout_p_exact36")))
}
numreads_text <- function(sym_) sprintf("%s Salmon NumReads (nf-core merged) group means %s", sym_, gm_str(sym_, "salmonNumReads"))
slx_sum <- slx_all[slx_all$feature == "SLX1A+SLX1B" & grepl("NumReads", slx_all$quant), levels(grp)]
p20min <- min(t1_summary$min_p_perm20)
fig3_r <- function(m) { x <- t1_mod[t1_mod$quant == "original" & t1_mod$context == m & t1_mod$in_module_of_this_context & t1_mod$gene_id != NKX_ID, ]
  sprintf("Fig 3 view (original, %s): r with NKX3-1 %.2f-%.2f across 6 samples, but exact permutation p >= %.1f for every gene (20 balanced relabelings)", m, min(x$r_pearson), max(x$r_pearson), p20min) }
rw <- function(q, col) rew_summary[rew_summary$quant == q, col]

named <- list(i = c("SLX1A", "BCL11B", "ZNF331"), ii = c("TPGS1", "SRMS"), iii = "BCL11B", iv = "LRIG1")
claim_ctx <- list(i = "basal", ii = "inflammation", iii = c("basal", "inflammation"), iv = "inflammation")
is_art <- function(sym_, ctxs) all(vapply(ctxs, function(cx) rsv(sym_, paste0("artifact_", cx)), logical(1)))
corr_support <- function(sym_, m) {   # supported in corrected data: NKX3.1 effect AND co-expression beyond group membership
  if (!rsv(sym_, "corrected_quantifiable")) return(FALSE)
  cn <- ifelse(m == "basal", "NKX_in_ctrl", "NKX_in_CM")
  isTRUE(rsv(sym_, paste0("corrected_", cn, "_padj")) < 0.05) &&
    any(c(t2v("corrected", m, sym_, "partial_padj_module"), t2v("corrected", m, sym_, "dose_padj_module"),
          t2v("corrected", m, sym_, "heldout_padj_module")) < 0.1, na.rm = TRUE)
}
# Verdict rule (fixed before inspecting outcomes): 'artifact of quantification' when >= half of the named genes
# meet the artifact rule; 'supported' when every named gene is supported in the corrected data (claim ii
# additionally needs a context-specific NKX3.1 effect that survives the selection-aware permutation);
# otherwise 'not supported'.
verdict_rule <- function(cl) {
  syms <- named[[cl]]; ctxs <- claim_ctx[[cl]]
  arts <- vapply(syms, is_art, logical(1), ctxs = ctxs)
  sup  <- vapply(syms, function(s) all(vapply(ctxs, function(cx) corr_support(s, cx), logical(1))), logical(1))
  if (cl == "ii" && !(isTRUE(rw("corrected", "p_selection_aware_T") < 0.05) && all(sup)))
    return(list(v = "not supported", arts = arts, sup = sup))
  list(v = if (mean(arts) >= 0.5) "artifact of quantification" else if (all(sup)) "supported" else "not supported",
       arts = arts, sup = sup)
}
claim_text <- c(
  i   = "(i) Basal: NKX3.1 co-expressed with a module incl. SLX1A, BCL11B, ZNF331 (Fig 3D)",
  ii  = "(ii) Inflammation: module shifts to a distinct one incl. TPGS1 and SRMS (Fig 3E)",
  iii = "(iii) BCL11B remains associated with NKX3.1 across conditions",
  iv  = "(iv) LRIG1 induced under inflammation")
dec_c <- list(
  i = c(gene_text("SLX1A", "basal", "corrected"),
        sprintf("SLX1A+SLX1B (near-identical paralogues) summed Salmon NumReads %s: no NKX3.1 change", paste(round(unlist(slx_sum)), collapse = "/")),
        gene_text("BCL11B", "basal", "corrected"), gene_text("ZNF331", "basal", "corrected"),
        mod_text("corrected", "basal"), wg_text("corrected")),
  ii = c(gene_text("TPGS1", "inflammation", "corrected"), gene_text("SRMS", "inflammation", "corrected"),
         mod_text("corrected", "inflammation"),
         sprintf("rewiring: %d/%d published module genes (excl. NKX3-1) with DESeq2 NKX3.1xCM interaction padj<0.1; Fig3E-only module-score interaction t=%.2f (p %s)",
                 sum(rw_gene$corrected_interaction_padj[rw_gene$gene_id != NKX_ID] < 0.1, na.rm = TRUE),
                 sum(!is.na(rw_gene$corrected_interaction_padj[rw_gene$gene_id != NKX_ID])),
                 rw_mod$interaction_t[rw_mod$quant == "corrected" & rw_mod$gene_set == "Fig3E_only"],
                 fmt_p(rw_mod$interaction_p[rw_mod$quant == "corrected" & rw_mod$gene_set == "Fig3E_only"])),
         sprintf("modules re-derived from corrected counts with the published procedure overlap the published boxes by Jaccard %.2f (basal) and %.2f (inflammation); their context-specificity T=%.2f exceeds the selection-aware null (p=%s), but the inflammation-only set (%s) is raised only in NKX3.1-inflammation libraries (%s), which also carry the most reads outside annotated exons (%s vs %s in the other nine libraries) - a library-level pattern, not the published genes",
                 rw("corrected", "jaccard_with_published_basal"), rw("corrected", "jaccard_with_published_inflammation"),
                 rw("corrected", "observed_T"), fmt_p(rw("corrected", "p_selection_aware_T")),
                 rw("corrected", "inflammation_only_n_noncoding_median_baseMean"), rw("corrected", "inflammation_only_score_group_means"),
                 nofeat_rng("NKX3.1-inflammation"), nofeat_other)),
  iii = c(gene_text("BCL11B", "basal", "corrected")),
  iv = c(gene_text("LRIG1", "inflammation", "corrected"),
         sprintf("LRIG1 NKX3.1xCM interaction (corrected) log2FC %.2f, padj %s; CM effect within NKX3.1 cells log2FC %.2f, padj %s",
                 rsv("LRIG1", "corrected_interaction_LFC_MLE"), fmt_p(rsv("LRIG1", "corrected_interaction_padj")),
                 rsv("LRIG1", "corrected_CM_in_NKX_LFC_MLE"), fmt_p(rsv("LRIG1", "corrected_CM_in_NKX_padj")))))
dec_o <- list(
  i = c(fig3_r("basal"), gene_text("SLX1A", "basal", "original"), gene_text("BCL11B", "basal", "original"),
        gene_text("ZNF331", "basal", "original"), mod_text("original", "basal"), wg_text("original")),
  ii = c(fig3_r("inflammation"), gene_text("TPGS1", "inflammation", "original"), gene_text("SRMS", "inflammation", "original"),
         sprintf("TPGS1 and SRMS are also inside the Fig 3D (basal) box; %d/%d Fig 3E genes are shared with Fig 3D (Jaccard %.2f; null median under random contexts %.2f, p=%s)",
                 length(shared_ids), length(modules$inflammation), rw("original", "observed_jaccard_basal_vs_inflammation"),
                 rw("original", "null_median_jaccard"), fmt_p(rw("original", "p_jaccard_as_low_or_lower"))),
         sprintf("module-score NKX3.1xCM interaction t=%.2f (Fig3E-only genes) and t=%.2f (Fig3D-only genes), but the identical selection applied to random context labels gives a context-specificity statistic at least as large (T=%.2f) in %.0f%% of 400 relabelings (selection-aware p=%s); %d/%d module genes (excl. NKX3-1) with DESeq2 interaction padj<0.1 (min padj %s)",
                 rw("original", "observed_t_inflammation_only"), rw("original", "observed_t_basal_only"), rw("original", "observed_T"),
                 100 * rw("original", "p_selection_aware_T"), fmt_p(rw("original", "p_selection_aware_T")),
                 sum(rw_gene$original_interaction_padj[rw_gene$gene_id != NKX_ID] < 0.1, na.rm = TRUE),
                 sum(!is.na(rw_gene$original_interaction_padj[rw_gene$gene_id != NKX_ID])),
                 fmt_p(min(rw_gene$original_interaction_padj[rw_gene$gene_id != NKX_ID], na.rm = TRUE))),
         mod_text("original", "inflammation")),
  iii = c(gene_text("BCL11B", "basal", "original"), gene_text("BCL11B", "inflammation", "original"), numreads_text("BCL11B")),
  iv = c(gene_text("LRIG1", "inflammation", "original"), numreads_text("LRIG1")))
ho <- function(q, m, st_, col) t3v(q, m, "heldout_context_groupcentred", st_, col)
summ <- c(
  i = sprintf("SLX1A and BCL11B have no unique-read support (at most %d and %d reads per sample; the SLX1A gain is offset by SLX1B, summed NumReads %s); ZNF331 is expressed but not NKX3.1-associated (log2FC %.2f, padj %s). In samples not used to pick the module, its genes do not co-vary with NKX3-1 beyond group membership (mean r %.2f, matched-null p = %s corrected; %.2f, p = %s original). In-sample 'co-expression' (6-sample r up to 0.99) cannot reach p < %.1f by label permutation.",
              rsv("SLX1A", "corrected_max_unique_reads_per_sample"), rsv("BCL11B", "corrected_max_unique_reads_per_sample"),
              paste(round(unlist(slx_sum)), collapse = "/"), rsv("ZNF331", "corrected_NKX_in_ctrl_LFC_MLE"), fmt_p(rsv("ZNF331", "corrected_NKX_in_ctrl_padj")),
              ho("corrected", "basal", "assoc", "observed"), fmt_p(ho("corrected", "basal", "assoc", "p_empirical_upper")),
              ho("original", "basal", "assoc", "observed"), fmt_p(ho("original", "basal", "assoc", "p_empirical_upper")), p20min),
  ii = sprintf("TPGS1 and SRMS are already inside the basal (Fig 3D) box and half of the Fig 3E genes are shared with it; their NKX3.1 effects are not supported by unique reads (TPGS1 log2FC %.2f, padj %s; SRMS at most %d reads per sample). Running the published selection on random context labels gives an equal or larger context switch %.0f%% of the time (selection-aware p = %s), and no module gene has a significant NKX3.1 x CM interaction in either quantification.",
               rsv("TPGS1", "corrected_NKX_in_CM_LFC_MLE"), fmt_p(rsv("TPGS1", "corrected_NKX_in_CM_padj")), rsv("SRMS", "corrected_max_unique_reads_per_sample"),
               100 * rw("original", "p_selection_aware_T"), fmt_p(rw("original", "p_selection_aware_T"))),
  iii = sprintf("BCL11B has at most %d unique reads per sample in every group (normalised means %s); the original %s pattern (log2FC %.1f basal, %.1f inflammation) exists only in Salmon read assignment (NumReads %s).",
                rsv("BCL11B", "corrected_max_unique_reads_per_sample"), gm_str("BCL11B", "corrected"), gm_str("BCL11B", "original"),
                rsv("BCL11B", "original_NKX_in_ctrl_LFC_MLE"), rsv("BCL11B", "original_NKX_in_CM_LFC_MLE"), gm_str("BCL11B", "salmonNumReads")),
  iv = sprintf("LRIG1 is well covered (normalised unique-read means %s) and unchanged by NKX3.1 under CM (log2FC %.2f, padj %s); the original log2FC %.2f (padj %s) arises from low Control-inflammation values in the manuscript's tximport lengthScaledTPM table (%s) that are absent from the Salmon NumReads (%s).",
               gm_str("LRIG1", "corrected"), rsv("LRIG1", "corrected_NKX_in_CM_LFC_MLE"), fmt_p(rsv("LRIG1", "corrected_NKX_in_CM_padj")),
               rsv("LRIG1", "original_NKX_in_CM_LFC_MLE"), fmt_p(rsv("LRIG1", "original_NKX_in_CM_padj")),
               gm_str("LRIG1", "original"), gm_str("LRIG1", "salmonNumReads")))
vt <- do.call(rbind, lapply(names(named), function(cl) {
  vr <- verdict_rule(cl)
  data.frame(claim = claim_text[[cl]], verdict = vr$v, summary = summ[[cl]],
             named_genes_artifact = paste(named[[cl]][vr$arts], collapse = ","),
             named_genes_supported_in_corrected = paste(named[[cl]][vr$sup], collapse = ","),
             group_order_of_means = "Control-basal/Control-inflammation/NKX3.1-basal/NKX3.1-inflammation (DESeq2-normalised counts)",
             decisive_corrected = paste(dec_c[[cl]], collapse = " | "),
             decisive_original  = paste(dec_o[[cl]], collapse = " | "))
}))
wtsv(vt, "verdict_table.tsv")

# =============================================================================
# 8. Save everything for the figure script
# =============================================================================
saveRDS(list(modules = modules, module_def = module_def, union_ids = union_ids, gene_class = gene_class,
             shared_ids = shared_ids, hm_store = hm_store, t1_mod = t1_mod, t1_summary = t1_summary,
             t2 = t2, t2_gw_summary = t2_gw_summary, t3 = t3, t3_null = t3_null, wg = wg, wg_summary = wg_summary,
             rs = rs, rw_gene = rw_gene, rw_mod = rw_mod, ctx_cor = ctx_cor, vt = vt, slx = slx_all,
             rew_null = rew_null, rew_summary = rew_summary, rew_mods = rew_mods, fc_sum = fc_sum,
             vst = list(corrected = Q$corrected$V[intersect(c(NKX_ID, union_ids), rownames(Q$corrected$V)), ],
                        corrected_allgenes = vst_fc_all[c(NKX_ID, union_ids), ],
                        original = Q$original$V[c(NKX_ID, union_ids), ]),
             resid = resid_store, cd = cd, traits = traits, NKX_ID = NKX_ID, N_PERM = N_PERM),
        file.path(out, "modules_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_14_modules.txt"))

message("\n==== Summary ====")
print(t1_summary); print(t2_gw_summary[, 1:9]); print(t3[, c(1:3, 5:10)]); print(t(wg_summary))
print(rw_mod); print(ctx_cor); print(rew_summary[, -(3:4)]); print(vt[, 1:4])
