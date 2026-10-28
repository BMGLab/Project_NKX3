#!/usr/bin/env Rscript
# =============================================================================
# 50_verify_fig3_03_loo_figure.R -- task verify_round2/fig3 (adversarial check of fig3_alternative)
#
#  1. Leave-one-library-out Hallmark GSEA for every library of each candidate panel (DESeq2 refit on 11
#     libraries, same design/ranking/GSEA settings as the corrected analysis; validated in script 02); all five
#     focal sets are recorded for each panel contrast (panel_set flags the sets drawn in that panel),
#     plus GC-corrected without #12 for panel C. Records NES, gene-permutation FDR and the median log2FC of
#     the displayed leading-edge genes.
#  2. Diagnostic figure (v03_fig3_verification_diagnostics): A panel-C score vs library GC-bias slope;
#     B MYC/OXPHOS (pink) module eigengene vs % reads outside exons; C NES of the panel sets as in Fig 2,
#     GC-corrected, and leaving out each library in turn.
# Needs the outputs of 50_verify_fig3_02_recompute.R. Inputs read-only.
# Run: /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/50_verify_fig3_03_loo_figure.R
# =============================================================================
if (!isTRUE(l10n_info()[["UTF-8"]])) for (loc in c("en_US.UTF-8", "C.UTF-8"))
  if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(msigdbr); library(clusterProfiler); library(BiocParallel)
  library(ggplot2); library(patchwork); library(ggrepel)
})
options(stringsAsFactors = FALSE, warn = 1)
register(SerialParam()); BPP <- MulticoreParam(16)
t0 <- Sys.time()
proj <- "/mnt/ssd1/Projects/NKX3_wd"; revd <- file.path(proj, "results/revision_2026"); out <- file.path(revd, "verify_round2/fig3")
FONT <- "Liberation Sans"
GRP_COL <- c("Control-basal" = "#909FC9", "Control-inflammation" = "#E28DBB", "NKX3.1-basal" = "#57C5D1", "NKX3.1-inflammation" = "#F5906A")
log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
wtsv <- function(x, f) fwrite(as.data.table(x), file.path(out, f), sep = "\t", quote = FALSE, na = "NA")
umin <- function(x) gsub("-", "−", format(x, trim = TRUE))

core <- readRDS(file.path(revd, "core/core.rds"))
cd <- core$cd; ann <- core$ann; S <- cd$sample_id; cd$short <- paste0("#", cd$num)
gtab <- fread(file.path(revd, "markers_qc/gene_gc_length.tsv")); GC <- setNames(gtab$gc, gtab$gene_id)
H_tab <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
T2G <- as.data.frame(unique(H_tab[!is.na(ncbi_gene) & ncbi_gene != "", .(term = gs_name, gene = as.character(ncbi_gene))]))
FOCAL <- c(MYC = "HALLMARK_MYC_TARGETS_V1", OXPHOS = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           E2F = "HALLMARK_E2F_TARGETS", G2M = "HALLMARK_G2M_CHECKPOINT", AJ = "HALLMARK_APICAL_JUNCTION")
PANELS <- list(A = list(contrast = "CM_in_vector", sets = FOCAL[c("MYC", "OXPHOS")], groups = c("vector_ctrl", "vector_CM")),
               B = list(contrast = "NKX_in_ctrl",  sets = FOCAL[c("E2F", "G2M")],   groups = c("vector_ctrl", "NKX3.1_ctrl")),
               C = list(contrast = "NKX_in_CM",    sets = FOCAL["AJ"],              groups = c("vector_CM", "NKX3.1_CM")))
le_tab <- fread(file.path(revd, "fig3_alternative/fig3alt_leading_edge_genes.tsv"))
disp <- lapply(setNames(nm = names(PANELS)), function(pn) le_tab[panel == pn & displayed == TRUE]$gene_id)

dds0 <- core$fc$dds_full; cnt <- counts(dds0)[, S]; sf_all <- sizeFactors(dds0)[S]
coldat <- data.frame(row.names = S, NKX = factor(as.character(cd$NKX), levels = c("vector", "NKX3.1")),
                     CM = factor(as.character(cd$CM), levels = c("ctrl", "CM")))
gc_nf <- function(samples) {                                   # per-library loess GC offsets (overall reference), as script 02
  m <- cnt[, samples]; sf <- estimateSizeFactorsForMatrix(m)
  g <- GC[rownames(m)]; g[is.na(g)] <- median(g, na.rm = TRUE)
  nm <- sweep(m, 2, sf, "/"); y <- log2(nm + 1); dev <- y - rowMeans(y)
  train <- rowMeans(nm) >= 50 & !is.na(GC[rownames(m)]); rng <- range(g[train])
  fit <- sapply(samples, function(s) predict(loess(dev[train, s] ~ g[train], span = 0.3, degree = 1, family = "symmetric"),
                                             newdata = pmin(pmax(g, rng[1]), rng[2])))
  nf <- sweep(2^fit, 2, sf, "*"); nf / exp(rowMeans(log(nf)))
}
fit_contrasts <- function(samples, nf = NULL) {
  d <- DESeqDataSetFromMatrix(cnt[, samples], coldat[samples, ], ~ NKX + CM + NKX:CM)
  if (!is.null(nf)) normalizationFactors(d) <- nf
  d <- DESeq(d, quiet = TRUE, parallel = TRUE, BPPARAM = BPP)
  list(CM_in_vector = results(d, name = "CM_CM_vs_ctrl"), NKX_in_ctrl = results(d, name = "NKX_NKX3.1_vs_vector"),
       NKX_in_CM = results(d, contrast = list(c("NKX_NKX3.1_vs_vector", "NKXNKX3.1.CMCM"))))
}
ranks_of <- function(r) {
  d <- data.table(gene_id = rownames(r), stat = r$stat)[, entrez := ann[gene_id, "entrez"]][!is.na(entrez) & is.finite(stat)]
  d <- d[order(-abs(stat))][!duplicated(entrez)]; sort(setNames(d$stat, d$entrez), decreasing = TRUE)
}
run_gsea <- function(r) {
  set.seed(42)
  as.data.table(as.data.frame(GSEA(geneList = ranks_of(r), TERM2GENE = T2G, exponent = 1, minGSSize = 15, maxGSSize = 500, eps = 0,
    pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, seed = FALSE, by = "fgsea")))[, .(ID, NES, p.adjust)]
}
rows <- list()
add_rows <- function(fits, analysis, dropped, panels) for (pn in panels) {
  p <- PANELS[[pn]]; r <- fits[[p$contrast]]; g <- run_gsea(r)
  for (s in unname(FOCAL)) rows[[length(rows) + 1]] <<- data.table(analysis = analysis, dropped = dropped, panel = pn, contrast = p$contrast,
    set = s, panel_set = s %in% p$sets,
    NES = round(g[ID == s]$NES, 3), FDR = signif(g[ID == s]$p.adjust, 3), n_sets_FDR_lt_0.05 = sum(g$p.adjust < 0.05),
    median_log2FC_displayed = round(median(r[intersect(disp[[pn]], rownames(r)), "log2FoldChange"], na.rm = TRUE), 3))
}
# VERIFY_FIG3_FIGURE_ONLY=1 redraws the figure from the saved table during development only; the recorded run recomputes.
if (Sys.getenv("VERIFY_FIG3_FIGURE_ONLY") == "1" && file.exists(file.path(out, "v03_gsea_leave_one_library_out.tsv"))) {
  loo <- fread(file.path(out, "v03_gsea_leave_one_library_out.tsv")); log_msg("DEV: leave-one-out table read from disk")
} else {
log_msg("full-data fits (plain, GC-corrected)")
add_rows(fit_contrasts(S), "as Fig 2 (corrected)", "none", names(PANELS))
add_rows(fit_contrasts(S, gc_nf(S)), "GC-corrected", "none", names(PANELS))
for (s in S) {
  lib <- cd$short[S == s]; grp <- as.character(cd$group[S == s])
  pans <- names(PANELS)[vapply(PANELS, function(p) grp %in% p$groups, TRUE)]
  log_msg("leave out %s (panels %s)", lib, paste(pans, collapse = ","))
  add_rows(fit_contrasts(setdiff(S, s)), "leave one library out", lib, pans)
}
s12 <- setdiff(S, "24004-12_S15")
add_rows(fit_contrasts(s12, gc_nf(s12)), "GC-corrected, without #12", "#12", "C")
loo <- rbindlist(rows)
ref <- loo[analysis == "as Fig 2 (corrected)", .(panel, set, NES_full = NES)]
loo <- merge(loo, ref, by = c("panel", "set"))
loo_sum <- loo[analysis == "leave one library out", .(contrast = contrast[1], panel_set = panel_set[1], NES_full = NES_full[1],
  NES_loo_min = min(NES), NES_loo_max = max(NES), library_giving_max = dropped[which.max(NES)], FDR_at_max = FDR[which.max(NES)],
  n_loo_FDR_lt_0.05 = sum(FDR < 0.05), n_loo = .N), by = .(panel, set)]
print(loo[, .(analysis, dropped, panel, set = sub("HALLMARK_", "", set), NES, FDR, median_log2FC_displayed)])
print(loo_sum)
wtsv(loo[order(panel, set, analysis, dropped)], "v03_gsea_leave_one_library_out.tsv")
wtsv(loo_sum, "v03_gsea_leave_one_library_out_summary.tsv")
}

# -----------------------------------------------------------------------------
# figure (panel sets only; the other focal sets of each contrast are in the tables)
# -----------------------------------------------------------------------------
log_msg("figure")
ps <- fread(file.path(out, "v02_C_panel_scores_by_library.tsv"))
pk <- fread(file.path(out, "v02_B_pink_module_by_library.tsv"))
mt <- fread(file.path(out, "v02_B_module_design_tests.tsv"))[network == "sensitivity (17,405 genes)" & module == "pink"]
gcr <- fread(file.path(out, "v02_C_panel_gc_confound.tsv"))[panel == "C"]
th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.text = element_text(colour = "black", size = 6.5), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), plot.title = element_text(size = 7.5, face = "bold"),
        plot.subtitle = element_text(size = 6.3, colour = "grey25"), plot.title.position = "plot",
        legend.text = element_text(size = 6.5), plot.background = element_rect(fill = "white", colour = NA))
ps[, group := factor(group, levels = names(GRP_COL))]; pk[, group := factor(group, levels = names(GRP_COL))]
pA <- ggplot(ps, aes(gc_slope, C, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_point(shape = 21, size = 2.3, stroke = 0.3, colour = "grey20") +
  geom_text_repel(aes(label = library), size = 2.1, family = FONT, min.segment.length = 0, segment.size = 0.2,
                  seed = 3, box.padding = 0.35, point.padding = 0.3, force = 2, max.overlaps = Inf) +
  scale_fill_manual(values = GRP_COL, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(shape = 21, size = 2.3))) +
  scale_x_continuous(labels = umin) + scale_y_continuous(labels = umin) +
  labs(title = "A  Apical-junction genes (panel C) vs library GC bias",
       subtitle = sprintf("r = %s across 12 libraries, %s within groups",
                          umin(sprintf("%.2f", gcr$r_score_vs_gc_slope_12lib)), umin(sprintf("%.2f", gcr$r_score_vs_gc_slope_within_group))),
       x = "Library GC-bias slope (log2 per 10% GC)", y = "Mean row-centred VST, 30 panel C genes") + th
pB <- ggplot(pk, aes(pct_reads_outside_exons, pink_eigengene, fill = group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_point(aes(shape = CM), size = 2.3, stroke = 0.3, colour = "grey20") +
  geom_text_repel(aes(label = library), size = 2.1, family = FONT, min.segment.length = 0, segment.size = 0.2,
                  seed = 3, box.padding = 0.25, point.padding = 0.15, max.overlaps = Inf) +
  scale_fill_manual(values = GRP_COL, name = NULL, drop = FALSE, guide = "none") +
  scale_shape_manual(values = c(ctrl = 21, CM = 24), guide = "none") +
  scale_y_continuous(labels = umin) +
  labs(title = "B  MYC/OXPHOS WGCNA module (pink, 1,018 genes)",
       subtitle = sprintf("r with CM %s (exact P = %s, the smallest possible with 12 libraries)\nr with %% reads outside exons %s across libraries, %s within groups",
                          umin(sprintf("%.2f", mt$r_CM)), sprintf("%.4f", mt$p_exact924_CM), umin(sprintf("%.2f", mt$r_pct_reads_outside_exons)),
                          umin(sprintf("%.2f", mt$r_within_group_pct_outside))),
       x = "Reads outside annotated exons (%)", y = "Module eigengene (triangles: CM)") + th
SET_LAB <- c(HALLMARK_MYC_TARGETS_V1 = "A: MYC targets V1", HALLMARK_OXIDATIVE_PHOSPHORYLATION = "A: OXPHOS",
             HALLMARK_E2F_TARGETS = "B: E2F targets", HALLMARK_G2M_CHECKPOINT = "B: G2M checkpoint", HALLMARK_APICAL_JUNCTION = "C: Apical junction")
fl <- copy(loo)[analysis != "GC-corrected, without #12" & panel_set == TRUE]
lev <- rev(SET_LAB)
fl[, an := factor(analysis, levels = c("as Fig 2 (corrected)", "GC-corrected", "leave one library out"))]
fl[, y := match(SET_LAB[set], lev) + c(0.22, 0, -0.22)[as.integer(an)]]
fl[, fdr_ok := factor(ifelse(FDR < 0.05, "FDR < 0.05", "FDR \u2265 0.05"), levels = c("FDR < 0.05", "FDR \u2265 0.05"))]
lab_pts <- fl[an == "leave one library out", .SD[which.max(abs(NES - NES_full))], by = set]
pC <- ggplot(fl, aes(NES, y)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_hline(yintercept = seq_along(lev) + 0.5, linewidth = 0.2, colour = "grey92") +
  geom_point(data = fl[an == "leave one library out"], aes(shape = fdr_ok), colour = "grey45", size = 1.5, stroke = 0.5) +
  geom_point(data = fl[an != "leave one library out"], aes(colour = an, shape = fdr_ok), size = 2.1, stroke = 0.7) +
  geom_text(data = lab_pts, aes(label = paste("without", dropped)), size = 2.1, family = FONT, colour = "grey25",
            vjust = 1.9) +
  scale_colour_manual(values = c("as Fig 2 (corrected)" = "#000000", "GC-corrected" = "#E69F00"), name = NULL,
                      labels = c("All 12 libraries (as Fig 2)", "GC-corrected")) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR \u2265 0.05" = 1), name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-3.05, 1), breaks = seq(-3, 1, 0.5), labels = umin) +
  scale_y_continuous(breaks = seq_along(lev), labels = unname(lev), limits = c(0.25, length(lev) + 0.5), expand = c(0, 0)) +
  labs(title = "C  Hallmark GSEA of the candidate-panel sets (gene-permutation FDR, as Fig 2)",
       subtitle = "grey: each library of the contrast left out in turn (DESeq2 refit on 11 libraries); label = largest change",
       x = "NES (treatment minus reference)", y = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)), shape = guide_legend(order = 2)) + th +
  theme(legend.position = "bottom", axis.ticks.y = element_blank(), axis.line.y = element_blank())
fig <- (pA + pB + plot_layout(guides = "collect") & theme(legend.position = "bottom")) / pC + plot_layout(heights = c(1, 0.95))
for (ext in c("pdf", "png")) {
  f <- file.path(out, paste0("v03_fig3_verification_diagnostics.", ext))
  if (ext == "pdf") cairo_pdf(f, width = 7, height = 6.4, family = FONT) else png(f, width = 7, height = 6.4, units = "in", res = 300, type = "cairo", family = FONT)
  print(fig); invisible(dev.off())
}
writeLines(c(sprintf("50_verify_fig3_03_loo_figure.R finished %s (%.1f min)", format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "mins"))),
             "", capture.output(sessionInfo())), file.path(out, "sessionInfo_50_verify_fig3_03.txt"))
log_msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
