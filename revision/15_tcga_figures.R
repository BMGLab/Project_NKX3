#!/usr/bin/env Rscript
# =============================================================================
# 15_tcga_figures.R -- figures for the independent-dataset check (TCGA-PRAD + public RWPE1 NKX3.1 data)
# Inputs : results/revision_2026/tcga/tcga_results.rds (15_tcga.R), rwpe1_results.rds and
#          tables/rwpe1_<GSE>_limma_all_genes.tsv.gz (15_tcga_geo_rwpe1.R)
# Outputs: results/revision_2026/tcga/figures/tcga_fig{A..F}_*.{pdf,png,svg} + two combined pages
# Run    : /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_figures.R
# Colour : reported direction uses the deck's up/down pair (green #1A7F1A / dark red #A50F15). That pair is not
#          separable for deuteranopes (validator dE 4.9), so direction is also carried by facet, triangle orientation
#          and text; colour is never the only channel.
# Sign   : treatment minus reference everywhere (NKX3.1 minus vector control); TCGA = correlation with NKX3-1.
# =============================================================================
# figure_qa (2026-10-21): labels need a UTF-8 character locale (plotmath minus signs, non-ASCII marks);
# a clean shell (env -i) starts in the C locale and renders them as "<U+2212>" or "..". Rendering only.
if (!isTRUE(l10n_info()[["UTF-8"]])) for (loc in c("en_US.UTF-8", "C.UTF-8"))
  if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(data.table); library(scales); library(svglite) })
set.seed(20261021)
out     <- "/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/tcga"
fig_dir <- file.path(out, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
R1 <- readRDS(file.path(out, "tcga_results.rds"))
R2 <- readRDS(file.path(out, "rwpe1_results.rds"))

FONT <- "Arial"                     # resolved to Liberation Sans (metric-identical) by fontconfig for PDF/PNG
COL_UP <- "#1A7F1A"; COL_DN <- "#A50F15"; INK <- "#222222"; MUTED <- "#666666"
MINUS <- "−"
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 2)),
          plot.subtitle = element_text(size = 6.2, colour = MUTED, hjust = 0, margin = margin(b = 3), lineheight = 1.05),
          plot.title.position = "plot",
          strip.background = element_blank(), strip.text = element_text(size = 6.6, face = "bold", colour = INK),
          axis.line = element_line(linewidth = 0.3, colour = INK), axis.ticks = element_line(linewidth = 0.3, colour = INK),
          axis.text = element_text(size = 6.2, colour = INK), axis.title = element_text(size = 6.6, colour = INK),
          legend.text = element_text(size = 6.2), legend.title = element_text(size = 6.2, face = "bold"),
          legend.key.size = unit(7, "pt"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(4, 6, 4, 4))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white",
         device = "png", type = "cairo")
  f_svg <- file.path(fig_dir, paste0(name, ".svg"))
  svglite(f_svg, width = w, height = h, bg = "white"); print(p); invisible(dev.off())
  # name Arial first so PowerPoint/Illustrator use Arial; Liberation Sans (metric-identical) as fallback
  s <- readLines(f_svg, warn = FALSE)
  writeLines(gsub('font-family: "Liberation Sans"', 'font-family: Arial, "Liberation Sans", Helvetica, sans-serif', s, fixed = TRUE), f_svg)
}
fmt_minus <- function(s) gsub("-", MINUS, s, fixed = TRUE)
num_lab <- function(x) fmt_minus(ifelse(abs(x) < 1e-9, "0", sub("^(-?)0\\.", "\\1.", sprintf("%.1f", x))))
int_lab <- function(x) fmt_minus(format(x, trim = TRUE))
n_tum <- R1$prov$n_primary_tumour_used; n_nor <- R1$prov$n_solid_normal_used; n_expr <- R1$prov$n_expressed_tumour

# ---- shared gene ordering / labels (panels A and F) ----
CT <- as.data.table(R1$CT)
CT[, tag := fifelse(in_module_3D & in_module_3E, " [D,E]", fifelse(in_module_3D, " [D]", fifelse(in_module_3E, " [E]", "")))]
CT[, glabel := paste0(symbol, fifelse(in_abstract, "*", ""), tag)]
CAT_LEV <- c("Reported activated", "Repressed", "Fig 3D/E module only")
CT[, cat_short := factor(fcase(category == "Activated (text/Fig 1)", CAT_LEV[1],
                               category == "Repressed (text/Fig 1)", CAT_LEV[2], default = CAT_LEV[3]), levels = CAT_LEV)]
DIR_LEV <- c("Reported up with NKX3.1", "Reported down with NKX3.1")
CT[, dir := factor(fifelse(expected_sign > 0, DIR_LEV[1], DIR_LEV[2]), levels = DIR_LEV)]
setorder(CT, cat_short, rho_raw)
CT[, glabel := factor(glabel, levels = unique(glabel))]

dir_scales <- list(
  scale_shape_manual(values = setNames(c(24, 25), DIR_LEV), name = NULL, drop = FALSE),
  scale_colour_manual(values = setNames(c(COL_UP, COL_DN), DIR_LEV), name = NULL, drop = FALSE),
  scale_fill_manual(values = c(setNames(c(COL_UP, COL_DN), DIR_LEV), none = "white"), breaks = c(DIR_LEV[1], "none"),
                    labels = c("Expressed / passes filter", "Below expression filter (open)"), name = NULL, drop = FALSE),
  guides(shape = guide_legend(order = 1, override.aes = list(fill = c(COL_UP, COL_DN), colour = c(COL_UP, COL_DN))),
         colour = "none",
         fill = guide_legend(order = 2, override.aes = list(shape = 24, colour = INK, fill = c(INK, "white")))))

gene_panel <- function(D, bands, xlab, title, subtitle) {
  ggplot(D, aes(y = glabel)) +
    geom_rect(data = bands, aes(xmin = q025, xmax = q975, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#EDEDED") +
    geom_rect(data = bands, aes(xmin = q25, xmax = q75, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#D2D2D2") +
    geom_vline(data = bands, aes(xintercept = q50), colour = "#8F8F8F", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
    { if ("lo" %in% names(D)) geom_linerange(aes(xmin = lo, xmax = hi, colour = dir), linewidth = 0.3, alpha = 0.7, na.rm = TRUE) } +
    geom_point(aes(x = x, shape = dir, colour = dir, fill = fillcol), size = 1.25, stroke = 0.45, na.rm = TRUE) +
    dir_scales +
    facet_grid(cat_short ~ col, scales = "free_y", space = "free_y") +
    labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
    theme_pub() +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"),
          strip.text.y = element_text(angle = -90, size = 6.2), panel.spacing.x = unit(7, "pt"),
          panel.spacing.y = unit(3, "pt"), axis.text.y = element_text(size = 5.9))
}

# -----------------------------------------------------------------------------
# A. reported genes vs the genome-wide null (TCGA-PRAD primary tumours)
# -----------------------------------------------------------------------------
ADJ_LAB <- c(raw = "Unadjusted", AR_PROL = "Adjusted for AR +\nproliferation scores")
A_long <- rbind(
  CT[, .(glabel, cat_short, dir, expressed, col = ADJ_LAB[["raw"]], x = rho_raw, lo = rho_raw_lo, hi = rho_raw_hi)],
  CT[, .(glabel, cat_short, dir, expressed, col = ADJ_LAB[["AR_PROL"]], x = rho_AR_PROL, lo = rho_AR_PROL_lo, hi = rho_AR_PROL_hi)])
A_long[, col := factor(col, levels = ADJ_LAB)]
A_long[, fillcol := factor(fifelse(expressed, as.character(dir), "none"), levels = c(DIR_LEV, "none"))]
NL <- as.data.table(R1$NULLS)[cohort == "Primary tumour" & adjustment %in% names(ADJ_LAB)]
NL[, col := factor(ADJ_LAB[adjustment], levels = ADJ_LAB)]
pA <- gene_panel(A_long, NL,
  xlab = sprintf("Spearman ρ with NKX3-1 (TCGA-PRAD primary tumours, n = %d)", n_tum),
  title = "Reported NKX3.1 targets vs NKX3-1 in TCGA-PRAD",
  subtitle = sprintf("Grey: 25-75%% and 2.5-97.5%% of the NKX3-1 correlations of all\n%s expressed genes; line = median. Bars: 95%% CI.\n* Abstract gene; [D]/[E] member of the Fig 3D/3E module.",
                     format(n_expr, big.mark = ","))) +
  scale_x_continuous(breaks = seq(-0.6, 0.6, 0.3), labels = num_lab)
save_fig(pA, "tcga_figA_claimed_genes_vs_genomewide_null", 3.45, 5.8)

# -----------------------------------------------------------------------------
# F. the same genes in NKX3.1-overexpressing RWPE1 cells (public GEO data; limma t, NKX3.1 WT minus vector)
# -----------------------------------------------------------------------------
RW_LAB <- c(GSE81439 = "RWPE1, GSE81439\n(6 vs 6)", GSE115335 = "RWPE1 + paraquat,\nGSE115335 (6 vs 6)")
XLIM_F <- c(-12, 18)
CLr <- rbind(as.data.table(R2$A$cl), as.data.table(R2$B$cl))
CLr <- merge(CLr, CT[, .(symbol, glabel, cat_short, dir)], by = "symbol")
CLr[, col := factor(RW_LAB[dataset], levels = RW_LAB)]
CLr[, fillcol := factor(fifelse(passes_filter %in% TRUE, as.character(dir), "none"), levels = c(DIR_LEV, "none"))]
F_long <- CLr[, .(glabel, cat_short, dir, fillcol, col, t, x = pmin(pmax(t, XLIM_F[1]), XLIM_F[2]))]
F_cap <- F_long[!is.na(t) & (t > XLIM_F[2] | t < XLIM_F[1])]
F_bands <- rbindlist(lapply(names(RW_LAB), function(a) {
  tt <- fread(file.path(out, "tables", sprintf("rwpe1_%s_limma_all_genes.tsv.gz", a)))
  v <- tt[passes_filter == TRUE & symbol != "NKX3-1", t]
  q <- quantile(v, c(.025, .25, .5, .75, .975))
  data.table(col = factor(RW_LAB[[a]], levels = RW_LAB), q025 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q975 = q[5])
}))
nkx_fc <- c(R2$A$nkx$logFC, R2$B$nkx$logFC)
pF <- gene_panel(F_long, F_bands,
  xlab = "Moderated t, NKX3.1 WT minus vector (limma)",
  title = "Same genes after NKX3.1 overexpression in RWPE1",
  subtitle = sprintf("Public data (GEO). Grey: 25-75%% and 2.5-97.5%% of t over all\ntested genes. NKX3-1 itself: log2FC +%.1f and +%.1f. No point: gene\nabsent from the GEO matrix. * Abstract gene; [D]/[E] Fig 3D/3E module.",
                     nkx_fc[1], nkx_fc[2])) +
  geom_text(data = F_cap, aes(x = x, y = glabel, label = sprintf("t = %.0f", t)), hjust = 1.25, size = 5.4 / .pt,
            family = FONT, colour = INK) +
  scale_x_continuous(limits = XLIM_F, breaks = c(-10, 0, 10), labels = int_lab)
save_fig(pF, "tcga_figF_claimed_genes_RWPE1_NKX3-1_overexpression", 3.45, 5.8)

# -----------------------------------------------------------------------------
# B. example scatters: Abstract genes vs NKX3-1 (+ AR score as a reference)
# -----------------------------------------------------------------------------
X <- R1$scatter$X; sym <- R1$scatter$sym; idn <- R1$scatter$id_nkx
ab <- c("BCL11B", "TPGS1", "SSC4D", "USP3", "MDM4", "NTNG1", "CSH2", "PIP4P2")
REF <- "AR score (reference)"
B_df <- rbindlist(lapply(ab, function(g) data.table(gene = g, x = X[idn, ], y = X[sym[[g]], ])))
B_df <- rbind(B_df, data.table(gene = REF, x = X[idn, ], y = R1$scatter$ar))
ann <- CT[symbol %in% ab, .(gene = symbol, rho = rho_raw, pctl = pctl_raw, expected_sign, median_TPM)]
ann <- rbind(ann, data.table(gene = REF, rho = R1$prov$context$rho_NKX3_1_vs_AR_score, pctl = NA_real_,
                             expected_sign = 0L, median_TPM = NA_real_))
ord_suffix <- function(p) { r <- round(p); ifelse(r %% 100 %in% 11:13, "th", c("th", "st", "nd", "rd", rep("th", 6))[r %% 10 + 1]) }
ann[, title := fcase(expected_sign > 0, paste0(gene, " (up)"), expected_sign < 0, paste0(gene, " (down)"), default = gene)]
ann[, txt := fmt_minus(fifelse(is.na(pctl), sprintf("ρ = %.2f", rho),
                     sprintf("ρ = %.2f (%.0f%s pct)\nmedian TPM %.1f", rho, pctl, ord_suffix(pctl), median_TPM)))]
lev <- ann$title[match(c(ab, REF), ann$gene)]
B_df[, title := factor(ann$title[match(gene, ann$gene)], levels = lev)]
ann[, title := factor(title, levels = lev)]
pB <- ggplot(B_df, aes(x, y)) +
  geom_point(size = 0.35, alpha = 0.35, colour = "#505050", stroke = 0) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt), hjust = -0.04, vjust = 1.15, size = 5.4 / .pt,
            family = FONT, colour = INK, lineheight = 0.9) +
  facet_wrap(~ title, ncol = 3, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.38)), labels = int_lab) +
  labs(x = "NKX3-1, log2(TPM + 1)", y = "Gene, log2(TPM + 1)   [AR score: mean z]",
       title = sprintf("Abstract genes vs NKX3-1 (TCGA-PRAD tumours, n = %d)", n_tum),
       subtitle = "(up)/(down): direction reported in the manuscript. pct = percentile of ρ\namong the NKX3-1 correlations of all expressed genes.") +
  theme_pub() + theme(strip.text = element_text(size = 6.2), axis.text = element_text(size = 5.8), panel.spacing = unit(4, "pt"))
save_fig(pB, "tcga_figB_example_scatters", 3.45, 3.5)

# -----------------------------------------------------------------------------
# C. Hallmark programme scores vs NKX3-1, against expression-matched random gene-set scores
# -----------------------------------------------------------------------------
PW <- as.data.table(R1$PW); PWC <- as.data.table(R1$PWC)
C_ADJ <- c(raw = "Unadjusted", AR = "AR-adjusted", AR_STROMA = "AR + stroma")
COH_LAB <- c("Primary tumour" = sprintf("Tumours\nn = %d", n_tum), "Solid tissue normal" = sprintf("Normal\nn = %d", n_nor))
PROG <- c("E2F targets", "G2M checkpoint", "Mitotic spindle", "MYC targets V1", "Apical junction")
cdat <- merge(PW[adjustment %in% names(C_ADJ), .(cohort, label, adjustment, rho, ci_lo, ci_hi)],
              PWC[adjustment %in% names(C_ADJ), .(cohort, label, adjustment, null_mean, null_q025, null_q975, percentile)],
              by = c("cohort", "label", "adjustment"))
cdat[, adj := factor(C_ADJ[adjustment], levels = C_ADJ)]
cdat[, coh := factor(COH_LAB[cohort], levels = COH_LAB)]
cdat[, label := factor(label, levels = rev(PROG))]
pC <- ggplot(cdat, aes(y = label)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_linerange(aes(xmin = null_q025, xmax = null_q975), linewidth = 3.2, colour = "#D6D6D6") +
  geom_point(aes(x = null_mean), shape = 124, size = 2.2, colour = "#8F8F8F") +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.4, colour = INK) +
  geom_point(aes(x = rho), size = 1.3, colour = INK) +
  facet_grid(coh ~ adj) +
  scale_x_continuous(limits = c(-0.4, 0.75), breaks = c(-0.3, 0, 0.3, 0.6), labels = num_lab, oob = scales::squish) +
  labs(x = "Spearman ρ, NKX3-1 vs programme score", y = NULL,
       title = "Programme scores vs NKX3-1 in TCGA-PRAD",
       subtitle = "Black: programme score (mean z of expressed Hallmark genes), 95% CI.\nGrey bar: 95% range (tick: mean) of 10,000 random gene sets of equal\nsize matched on expression decile.") +
  theme_pub() + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6.2), panel.spacing.x = unit(6, "pt"),
                      panel.spacing.y = unit(5, "pt"))
save_fig(pC, "tcga_figC_hallmark_scores_vs_NKX3-1", 3.45, 2.65)

# -----------------------------------------------------------------------------
# D. module coherence vs expression-matched random gene sets (primary tumours)
# -----------------------------------------------------------------------------
COH <- as.data.table(R1$COH)[cohort == "Primary tumour" & adjustment == "raw" & scope == "expressed"]
CNL <- as.data.table(R1$COH_NULL)[cohort == "Primary tumour"]
D_sets <- c(module_3D = "Fig 3D\nmodule", module_3E = "Fig 3E\nmodule",
            posctrl_AR_score_genes = "AR genes\n(control)", posctrl_proliferation_genes = "Prolif. genes\n(control)")
CNL <- CNL[gene_set %in% names(D_sets)]; COH <- COH[gene_set %in% names(D_sets)]
lab_of <- function(gs) paste0(D_sets[gs], "\nk = ", COH$n_genes[match(gs, COH$gene_set)])
CNL[, set := factor(lab_of(gene_set), levels = lab_of(names(D_sets)))]
COH[, set := factor(lab_of(gene_set), levels = lab_of(names(D_sets)))]
fmt_p <- function(p) ifelse(p < 0.001, "p < .001", ifelse(p < 0.01, sub("0\\.", ".", sprintf("p = %.3f", p)),
                                                          sub("0\\.", ".", sprintf("p = %.2f", p))))
COH[, lab := sprintf("obs. %s\n%s", sub("^0\\.", ".", sprintf("%.2f", obs_mean_pairwise_rho)), fmt_p(p_perm))]
COH[, hj := fifelse(gene_set %in% c("module_3D", "module_3E"), -0.06, 1.06)]
pD <- ggplot(CNL, aes(x = null)) +
  geom_histogram(bins = 40, fill = "#BDBDBD", colour = NA) +
  geom_vline(data = COH, aes(xintercept = obs_mean_pairwise_rho), colour = INK, linewidth = 0.5) +
  geom_label(data = COH, aes(x = obs_mean_pairwise_rho, y = Inf, label = lab, hjust = hj), vjust = 1.1,
             size = 5.3 / .pt, family = FONT, colour = INK, lineheight = 0.9, fill = "white", linewidth = 0,
             label.padding = unit(1, "pt")) +
  facet_wrap(~ set, nrow = 1, scales = "free") +
  scale_x_continuous(breaks = pretty_breaks(3), labels = function(x) fmt_minus(sub("^(-?)0\\.", "\\1.", format(x, trim = TRUE, drop0trailing = TRUE)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Mean pairwise Spearman ρ among set genes", y = "Random sets",
       title = "Module coherence in TCGA-PRAD primary tumours",
       subtitle = "Grey: 10,000 random sets of equal size matched on expression\ndecile; line: observed set. k = expressed genes in the set.") +
  theme_pub() + theme(strip.text = element_text(size = 6.0), axis.text.y = element_text(size = 5.5),
                      axis.text.x = element_text(size = 5.6), panel.spacing.x = unit(5, "pt"))
save_fig(pD, "tcga_figD_module_coherence", 3.45, 2.0)

# -----------------------------------------------------------------------------
# E. Direction of the five programmes across datasets (preranked fgsea NES)
# -----------------------------------------------------------------------------
HM5 <- c(HALLMARK_E2F_TARGETS = "E2F", HALLMARK_G2M_CHECKPOINT = "G2M", HALLMARK_MITOTIC_SPINDLE = "Mitotic\nspindle",
         HALLMARK_MYC_TARGETS_V1 = "MYC\nV1", HALLMARK_APICAL_JUNCTION = "Apical\njunction")
GS <- as.data.table(R1$GS); FG <- rbind(as.data.table(R2$A$fg), as.data.table(R2$B$fg))
ROWS <- c("DU145 NKX3.1 vs control (basal)" = "Basal",
          "DU145 NKX3.1 vs control (basal), expression trend removed" = "Basal, detrended",
          "DU145 NKX3.1 vs control (inflammation)" = "Inflammation",
          "DU145 NKX3.1 vs control (inflammation), expression trend removed" = "Inflam., detrended",
          "RWPE1 NKX3.1 WT vs vector (GSE81439)" = "GSE81439",
          "RWPE1 NKX3.1 WT vs vector (GSE81439), expression trend removed" = "GSE81439, detrended",
          "RWPE1 NKX3.1 WT vs vector, paraquat (GSE115335)" = "GSE115335 (+PQ)",
          "RWPE1 NKX3.1 WT vs vector, paraquat (GSE115335), expression trend removed" = "GSE115335, detrended",
          "TCGA tumour: Spearman rho with NKX3-1" = "ρ unadjusted",
          "TCGA tumour: partial rho | AR score" = "ρ AR-adjusted",
          "TCGA tumour: partial rho | AR + proliferation + stroma" = "ρ AR+prolif.+stroma")
GRP <- function(r) fifelse(grepl("^DU145", r), "DU145", fifelse(grepl("^RWPE1", r), "RWPE1", "TCGA"))
E_df <- rbind(GS[, .(ranking, pathway, NES, padj)], FG[, .(ranking, pathway, NES, padj)])[pathway %in% names(HM5) & ranking %in% names(ROWS)]
E_df[, row := factor(ROWS[ranking], levels = rev(ROWS))]
E_df[, col := factor(HM5[pathway], levels = HM5)]
E_df[, grp := factor(GRP(ranking), levels = c("DU145", "RWPE1", "TCGA"))]
E_df[, txt := sprintf("%s%.1f%s", ifelse(NES < 0, MINUS, ""), abs(NES), fifelse(padj < 0.05, "*", ""))]
lim <- ceiling(max(abs(E_df$NES)) * 2) / 2
pE <- ggplot(E_df, aes(x = col, y = row, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = txt, colour = abs(NES) > 0.55 * lim), size = 5.6 / .pt, family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient2(low = COL_DN, mid = "#F4F4F4", high = COL_UP, midpoint = 0, limits = c(-lim, lim),
                       name = "NES", breaks = c(-2, -1, 0, 1, 2), labels = int_lab) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Five programmes across datasets (fgsea NES)",
       # figure_qa (2026-10-21): DU145 rows now say which contrast they are (NKX3.1 minus vector); presentation only
       subtitle = "DU145: Wald z (NKX3.1 minus vector), corrected counts; RWPE1:\nlimma t (NKX3.1 WT minus vector); TCGA: correlation with NKX3-1.\nDetrended: expression trend removed by robust loess. * BH q < 0.05.") +
  theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(size = 6, lineheight = 0.9),
        axis.text.y = element_text(size = 6), legend.position = "bottom", legend.key.height = unit(6, "pt"),
        legend.key.width = unit(18, "pt"), legend.title = element_text(vjust = 0.9), strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 6.4), panel.spacing.y = unit(4, "pt"))
save_fig(pE, "tcga_figE_programme_direction_across_datasets", 3.45, 3.3)

# -----------------------------------------------------------------------------
# Combined pages (7.1 in wide): page 1 = A | F ; page 2 = (B / C) | (D / E)
# -----------------------------------------------------------------------------
tagt <- theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
p1 <- (pA + labs(tag = "A") + tagt) | (pF + labs(tag = "F") + tagt)
save_fig(p1, "tcga_figure_combined_page1_genes", 7.1, 5.95)
wr <- function(p, tag) wrap_elements(full = p + labs(tag = tag) + tagt)   # no axis alignment across panels
left  <- wr(pB, "B") / wr(pC, "C") + plot_layout(heights = c(3.5, 2.65))
right <- wr(pD, "D") / wr(pE, "E") + plot_layout(heights = c(2.0, 3.3))
p2 <- wrap_elements(full = left) | wrap_elements(full = right)
save_fig(p2, "tcga_figure_combined_page2_programmes_modules", 7.1, 6.4)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_15_tcga_figures.txt"))
message("figures done")
