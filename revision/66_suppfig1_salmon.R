#!/usr/bin/env Rscript
# 66_suppfig1_salmon.R - Supplementary Fig. 1C-D (PCA of the top-500 variable genes; Pearson sample correlation) from
# the ORIGINAL Salmon/tximport counts (core$sal$vst_blind, blind VST of the DESeq2 fit used throughout),
# without the plasmid-flagged genes (salmon_v3/plasmid_flagged_genes.tsv).
# Plot code copied unchanged from 17_markers_qc.R (SuppFig1C/1D candidates); only the input matrix differs.
# Outputs: results/revision_2026/salmon_v3/suppfig1/
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(ggrepel); library(dplyr); library(patchwork)
  library(ComplexHeatmap); library(circlize); library(grid)})
proj <- "/mnt/ssd1/Projects/NKX3_wd"
out <- file.path(proj, "results/revision_2026/salmon_v3/suppfig1"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
SEED <- 17; set.seed(SEED)
FONT <- "Arial"   # resolved by fontconfig to Liberation Sans (metric-compatible) for PDF/PNG; SVG keeps "Arial"
# register "Arial" (Helvetica metrics) for the PostScript-based null devices some packages use to measure text
if (!"Arial" %in% names(grDevices::pdfFonts())) grDevices::pdfFonts(Arial = grDevices::pdfFonts()$Helvetica)
if (!"Arial" %in% names(grDevices::postscriptFonts())) grDevices::postscriptFonts(Arial = grDevices::postscriptFonts()$Helvetica)
grp_lv  <- c("vector_ctrl", "vector_CM", "NKX3.1_ctrl", "NKX3.1_CM")
grp_lab <- c(vector_ctrl = "Control-basal", vector_CM = "Control-inflammation",
             NKX3.1_ctrl = "NKX3.1-basal", NKX3.1_CM = "NKX3.1-inflammation")
grp_col <- c("Control-basal" = "#909FC9", "Control-inflammation" = "#E28DBB",
             "NKX3.1-basal" = "#57C5D1", "NKX3.1-inflammation" = "#F5906A")
# the prescribed group colours fail CVD separation (validated), so every scatter also encodes group by shape
grp_shape <- c("Control-basal" = 21, "Control-inflammation" = 22, "NKX3.1-basal" = 24, "NKX3.1-inflammation" = 23)
col_up <- "#1A7F1A"; col_dn <- "#A50F15"                       # deck colours for up / down marks
div_lo <- "#C51B7D"; div_mid <- "#F0EFEC"; div_hi <- "#1A7F1A"  # diverging scale: CVD-validated green/magenta pair
U_MINUS <- "\u2212"; U_TIMES <- "\u00d7"; U_GE <- "\u2265"; U_LE <- "\u2264"

theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(line = element_line(linewidth = 0.3), axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3), axis.text = element_text(size = base - 0.5, colour = "black"),
          axis.title = element_text(size = base), plot.title = element_text(size = base + 0.5, face = "bold"),
          plot.subtitle = element_text(size = base - 1), plot.caption = element_text(size = base - 1.5, hjust = 0),
          legend.title = element_text(size = base - 0.5), legend.text = element_text(size = base - 1),
          legend.key.size = unit(2.8, "mm"), strip.background = element_blank(),
          strip.text = element_text(size = base - 0.5, face = "bold"),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = NA, colour = NA))
}

FIGS <- character(0); TABS <- character(0)
save_gg <- function(p, name, w, h) {
  f <- file.path(out, name)
  ggsave(paste0(f, ".pdf"), p, width = w, height = h, units = "in", device = cairo_pdf, bg = "white")
  ggsave(paste0(f, ".png"), p, width = w, height = h, units = "in", dpi = 300, device = "png", type = "cairo", bg = "white")
  ggsave(paste0(f, ".svg"), p, width = w, height = h, units = "in", device = svglite::svglite, bg = "white")
  FIGS <<- c(FIGS, paste0(f, ".pdf"))
  invisible(f)
}
save_grid <- function(draw_fun, name, w, h) {
  f <- file.path(out, name)
  cairo_pdf(paste0(f, ".pdf"), width = w, height = h, family = FONT); draw_fun(); invisible(dev.off())
  png(paste0(f, ".png"), width = w, height = h, units = "in", res = 300, type = "cairo", family = FONT, bg = "white")
  draw_fun(); invisible(dev.off())
  svglite::svglite(paste0(f, ".svg"), width = w, height = h, bg = "white"); draw_fun(); invisible(dev.off())
  FIGS <<- c(FIGS, paste0(f, ".pdf"))
  invisible(f)
}
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
cd <- core$cd; samples <- cd$sample_id
cd$label <- factor(cd$label, levels = unname(grp_lab)); cd$short <- paste0("#", cd$num)
pca_top <- function(mat, ntop = 500) {
  rv <- matrixStats::rowVars(mat)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]
  p <- prcomp(t(mat[sel, , drop = FALSE]), center = TRUE, scale. = FALSE)
  list(p = p, sel = rownames(mat)[sel], ve = p$sdev^2 / sum(p$sdev^2))
}
pdf_scores <- function(P, k1 = 1, k2 = 2) {
  x <- P$p$x; sid <- rownames(x)
  data.frame(sample_id = sid, a = x[, k1], b = x[, k2], group = factor(cd[sid, "label"], levels = unname(grp_lab)),
             short = cd[sid, "short"])
}
# identical scale + guide specs in every panel so that patchwork merges the group legends
grp_scales <- function(rows = 1) list(
  scale_fill_manual(values = grp_col, name = NULL, drop = FALSE),
  scale_shape_manual(values = grp_shape, name = NULL, drop = FALSE),
  guides(fill = guide_legend(nrow = rows), shape = guide_legend(nrow = rows)))
pca_plot <- function(P, title, subtitle = NULL, k1 = 1, k2 = 2, legend = TRUE, hl = NULL, rows = 1) {
  d <- pdf_scores(P, k1, k2)
  p <- ggplot(d, aes(a, b)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(fill = group, shape = group), size = 2.1, stroke = 0.35) +
    geom_text_repel(aes(label = short), size = 5.4 / .pt, family = FONT, seed = SEED, box.padding = 0.18,
                    point.padding = 0.12, min.segment.length = 0.2, segment.size = 0.15, max.overlaps = 50) +
    grp_scales(rows) +
    labs(x = sprintf("PC%d (%.1f%%)", k1, 100 * P$ve[k1]), y = sprintf("PC%d (%.1f%%)", k2, 100 * P$ve[k2]),
         title = title, subtitle = subtitle) +
    theme_pub(7) + theme(legend.position = if (legend) "bottom" else "none", plot.title.position = "plot")
  if (!is.null(hl)) p <- p + geom_point(data = d[d$sample_id %in% hl, ], shape = 1, size = 3.6, stroke = 0.35, colour = "grey20")
  p
}
FLAG <- read.delim(file.path(proj, "results/revision_2026/salmon_v3/plasmid_flagged_genes.tsv"))$gene_id
vb <- core$sal$vst_blind[, samples]
vb <- vb[!rownames(vb) %in% FLAG, ]   # plasmid-flagged genes excluded, as in the gene-set analyses
P5 <- pca_top(vb, 500)
p_s1c <- pca_plot(P5, NULL, NULL, legend = TRUE, rows = 2) +
  labs(title = NULL, subtitle = NULL) + theme(legend.text = element_text(size = 6.2))
save_gg(p_s1c, "SuppFig1C_PCA_top500_salmon", 3.4, 3.4)
cm5 <- cor(vb[P5$sel, samples])
lab_hm <- cd[samples, "short"]
dimnames(cm5) <- list(lab_hm, lab_hm)
hc <- hclust(as.dist(1 - cm5), method = "complete")
write.table(data.frame(library = lab_hm, sample_id = samples, group = cd[samples, "label"], round(cm5, 4),
                     cluster_order = match(lab_hm, lab_hm[hc$order]), check.names = FALSE), file.path(out, "sample_correlation_top500_salmon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
col_fun <- circlize::colorRamp2(c(min(cm5), (min(cm5) + 1) / 2, 1), c("#F7FBFF", "#6BAED6", "#08306B"))
gp_txt <- function(sz, face = "plain") gpar(fontfamily = FONT, fontsize = sz, fontface = face)
grp_vec <- as.character(cd[samples, "label"])
draw_hm <- function() {
  ha_top <- HeatmapAnnotation(Group = grp_vec, col = list(Group = grp_col), show_annotation_name = FALSE,
                              simple_anno_size = unit(1.8, "mm"),
                              annotation_legend_param = list(Group = list(title = "Group", title_gp = gp_txt(6.5, "bold"),
                                                                          labels_gp = gp_txt(6), at = unname(grp_lab), nrow = 2,
                                                                          grid_height = unit(2.6, "mm"), grid_width = unit(2.6, "mm"))))
  ha_left <- rowAnnotation(Group = grp_vec, col = list(Group = grp_col), show_annotation_name = FALSE,
                           simple_anno_size = unit(1.8, "mm"), show_legend = FALSE)
  hm <- Heatmap(cm5, name = "Pearson r", col = col_fun, cluster_rows = hc, cluster_columns = hc,
                row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
                row_names_gp = gp_txt(6), column_names_gp = gp_txt(6), column_names_rot = 90,
                top_annotation = ha_top, left_annotation = ha_left, rect_gp = gpar(col = "white", lwd = 0.4),
                heatmap_legend_param = list(title_gp = gp_txt(6.5, "bold"), labels_gp = gp_txt(6), direction = "horizontal",
                                            legend_width = unit(22, "mm"), grid_height = unit(2.4, "mm"),
                                            at = round(seq(min(cm5), 1, length.out = 3), 2)),
                width = unit(12 * 3.4, "mm"), height = unit(12 * 3.4, "mm"))
  draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE,
       legend_gap = unit(4, "mm"), padding = unit(c(2, 3, 2, 3), "mm"))
}
save_grid(draw_hm, "SuppFig1D_correlation_salmon", 3.4, 3.3)
ve <- data.frame(PC = paste0("PC", 1:11), pct_variance = round(100 * P5$ve[1:11], 2))
write.table(ve, file.path(out, "pca_variance_explained_salmon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sc <- data.frame(sample_id = rownames(P5$p$x), group = cd[rownames(P5$p$x), "label"], round(P5$p$x[, 1:3], 3))
write.table(sc, file.path(out, "pca_scores_salmon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(ve[1:4, ]); print(sc)
cat("min / max off-diagonal r:", round(range(cm5[upper.tri(cm5)]), 4), "\n")
