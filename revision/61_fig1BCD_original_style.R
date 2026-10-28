#!/usr/bin/env Rscript
# 61_fig1BCD_original_style.R - Fig. 1B-D regenerated with the ORIGINAL code and data
#   data : the submitted Salmon/tximport VST (results/rna_downstreams_out_v15/vsd.rds), unchanged
#   code : statistics and plot copied from scripts/rna_downstreams_salmon.R (run_comparison, PLOT 2),
#          same thresholds (|LFC| >= 0.25, |z| >= 2, mean VST >= 2), colours, size and top-20 labelling
#   changes (Reviewer 1): panel B is computed as CM - control (it was control - CM, opposite to C and D);
#          labels of genes whose Salmon counts come mainly from plasmid reads carry a dagger
#          (salmon_v3/plasmid_flagged_genes.tsv; see 60_plasmid_flags.py); ggrepel seed fixed.
#   The plots show a signed Welch z statistic, so the legend calls them signed z-statistic plots, not volcano plots.
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(ggrepel); library(readr)})
proj <- "/mnt/ssd1/Projects/NKX3_wd"
out_dir <- file.path(proj, "results/revision_2026/salmon_v3/fig1BCD"); dir.create(out_dir, showWarnings = FALSE)
vsd <- readRDS(file.path(proj, "results/rna_downstreams_out_v15/vsd.rds"))
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
flag <- read_tsv(file.path(proj, "results/revision_2026/salmon_v3/plasmid_flagged_genes.tsv"), show_col_types = FALSE)
vsd_mat <- assay(vsd)
sym <- core$ann[rownames(vsd_mat), "symbol"]
sym[is.na(sym) | sym == ""] <- rownames(vsd_mat)[is.na(sym) | sym == ""]
cond_vec <- factor(as.character(colData(vsd)$condition))

min_expr_any <- 2.0; min_abs_lfc <- 0.25; min_abs_z <- 2.0; top_label <- 20
panels <- list(B = c("HM-vector-CM", "HM-vector-C"),     # was c("HM-vector-C", "HM-vector-CM")
               C = c("HM-NKX3.1-C", "HM-vector-C"),
               D = c("HM-NKX3.1-CM", "HM-vector-CM"))
summ <- list()
for (pn in names(panels)) {
  g1n <- panels[[pn]][1]; g2n <- panels[[pn]][2]
  i1 <- which(cond_vec == g1n); i2 <- which(cond_vec == g2n)
  mu1 <- rowMeans(vsd_mat[, i1, drop = FALSE]); mu2 <- rowMeans(vsd_mat[, i2, drop = FALSE])
  v1 <- apply(vsd_mat[, i1, drop = FALSE], 1, var); v2 <- apply(vsd_mat[, i2, drop = FALSE], 1, var)
  n1 <- length(i1); n2 <- length(i2)
  lfc <- mu1 - mu2
  se <- sqrt(pmax(v1, 1e-6)/n1 + pmax(v2, 1e-6)/n2)
  z <- lfc / se
  p <- 2 * pnorm(-abs(z)); q <- p.adjust(p, method = "BH")
  expressed_ok <- pmax(mu1, mu2) >= min_expr_any
  sig <- expressed_ok & (abs(lfc) >= min_abs_lfc) & (abs(z) >= min_abs_z)
  res <- data.frame(gene_id = rownames(vsd_mat), gene = sym, mean_g1 = mu1, mean_g2 = mu2, lfc_vst = lfc, se = se,
                    z = z, pval = p, padj = q, expressed_ok = expressed_ok, significant = sig, stringsAsFactors = FALSE)
  res <- res[order(-abs(res$z), -abs(res$lfc_vst)), ]
  res$plasmid_flag <- res$gene_id %in% flag$gene_id
  res$label <- ifelse(res$plasmid_flag, paste0(res$gene, "†"), res$gene)
  base <- sprintf("Fig1%s_%s_vs_%s_signed_z", pn, g1n, g2n)
  write_tsv(res, file.path(out_dir, paste0(base, "_results.tsv")))
  lab_idx <- head(seq_len(nrow(res)), top_label)
  q <- ggplot(res, aes(x = lfc_vst, y = z, color = significant)) +
    geom_hline(yintercept = c(-min_abs_z, min_abs_z), linetype = 3) +
    geom_vline(xintercept = c(-min_abs_lfc, min_abs_lfc), linetype = 3) +
    geom_point(size = 1) +
    scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70")) +
    labs(title = sprintf("z vs LFC (VST): %s vs %s", g1n, g2n),
         x = sprintf("LFC (VST): %s − %s", g1n, g2n),
         y = "z (Welch approx)") +
    theme_bw(base_size = 12) +
    ggrepel::geom_text_repel(data = res[lab_idx, ], aes(label = label), size = 3, max.overlaps = Inf,
                             box.padding = 0.3, show.legend = FALSE, seed = 42)
  # the submitted Fig. 1 hid the plot title behind white boxes; "_notitle" is the same plot without it, for the figure
  for (v in c("", "_notitle")) for (dev in c("pdf", "png")) {
    f <- file.path(out_dir, paste0(base, v, ".", dev))
    if (dev == "pdf") cairo_pdf(f, width = 7.5, height = 6) else png(f, width = 7.5, height = 6, units = "in", res = 300, type = "cairo")
    print(if (v == "") q else q + labs(title = NULL)); dev.off()
  }
  summ[[pn]] <- data.frame(panel = pn, comparison = paste(g1n, "-", g2n), n_highlighted = sum(sig),
                           labelled = paste(res$label[lab_idx], collapse = ", "))
}
s <- do.call(rbind, summ); write_tsv(s, file.path(out_dir, "Fig1BCD_summary.tsv")); print(s)
# check against the submitted tables: identical statistics (B with sign reversed)
v15 <- c(B = "cmp_HM-vector-C_vs_HM-vector-CM_results.tsv", C = "cmp_HM-NKX3.1-C_vs_HM-vector-C_results.tsv",
         D = "cmp_HM-NKX3.1-CM_vs_HM-vector-CM_results.tsv")
for (pn in names(v15)) {
  o <- read_tsv(file.path(proj, "results/rna_downstreams_out_v15", v15[[pn]]), show_col_types = FALSE)
  m <- read_tsv(file.path(out_dir, sprintf("Fig1%s_%s_vs_%s_signed_z_results.tsv", pn, panels[[pn]][1], panels[[pn]][2])), show_col_types = FALSE)
  sgn <- if (pn == "B") -1 else 1
  # rows are in the same order (|z|, then |LFC|); symbols are not unique (7 duplicates), so compare by position
  cat(pn, ": same gene order:", identical(o$gene, m$gene), " max |z diff| =", max(abs(sgn * o$z - m$z)),
      " same highlighted set:", identical(o$significant, m$significant), "\n")
}
