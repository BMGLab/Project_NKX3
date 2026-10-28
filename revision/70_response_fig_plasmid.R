#!/usr/bin/env Rscript
# 70_response_fig_plasmid.R - Response Figure R1 (Reviewer 1, Major 1): plasmid-derived reads and the genes they reach.
#   A: plasmid sequences (elements summed per plasmid) in the raw reads of all 12 libraries (rootcause/tables/12_plasmid_screen_all_samples.tsv;
#      first 2 M read pairs per library, reads per million screened reads)
#   B: Salmon counts (DESeq2-normalized, original counts = core$sal$norm) of the plasmid-flagged genes that pass the
#      expression filter, per library, with NKX3-1 (the transfected gene) for reference
# Output: results/revision_2026/salmon_v3/response/R1_plasmid_reads_and_flagged_genes.{pdf,png}
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
proj <- "/mnt/ssd1/Projects/NKX3_wd"; rv <- file.path(proj, "results/revision_2026")
out <- file.path(rv, "salmon_v3/response"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
core <- readRDS(file.path(rv, "core/core.rds")); cd <- core$cd; ann <- core$ann
lab_lv <- c("Control-basal", "Control-inflammation", "NKX3.1-basal", "NKX3.1-inflammation")
grp_col <- c("Control-basal" = "#909FC9", "Control-inflammation" = "#E28DBB", "NKX3.1-basal" = "#57C5D1", "NKX3.1-inflammation" = "#F5906A")
lib_lab <- setNames(paste0("#", cd$num), cd$sample_id)

scr <- read.delim(file.path(rv, "rootcause/tables/12_plasmid_screen_all_samples.tsv"), check.names = FALSE)
el <- setdiff(colnames(scr), c("sample", "group", "reads_screened"))
# elements are summed per plasmid (the element-level table is rootcause/tables/12_plasmid_screen_all_samples.tsv)
nkx_el <- c("CMV/T7_5UTR_to_attB1", "attB2_to_FP_ATG", "FP_mid(GFP-family)", "FP_aa203_Y(YFP-type)", "BlastR_bsd_ORF", "NKX3-1_CDS_start")
ctl_el <- setdiff(el, nkx_el)
a <- data.frame(sample = rep(scr$sample, 2),
                element = rep(c("NKX3.1 plasmid", "control plasmid"), each = nrow(scr)),
                rpm = c(rowSums(scr[, nkx_el]), rowSums(scr[, ctl_el]))) %>%
  mutate(lib = factor(lib_lab[sample], levels = lib_lab[cd$sample_id[order(factor(cd$label, levels = lab_lv), cd$num)]]),
         element = factor(element, levels = rev(unique(element))))
pA <- ggplot(a, aes(lib, element, fill = rpm)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(rpm > 0, format(round(rpm, 1), nsmall = 1), "0")), size = 2.3, colour = ifelse(a$rpm > 150, "white", "grey15")) +
  scale_fill_gradient(low = "#F4F4F2", high = "#08306B", name = "reads per million\nscreened reads", trans = "sqrt") +
  labs(title = "A   Plasmid sequences in the raw reads", x = "library (grouped as in Fig. 1)", y = NULL) +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

flag <- read.delim(file.path(rv, "salmon_v3/plasmid_flagged_genes.tsv"))
ids <- c("ENSG00000167034", intersect(flag$gene_id, rownames(core$sal$norm)))
b <- data.frame(gene_id = rep(ids, each = nrow(cd)), sample_id = rep(cd$sample_id, length(ids))) %>%
  mutate(count = core$sal$norm[cbind(gene_id, sample_id)], symbol = ann[gene_id, "symbol"],
         symbol = ifelse(gene_id %in% flag$gene_id, paste0(symbol, "†"), symbol),
         symbol = factor(symbol, levels = unique(symbol)), group = factor(cd[sample_id, "label"], levels = lab_lv))
pB <- ggplot(b, aes(group, count + 1, fill = group)) +
  geom_point(shape = 21, size = 1.9, stroke = 0.3, position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  facet_wrap(~ symbol, nrow = 2, scales = "free_y") + scale_y_log10() +
  scale_fill_manual(values = grp_col, name = NULL) +
  labs(title = "B   Salmon counts of NKX3-1 and the plasmid-flagged genes (†)", x = NULL, y = "normalized count + 1 (log scale)") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom",
                                  plot.title = element_text(face = "bold"), strip.background = element_blank())
fig <- pA / pB + plot_layout(heights = c(0.28, 1.15))
for (dev in c("pdf", "png")) {
  f <- file.path(out, paste0("R1_plasmid_reads_and_flagged_genes.", dev))
  if (dev == "pdf") cairo_pdf(f, width = 7.5, height = 6.6) else png(f, width = 7.5, height = 6.6, units = "in", res = 250, type = "cairo")
  print(fig); dev.off()
}
cat("written", out, "\n")
