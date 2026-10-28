#!/usr/bin/env Rscript
# 62_salmon_named_genes.R - Reviewer 1, Major 1: baseMean, normalised counts and apeglm-shrunken fold changes
# for every gene named in Fig. 1B-D (the 20 labelled genes per panel) or in the text, from the ORIGINAL
# Salmon/tximport counts (core$sal: DESeq2 ~ NKX + CM + NKX:CM on all 12 libraries, apeglm shrinkage).
#   * plasmid flag: salmon_v3/plasmid_flagged_genes.tsv (60_plasmid_flags.py)
#   * paralogue / readthrough families (defined below): member counts summed per
#     library and tested as one feature (rows appended to the same model; size factors of the gene-level fit)
# Outputs (results/revision_2026/salmon_v3/named_genes/):
#   named_genes.tsv, named_genes_normalized_counts.tsv, family_totals.tsv,
#   SuppFig_named_genes_shrunken_FC.{pdf,png}
suppressPackageStartupMessages({library(DESeq2); library(BiocParallel); library(ggplot2); library(readr)
  library(dplyr); library(tidyr); library(patchwork)})
register(MulticoreParam(12)); set.seed(42)
proj <- "/mnt/ssd1/Projects/NKX3_wd"; rv <- file.path(proj, "results/revision_2026")
out <- file.path(rv, "salmon_v3/named_genes"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
core <- readRDS(file.path(rv, "core/core.rds")); ann <- core$ann; cd <- core$cd
flag <- read_tsv(file.path(rv, "salmon_v3/plasmid_flagged_genes.tsv"), show_col_types = FALSE)

# ---- genes named in Fig 1B-D (top-20 labels per panel) and in the text
panels <- c(B = "Fig1B_HM-vector-CM_vs_HM-vector-C", C = "Fig1C_HM-NKX3.1-C_vs_HM-vector-C",
            D = "Fig1D_HM-NKX3.1-CM_vs_HM-vector-CM")
contrast_of <- c(B = "CM_in_vector", C = "NKX_in_ctrl", D = "NKX_in_CM")
lab <- bind_rows(lapply(names(panels), function(p) {
  r <- read_tsv(file.path(rv, "salmon_v3/fig1BCD", paste0(panels[[p]], "_signed_z_results.tsv")), show_col_types = FALSE)
  data.frame(panel = p, rank = 1:20, gene_id = r$gene_id[1:20])
}))
text_genes <- c("ENSG00000167034")                                  # NKX3-1
# genes reported as RNA-seq results in the submitted manuscript (text, Fig. 1E-G; tier A of reliability/named_genes_curated.tsv),
# without the protein-marker genes (CD44, PROM1, CDH2), which are not RNA-seq results
tierA <- read.delim(file.path(rv, "reliability/named_genes_curated.tsv"))
tierA <- tierA$gene_id[startsWith(tierA$tier, "A") & !tierA$symbol %in% c("CD44", "PROM1", "CDH2")]
text_genes <- c(text_genes, tierA)
named <- unique(c(lab$gene_id, text_genes, flag$gene_id))
res <- core$sal$res
tested <- res$CM_in_vector$gene_id
cat("named genes:", length(named), " not in the tested set:", paste(ann[setdiff(named, tested), "symbol"], collapse = ", "), "\n")
named <- intersect(named, tested)

# ---- per-gene table
norm <- core$sal$norm
grp <- cd[colnames(norm), "label"]
wide <- data.frame(gene_id = named, symbol = ann[named, "symbol"], biotype = ann[named, "biotype"])
wide$labelled_in <- sapply(named, function(g) paste(lab$panel[lab$gene_id == g], collapse = ","))
wide$named_in_submitted_text <- ifelse(named %in% tierA, "yes", "")
wide$plasmid_derived <- ifelse(named %in% flag$gene_id, "yes", "")
# paralogue / readthrough groups named by Reviewer 1 (and the GH/CSH cluster, the ASDURF/ASNSD1 upstream ORF)
fam <- list("SLX1A/SLX1B paralogues + SLX1A-SULT1A3, SLX1B-SULT1A4 readthroughs" = c("SLX1A", "SLX1B", "SLX1A-SULT1A3", "SLX1B-SULT1A4"),
            "RGPD1-RGPD8 (RGPD7 not annotated in Ensembl 99)" = c("RGPD1", "RGPD2", "RGPD3", "RGPD4", "RGPD5", "RGPD6", "RGPD8"),
            "COMMD3 + BMI1 + COMMD3-BMI1 readthrough" = c("COMMD3", "BMI1", "COMMD3-BMI1"),
            "GH1/GH2/CSH1/CSH2/CSHL1 cluster" = c("GH1", "GH2", "CSH1", "CSH2", "CSHL1"),
            "AC131160.1 locus (AC131160.1 + MAP6D1 + PARL)" = c("AC131160.1", "MAP6D1", "PARL"),
            "ASDURF + ASNSD1 (upstream ORF of ASNSD1)" = c("ASDURF", "ASNSD1"))
wide$family_or_locus <- sapply(wide$symbol, function(s) paste(names(fam)[sapply(fam, function(m) s %in% m)], collapse = "; "))
wide$baseMean <- res$CM_in_vector$baseMean[match(named, res$CM_in_vector$gene_id)]
for (cn in c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM", "interaction")) {
  r <- res[[cn]][match(named, res[[cn]]$gene_id), ]
  wide[[paste0(cn, "_log2FC_MLE")]] <- r$log2FC_MLE;    wide[[paste0(cn, "_SE_MLE")]] <- r$lfcSE_MLE
  wide[[paste0(cn, "_log2FC_apeglm")]] <- r$log2FC_shrunk; wide[[paste0(cn, "_SE_apeglm")]] <- r$lfcSE_shrunk
  wide[[paste0(cn, "_padj")]] <- r$padj
}
for (l in unique(cd$label)) {
  v <- norm[named, grp == l, drop = FALSE]
  wide[[paste0("mean_", l)]] <- rowMeans(v); wide[[paste0("sd_", l)]] <- apply(v, 1, sd)
}
wide <- wide[order(wide$labelled_in == "", -wide$baseMean), ]
write_tsv(wide, file.path(out, "named_genes.tsv"))
nc <- data.frame(gene_id = named, symbol = ann[named, "symbol"], norm[named, cd$sample_id], check.names = FALSE)
write_tsv(nc, file.path(out, "named_genes_normalized_counts.tsv"))

# ---- paralogue / readthrough family totals (same model, same size factors)
cs <- core$counts_sal
sym_all <- ann[rownames(cs), "symbol"]
famrows <- t(sapply(fam, function(m) colSums(cs[sym_all %in% m, cd$sample_id, drop = FALSE])))
rownames(famrows) <- paste0("FAMILY:", names(fam))
m <- round(rbind(cs[core$sal$keep, cd$sample_id], famrows))
specs <- list(CM_in_vector = list(NKX = c("vector", "NKX3.1"), CM = c("ctrl", "CM"), coef = "CM_CM_vs_ctrl"),
              NKX_in_ctrl  = list(NKX = c("vector", "NKX3.1"), CM = c("ctrl", "CM"), coef = "NKX_NKX3.1_vs_vector"),
              NKX_in_CM    = list(NKX = c("vector", "NKX3.1"), CM = c("CM", "ctrl"), coef = "NKX_NKX3.1_vs_vector"))
sf <- sizeFactors(core$sal$dds_full)[cd$sample_id]
famres <- list()
for (cn in names(specs)) {
  sp <- specs[[cn]]; cdx <- cd
  cdx$NKX <- factor(as.character(cdx$NKX), levels = sp$NKX); cdx$CM <- factor(as.character(cdx$CM), levels = sp$CM)
  d <- DESeqDataSetFromMatrix(m, cdx, ~ NKX + CM + NKX:CM); sizeFactors(d) <- sf
  d <- DESeq(d, quiet = TRUE, parallel = TRUE)
  r <- results(d, name = sp$coef, alpha = 0.05)
  s <- lfcShrink(d, coef = sp$coef, type = "apeglm", res = r, quiet = TRUE, parallel = TRUE)
  k <- grep("^FAMILY:", rownames(r))
  famres[[cn]] <- data.frame(family = sub("^FAMILY:", "", rownames(r)[k]), contrast = cn, baseMean = r$baseMean[k],
                             log2FC_MLE = r$log2FoldChange[k], SE_MLE = r$lfcSE[k],
                             log2FC_apeglm = s$log2FoldChange[k], SE_apeglm = s$lfcSE[k], padj = r$padj[k])
}
famtab <- bind_rows(famres) %>%
  mutate(members = sapply(family, function(f) paste(fam[[f]], collapse = "+")))
famnorm <- sweep(famrows, 2, sf, "/")
famtab <- famtab %>% left_join(data.frame(family = sub("^FAMILY:", "", rownames(famnorm)), round(famnorm, 1),
                                          check.names = FALSE), by = "family")
write_tsv(famtab, file.path(out, "family_totals.tsv"))

# ---- supplementary figure: MLE vs apeglm-shrunken log2FC, with baseMean, for the labelled genes of each panel
ptitle <- c(B = "A   CM − control, vector cells (Fig. 1B)", C = "B   NKX3.1 − vector, basal (Fig. 1C)",
            D = "C   NKX3.1 − vector, CM (Fig. 1D)")
plots <- lapply(names(panels), function(p) {
  cn <- contrast_of[[p]]; r <- res[[cn]]
  g <- lab$gene_id[lab$panel == p]
  d <- r[match(g, r$gene_id), ] %>%
    mutate(symbol = ann[gene_id, "symbol"],
           name = ifelse(gene_id %in% flag$gene_id, paste0(symbol, "†"), symbol),
           name = make.unique(name, sep = " "),
           bm = formatC(baseMean, format = "f", digits = 0, big.mark = ","))
  d$name <- factor(d$name, levels = d$name[order(d$log2FC_MLE)])
  long <- bind_rows(
    data.frame(name = d$name, est = d$log2FC_MLE, se = d$lfcSE_MLE, type = "unshrunken (MLE)"),
    data.frame(name = d$name, est = d$log2FC_shrunk, se = d$lfcSE_shrunk, type = "apeglm-shrunken"))
  long$type <- factor(long$type, levels = c("unshrunken (MLE)", "apeglm-shrunken"))
  xr <- range(c(long$est - 1.96 * long$se, long$est + 1.96 * long$se), na.rm = TRUE)
  ggplot(long, aes(x = est, y = name, colour = type, shape = type)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0, linewidth = 0.4,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 1.8, position = position_dodge(width = 0.6)) +
    geom_text(data = d, aes(x = xr[2] + 0.08 * diff(xr), y = name, label = bm), inherit.aes = FALSE,
              hjust = 0, size = 2.8, colour = "grey20") +
    annotate("text", x = xr[2] + 0.08 * diff(xr), y = length(g) + 0.9, label = "baseMean", hjust = 0, size = 2.8,
             fontface = "bold") +
    scale_colour_manual(values = c("unshrunken (MLE)" = "grey60", "apeglm-shrunken" = "#D55E00"), name = NULL) +
    scale_shape_manual(values = c("unshrunken (MLE)" = 1, "apeglm-shrunken" = 16), name = NULL) +
    scale_x_continuous(limits = c(xr[1], xr[2] + 0.32 * diff(xr)), breaks = pretty(xr), expand = c(0.02, 0)) +
    coord_cartesian(clip = "off", ylim = c(1, length(g) + 0.5)) +
    labs(title = ptitle[[p]], x = expression(log[2] ~ "fold change (± 95% CI)"), y = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          plot.margin = margin(14, 6, 4, 4))
})
fig <- wrap_plots(plots, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
for (dev in c("pdf", "png")) {
  f <- file.path(out, paste0("SuppFig_named_genes_shrunken_FC.", dev))
  if (dev == "pdf") cairo_pdf(f, width = 13, height = 5.8) else png(f, width = 13, height = 5.8, units = "in", res = 300, type = "cairo")
  print(fig); dev.off()
}
cat("written:", out, "\n")
print(as.data.frame(famtab[, c("family", "contrast", "baseMean", "log2FC_MLE", "log2FC_apeglm", "padj")]), digits = 3)
