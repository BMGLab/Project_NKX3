#!/usr/bin/env Rscript
# =============================================================================
# 14_modules_figures.R -- supplementary figure for Reviewer 1, major point 5
# (formal tests of the "NKX3.1 co-expression module" claims). Reads the results
# written by 14_modules.R; run that script first.
#   A  correlation of every Fig 3D/3E module gene with NKX3-1: Fig-3 view (6 samples)
#      versus co-variation beyond group membership (within-group / dose-response)
#   B  module coherence and association with NKX3-1 versus expression-matched
#      random gene sets, in-sample versus held-out context
#   C  dose-response among the 6 NKX3.1-transfected samples for six named genes
#   D  WGCNA module eigengene - design correlations
#   E  selection-aware permutation null for "context-dependent rewiring"
# Run: /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/14_modules_figures.R
# =============================================================================
# figure_qa (2026-10-21): labels need a UTF-8 character locale (plotmath minus signs, non-ASCII marks);
# a clean shell (env -i) starts in the C locale and renders them as "<U+2212>" or "..". Rendering only.
if (!isTRUE(l10n_info()[["UTF-8"]])) for (loc in c("en_US.UTF-8", "C.UTF-8"))
  if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
proj <- "/mnt/ssd1/Projects/NKX3_wd"
out  <- file.path(proj, "results/revision_2026/modules")
set.seed(42)
r  <- readRDS(file.path(out, "modules_results.rds"))
cd <- r$cd; NKX_ID <- r$NKX_ID
is_nkx <- cd$NKX == "NKX3.1"; is_cm <- cd$CM == "CM"

FONT <- "Helvetica"
grp_cols  <- c("Control-basal" = "#909FC9", "Control-inflammation" = "#E28DBB",
               "NKX3.1-basal" = "#57C5D1", "NKX3.1-inflammation" = "#F5906A")
quant_lab <- c(original = "Original (Salmon)", corrected = "Corrected (unique reads)")
th <- theme_classic(base_size = 7, base_family = FONT) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 3)),
        axis.text = element_text(size = 6.5, colour = "black"), axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        legend.key.size = unit(8, "pt"), strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold"), strip.clip = "off", plot.margin = margin(4, 4, 4, 4),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))
save_fig <- function(p, name, w, h) {
  cairo_pdf(file.path(out, paste0(name, ".pdf")), width = w, height = h, family = FONT)   # embeds fonts, keeps Unicode marks
  print(p); invisible(dev.off())
  png(file.path(out, paste0(name, ".png")), width = w, height = h, units = "in", res = 300, type = "cairo", family = FONT)
  print(p); invisible(dev.off())
}
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(signif(p, 2), format = "fg", digits = 2))
perms <- function(v) if (length(v) <= 1) matrix(v, 1) else
  do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], perms(v[-i]))))
within_group_perms <- function(g) {
  idx <- split(seq_along(g), droplevels(factor(g))); P <- lapply(idx, perms)
  grid <- expand.grid(lapply(P, function(m) seq_len(nrow(m))))
  sapply(seq_len(nrow(grid)), function(k) { o <- seq_along(g); for (j in seq_along(idx)) o[idx[[j]]] <- P[[j]][grid[k, j], ]; o })
}

# ---------------------------------------------------------------- gene order / classes
md <- r$module_def
ob <- md[md$module == "basal", ];        ob <- ob[order(ob$heatmap_row), ]
oi <- md[md$module == "inflammation", ]; oi <- oi[order(oi$heatmap_row), ]
g_order <- setdiff(unique(c(ob$gene_id, oi$gene_id)), NKX_ID)
sym <- setNames(md$symbol, md$gene_id)[g_order]; sym <- sym[!duplicated(names(sym))]
class_lab <- c(Fig3D_only = "Fig 3D only", both = "Fig 3D and 3E", Fig3E_only = "Fig 3E only")
gclass <- factor(class_lab[r$gene_class[g_order]], levels = class_lab)
quantifiable_fc <- setNames(r$rs$corrected_quantifiable, r$rs$gene_id)

# =============================================================================
# A. correlation with NKX3-1: Fig-3 view vs beyond-group statistics
# =============================================================================
resid_ctx <- function(V, ctx) {     # within-group residual r with NKX3-1 in one context (6 samples), exact p (36 shuffles)
  s <- which(cd$CM == ifelse(ctx == "basal", "ctrl", "CM")); g <- droplevels(cd$group[s])
  X <- V[, s]; G <- model.matrix(~ 0 + g); R <- X - (X %*% G %*% diag(1 / colSums(G))) %*% t(G)
  R <- R[rowSums(R^2) > 1e-12, , drop = FALSE]; Rz <- R / sqrt(rowSums(R^2)); nz <- Rz[NKX_ID, ]
  rr <- as.vector(Rz %*% nz); PI <- within_group_perms(g)
  pp <- rowMeans(abs(Rz %*% matrix(nz[PI], nrow = length(s))) >= abs(rr) - 1e-9)
  data.frame(gene_id = rownames(R), r = rr, p = pp)
}
colsA <- c(r6_basal = "r, basal\n(Fig 3D view)", r6_infl = "r, inflam.\n(Fig 3E view)",
           res_basal = "within-group r,\nbasal", res_infl = "within-group r,\ninflam.",
           dose = "dose-response r,\nNKX3.1 cells")
A <- do.call(rbind, lapply(c("original", "corrected"), function(q) {
  V <- r$vst[[q]]
  t1 <- r$t1_mod[r$t1_mod$quant == q, ]
  t2 <- r$t2[r$t2$quant == q, ]; t2 <- t2[!duplicated(t2$gene_id), ]
  rb <- resid_ctx(V, "basal"); ri <- resid_ctx(V, "inflammation")
  one <- function(stat, rv, pv) data.frame(quant = q, gene_id = g_order, stat = stat, r = rv, p = pv)
  d <- rbind(
    one("r6_basal", t1$r_pearson[t1$context == "basal"][match(g_order, t1$gene_id[t1$context == "basal"])], NA),
    one("r6_infl",  t1$r_pearson[t1$context == "inflammation"][match(g_order, t1$gene_id[t1$context == "inflammation"])], NA),
    one("res_basal", rb$r[match(g_order, rb$gene_id)], rb$p[match(g_order, rb$gene_id)]),
    one("res_infl",  ri$r[match(g_order, ri$gene_id)], ri$p[match(g_order, ri$gene_id)]),
    one("dose", t2$dose_partial_r[match(g_order, t2$gene_id)], t2$dose_p[match(g_order, t2$gene_id)]))
  if (q == "corrected") d$r[!quantifiable_fc[d$gene_id]] <- NA
  d
}))
A$symbol <- factor(sym[A$gene_id], levels = rev(unname(sym)))
A$class  <- gclass[match(A$gene_id, g_order)]
A$stat   <- factor(colsA[A$stat], levels = colsA)
A$quant  <- factor(quant_lab[A$quant], levels = quant_lab)
A$mark   <- ifelse(is.na(A$r), "×", ifelse(!is.na(A$p) & A$p < 0.05, "•", ""))
pA <- ggplot(A, aes(stat, symbol, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = mark), size = 2.4, family = FONT, colour = "grey15", vjust = 0.45) +
  facet_grid(class ~ quant, scales = "free", space = "free", switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), na.value = "grey85", name = "Pearson r with NKX3-1") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  scale_x_discrete(position = "top") +
  # figure_qa (2026-10-21): caption explains the cell marks (presentation only)
  labs(x = NULL, y = NULL, title = "Correlation of module genes with NKX3-1",
       caption = "\u00d7 not quantifiable with corrected counts (grey)\n\u2022 nominal p < 0.05 (within-group or dose-response r)") + th +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
        axis.text.y = element_text(face = "italic", size = 6), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = 6.5, face = "plain"),
        strip.text.x = element_text(size = 6.5), panel.spacing = unit(2, "pt"),
        legend.position = "bottom", legend.key.width = unit(18, "pt"), legend.key.height = unit(5, "pt"),
        legend.title = element_text(size = 6.5), legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 6, colour = "grey25", hjust = 0, lineheight = 1.05),
        plot.caption.position = "plot")
save_fig(pA, "modules_panelA_correlation_heatmap", 3.5, 5.0)

# =============================================================================
# B. coherence / association with NKX3-1 vs expression-matched random gene sets
# =============================================================================
set_lab <- c(all12_groupcentred = "in-sample", heldout_context_groupcentred = "held-out")
mod_lab <- c(basal = "Fig 3D module", inflammation = "Fig 3E module")
stat_lab <- c(coherence = "Coherence\n(mean pairwise r)", assoc_NKX3_1 = "Association with NKX3-1\n(mean r)")
B <- do.call(rbind, lapply(names(r$t3_null), function(k) {
  parts <- strsplit(k, " ")[[1]]; if (!parts[3] %in% names(set_lab)) return(NULL)
  x <- r$t3_null[[k]]
  do.call(rbind, lapply(rownames(x$null), function(st) {
    nl <- x$null[st, ]
    data.frame(quant = parts[1], module = parts[2], set = parts[3], stat = st, obs = x$obs[[st]],
               lo = unname(quantile(nl, 0.025)), med = median(nl), hi = unname(quantile(nl, 0.975)),
               p = (1 + sum(nl >= x$obs[[st]])) / (length(nl) + 1), n = x$n)
  }))
}))
row_levels <- c("Fig 3D module, in-sample", "Fig 3D module, held-out", "Fig 3E module, in-sample", "Fig 3E module, held-out")
B$row   <- factor(paste(mod_lab[B$module], set_lab[B$set], sep = ", "), levels = rev(row_levels))
B$quant <- factor(c(original = "Original", corrected = "Corrected")[B$quant], levels = c("Original", "Corrected"))
B$stat  <- factor(stat_lab[B$stat], levels = stat_lab)
B$sig   <- B$p < 0.05
B$lab   <- paste0("p = ", fmt_p(B$p))
pB <- ggplot(B, aes(y = row)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 2.4, colour = "grey85", orientation = "y") +
  geom_point(aes(x = med), shape = 124, size = 2.2, colour = "grey40") +
  geom_point(aes(x = obs, fill = sig), shape = 21, size = 1.7, stroke = 0.4, colour = "black") +
  geom_text(aes(x = Inf, label = lab), hjust = 1.05, size = 2.1, family = FONT) +
  facet_grid(quant ~ stat) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4), labels = c("0", "0.2", "0.4")) +
  coord_cartesian(xlim = c(min(B$lo, B$obs) - 0.03, 0.95)) +
  labs(x = "Value (bar: central 95% of 10,000 expression-matched\nrandom gene sets; tick: their median)", y = NULL,
       title = "Coherence beyond group membership",
       subtitle = "in-sample: 12 samples, group-centred\nheld-out: 6 samples of the context not used to pick the module\nfilled circle: p < 0.05 (one-sided, vs random sets)") + th +
  theme(axis.text.y = element_text(size = 6.2), strip.text.x = element_text(size = 6.5), strip.text.y = element_text(size = 6.5),
        panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.3), axis.line = element_blank(),
        axis.title.x = element_text(size = 6.2), plot.subtitle = element_text(size = 6.2, colour = "grey25", margin = margin(b = 3)))
save_fig(pB, "modules_panelB_coherence", 3.5, 3.2)

# =============================================================================
# C. dose-response among the 6 NKX3.1-transfected samples
# =============================================================================
genesC <- c("BCL11B", "SSC4D", "TPGS1", "SRMS", "ZNF331", "LRIG1")
idC <- setNames(r$rs$gene_id, r$rs$symbol)[genesC]
mkC <- function(q, g, first_row, first_col, last_row) {
  V <- if (q == "corrected") r$vst$corrected_allgenes else r$vst$original
  id <- idC[[g]]
  d <- data.frame(x = V[NKX_ID, is_nkx], y = V[id, is_nkx], cm = is_cm[is_nkx], label = cd$label[is_nkx])
  ref <- data.frame(y = c(mean(V[id, !is_nkx & !is_cm]), mean(V[id, !is_nkx & is_cm])),
                    label = c("Control-basal", "Control-inflammation"))
  t2 <- r$t2[r$t2$quant == q & r$t2$symbol == g, ][1, ]
  ok <- !(q == "corrected" && !quantifiable_fc[[id]])
  p <- ggplot(d, aes(x, y)) +
    geom_hline(data = ref, aes(yintercept = y, colour = label), linetype = "22", linewidth = 0.4) +
    geom_blank(data = ref, aes(x = mean(d$x), y = y))
  if (ok) {
    f <- lm(y ~ x + cm, data = d)
    segs <- do.call(rbind, lapply(c(FALSE, TRUE), function(cmv) {
      xx <- range(d$x[d$cm == cmv]); data.frame(x = xx[1], xend = xx[2],
        y = predict(f, data.frame(x = xx[1], cm = cmv)), yend = predict(f, data.frame(x = xx[2], cm = cmv)),
        label = ifelse(cmv, "NKX3.1-inflammation", "NKX3.1-basal"))
    }))
    p <- p + geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend, colour = label), linewidth = 0.5)
    ann <- sprintf("slope %.2f\np = %s", t2$dose_slope, fmt_p(t2$dose_p))
  } else ann <- sprintf("not quantifiable\n(≤%d reads/sample)", r$rs$corrected_max_unique_reads_per_sample[r$rs$symbol == g])
  yr <- range(c(d$y, ref$y)); pad <- max(0.25, diff(yr) * 0.12)
  p + geom_point(aes(fill = label), shape = 21, size = 1.6, stroke = 0.3, colour = "black") +
    annotate("text", x = -Inf, y = Inf, label = ann, hjust = -0.06, vjust = 1.15, size = 2.0, family = FONT, lineheight = 0.9) +
    scale_colour_manual(values = grp_cols, limits = names(grp_cols), name = NULL, drop = FALSE) +
    scale_fill_manual(values = grp_cols, limits = names(grp_cols), name = NULL, drop = FALSE) +
    scale_y_continuous(limits = c(yr[1] - 0.15 * pad, yr[2] + 1.9 * pad), breaks = scales::breaks_pretty(3)) +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    labs(title = if (first_row) g else NULL, x = if (last_row) "NKX3-1 (VST)" else NULL,
         y = if (first_col) paste0(ifelse(q == "original", "Original", "Corrected"), "\ngene (VST)") else NULL) +
    th + theme(plot.title = element_text(face = "bold.italic", size = 7, hjust = 0.5), plot.margin = margin(2, 3, 2, 3),
               axis.text = element_text(size = 6))
}
plC <- list()
for (q in c("original", "corrected")) for (g in genesC)
  plC[[paste(q, g)]] <- mkC(q, g, q == "original", g == genesC[1], q == "corrected")
leg_df <- data.frame(x = rep(1:2, 4), y = rep(1:4, each = 2), label = factor(rep(names(grp_cols), each = 2), levels = names(grp_cols)))
leg_lab <- c("Control-basal (group mean)", "Control-inflammation (group mean)", "NKX3.1-basal (samples, fit)", "NKX3.1-inflammation (samples, fit)")
dummy <- ggplot(leg_df, aes(x, y, colour = label, fill = label, linetype = label, shape = label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.7, stroke = 0.3, na.rm = TRUE) +
  scale_colour_manual(values = grp_cols, labels = leg_lab, name = NULL) + scale_fill_manual(values = grp_cols, labels = leg_lab, name = NULL) +
  scale_linetype_manual(values = c("22", "22", "solid", "solid"), labels = leg_lab, name = NULL) +
  scale_shape_manual(values = c(NA, NA, 21, 21), labels = leg_lab, name = NULL) +
  th + theme(legend.position = "bottom", legend.key.width = unit(16, "pt"), legend.text = element_text(size = 6.5))
leg <- cowplot::get_plot_component(dummy, "guide-box-bottom", return_all = TRUE)
if (is.list(leg) && !inherits(leg, "grob")) leg <- leg[[1]]
plC <- lapply(plC, function(p) p + theme(legend.position = "none"))
pC <- (wrap_plots(plC, nrow = 2) / wrap_elements(full = leg) + plot_layout(heights = c(1, 0.07))) +
  plot_annotation(title = "Dose-response among NKX3.1-transfected samples: lm(gene ~ NKX3-1 + CM), 3 residual df",
                  theme = theme(plot.title = element_text(size = 8, face = "bold", family = FONT, hjust = 0.5)))
save_fig(pC, "modules_panelC_dose_response", 7, 3.1)

# =============================================================================
# D. WGCNA module eigengene - design correlations
# =============================================================================
trait_lab <- c(NKX = "NKX3.1", CM = "CM", NKXxCM = "NKX3.1×CM")
D <- do.call(rbind, lapply(c("original", "corrected"), function(q) {
  w <- r$wg[[q]]; mt <- w$module_trait
  ge <- w$genes[w$genes$gene_id != NKX_ID & w$genes$in_top5000, ]
  ncl <- table(factor(ge$module, levels = unique(mt$module)))
  ordm <- mt$module[mt$trait == "NKX"][order(mt$r[mt$trait == "NKX"])]
  mt$lab <- sprintf("%s (%d)%s", mt$module, mt$size, ifelse(mt$module == w$nkx_module, " ◀ NKX3-1", ""))
  mt$lab <- factor(mt$lab, levels = unique(mt$lab[match(ordm, mt$module)]))
  mt$x <- match(mt$trait, names(trait_lab)); mt$quant <- quant_lab[[q]]
  mt$ncl <- as.integer(ncl[mt$module])
  lib <- mt[mt$trait == "NKX", ]                              # library covariate: one extra column per module
  lib$trait <- "nonexonic"; lib$x <- 4.2; lib$r <- lib$r_with_pct_reads_outside_exons; lib$padj_perm <- NA
  rbind(mt, lib)
}))
D$quant <- factor(D$quant, levels = quant_lab)
D$txt <- paste0(sprintf("%.2f", D$r), ifelse(!is.na(D$padj_perm) & D$padj_perm < 0.05, "*", ""))
Dsw <- unique(D[, c("quant", "lab", "module", "ncl")])
pD <- ggplot(D, aes(x, lab)) +
  geom_tile(aes(fill = r), colour = "white", linewidth = 0.3) +
  geom_text(aes(label = txt), size = 2.1, family = FONT) +
  geom_point(data = Dsw, aes(x = 0.38, y = lab, colour = module), shape = 15, size = 2.3) +
  geom_text(data = Dsw, aes(x = 5.05, y = lab, label = ifelse(ncl > 0, ncl, "–")), size = 2.1, family = FONT) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Eigengene r") +
  scale_x_continuous(breaks = c(1:3, 4.2, 5.05), labels = c(trait_lab, "% reads outside exons", "claimed genes"),
                     expand = expansion(add = c(0.25, 0.25))) +
  facet_wrap(~ quant, scales = "free_y") +
  labs(x = NULL, y = NULL, title = "WGCNA (signed, bicor; top 5,000 variable genes)") + th +
  guides(fill = guide_colourbar(title.position = "left", title.vjust = 1)) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6, angle = 40, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 6), legend.position = "bottom", legend.key.height = unit(5, "pt"),
        legend.key.width = unit(18, "pt"), legend.margin = margin(0, 0, 0, 0), panel.spacing = unit(12, "pt"))
save_fig(pD, "modules_panelD_wgcna_module_trait", 4.6, 2.9)

# =============================================================================
# E. selection-aware permutation null for "rewiring"
# =============================================================================
E <- do.call(rbind, lapply(c("original", "corrected"), function(q) data.frame(quant = quant_lab[[q]], T = r$rew_null[[q]]$T)))
Eo <- data.frame(quant = quant_lab[c("original", "corrected")],
                 T = r$rew_summary$observed_T[match(c("original", "corrected"), r$rew_summary$quant)],
                 p = r$rew_summary$p_selection_aware_T[match(c("original", "corrected"), r$rew_summary$quant)],
                 note = c("published Fig 3D/E modules", "modules re-derived\nfrom corrected counts"))
E$quant <- factor(E$quant, levels = quant_lab); Eo$quant <- factor(Eo$quant, levels = quant_lab)
pE <- ggplot(E, aes(T)) +
  geom_histogram(binwidth = 0.5, fill = "grey80", colour = "white", linewidth = 0.15, boundary = 0) +
  geom_vline(data = Eo, aes(xintercept = T), colour = "#A50F15", linewidth = 0.5) +
  geom_text(data = Eo, aes(x = Inf, y = Inf, label = paste0("observed T = ", sprintf("%.1f", T), "\np = ", fmt_p(p), "\n", note)),
            hjust = 1.02, vjust = 1.1, size = 2.0, family = FONT, lineheight = 0.9) +
  facet_wrap(~ quant, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.75))) +
  labs(x = "Context-specificity T\n(t inflammation-only − t basal-only)", y = "Relabelings (of 400)",
       title = "Rewiring vs selection null") + th
save_fig(pE, "modules_panelE_rewiring_null", 2.4, 2.7)

# =============================================================================
# combined supplementary figure + WGCNA soft-threshold QC
# =============================================================================
fig <- (wrap_elements(full = pA) | wrap_elements(full = pB)) /
  wrap_elements(full = pC) /
  ((wrap_elements(full = pD) | wrap_elements(full = pE)) + plot_layout(widths = c(1.9, 1))) +
  plot_layout(heights = c(5.0, 3.1, 2.9)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 10, face = "bold", family = FONT))
save_fig(fig, "modules_suppfig_combined", 7, 11.0)

sft <- do.call(rbind, lapply(names(r$wg), function(q) data.frame(r$wg[[q]]$sft, qlab = quant_lab[[q]])))
sft$qlab <- factor(sft$qlab, levels = quant_lab)
ps1 <- ggplot(sft, aes(Power, signed_R2, colour = qlab)) +
  geom_hline(yintercept = 0.8, linetype = "22", linewidth = 0.3, colour = "grey50") +
  geom_line(linewidth = 0.4) + geom_point(size = 0.9) +
  geom_point(data = sft[sft$chosen, ], shape = 21, size = 2.4, stroke = 0.6, fill = NA, show.legend = FALSE) +
  scale_colour_manual(values = c("#5E5E5E", "#1A7F1A"), name = NULL) +
  labs(x = "Soft-threshold power", y = "Signed scale-free R²") + th + theme(legend.position = "bottom")
ps2 <- ggplot(sft, aes(Power, mean.k., colour = qlab)) + geom_line(linewidth = 0.4) + geom_point(size = 0.9) +
  scale_y_log10() + scale_colour_manual(values = c("#5E5E5E", "#1A7F1A"), name = NULL) +
  labs(x = "Soft-threshold power", y = "Mean connectivity") + th + guides(colour = "none")
save_fig((ps1 | ps2) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0)),
         "modules_wgcna_soft_threshold", 3.5, 2.1)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_14_modules_figures.txt"))
message("figures written to ", out)
