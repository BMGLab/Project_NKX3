#!/usr/bin/env Rscript
# 65_pub_figures_salmon_compact_GO.R - generated from 41b_pub_figures_compact_GO.R: identical figure code; input = results/revision_2026/salmon_v3/gsea (Salmon counts,
# plasmid-flagged genes removed; quant "salmon", trend-adjusted "salmonadj"); output = salmon_v3/pub_figures/compact_GO
# 41b: compact variant of 41_pub_figures.R for the GO:BP figure only (up to 5 sets per direction), so that
# the GO panels and the flow-cytometry panels fit on one journal page. Generated from 41_pub_figures.R by sed.
# =============================================================================
# 41_pub_figures.R  -  publication versions of the GSEA NES dot plots
#   Figure 2 (MSigDB Hallmark), GO:BP figure (submitted Fig 4A-C), Supplementary C2 figure
#
# Co-author request: "remove the writings underneath the plots of Figure 2A-B-C and the GO
# figure A-B-C; make them more ready to publication".
#
# CONTENT IS UNCHANGED from round 1 (12_gsea_figures.R):
#   - corrected GSEA tables  results/revision_2026/gsea/tables/gsea_salmon_<coll>_<contrast>.tsv
#     (featureCounts, genes ranked by the DESeq2 Wald statistic, treatment - reference)
#   - same selection: up to 10 sets per direction with padj < 0.05, smallest padj first
#     (ties: larger |NES|), displayed ordered by NES
#   - same marks:  double dagger = not significant (padj >= 0.05, untested, or sign reversed) in the
#                  expression-trend-adjusted run (gsea_fcadj_*);
#                  dagger = significant only because NKX3-1 itself is in the set (verdict re-derived
#                  from tables/circularity_minus_NKX3-1/gsea_salmon_*_minusNKX3-1.tsv with the 12_gsea.R rule)
#   Both the selection and the marks are cross-checked against the round-1 objects
#   (gsea_results_all.rds, circularity_NKX3-1_removed.tsv); the script stops on any difference.
#
# PRESENTATION CHANGES:
#   - no caption / footnote text in the image (it goes into the figure legend)
#   - one shared legend per figure (fill = -log10 adjusted P, size = set size)
#   - Liberation Sans (Arial metrics), 8 pt titles, 7 pt everything else, at final print size
#   - all panels of a figure share the label-column width, the plot width and the x range
#   - set names wrapped at word boundaries with balanced line lengths; no name is abbreviated;
#     row heights follow the number of label lines, so wrapped names never touch their neighbours
#   - single panels 3.45 in wide (no legend; the legend is a separate file), combined figure
#     7.0 in wide (2 x 3.45 in + 0.10 in gutter), height <= 9.5 in
#
# Outputs: results/revision_2026/pub_figures/
#   <fig>_{A,B,C}.{pdf,png}, <fig>_combined.{pdf,png}, <fig>_legend.{pdf,png}
#   with <fig> in Fig2_Hallmark, GO_BP, SuppC2;  pub_figures_sets.tsv (every plotted set),
#   pub_figures_panel_counts.tsv (numbers for the figure legends), pub_figures_label_edits.tsv,
#   pub_figures_crosscheck.tsv, pub_figures_layout.tsv, pub_figures_files.tsv,
#   sessionInfo_65_pub_figures_salmon_compact_GO.txt
#
# Run from a clean shell:
#   /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/41_pub_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(stringr); library(scales)
  library(grid); library(gtable)
})
set.seed(42)
# plotmath (legend title) needs a UTF-8 character locale; a bare shell (env -i) starts in "C"
if (!l10n_info()[["UTF-8"]]) {
  for (loc in c("C.UTF-8", "C.utf8", "en_US.UTF-8", "en_US.utf8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  if (!l10n_info()[["UTF-8"]]) stop("no UTF-8 locale available for LC_CTYPE")
}

proj  <- "/mnt/ssd1/Projects/NKX3_wd"
gdir  <- file.path(proj, "results/revision_2026/salmon_v3/gsea")
tabd  <- file.path(gdir, "tables")
circd <- file.path(tabd, "circularity_minus_NKX3-1")
out   <- file.path(proj, "results/revision_2026/salmon_v3/pub_figures/compact_GO")
work  <- "/mnt/ssd1/Projects/NKX3_wd/work/pub_figures_compact_salmon"
for (d in c(out, work)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# text is measured on the device that will draw it: keep a cairo device open for measuring
MEASURE_DEV <- file.path(work, "measure_device.pdf")
FONT <- "Liberation Sans"                 # Arial-metric; Arial itself is not installed here
cairo_pdf(MEASURE_DEV, width = 7, height = 10, family = FONT)
MEASURE_DEV_ID <- dev.cur()

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
PADJ    <- 0.05
N_SIDE  <- 5   # compact variant: up to 5 sets per direction
FS_TITLE <- 8; FS <- 7                    # pt at final print size
LH      <- 0.90                           # line height of wrapped set names
# R draws the lines of a multi-line string 1.2 x lineheight x fontsize apart (cra[1] = 1.2 ps),
# so each extra label line adds PITCH_IN to its row; a one-line row is FIGS[[coll]]$row1 tall (text + gap)
PITCH_IN <- 1.2 * LH * FS / 72            # 7.56 pt
TOP_IN  <- 0.30                           # header strip (direction labels) inside the panel
PAD_IN  <- 0.04                           # space under the last row
PANEL_W <- 3.45                           # single panel = one journal column
FIG_W   <- 7.00                           # combined figure = two columns
GUTTER  <- FIG_W - 2 * PANEL_W            # 0.10 in between the two columns
ROW_GAP <- 0.14                           # vertical gap between the two panel rows
FIG_HMAX <- 9.5
COL_UP  <- "#1A7F1A"; COL_DN <- "#A50F15" # round-1 green / red
SEQ     <- c("#cde2fb", "#6da7ec", "#256abf", "#0d366b")   # round-1 one-hue ramp (light -> dark)
SIZE_RANGE <- c(1.0, 3.5)                 # ggplot size range for 10-500 genes (size 3.5 = ~2.8 mm dot)
DAGGER  <- intToUtf8(8224); DDAGGER <- intToUtf8(8225)
GEQ     <- intToUtf8(8805); MINUS <- intToUtf8(8722)
LARR    <- intToUtf8(8592); RARR <- intToUtf8(8594)

CONTRASTS <- c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM")
TAGS      <- c("A", "B", "C")
# per figure: file stem, label wrap width (characters) and the height of a one-line row (inches;
# 0.13 in = 9.4 pt for 7-pt text). C2 wraps at 32 characters (round 1: 30) so that the combined
# C2 figure stays under 9.5 in without shortening any set name.
FIGS <- list(

  gobp     = list(stem = "GO_BP_top5",    wrap = 44, row1 = 0.13, name = "GO Biological Process")

)
# set names shortened for display (none needed so far); ID -> displayed text
ABBREV <- c()

CINFO <- fread(file.path(gdir, "contrast_info.tsv"))
DEP_VERDICTS <- c("significant only because of NKX3-1", "sign reverses without NKX3-1",
                  "significant with NKX3-1; too small to test without it")

# ---------------------------------------------------------------------------
# Tables, selection and marks (round-1 rules)
# ---------------------------------------------------------------------------
read_tab <- function(q, coll, cn) {
  f <- file.path(tabd, sprintf("gsea_%s_%s_%s.tsv", q, coll, cn))
  if (!file.exists(f)) stop("missing table: ", f)
  fread(f, sep = "\t", quote = "")
}
select_sets <- function(tab) {             # identical to 12_gsea_figures.R::select_sets
  sig <- tab[p.adjust < PADJ]
  up <- head(sig[NES > 0][order(p.adjust, -NES)], N_SIDE)
  dn <- head(sig[NES < 0][order(p.adjust, NES)], N_SIDE)
  list(sel = rbind(up, dn), n_up = sig[NES > 0, .N], n_dn = sig[NES < 0, .N], n_tested = nrow(tab))
}
# NKX3-1 circularity: re-derive the 12_gsea.R verdict from the minus-NKX3-1 re-run
circ_verdicts <- function(coll, cn, tab) {
  ids <- tab[contains_NKX3_1 == TRUE & p.adjust < PADJ, ID]
  if (!length(ids)) return(data.table(ID = character(), NES_without = numeric(), padj_without = numeric(),
                                      verdict = character()))
  f <- file.path(circd, sprintf("gsea_salmon_%s_%s_minusNKX3-1.tsv", coll, cn))
  if (!file.exists(f)) stop("significant NKX3-1-containing sets but no circularity table: ", f)
  wo <- fread(f, sep = "\t", quote = "")
  d <- merge(tab[ID %in% ids, .(ID, NES_with = NES, padj_with = p.adjust)],
             wo[, .(ID, NES_without = NES, padj_without = p.adjust)], by = "ID", all.x = TRUE)
  d[, verdict := fifelse(is.na(padj_without), "significant with NKX3-1; too small to test without it",
                 fifelse(padj_without >= PADJ, "significant only because of NKX3-1",
                 fifelse(sign(NES_with) == sign(NES_without), "robust to removing NKX3-1",
                         "sign reverses without NKX3-1")))]
  d[, .(ID, NES_without, padj_without, verdict)]
}

sentence_case <- function(x) ifelse(grepl("^[a-z]+[A-Z]", x), x,          # keep rRNA, mRNA, ncRNA ...
                                    paste0(toupper(substr(x, 1, 1)), substring(x, 2)))
base_label <- function(id, desc, coll) {
  x <- switch(coll,
              hallmark = gsub("_", " ", sub("^HALLMARK_", "", id)),
              c2       = gsub("_", " ", id),
              gobp     = sentence_case(desc))
  hit <- id %in% names(ABBREV)
  x[hit] <- ABBREV[id[hit]]
  x
}

# rendered width (inches) of single-line strings at the label font size
text_w <- function(s) vapply(s, function(z) convertWidth(grobWidth(textGrob(z, gp = gpar(fontsize = FS, fontfamily = FONT))),
                                                         "in", valueOnly = TRUE), 1, USE.NAMES = FALSE)
# Balanced wrapping: the number of lines is that of a greedy wrap at `width` characters (as round 1);
# the break points are then chosen to minimise the widest rendered line (mark included on the last
# line), so names break evenly ("MEBARKI HCC PROGENITOR / FZD8CRD DN", not "... FZD8CRD / DN").
balanced_wrap <- function(x, width, mark = "") {
  sfx <- if (nzchar(mark)) paste0(" ", mark) else ""
  k <- str_count(str_wrap(x, width = width), "\n") + 1L
  if (k == 1L) return(paste0(x, sfx))
  w <- strsplit(x, " ", fixed = TRUE)[[1]]; n <- length(w)
  best <- NULL; best_score <- c(Inf, Inf)
  for (cc in combn(n - 1L, k - 1L, simplify = FALSE)) {
    idx <- c(0L, cc, n)
    lines <- vapply(seq_len(k), function(j) paste(w[(idx[j] + 1L):idx[j + 1L]], collapse = " "), "")
    lines[k] <- paste0(lines[k], sfx)
    lw <- text_w(lines)
    score <- c(round(max(lw), 3), round(sum((max(lw) - lw)^2), 4))   # widest line, then raggedness
    if (score[1] < best_score[1] || (score[1] == best_score[1] && score[2] < best_score[2])) {
      best <- lines; best_score <- score
    }
  }
  paste(best, collapse = "\n")
}

panel_data <- function(coll, cn) {
  tab <- read_tab("salmon", coll, cn)
  s   <- select_sets(tab)
  d   <- copy(s$sel)
  adj <- read_tab("salmonadj", coll, cn)[, .(ID, NES_trend_adj = NES, padj_trend_adj = p.adjust)]
  d   <- merge(d, adj, by = "ID", all.x = TRUE, sort = FALSE)
  d[, ddagger := is.na(padj_trend_adj) | padj_trend_adj >= PADJ | sign(NES_trend_adj) != sign(NES)]
  cv  <- circ_verdicts(coll, cn, tab)
  d   <- merge(d, cv, by = "ID", all.x = TRUE, sort = FALSE)
  d[, dagger := !is.na(verdict) & verdict %in% DEP_VERDICTS]
  d[, name := base_label(ID, Description, coll)]
  d[, mark := paste0(fifelse(dagger, DAGGER, ""), fifelse(ddagger, DDAGGER, ""))]
  d[, label := mapply(balanced_wrap, name, FIGS[[coll]]$wrap, mark, USE.NAMES = FALSE)]
  d <- d[order(NES)]
  d[, nl := str_count(label, "\n") + 1L]
  d[, h := FIGS[[coll]]$row1 + (nl - 1L) * PITCH_IN]
  # QA: rendered label height must fit inside its row (no overlap with the neighbouring rows)
  d[, text_h := vapply(label, function(z) convertHeight(grobHeight(textGrob(z, gp = gpar(fontsize = FS,
                         fontfamily = FONT, lineheight = LH))), "in", valueOnly = TRUE), 1)]
  d[, clearance := h - text_h]
  if (any(d$clearance <= 0)) stop(sprintf("label taller than its row in %s %s", coll, cn))
  d[, y := cumsum(h) - h / 2 + PAD_IN]
  d[, dirn := factor(ifelse(NES > 0, "up", "down"), levels = c("up", "down"))]
  d[, `:=`(collection = coll, contrast = cn)]
  attr(d, "H") <- sum(d$h) + PAD_IN + TOP_IN
  attr(d, "counts") <- data.table(collection = coll, contrast = cn, n_tested = s$n_tested,
                                  n_sig = s$n_up + s$n_dn, n_sig_up = s$n_up, n_sig_down = s$n_dn,
                                  n_shown_up = d[NES > 0, .N], n_shown_down = d[NES < 0, .N],
                                  n_dagger = sum(d$dagger), n_ddagger = sum(d$ddagger),
                                  n_nkx_sets_significant = tab[contains_NKX3_1 == TRUE & p.adjust < PADJ, .N],
                                  min_label_clearance_in = round(min(d$clearance), 3))
  d
}

# ---------------------------------------------------------------------------
# Cross-check against round 1 (selection and marks must be identical)
# ---------------------------------------------------------------------------
G1 <- readRDS(file.path(gdir, "gsea_results_all.rds"))
CIRC1 <- fread(file.path(gdir, "circularity_NKX3-1_removed.tsv"), sep = "\t", quote = "")
round1_marks <- function(coll, cn) {
  tab <- G1$ALL[[paste("salmon", coll, cn, sep = "_")]]
  sel <- copy(select_sets(tab)$sel)
  a <- G1$ALL[[paste("salmonadj", coll, cn, sep = "_")]]
  m <- merge(sel[, .(ID, NES)], a[, .(ID, NES_adj = NES, padj_adj = p.adjust)], by = "ID", all.x = TRUE)
  nonrob <- m[is.na(padj_adj) | padj_adj >= PADJ | sign(NES_adj) != sign(NES), ID]
  dep <- G1$CIRC[quant == "salmon" & collection == coll & contrast == cn & verdict %in% DEP_VERDICTS, ID]
  sel[, .(ID, NES, p.adjust, dagger1 = ID %in% dep, ddagger1 = ID %in% nonrob)]
}

PD <- list(); CHECK <- list()
for (coll in names(FIGS)) for (cn in CONTRASTS) {
  d  <- panel_data(coll, cn)
  r1 <- round1_marks(coll, cn)
  m  <- merge(d[, .(ID, NES, p.adjust, dagger, ddagger, verdict)], r1, by = "ID", all = TRUE,
              suffixes = c("", ".r1"))
  cv1 <- CIRC1[quant == "salmon" & collection == coll & contrast == cn, .(ID, verdict_r1 = verdict)]
  m  <- merge(m, cv1, by = "ID", all.x = TRUE)
  ok <- nrow(m) == nrow(d) && nrow(m) == nrow(r1) && !anyNA(m$NES) && !anyNA(m$NES.r1) &&
        isTRUE(all.equal(m$NES, m$NES.r1, tolerance = 1e-9)) &&
        isTRUE(all.equal(m$p.adjust, m$p.adjust.r1, tolerance = 1e-9)) &&
        all(m$dagger == m$dagger1) && all(m$ddagger == m$ddagger1) &&
        all(is.na(m$verdict) == is.na(m$verdict_r1)) && all(m$verdict == m$verdict_r1, na.rm = TRUE)
  CHECK[[length(CHECK) + 1]] <- data.table(collection = coll, contrast = cn, n_sets = nrow(d),
                                           n_sets_round1 = nrow(r1), same_sets_values_marks = ok)
  if (!ok) { print(m); stop(sprintf("selection/marks differ from round 1: %s %s", coll, cn)) }
  PD[[paste(coll, cn, sep = "_")]] <- d
}
CHECK <- rbindlist(CHECK)
print(CHECK)
rm(G1); invisible(gc())

# ---------------------------------------------------------------------------
# Plot pieces
# ---------------------------------------------------------------------------
theme_pub <- function() {
  theme_classic(base_size = FS, base_family = FONT) +
    theme(plot.title = element_text(size = FS_TITLE, face = "bold", margin = margin(b = 1.5)),
          plot.subtitle = element_text(size = FS, colour = "grey30", margin = margin(b = 3)),
          plot.title.position = "plot",
          axis.text = element_text(size = FS, colour = "black"),
          axis.text.y = element_text(lineheight = LH, margin = margin(r = 2.5)),
          axis.text.x = element_text(margin = margin(t = 1.5)),
          axis.title = element_text(size = FS),
          axis.title.x = element_text(margin = margin(t = 3)),
          axis.line = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.x = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(2, "pt"),
          panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.25),
          legend.position = "none",
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(3, 6, 3, 3))
}
minus_lab <- function(x) ifelse(x < 0, paste0(MINUS, abs(x)), as.character(x))

fig_limits <- function(ds) {               # as 12_gsea_figures.R::fig_limits
  sel <- rbindlist(ds, fill = TRUE)
  raw <- max(-log10(sel$p.adjust)); fm <- if (raw > 10) 10 else max(2, ceiling(raw))
  list(fill_max = fm, capped = raw > fm, xmax = max(2.5, ceiling(max(abs(sel$NES)) * 2) / 2 + 0.25))
}
fill_scale <- function(lim) {
  fb <- pretty(c(-log10(PADJ), lim$fill_max), n = 4); fb <- fb[fb >= -log10(PADJ) & fb <= lim$fill_max]
  fl <- as.character(fb); if (lim$capped) fl[length(fl)] <- paste0(GEQ, fl[length(fl)])
  scale_fill_gradientn(colours = SEQ, limits = c(-log10(PADJ), lim$fill_max), breaks = fb, labels = fl,
                       oob = squish, name = bquote(.(MINUS) * log[10] ~ "adjusted" ~ italic(P)),
                       guide = guide_colourbar(order = 1,
                                               theme = theme(legend.key.height = unit(0.95, "in"),
                                                             legend.key.width = unit(0.11, "in"),
                                                             legend.ticks = element_blank())))
}
size_scale <- function() scale_size(range = SIZE_RANGE, limits = c(10, 500), breaks = c(25, 100, 250, 500),
                                    name = "Set size",
                                    guide = guide_legend(order = 2, override.aes = list(fill = "grey70")))

dir_label <- function(lbl, side) {         # "lower with CM" -> "<- lower" / "with CM"
  w <- strsplit(lbl, " ", fixed = TRUE)[[1]]
  first <- w[1]; rest <- paste(w[-1], collapse = " ")
  if (side == "down") paste0(LARR, " ", first, "\n", rest) else paste0(first, " ", RARR, "\n", rest)
}

nes_panel <- function(d, coll, cn, tag, lim, pw = NULL, xt_hjust = 0.5) {
  ci <- CINFO[contrast == cn]
  H  <- attr(d, "H"); xmax <- lim$xmax
  p <- ggplot(d, aes(x = NES, y = y)) +
    geom_vline(xintercept = 0, linetype = "22", linewidth = 0.3, colour = "grey40") +
    geom_segment(aes(x = 0, xend = NES, yend = y, colour = dirn), linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(size = setSize, fill = pmin(-log10(p.adjust), lim$fill_max)),
               shape = 21, colour = "grey15", stroke = 0.25) +
    scale_colour_manual(values = c(up = COL_UP, down = COL_DN), guide = "none") +
    fill_scale(lim) + size_scale() +
    scale_x_continuous(limits = c(-xmax, xmax), breaks = seq(-floor(xmax), floor(xmax)), labels = minus_lab,
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = d$y, labels = d$label, limits = c(0, H), expand = c(0, 0)) +
    annotate("text", x = -xmax, y = H - TOP_IN / 2, label = dir_label(ci$down_label, "down"), hjust = 0,
             vjust = 0.5, size = FS / .pt, colour = COL_DN, family = FONT, fontface = "bold", lineheight = LH) +
    annotate("text", x = xmax, y = H - TOP_IN / 2, label = dir_label(ci$up_label, "up"), hjust = 1,
             vjust = 0.5, size = FS / .pt, colour = COL_UP, family = FONT, fontface = "bold", lineheight = LH) +
    labs(title = paste0(tag, "    ", ci$title), subtitle = ci$subtitle,
         x = "Normalized enrichment score (NES)", y = NULL) +
    theme_pub() +
    theme(panel.heights = unit(H, "in"), axis.title.x = element_text(hjust = xt_hjust, margin = margin(t = 3)))
  if (!is.null(pw)) p <- p + theme(panel.widths = unit(pw, "in"))
  p
}

in_w <- function(u) convertWidth(u, "in", valueOnly = TRUE)
in_h <- function(u) convertHeight(u, "in", valueOnly = TRUE)
grob_w <- function(g) sum(in_w(g$widths))
grob_h <- function(g) sum(in_h(g$heights))
col_of <- function(g, nm) unique(g$layout$l[g$layout$name == nm])
row_of <- function(g, nm) unique(g$layout$t[g$layout$name == nm])

# build the three panels of one figure with one label-column width and one plot width
build_panels <- function(coll) {
  ds  <- lapply(CONTRASTS, function(cn) PD[[paste(coll, cn, sep = "_")]])
  lim <- fig_limits(ds)
  g0  <- lapply(seq_along(ds), function(i) ggplotGrob(nes_panel(ds[[i]], coll, CONTRASTS[i], TAGS[i], lim)))
  ax  <- col_of(g0[[1]], "axis-l"); pc <- col_of(g0[[1]], "panel")
  AW  <- max(vapply(g0, function(g) in_w(g$widths[ax]), 1))
  fixed <- sum(in_w(g0[[1]]$widths[-c(ax, pc)])) + AW
  PW  <- PANEL_W - fixed
  # x-axis title: centred under the plot if it fits within the plot plus the right margin
  tw  <- in_w(grobWidth(textGrob("Normalized enrichment score (NES)", gp = gpar(fontsize = FS, fontfamily = FONT))))
  rmargin <- in_w(unit(6, "pt"))
  hj  <- if (tw <= PW + 2 * rmargin - 0.01) 0.5 else 1
  gs  <- lapply(seq_along(ds), function(i) {
    g <- ggplotGrob(nes_panel(ds[[i]], coll, CONTRASTS[i], TAGS[i], lim, pw = PW, xt_hjust = hj))
    g$widths[ax] <- unit(AW, "in")
    g
  })
  ws <- vapply(gs, grob_w, 1)
  if (any(abs(ws - PANEL_W) > 1e-6)) stop("panel width is not ", PANEL_W, " in: ", paste(round(ws, 4), collapse = ", "))
  list(gs = gs, ds = ds, lim = lim, PW = PW, AW = AW, xtitle_w = tw, xtitle_hjust = hj,
       top_off = vapply(gs, function(g) sum(in_h(g$heights[seq_len(row_of(g, "panel") - 1)])), 1))
}

# shared legend: colour bar and size key side by side, titles on top
legend_grob <- function(coll, lim) {
  d <- rbindlist(lapply(CONTRASTS, function(cn) PD[[paste(coll, cn, sep = "_")]]))
  p <- ggplot(d, aes(NES, y)) +
    geom_point(aes(size = setSize, fill = pmin(-log10(p.adjust), lim$fill_max)), shape = 21,
               colour = "grey15", stroke = 0.25) +
    fill_scale(lim) + size_scale() + theme_pub() +
    theme(legend.position = "right", legend.direction = "vertical", legend.box = "horizontal",
          legend.box.just = "top", legend.justification = c(0, 1), legend.title.position = "top",
          legend.title = element_text(size = FS, margin = margin(b = 3)),
          legend.text = element_text(size = FS), legend.spacing.x = unit(0.22, "in"),
          legend.key.spacing.y = unit(1, "pt"),
          legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
          legend.background = element_blank(), legend.box.background = element_blank())
  g <- ggplotGrob(p)
  leg <- g$grobs[[which(g$layout$name == "guide-box-right")]]
  list(full = leg, w = sum(in_w(leg$widths)), h = sum(in_h(leg$heights)))
}

# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
FILES <- character(0)
save_draw <- function(draw_fun, base, w, h) {
  f <- file.path(out, paste0(base, c(".pdf", ".png")))
  cairo_pdf(f[1], width = w, height = h, family = FONT); draw_fun(); invisible(dev.off())
  png(f[2], width = w, height = h, units = "in", res = 300, type = "cairo", bg = "white", family = FONT)
  draw_fun(); invisible(dev.off())
  if (MEASURE_DEV_ID %in% dev.list()) dev.set(MEASURE_DEV_ID)
  FILES <<- c(FILES, f)
  invisible(f)
}
place <- function(g, x, y_top, w, h, total_h) {
  pushViewport(viewport(x = unit(x, "in"), y = unit(total_h - y_top, "in"), width = unit(w, "in"),
                        height = unit(h, "in"), just = c("left", "top")))
  grid.draw(g); popViewport()
}

LAYOUT <- list()
for (coll in names(FIGS)) {
  stem <- FIGS[[coll]]$stem
  B <- build_panels(coll)
  hs <- vapply(B$gs, grob_h, 1)
  # single panels
  for (i in seq_along(B$gs)) {
    g <- B$gs[[i]]; h <- hs[i]
    save_draw(function() { grid.newpage(); grid.draw(g) }, sprintf("%s_%s", stem, TAGS[i]), PANEL_W, h)
  }
  # legend (also saved on its own, for assembling single panels)
  L <- legend_grob(coll, B$lim)
  pad <- 0.06
  save_draw(function() {
    grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
    place(L$full, pad, pad, L$w, L$h, L$h + 2 * pad)
  }, sprintf("%s_legend", stem), L$w + 2 * pad, L$h + 2 * pad)
  # combined: A | B over C | legend
  top_h <- max(hs[1:2]); leg_x <- PANEL_W + GUTTER + 0.30; leg_y <- top_h + ROW_GAP + B$top_off[3]
  bottom_h <- max(hs[3], B$top_off[3] + L$h)
  total_h <- top_h + ROW_GAP + bottom_h
  LAYOUT[[coll]] <- data.table(figure = stem, panel_width_in = PANEL_W, plot_area_width_in = round(B$PW, 3),
                               label_column_in = round(B$AW, 3), xtitle_width_in = round(B$xtitle_w, 3),
                               xtitle_hjust = B$xtitle_hjust, height_A = round(hs[1], 3), height_B = round(hs[2], 3),
                               height_C = round(hs[3], 3), legend_w = round(L$w, 3), legend_h = round(L$h, 3),
                               combined_w = FIG_W, combined_h = round(total_h, 3),
                               fill_cap = B$lim$fill_max, fill_capped = B$lim$capped, x_limit = B$lim$xmax)
  if (total_h > FIG_HMAX) warning(sprintf("%s combined figure is %.2f in tall (> %.1f in)", stem, total_h, FIG_HMAX))
  gA <- B$gs[[1]]; gB <- B$gs[[2]]; gC <- B$gs[[3]]
  save_draw(function() {
    grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
    place(gA, 0, 0, PANEL_W, hs[1], total_h)
    place(gB, PANEL_W + GUTTER, 0, PANEL_W, hs[2], total_h)
    place(gC, 0, top_h + ROW_GAP, PANEL_W, hs[3], total_h)
    place(L$full, leg_x, leg_y, L$w, L$h, total_h)
  }, sprintf("%s_combined", stem), FIG_W, total_h)
}
LAYOUT <- rbindlist(LAYOUT)
print(LAYOUT)

# ---------------------------------------------------------------------------
# Tables for the legends and the record
# ---------------------------------------------------------------------------
SETS <- rbindlist(lapply(names(FIGS), function(coll) rbindlist(lapply(seq_along(CONTRASTS), function(i) {
  d <- copy(PD[[paste(coll, CONTRASTS[i], sep = "_")]])[order(-NES)]
  d[, .(figure = FIGS[[coll]]$stem, panel = TAGS[i], collection, contrast, row_top_to_bottom = seq_len(.N), ID,
        Description, label_shown = gsub("\n", " ", label, fixed = TRUE),
        label_line_breaks = gsub("\n", " | ", label, fixed = TRUE), label_lines = nl,
        NES = round(NES, 4), padj = signif(p.adjust, 4), setSize,
        dagger_NKX3_1_dependent = dagger, ddagger_trend_sensitive = ddagger,
        NES_trend_adj = round(NES_trend_adj, 4), padj_trend_adj = signif(padj_trend_adj, 4),
        circularity_verdict = verdict, NES_without_NKX3_1 = round(NES_without, 4),
        padj_without_NKX3_1 = signif(padj_without, 4))]
}))))
fwrite(SETS, file.path(out, "pub_figures_sets.tsv"), sep = "\t")
COUNTS <- rbindlist(lapply(names(PD), function(k) attr(PD[[k]], "counts")))
COUNTS[, figure := vapply(collection, function(co) FIGS[[co]]$stem, "")]
COUNTS[, panel := TAGS[match(contrast, CONTRASTS)]]
setcolorder(COUNTS, c("figure", "panel"))
fwrite(COUNTS, file.path(out, "pub_figures_panel_counts.tsv"), sep = "\t")
print(COUNTS)
EDITS <- SETS[, .(figure, panel, ID, label_shown)][ID %in% names(ABBREV)]
if (!nrow(EDITS)) EDITS <- data.table(note = paste(
  "No set name was abbreviated. Displayed names = MSigDB/GO names (Hallmark: HALLMARK_ prefix removed;",
  "Hallmark/C2: underscores replaced by spaces; GO: first letter capitalised), wrapped at word boundaries."))
fwrite(EDITS, file.path(out, "pub_figures_label_edits.tsv"), sep = "\t")
fwrite(CHECK, file.path(out, "pub_figures_crosscheck.tsv"), sep = "\t")
fwrite(LAYOUT, file.path(out, "pub_figures_layout.tsv"), sep = "\t")

if (MEASURE_DEV_ID %in% dev.list()) invisible(dev.off(MEASURE_DEV_ID))
unlink(MEASURE_DEV)
fwrite(data.table(file = FILES, bytes = file.size(FILES)), file.path(out, "pub_figures_files.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_65_pub_figures_salmon_compact_GO.txt"))
message(sprintf("done: %d figure files", length(FILES)))
