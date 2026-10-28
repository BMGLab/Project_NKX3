#!/usr/bin/env Rscript
# =============================================================================
# 63_gsea_salmon.R - generated from 12_gsea.R (same methods, parameters and seeds); the main quantification is now
# the ORIGINAL Salmon/tximport counts (core$sal: DESeq2 ~ NKX + CM + NKX:CM, all 12 libraries), with the genes whose
# Salmon counts come mainly from plasmid reads (salmon_v3/plasmid_flagged_genes.tsv) removed from every ranked list.
#   quant "salmon"    = main analysis (flagged genes removed)
#   quant "salmonadj" = expression-trend sensitivity (Wald statistic centred within 20 expression bins)
#   quant "sal"       = Salmon counts without removing the flagged genes (shows the removal does not matter)
#   quant "v15rerun"  = original ranking metric, as in 12_gsea.R
# Outputs: results/revision_2026/salmon_v3/gsea/ (same file layout as results/revision_2026/gsea/)
# The 12_gsea.R header follows.
# =============================================================================
# =============================================================================
# 12_gsea.R  -  Revision task "gsea" (Reviewer 1, major point 2)
#
# Rebuilds every GSEA in the manuscript (Fig 2A-C Hallmark, Fig 4A-C GO:BP,
# Suppl Fig 3B-D MSigDB C2) as preranked GSEA reported with NES and with ONE
# sign convention: treatment minus reference (CM - control; NKX3.1 - vector).
#
#   ranking   : DESeq2 Wald statistic of the factorial model (core$<quant>$res[[contrast]]$stat)
#   IDs       : Ensembl -> Entrez via core$ann$entrez; NA dropped; duplicated Entrez -> max |stat|
#   method    : clusterProfiler::GSEA (Hallmark, whole C2; msigdbr) and clusterProfiler::gseGO
#               (BP, org.Hs.eg.db, ENTREZID); fgsea multilevel, eps = 0, minGSSize = 15,
#               maxGSSize = 500, BH, pvalueCutoff = 1 (every tested set is reported),
#               set.seed(42) before every call (results verified to be reproducible)
#   quant     : fc  = featureCounts recount (main, corrected analysis)
#               sal = original Salmon counts, same model/ranking (sensitivity)
#               fcadj = fc Wald statistic centred within 20 equal-frequency bins of log10(baseMean)
#                       (sensitivity: removes the expression-dependent trend of the statistic)
#               v15rerun = original ranking metric (Welch VST mean difference on the
#                          published blind VST), sign-corrected, original clusterProfiler
#                          defaults (minGSSize 10, eps 1e-10), current gene-set release
#   also      : published v15 tables sign-corrected (NES x -1 for the CM comparison),
#               concordance tables, and a circularity check (NKX3-1 removed from the ranking)
#
# Downstream (sourced at the end): 12_gsea_figures.R (NES dot plots) and
# 12_gsea_claims.R (claim-by-claim direction table).
#
# Run from a clean shell:
#   /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/12_gsea.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)            # class definitions for objects stored in core.rds
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(fgsea)
  library(BiocParallel)
  library(data.table)
  library(parallel)
})
options(stringsAsFactors = FALSE)

proj   <- "/mnt/ssd1/Projects/NKX3_wd"
script_dir <- file.path(proj, "scripts/revision")
out    <- file.path(proj, "results/revision_2026/salmon_v3/gsea")
tabd   <- file.path(out, "tables")
circd  <- file.path(tabd, "circularity_minus_NKX3-1")
work   <- "/mnt/ssd1/Projects/NKX3_wd/work/gsea_salmon"
v15dir <- file.path(proj, "results/rna_downstreams_out_v15")
for (d in c(out, tabd, circd, work)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
register(SerialParam())   # fgsea(nproc = 0) inside clusterProfiler uses bpparam(); we parallelise across runs
N_WORKERS  <- 12          # forked R processes (task limit: <= 16 threads)
PADJ       <- 0.05
GS_MIN     <- 15
GS_MAX     <- 500
NKX_ENTREZ <- "4824"      # NKX3-1
NKX_ENSG   <- "ENSG00000167034"

log_msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))

# ----------------------------------------------------------------------------
# Contrast metadata (labels used in tables and figures)
# ----------------------------------------------------------------------------
CINFO <- data.table(
  contrast  = c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM", "CM_in_NKX", "NKX_main", "CM_main", "interaction"),
  title     = c("CM vs control (vector cells)", "NKX3.1 vs control (basal)", "NKX3.1 vs control (inflammation)",
                "CM vs control (NKX3.1 cells)", "NKX3.1 vs control (main effect)", "CM vs control (main effect)",
                "NKX3.1 × CM interaction"),
  subtitle  = c("Control-inflammation vs Control-basal", "NKX3.1-basal vs Control-basal",
                "NKX3.1-inflammation vs Control-inflammation", "NKX3.1-inflammation vs NKX3.1-basal",
                "additive model ~ NKX + CM, all 12 samples", "additive model ~ NKX + CM, all 12 samples",
                "NKX3.1 effect under CM minus NKX3.1 effect under basal"),
  treatment = c("Control-inflammation", "NKX3.1-basal", "NKX3.1-inflammation", "NKX3.1-inflammation",
                "NKX3.1", "CM", "NKX3.1 under CM"),
  up_label  = c("higher with CM", "higher in NKX3.1", "higher in NKX3.1", "higher with CM",
                "higher in NKX3.1", "higher with CM", "NKX3.1 effect more positive under CM"),
  down_label= c("lower with CM", "lower in NKX3.1", "lower in NKX3.1", "lower with CM",
                "lower in NKX3.1", "lower with CM", "NKX3.1 effect more negative under CM"),
  v15_cmp   = c("HM-vector-C_vs_HM-vector-CM", "HM-NKX3.1-C_vs_HM-vector-C", "HM-NKX3.1-CM_vs_HM-vector-CM",
                NA, NA, NA, NA),
  v15_sign  = c(-1, 1, 1, NA, NA, NA, NA),        # published CM comparison was ranked control - CM
  fig_hallmark = c("Fig 2A", "Fig 2B", "Fig 2C", NA, NA, NA, NA),
  fig_gobp     = c("Fig 4A", "Fig 4B", "Fig 4C", NA, NA, NA, NA),
  fig_c2       = c("Suppl Fig 3B", "Suppl Fig 3C", "Suppl Fig 3D", NA, NA, NA, NA)
)
MAIN3      <- CINFO$contrast[1:3]
ALL7       <- CINFO$contrast
COLLS      <- c("hallmark", "c2", "gobp")
COLL_LABEL <- c(hallmark = "MSigDB Hallmark", c2 = "MSigDB C2 (all)", gobp = "GO:BP")
V15_FILE   <- c(hallmark = "MSIG_Hallmark_GSEA", c2 = "MSIG_C2_GSEA", gobp = "GO_BP_GSEA")
fwrite(CINFO, file.path(out, "contrast_info.tsv"), sep = "\t")

# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
log_msg("reading core.rds")
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
ann  <- core$ann
cd   <- core$cd
FLAG <- fread(file.path(proj, "results/revision_2026/salmon_v3/plasmid_flagged_genes.tsv"))$gene_id
core$salmon <- core$sal
core$salmon$res <- lapply(core$sal$res, function(r) r[!r$gene_id %in% FLAG, ])
log_msg("plasmid-flagged genes removed from the Salmon ranked lists: %d of %d present",
        sum(FLAG %in% core$sal$res$CM_in_vector$gene_id), length(FLAG))

# gene sets: MSigDB via msigdbr (Entrez = ncbi_gene); GO:BP is built inside gseGO from org.Hs.eg.db
msig <- function(coll) as.data.table(msigdbr(species = "Homo sapiens", collection = coll))
H_tab  <- msig("H")
C2_tab <- msig("C2")
T2G <- list(
  hallmark = as.data.frame(unique(H_tab[!is.na(ncbi_gene) & ncbi_gene != "", .(term = gs_name, gene = as.character(ncbi_gene))])),
  c2       = as.data.frame(unique(C2_tab[!is.na(ncbi_gene) & ncbi_gene != "", .(term = gs_name, gene = as.character(ncbi_gene))]))
)
SET_META <- unique(rbind(H_tab[, .(ID = gs_name, subcollection = "H", gs_description)],
                         C2_tab[, .(ID = gs_name, subcollection = gs_subcollection, gs_description)]))
SET_META <- SET_META[!duplicated(ID)]
MSIGDB_VERSION <- unique(H_tab$db_version)
log_msg("msigdbr %s (MSigDB %s): Hallmark %d sets, C2 %d sets",
        as.character(packageVersion("msigdbr")), MSIGDB_VERSION,
        uniqueN(T2G$hallmark$term), uniqueN(T2G$c2$term))

# ----------------------------------------------------------------------------
# Ranked lists
# ----------------------------------------------------------------------------
# corrected pipeline: DESeq2 Wald statistic, treatment - reference
ranks_from_res <- function(res) {
  d <- data.table(gene_id = res$gene_id, stat = res$stat)
  d[, entrez := ann[gene_id, "entrez"]]
  d <- d[!is.na(entrez) & is.finite(stat)]
  d <- d[order(-abs(stat))][!duplicated(entrez)]          # duplicated Entrez: keep max |stat|
  sort(setNames(d$stat, d$entrez), decreasing = TRUE)
}
# original v15 metric: Welch VST mean difference on the published blind VST, now treatment - reference;
# original ID handling: bitr(ENSEMBL -> ENTREZID), duplicated Entrez collapsed by max |difference|
vst15 <- core$vst_sal_v15
grp15 <- as.character(cd[colnames(vst15), "group"])
map15 <- as.data.table(suppressMessages(suppressWarnings(
  bitr(rownames(vst15), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))))
ranks_v15 <- function(g_treat, g_ref) {
  lfc <- rowMeans(vst15[, grp15 == g_treat, drop = FALSE]) - rowMeans(vst15[, grp15 == g_ref, drop = FALSE])
  d <- merge(map15, data.table(ENSEMBL = names(lfc), lfc = lfc), by = "ENSEMBL")
  d <- d[!is.na(ENTREZID) & is.finite(lfc)][order(-abs(lfc))][!duplicated(ENTREZID)]
  sort(setNames(d$lfc, d$ENTREZID), decreasing = TRUE)
}
V15_GROUPS <- list(CM_in_vector = c("vector_CM", "vector_ctrl"),
                   NKX_in_ctrl  = c("NKX3.1_ctrl", "vector_ctrl"),
                   NKX_in_CM    = c("NKX3.1_CM", "vector_CM"))

# expression-trend sensitivity: centre the Wald statistic within 20 equal-frequency bins of log10(baseMean)
N_BINS <- 20
bin_centre <- function(res) {
  lbm <- log10(res$baseMean)
  bin <- cut(lbm, unique(quantile(lbm, seq(0, 1, length.out = N_BINS + 1))), include.lowest = TRUE)
  res$stat <- res$stat - ave(res$stat, bin, FUN = median)
  res
}
RANKS <- list()
for (q in c("salmon", "sal")) for (cn in ALL7) RANKS[[paste(q, cn, sep = ":")]] <- ranks_from_res(core[[q]]$res[[cn]])
for (cn in ALL7) RANKS[[paste("salmonadj", cn, sep = ":")]] <- ranks_from_res(bin_centre(core$salmon$res[[cn]]))

# how the Wald statistic depends on expression level (deciles of log10 baseMean)
TREND <- rbindlist(lapply(c("salmon", "sal"), function(q) rbindlist(lapply(ALL7, function(cn) {
  r <- as.data.table(core[[q]]$res[[cn]])
  r[, decile := cut(log10(baseMean), quantile(log10(baseMean), 0:10 / 10), include.lowest = TRUE, labels = FALSE)]
  r[, .(quant = q, contrast = cn, n_genes = .N, median_baseMean = round(median(baseMean), 1),
        mean_stat = round(mean(stat), 3), median_stat = round(median(stat), 3), frac_stat_negative = round(mean(stat < 0), 3),
        spearman_stat_vs_log10baseMean = round(cor(r$stat, log10(r$baseMean), method = "spearman"), 3)), by = decile][order(decile)]
}))))
fwrite(TREND, file.path(out, "stat_expression_trend.tsv"), sep = "\t")
for (cn in MAIN3) RANKS[[paste("v15rerun", cn, sep = ":")]] <- ranks_v15(V15_GROUPS[[cn]][1], V15_GROUPS[[cn]][2])

rank_summary <- rbindlist(lapply(names(RANKS), function(k) {
  r <- RANKS[[k]]; p <- strsplit(k, ":")[[1]]
  data.table(quant = p[1], contrast = p[2], n_ranked_genes = length(r),
             NKX3_1_rank = match(NKX_ENTREZ, names(r)), NKX3_1_metric = unname(r[NKX_ENTREZ]),
             metric_max = max(r), metric_min = min(r))
}))
fwrite(rank_summary, file.path(out, "ranked_list_summary.tsv"), sep = "\t")
print(rank_summary)

# symbols for leading edges (org.Hs.eg.db; same annotation release as the original analysis)
all_entrez <- unique(c(unlist(lapply(RANKS, names)), T2G$hallmark$gene))
ENTREZ2SYM <- suppressMessages(mapIds(org.Hs.eg.db, keys = all_entrez, column = "SYMBOL",
                                      keytype = "ENTREZID", multiVals = "first"))
to_symbols <- function(x) vapply(strsplit(x, "/", fixed = TRUE), function(e) {
  s <- ENTREZ2SYM[e]; s[is.na(s)] <- e[is.na(s)]; paste(s, collapse = "/") }, character(1))

# ----------------------------------------------------------------------------
# One GSEA run
# ----------------------------------------------------------------------------
run_one <- function(job) {
  ranks <- RANKS[[paste(job$quant, job$contrast, sep = ":")]]
  if (isTRUE(job$drop_nkx)) ranks <- ranks[names(ranks) != NKX_ENTREZ]
  # corrected pipeline: minGSSize 15, eps 0. The v15 re-run keeps the ORIGINAL clusterProfiler defaults
  # (minGSSize 10, eps 1e-10) so that it reproduces the published tables, incl. sets of 10-14 genes.
  gs_min <- if (job$quant == "v15rerun") 10 else GS_MIN
  eps    <- if (job$quant == "v15rerun") 1e-10 else 0
  set.seed(42)
  g <- if (job$collection == "gobp") {
    gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", exponent = 1,
          minGSSize = gs_min, maxGSSize = GS_MAX, eps = eps, pvalueCutoff = 1, pAdjustMethod = "BH",
          verbose = FALSE, seed = FALSE, by = "fgsea")
  } else {
    GSEA(geneList = ranks, TERM2GENE = T2G[[job$collection]], exponent = 1,
         minGSSize = gs_min, maxGSSize = GS_MAX, eps = eps, pvalueCutoff = 1, pAdjustMethod = "BH",
         verbose = FALSE, seed = FALSE, by = "fgsea")
  }
  df <- as.data.table(as.data.frame(g))
  gs <- g@geneSets
  nkx_sets <- names(gs)[vapply(gs, function(x) NKX_ENTREZ %in% x, logical(1))]
  df[, contains_NKX3_1 := ID %in% nkx_sets]
  df[, NKX3_1_in_core := vapply(strsplit(core_enrichment, "/", fixed = TRUE),
                                function(x) NKX_ENTREZ %in% x, logical(1))]
  list(job = job, table = df, n_ranked = length(ranks), nkx_sets = nkx_sets)
}

decorate <- function(df, job) {
  ci <- CINFO[contrast == job$contrast]
  df <- copy(df)
  df[, `:=`(quant = job$quant, collection = job$collection, contrast = job$contrast,
            NKX3_1_removed = isTRUE(job$drop_nkx),
            min_set_size = if (job$quant == "v15rerun") 10L else GS_MIN,
            ranking_metric = switch(job$quant,
                                    v15rerun = "Welch VST mean difference (published v15 VST), treatment - reference",
                                    salmonadj = "DESeq2 Wald statistic centred within 20 expression bins (log10 baseMean), treatment - reference",
                                    "DESeq2 Wald statistic, treatment - reference"))]
  df[, direction := fifelse(NES > 0, ci$up_label, ci$down_label)]
  df[, significant := p.adjust < PADJ]
  df[, n_core := lengths(strsplit(core_enrichment, "/", fixed = TRUE))]
  df[, core_enrichment_symbols := to_symbols(core_enrichment)]
  if (job$collection != "gobp") df <- merge(df, SET_META[, .(ID, subcollection)], by = "ID", all.x = TRUE, sort = FALSE)
  else df[, subcollection := "GO:BP"]
  setorder(df, p.adjust, pvalue)
  setcolorder(df, c("quant", "collection", "contrast", "NKX3_1_removed", "ID", "Description", "subcollection",
                    "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "significant", "direction",
                    "rank", "leading_edge", "n_core", "core_enrichment", "core_enrichment_symbols",
                    "contains_NKX3_1", "NKX3_1_in_core", "ranking_metric"))
  df
}

run_jobs <- function(jobs, label) {
  log_msg("%s: %d GSEA runs on %d workers", label, length(jobs), N_WORKERS)
  # GO:BP / C2 first so the long runs start early
  ord <- order(match(vapply(jobs, `[[`, "", "collection"), c("gobp", "c2", "hallmark")))
  res <- mclapply(jobs[ord], function(j) tryCatch(run_one(j), error = function(e) e),
                  mc.cores = N_WORKERS, mc.preschedule = FALSE)
  res <- res[order(ord)]
  bad <- vapply(res, function(x) inherits(x, "error") || is.null(x), logical(1))
  if (any(bad)) stop("GSEA failed for: ", paste(vapply(jobs[bad], function(j) paste(j, collapse = "/"), ""), collapse = "; "),
                     "\n", paste(vapply(res[bad], function(e) conditionMessage(e), ""), collapse = "\n"))
  res
}

# ----------------------------------------------------------------------------
# Stage 1: corrected (fc), Salmon sensitivity (sal), v15 ranking re-run
# ----------------------------------------------------------------------------
jobs1 <- list()
for (q in c("salmon", "sal", "salmonadj")) for (cn in ALL7) for (co in COLLS)
  jobs1[[length(jobs1) + 1]] <- list(quant = q, contrast = cn, collection = co, drop_nkx = FALSE)
for (cn in MAIN3) for (co in COLLS)
  jobs1[[length(jobs1) + 1]] <- list(quant = "v15rerun", contrast = cn, collection = co, drop_nkx = FALSE)

t0 <- Sys.time()
res1 <- run_jobs(jobs1, "stage 1")
log_msg("stage 1 done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

ALL <- list(); NKX_SETS <- list()
for (r in res1) {
  j <- r$job; key <- paste(j$quant, j$collection, j$contrast, sep = "_")
  ALL[[key]] <- decorate(r$table, j)
  NKX_SETS[[j$collection]] <- union(NKX_SETS[[j$collection]], r$nkx_sets)
  fwrite(ALL[[key]], file.path(tabd, sprintf("gsea_%s_%s_%s.tsv", j$quant, j$collection, j$contrast)),
         sep = "\t", quote = FALSE)
}

# ----------------------------------------------------------------------------
# Published v15 tables, sign-corrected to treatment - reference
# ----------------------------------------------------------------------------
PUB <- list()
for (cn in MAIN3) for (co in COLLS) {
  ci <- CINFO[contrast == cn]
  f <- file.path(v15dir, sprintf("dgea_%s_%s.tsv", ci$v15_cmp, V15_FILE[[co]]))
  p <- fread(f, sep = "\t", colClasses = list(character = c("ID", "core_enrichment")))
  p[, `:=`(quant = "v15published", collection = co, contrast = cn,
           sign_factor = ci$v15_sign, NES_as_published = NES, ES_as_published = enrichmentScore)]
  p[, `:=`(NES = NES * ci$v15_sign, enrichmentScore = enrichmentScore * ci$v15_sign)]
  p[, direction := fifelse(NES > 0, ci$up_label, ci$down_label)]
  p[, significant := p.adjust < PADJ]
  p[, core_enrichment_symbols := to_symbols(core_enrichment)]
  p[, contains_NKX3_1 := ID %in% NKX_SETS[[co]]]
  p[, NKX3_1_in_core := vapply(strsplit(core_enrichment, "/", fixed = TRUE), function(x) NKX_ENTREZ %in% x, logical(1))]
  p[, published_ranking := sprintf("Welch VST mean difference, %s (as published: %s)",
                                   ifelse(ci$v15_sign < 0, "sign flipped here to CM - control", "NKX3.1 - control"),
                                   gsub("_vs_", " minus ", ci$v15_cmp))]
  p[, published_figure := c(hallmark = ci$fig_hallmark, gobp = ci$fig_gobp, c2 = ci$fig_c2)[[co]]]
  setcolorder(p, c("quant", "collection", "contrast", "published_figure", "ID", "Description", "setSize",
                   "NES", "NES_as_published", "sign_factor", "enrichmentScore", "pvalue", "p.adjust",
                   "significant", "direction"))
  PUB[[paste(co, cn, sep = "_")]] <- p
  fwrite(p, file.path(tabd, sprintf("gsea_v15published_%s_%s.tsv", co, cn)), sep = "\t", quote = FALSE)
}

# re-run of the original ranking metric must reproduce the published (sign-corrected) NES
repro <- rbindlist(lapply(MAIN3, function(cn) rbindlist(lapply(COLLS, function(co) {
  p <- PUB[[paste(co, cn, sep = "_")]]; r <- ALL[[paste("v15rerun", co, cn, sep = "_")]]
  m <- merge(p[, .(ID, NES_pub = NES, padj_pub = p.adjust)], r[, .(ID, NES_rerun = NES, padj_rerun = p.adjust)],
             by = "ID", all.x = TRUE)
  data.table(contrast = cn, collection = co, n_published = nrow(p), n_matched = sum(!is.na(m$NES_rerun)),
             same_sign = sum(sign(m$NES_pub) == sign(m$NES_rerun), na.rm = TRUE),
             pearson_NES = if (sum(!is.na(m$NES_rerun)) > 2) round(cor(m$NES_pub, m$NES_rerun, use = "complete.obs"), 3) else NA_real_,
             median_abs_diff_NES = round(median(abs(m$NES_pub - m$NES_rerun), na.rm = TRUE), 3),
             still_padj_lt_0.05 = sum(m$padj_rerun < PADJ, na.rm = TRUE))
}))))
fwrite(repro, file.path(out, "v15_rerun_vs_published.tsv"), sep = "\t")
print(repro)

# ----------------------------------------------------------------------------
# Concordance of NES signs
# ----------------------------------------------------------------------------
get_tab <- function(quant, co, cn) {
  if (quant == "v15published") {
    # published values where listed; sets absent from the published table (n.s. there) filled from the re-run
    p <- PUB[[paste(co, cn, sep = "_")]][, .(ID, NES, padj = p.adjust, source = "published")]
    r <- ALL[[paste("v15rerun", co, cn, sep = "_")]][!ID %in% p$ID, .(ID, NES, padj = p.adjust, source = "v15 re-run (not in published table)")]
    return(rbind(p, r))
  }
  ALL[[paste(quant, co, cn, sep = "_")]][, .(ID, NES, padj = p.adjust, source = quant)]
}
concord_one <- function(qa, qb, co, cn) {
  a <- get_tab(qa, co, cn); b <- get_tab(qb, co, cn)
  sig_a <- if (qa == "v15published") a[source == "published" & padj < PADJ, ID] else a[padj < PADJ, ID]
  sig_b <- if (qb == "v15published") b[source == "published" & padj < PADJ, ID] else b[padj < PADJ, ID]
  u <- union(sig_a, sig_b)
  m <- merge(data.table(ID = u), a[, .(ID, NES_a = NES, padj_a = padj, src_a = source)], by = "ID", all.x = TRUE)
  m <- merge(m, b[, .(ID, NES_b = NES, padj_b = padj, src_b = source)], by = "ID", all.x = TRUE)
  m[, `:=`(sig_a = ID %in% sig_a, sig_b = ID %in% sig_b)]
  cmp <- m[!is.na(NES_a) & !is.na(NES_b)]
  both <- cmp[sig_a & sig_b]
  allc <- merge(a[, .(ID, NES_a = NES)], b[, .(ID, NES_b = NES)], by = "ID")
  list(summary = data.table(
    comparison = sprintf("%s vs %s", qa, qb), collection = co, contrast = cn,
    n_sig_a = length(sig_a), n_sig_b = length(sig_b), n_sig_either = length(u),
    n_comparable = nrow(cmp), n_not_tested_in_one = length(u) - nrow(cmp),
    n_same_sign = sum(sign(cmp$NES_a) == sign(cmp$NES_b)),
    pct_same_sign = if (nrow(cmp)) round(100 * mean(sign(cmp$NES_a) == sign(cmp$NES_b)), 1) else NA_real_,
    n_sig_both = nrow(both), n_sig_both_same_sign = sum(sign(both$NES_a) == sign(both$NES_b)),
    n_sig_both_opposite = sum(sign(both$NES_a) != sign(both$NES_b)),
    n_only_a = sum(m$sig_a & !m$sig_b), n_only_b = sum(!m$sig_a & m$sig_b),
    spearman_NES_all_sets = if (nrow(allc) > 2) round(cor(allc$NES_a, allc$NES_b, method = "spearman"), 3) else NA_real_,
    n_sets_all = nrow(allc)),
    detail = m[, `:=`(comparison = sprintf("%s vs %s", qa, qb), collection = co, contrast = cn)])
}
CONC <- list(); CONC_DET <- list()
for (cn in MAIN3) for (co in COLLS) for (pr in list(c("salmon", "v15published"), c("sal", "v15published"), c("salmon", "sal"), c("salmon", "salmonadj"))) {
  x <- concord_one(pr[1], pr[2], co, cn); CONC[[length(CONC) + 1]] <- x$summary; CONC_DET[[length(CONC_DET) + 1]] <- x$detail
}
for (cn in setdiff(ALL7, MAIN3)) for (co in COLLS) for (pr in list(c("salmon", "sal"), c("salmon", "salmonadj"))) {
  x <- concord_one(pr[1], pr[2], co, cn); CONC[[length(CONC) + 1]] <- x$summary; CONC_DET[[length(CONC_DET) + 1]] <- x$detail
}
CONC <- rbindlist(CONC); CONC_DET <- rbindlist(CONC_DET, fill = TRUE)
setcolorder(CONC_DET, c("comparison", "collection", "contrast"))
fwrite(CONC, file.path(out, "concordance_summary.tsv"), sep = "\t")
fwrite(CONC_DET, file.path(out, "concordance_sets.tsv"), sep = "\t")
print(CONC[, .(comparison, collection, contrast, n_sig_a, n_sig_b, n_comparable, pct_same_sign, n_sig_both, n_sig_both_opposite, spearman_NES_all_sets)])

# ----------------------------------------------------------------------------
# Stage 2: circularity - rerun with NKX3-1 removed from the ranked list wherever a
# significant set (corrected run, or published v15 table) contains NKX3-1
# ----------------------------------------------------------------------------
jobs2 <- list()
for (r in res1) {
  j <- r$job; tab <- ALL[[paste(j$quant, j$collection, j$contrast, sep = "_")]]
  need <- any(tab$significant & tab$contains_NKX3_1)
  if (j$quant == "v15rerun") {
    p <- PUB[[paste(j$collection, j$contrast, sep = "_")]]
    need <- need || any(p$significant & p$contains_NKX3_1)
  }
  if (need) { j$drop_nkx <- TRUE; jobs2[[length(jobs2) + 1]] <- j }
}
res2 <- run_jobs(jobs2, "stage 2 (NKX3-1 removed)")
CIRC <- list()
for (r in res2) {
  j <- r$job; key <- paste(j$quant, j$collection, j$contrast, sep = "_")
  w  <- ALL[[key]]; wo <- decorate(r$table, j)
  fwrite(wo, file.path(circd, sprintf("gsea_%s_%s_%s_minusNKX3-1.tsv", j$quant, j$collection, j$contrast)),
         sep = "\t", quote = FALSE)
  nkx_rank <- rank_summary[quant == j$quant & contrast == j$contrast, NKX3_1_rank]
  ids <- w[contains_NKX3_1 & significant, ID]
  pub <- NULL
  if (j$quant == "v15rerun") {
    pub <- PUB[[paste(j$collection, j$contrast, sep = "_")]]
    ids <- union(ids, pub[contains_NKX3_1 & significant, ID])
  }
  if (!length(ids)) next
  d <- merge(data.table(ID = ids),
             w[, .(ID, Description, setSize, NES_with = NES, padj_with = p.adjust, NKX3_1_in_core,
                   n_core, core_enrichment_symbols)], by = "ID", all.x = TRUE)
  d <- merge(d, wo[, .(ID, setSize_without = setSize, NES_without = NES, padj_without = p.adjust)],
             by = "ID", all.x = TRUE)
  if (!is.null(pub)) d <- merge(d, pub[, .(ID, NES_published_signcorrected = NES, padj_published = p.adjust,
                                           Description_pub = Description)], by = "ID", all.x = TRUE)
  else d[, `:=`(NES_published_signcorrected = NA_real_, padj_published = NA_real_, Description_pub = NA_character_)]
  d[is.na(Description), Description := Description_pub][, Description_pub := NULL]
  d[, `:=`(quant = j$quant, collection = j$collection, contrast = j$contrast, NKX3_1_rank_in_list = nkx_rank)]
  CIRC[[key]] <- d
}
CIRC <- rbindlist(CIRC, fill = TRUE)
CIRC[, verdict := fifelse(is.na(padj_with), "not tested in this run (set size outside the size window)",
                   fifelse(padj_with < PADJ & is.na(padj_without),
                           "significant with NKX3-1; too small to test without it",
                   fifelse(padj_with < PADJ & padj_without >= PADJ,
                           "significant only because of NKX3-1",
                   fifelse(padj_with < PADJ & padj_without < PADJ & sign(NES_with) == sign(NES_without),
                           "robust to removing NKX3-1",
                   fifelse(padj_with >= PADJ, "n.s. with NKX3-1 (significant only in published table)",
                           "sign reverses without NKX3-1")))))]
CIRC[, delta_NES := NES_without - NES_with]
setcolorder(CIRC, c("quant", "collection", "contrast", "ID", "Description", "setSize", "NKX3_1_rank_in_list",
                    "NKX3_1_in_core", "NES_with", "padj_with", "NES_without", "padj_without", "delta_NES", "verdict",
                    "NES_published_signcorrected", "padj_published"))
setorder(CIRC, quant, contrast, collection, padj_with)
fwrite(CIRC, file.path(out, "circularity_NKX3-1_removed.tsv"), sep = "\t")
print(CIRC[, .N, by = .(quant, contrast, collection, verdict)])

# ----------------------------------------------------------------------------
# Counts of significant sets per run
# ----------------------------------------------------------------------------
NSIG <- rbindlist(lapply(names(ALL), function(k) {
  t <- ALL[[k]]
  data.table(quant = t$quant[1], collection = t$collection[1], contrast = t$contrast[1], n_tested = nrow(t),
             n_sig = sum(t$significant), n_sig_up = sum(t$significant & t$NES > 0),
             n_sig_down = sum(t$significant & t$NES < 0))
}))
NSIG <- rbind(NSIG, rbindlist(lapply(names(PUB), function(k) {
  t <- PUB[[k]]
  data.table(quant = "v15published", collection = t$collection[1], contrast = t$contrast[1], n_tested = NA_integer_,
             n_sig = sum(t$significant), n_sig_up = sum(t$significant & t$NES > 0),
             n_sig_down = sum(t$significant & t$NES < 0))
})))
NSIG[, contrast := factor(contrast, levels = ALL7)]
setorder(NSIG, quant, collection, contrast)
fwrite(NSIG, file.path(out, "n_significant_sets.tsv"), sep = "\t")
print(NSIG)

saveRDS(list(ALL = ALL, PUB = PUB, CINFO = CINFO, RANKS = RANKS, rank_summary = rank_summary,
             CONC = CONC, CIRC = CIRC, NSIG = NSIG, SET_META = SET_META, NKX_SETS = NKX_SETS, TREND = TREND,
             params = list(PADJ = PADJ, GS_MIN = GS_MIN, GS_MAX = GS_MAX, eps = 0, seed = 42,
                           msigdbr = as.character(packageVersion("msigdbr")), msigdb = MSIGDB_VERSION,
                           clusterProfiler = as.character(packageVersion("clusterProfiler")),
                           fgsea = as.character(packageVersion("fgsea")),
                           org.Hs.eg.db = as.character(packageVersion("org.Hs.eg.db")))),
        file.path(out, "gsea_results_all.rds"))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_63_gsea_salmon.txt"))
log_msg("compute stage finished")

# ----------------------------------------------------------------------------
# Figures and claim table
# ----------------------------------------------------------------------------
log_msg("all done")
