#!/usr/bin/env Rscript
# =============================================================================
# 15_tcga_geo_rwpe1.R -- public NKX3.1 perturbation data in a human prostate epithelial line (RWPE1)
#
# The time-boxed GEO search (15_tcga_geo.py) found two series with a ready processed matrix
# (DESeq variance-stabilised values, gene symbols x samples; GRCh37 / TopHat-Cufflinks-DESeq):
#   GSE81439  RWPE1 control vector (CAD1-6) vs NKX3.1 WT (CAD13-18)            (Dutta et al., Science 2016)
#   GSE115335 RWPE1 control vector + paraquat (CAD7-12) vs NKX3.1 WT + paraquat (CAD25-30)
#             (Papachristodoulou et al., Cancer Discov 2021)
# Mutant arms (T164A, R52C) are not used. limma-trend on the VST values, NKX3.1 WT minus vector.
# Run after 15_tcga_geo.py (downloads) and 15_tcga.R (TCGA tables used for concordance):
#   /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_geo_rwpe1.R
# =============================================================================
suppressPackageStartupMessages({
  library(limma); library(data.table); library(msigdbr); library(fgsea); library(org.Hs.eg.db)
})
SEED <- 20261021; set.seed(SEED)
proj    <- "/mnt/ssd1/Projects/NKX3_wd"
out     <- file.path(proj, "results/revision_2026/tcga")
tab_dir <- file.path(out, "tables")
work    <- "/mnt/ssd1/Projects/NKX3_wd/work/tcga"
dl      <- file.path(work, "geo_download")

CLAIMS <- list(
  activated = c("BCL11B", "TPGS1", "SSC4D", "USP3", "MDM4", "SRMS", "IGLL5", "SLX1A", "ZNF331", "LRIG1", "EHD1",
                "DIPK1B", "CAPN12", "CHMP4A", "ZNF548", "TRIM66", "PLEKHA3", "CCDC66", "FRMD4A"),
  repressed = c("NTNG1", "CSH2", "PIP4P2", "C3orf20", "SLX1B", "AC131160.1"),
  module_3D = c("EHD1", "SLX1A", "BCL11B", "SSC4D", "ZNF331", "NPHP4", "KIAA0895", "AGAP6", "TRPM4", "D2HGDH",
                "ZC3H12A", "KLHL3", "ULK4", "TPGS1", "FANCA", "TRPV3", "DIPK1B", "MCTS1", "IGLL5", "SRMS", "ARHGEF19"),
  module_3E = c("DIPK1B", "SSC4D", "FRMD4A", "IGLL5", "CNIH4", "LRIG1", "SRMS", "CCN1", "BCL11B", "TPGS1",
                "PPIP5K1", "HSPA6"))
claimed <- unique(unlist(CLAIMS))
HALLMARKS <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE",
               "HALLMARK_MYC_TARGETS_V1", "HALLMARK_APICAL_JUNCTION")

# ---- curated GEO table (manual curation of the 75 search hits; see geo_nkx3-1_search_all_hits.tsv) ----
geo_cur <- data.frame(
  accession = c("GSE81439", "GSE115335", "GSE47030", "GSE28596", "GSE28264 / GSE40269", "GSE103509 / GSE103536",
                "GSE81424 / GSE81432 / GSE115337 / GSE15943 / GSE92999", "(none) PMID 19462257"),
  model = c("RWPE1 (HPV18-immortalised normal prostate epithelium; AR-low)", "RWPE1, all arms paraquat-treated",
            "PrEC LH (immortalised prostate epithelial cells)", "LNCaP", "LNCaP", "human/mouse fibroblast reprogramming",
            "mouse prostate (Nkx3.1-/-, Nkx3.1;Pten, tissue recombinants)", "PC-3 stable pcDNA3.1-NKX3.1 transfectants"),
  design = c("control vector (n=6) vs NKX3.1 WT (n=6) vs NKX3.1 T164A (n=4); stable, inducible pTRIPZ",
             "vector (n=6) vs NKX3.1 WT (n=6) vs R52C (n=6) vs T164A (n=4), + paraquat",
             "GFP vs NKX3.1, 7 h and 10 h, n=2 each", "DHT vs ethanol time course (no NKX3.1 perturbation)",
             "NKX3-1 / AR / FOXA1 ChIP-seq (binding, not expression)", "NKX3-1 in iPSC reprogramming (not prostate)",
             "in vivo mouse expression profiles", "cDNA microarray (22k genes), 1,953 genes >2-fold"),
  platform = c("RNA-seq, Illumina HiSeq 2000 (GPL11154)", "RNA-seq, Illumina HiSeq 2000 (GPL11154)",
               "expression array (GPL6171)", "expression array (GPL6255)", "ChIP-seq (GPL9115)", "RNA-seq",
               "arrays / RNA-seq", "not deposited (no GEO link from PubMed)"),
  processed_matrix = c("yes: DESeq VST values (txt)", "yes: DESeq VST values (xlsx)", "non-normalised array data only",
                       "series matrix", "peaks (BED)", "-", "-", "no"),
  reference = c("Dutta et al., Science 2016 (PMID 27339988)", "Papachristodoulou et al., Cancer Discov 2021 (PMID 33893149)",
                "GEO record (no PMID listed)", "Tan et al., Mol Cell Biol 2012 (PMID 22083957)",
                "Tan et al., Mol Cell Biol 2012; GEO", "Mai et al., Nat Cell Biol 2018", "various",
                "Zhang et al., Mol Biol Rep 2010;37:1505 (PMID 19462257)"),
  analysed_here = c("yes (WT vs vector)", "yes (WT vs vector, under paraquat)", "no (n=2, 7-10 h, raw array data)",
                    "no (not an NKX3.1 perturbation)", "no (not expression)", "no (not prostate)", "no (mouse in vivo)",
                    "no (data not public)"),
  stringsAsFactors = FALSE)
fwrite(geo_cur, file.path(tab_dir, "geo_nkx3-1_perturbation_datasets_curated.tsv"), sep = "\t", quote = FALSE)

# ---- helpers ----
read_vst <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE)
  m <- as.matrix(dt[, -1]); storage.mode(m) <- "double"; rownames(m) <- dt[[1]]
  m[!is.na(rownames(m)) & rownames(m) != "", , drop = FALSE]
}
update_symbols <- function(m) {      # GRCh37-era symbols -> current HGNC symbols; duplicates: keep highest mean
  cur <- suppressWarnings(alias2SymbolTable(rownames(m), species = "Hs"))
  cur[is.na(cur)] <- rownames(m)[is.na(cur)]
  o <- order(-rowMeans(m)); m <- m[o, ]; cur <- cur[o]
  keep <- !duplicated(cur); m <- m[keep, ]; rownames(m) <- cur[keep]; m
}
hm <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))
HM_SYM <- lapply(split(hm$gene_symbol, hm$gs_name), unique)

run_one <- function(acc, file, ctrl, trt, label, min_mean = 4) {
  m <- update_symbols(read_vst(file))
  m <- m[, c(ctrl, trt)]
  grp <- factor(rep(c("vector", "NKX3.1"), c(length(ctrl), length(trt))), levels = c("vector", "NKX3.1"))
  keep <- rowMeans(m) >= min_mean
  # claimed genes are always reported, even below the filter
  mk <- m[keep | rownames(m) %in% c(claimed, "NKX3-1"), ]
  design <- model.matrix(~ grp)
  fit <- eBayes(lmFit(mk, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpNKX3.1", number = Inf, sort.by = "none")
  tt$symbol <- rownames(tt); tt$passes_filter <- keep[rownames(tt)]
  # BH within the filtered genome
  tt$adj.P.Val <- NA_real_; tt$adj.P.Val[tt$passes_filter] <- p.adjust(tt$P.Value[tt$passes_filter], "BH")
  ft <- sort(tt$t[tt$passes_filter]); tt$pctl_t <- 100 * findInterval(tt$t, ft) / length(ft)   # vs filtered genes
  # claimed genes
  cl <- data.frame(symbol = claimed, stringsAsFactors = FALSE)
  for (k in names(CLAIMS)) cl[[paste0("in_", k)]] <- cl$symbol %in% CLAIMS[[k]]
  cl$expected_sign <- ifelse(cl$in_repressed, -1L, 1L)
  cl <- cbind(cl, tt[match(cl$symbol, tt$symbol), c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val", "pctl_t", "passes_filter")])
  cl$present <- !is.na(cl$logFC)
  cl$agree <- sign(cl$logFC) == cl$expected_sign
  cl$dataset <- acc; cl$comparison <- label
  # sign tests (genes present and passing the filter)
  base_up <- mean(tt$t[tt$passes_filter] > 0)
  st <- rbindlist(lapply(c("activated", "repressed", "module_3D", "module_3E", "all"), function(s) {
    sel <- (if (s == "all") rep(TRUE, nrow(cl)) else cl[[paste0("in_", s)]]) & cl$present & cl$passes_filter %in% TRUE
    if (!any(sel)) return(NULL)
    a <- sum(cl$agree[sel]); n <- sum(sel)
    p0 <- mean(ifelse(cl$expected_sign[sel] > 0, base_up, 1 - base_up))
    data.frame(dataset = acc, gene_set = s, n_genes = n, n_agree = a,
               p_binom_vs_0.5 = binom.test(a, n, 0.5, alternative = "greater")$p.value,
               base_rate = p0, p_binom_vs_base = binom.test(a, n, p0, alternative = "greater")$p.value,
               n_agree_padj05 = sum(cl$agree[sel] & cl$adj.P.Val[sel] < 0.05, na.rm = TRUE),
               n_opposite_padj05 = sum(!cl$agree[sel] & cl$adj.P.Val[sel] < 0.05, na.rm = TRUE))
  }))
  # competitive tests: camera on claimed sets and on the five Hallmarks; fgsea on moderated t for all Hallmarks
  idx <- ids2indices(c(lapply(CLAIMS, intersect, rownames(fit)), HM_SYM[HALLMARKS]), rownames(fit))
  cam <- camera(mk, idx, design, contrast = "grpNKX3.1", trend.var = TRUE)
  cam$set <- rownames(cam); cam$dataset <- acc
  tf <- tt[tt$passes_filter & tt$symbol != "NKX3-1", ]
  # expression-dependence check (as for DU145) and a detrended ranking: robust loess of t on AveExpr
  lo <- loess(t ~ AveExpr, data = tf, span = 0.3, family = "symmetric", degree = 1)
  tf$t_detrended <- tf$t - predict(lo, tf)
  tf$dec <- cut(tf$AveExpr, quantile(tf$AveExpr, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  trend <- data.frame(dataset = acc, AveExpr_decile = 1:10, median_AveExpr = tapply(tf$AveExpr, tf$dec, median),
                      median_t = tapply(tf$t, tf$dec, median), frac_negative = tapply(tf$t < 0, tf$dec, mean))
  fg <- NULL
  for (v in c("t", "t_detrended")) {
    set.seed(SEED)
    f <- suppressWarnings(fgsea(HM_SYM, setNames(tf[[v]], tf$symbol), minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000))
    fg <- rbind(fg, data.frame(dataset = acc, ranking = paste0("RWPE1 ", label, " (", acc, ")",
                                                               if (v == "t_detrended") ", expression trend removed" else ""),
                               pathway = f$pathway, size = f$size, NES = f$NES, pval = f$pval, padj = f$padj))
  }
  list(tt = tt, tf = tf, cl = cl, st = st, cam = cam, fg = fg, trend = trend, n_tested = sum(keep), n_total = nrow(m),
       nkx = tt[tt$symbol == "NKX3-1", c("AveExpr", "logFC", "t", "P.Value")])
}

A <- run_one("GSE81439", file.path(dl, "GSE81439_RNAseq_human_RWPE1.txt.gz"),
             paste0("CAD", 1:6), paste0("CAD", 13:18), "NKX3.1 WT vs vector")
B <- run_one("GSE115335", file.path(dl, "GSE115335_RNAseq_human_RWPE1.tsv.gz"),
             paste0("CAD", 7:12), paste0("CAD", 25:30), "NKX3.1 WT vs vector, paraquat")

CL <- rbind(A$cl, B$cl); ST <- rbind(A$st, B$st); CAM <- rbind(A$cam, B$cam); FG <- rbind(A$fg, B$fg)
fwrite(CL, file.path(tab_dir, "rwpe1_claimed_genes_limma.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(ST, file.path(tab_dir, "rwpe1_direction_sign_tests.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(CAM, file.path(tab_dir, "rwpe1_camera_claimed_sets_and_hallmarks.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(FG, file.path(tab_dir, "rwpe1_hallmark_fgsea.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(rbind(A$trend, B$trend), file.path(tab_dir, "rwpe1_limma_t_by_expression_decile.tsv"), sep = "\t", quote = FALSE)
for (x in list(A, B)) {
  f <- file.path(tab_dir, sprintf("rwpe1_%s_limma_all_genes.tsv.gz", x$cl$dataset[1]))
  fwrite(x$tt[order(x$tt$P.Value), ], f, sep = "\t", quote = FALSE, na = "NA")
}

# ---- global concordance with DU145 (corrected counts) and TCGA (NKX3-1 correlations) ----
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
du_tab <- function(cn) {        # DU145 Wald z, raw and with the expression trend removed (as in 15_tcga.R)
  r <- core$fc$res[[cn]]; r <- r[!is.na(r$stat) & r$gene_id != "ENSG00000167034", ]
  r$symbol <- core$ann$symbol[match(r$gene_id, core$ann$gene_id)]; r$lbm <- log10(r$baseMean)
  lo <- loess(stat ~ lbm, data = r, span = 0.3, family = "symmetric", degree = 1)
  r$stat_detrended <- r$stat - predict(lo, r); r
}
du <- du_tab("NKX_in_ctrl"); du2 <- du_tab("NKX_in_CM")
tc <- fread(file.path(tab_dir, "tcga_genomewide_NKX3-1_spearman_primary_tumours.tsv.gz"), data.table = FALSE)
tc <- tc[tc$expressed, ]
conc <- rbindlist(lapply(list(A, B), function(x) {
  t1 <- x$tf
  f <- function(a, b) { s <- intersect(a$symbol[!duplicated(a$symbol)], b$symbol[!duplicated(b$symbol)])
    c(n = length(s), rho = cor(a$v[match(s, a$symbol)], b$v[match(s, b$symbol)], method = "spearman", use = "complete.obs")) }
  R <- rbind(
    f(data.frame(symbol = t1$symbol, v = t1$t), data.frame(symbol = du$symbol, v = du$stat)),
    f(data.frame(symbol = t1$symbol, v = t1$t), data.frame(symbol = du2$symbol, v = du2$stat)),
    f(data.frame(symbol = t1$symbol, v = t1$t_detrended), data.frame(symbol = du$symbol, v = du$stat_detrended)),
    f(data.frame(symbol = t1$symbol, v = t1$t_detrended), data.frame(symbol = du2$symbol, v = du2$stat_detrended)),
    f(data.frame(symbol = t1$symbol, v = t1$t), data.frame(symbol = tc$symbol, v = tc$rho_raw)),
    f(data.frame(symbol = t1$symbol, v = t1$t), data.frame(symbol = tc$symbol, v = tc$rho_AR_PROL)),
    f(data.frame(symbol = t1$symbol, v = t1$t_detrended), data.frame(symbol = tc$symbol, v = tc$rho_AR_PROL)))
  data.frame(dataset = x$cl$dataset[1],
             comparison = c("DU145 NKX3.1 vs control (basal), Wald z", "DU145 NKX3.1 vs control (inflammation), Wald z",
                            "DU145 basal, both statistics expression-detrended", "DU145 inflammation, both statistics expression-detrended",
                            "TCGA-PRAD Spearman rho with NKX3-1", "TCGA-PRAD partial rho | AR + proliferation",
                            "TCGA-PRAD partial rho | AR + proliferation vs detrended RWPE1 t"),
             n_genes = R[, "n"], spearman_of_gene_statistics = R[, "rho"])
}))
fwrite(conc, file.path(tab_dir, "rwpe1_global_concordance.tsv"), sep = "\t", quote = FALSE)

summ <- list(GSE81439 = list(n_genes_total = A$n_total, n_tested = A$n_tested, NKX3_1 = as.list(A$nkx)),
             GSE115335 = list(n_genes_total = B$n_total, n_tested = B$n_tested, NKX3_1 = as.list(B$nkx)),
             filter = "mean VST >= 4 across the 12 samples (claimed genes always reported)",
             model = "limma-trend, robust eBayes, ~ group (NKX3.1 WT minus vector)", seed = SEED)
writeLines(jsonlite::toJSON(summ, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out, "rwpe1_summary.json"))
saveRDS(list(A = A[c("cl", "st", "cam", "fg", "nkx", "trend", "n_tested")], B = B[c("cl", "st", "cam", "fg", "nkx", "trend", "n_tested")], conc = conc,
             geo = geo_cur), file.path(out, "rwpe1_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_15_tcga_geo_rwpe1.txt"))
print(summ); print(ST); print(CAM[CAM$set %in% c(names(CLAIMS), HALLMARKS), ]); print(conc); print(rbind(A$trend, B$trend))
print(FG[FG$pathway %in% HALLMARKS, ])
