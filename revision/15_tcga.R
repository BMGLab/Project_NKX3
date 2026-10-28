#!/usr/bin/env Rscript
# =============================================================================
# 15_tcga.R -- Independent-dataset check in TCGA-PRAD (Reviewer 1, major point 5)
#
# Questions
#   1. Do the genes the manuscript reports as NKX3.1-"activated" / "repressed" (Abstract, Results,
#      Fig 1) and the Fig 3D / 3E "NKX3.1 co-expression modules" co-vary with NKX3-1 across prostate
#      tumours more than a random expressed gene does?
#   2. Does that hold after adjusting for AR activity and proliferation (and, as a sensitivity,
#      stromal content)?
#   3. Are the Fig 3D / 3E modules co-expressed (coherent) in tumours, compared with random
#      expression-matched gene sets?
#   4. Do the gene-set-level effects seen in the corrected DU145 data (NKX3.1 lowers E2F / G2M /
#      mitotic spindle / MYC / apical junction) have a correlational counterpart in prostate tumours?
#      The DU145 rankings are also re-run after removing a global expression-dependent trend in the
#      Wald statistics (found during this analysis; see README notes in the task summary).
#
# Data   TCGA-PRAD, UCSC Xena GDC hub, "STAR - TPM" (values log2(TPM+1)), GENCODE v36.
#        Sensitivity: legacy Xena TCGA hub HiSeqV2 (RSEM normalised counts, log2(x+1)).
#        If the GDC file cannot be obtained, the legacy file becomes the primary source.
# Caveat TCGA-PRAD = AR-positive primary tumours with endogenous NKX3-1 (frequent 8p21 copy-number
#        loss). DU145 = AR-negative cell line with ectopic NKX3.1. Correlation is not regulation.
#
# Order  15_tcga_geo.py -> 15_tcga.R -> 15_tcga_geo_rwpe1.R -> 15_tcga_figures.R  (15_tcga_run_all.sh)
# Run    /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(msigdbr); library(fgsea); library(ppcor)
  library(jsonlite); library(AnnotationDbi); library(org.Hs.eg.db)
})

SEED <- 20261021
set.seed(SEED)
B_PERM <- 10000                       # permutations for matched-null tests

proj    <- "/mnt/ssd1/Projects/NKX3_wd"
out     <- file.path(proj, "results/revision_2026/tcga")
tab_dir <- file.path(out, "tables")
work    <- "/mnt/ssd1/Projects/NKX3_wd/work/tcga"
dl_dir  <- file.path(work, "download")
for (d in c(out, tab_dir, dl_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...))

# ----------------------------------------------------------------------------
# 0. Downloads + provenance
# ----------------------------------------------------------------------------
SRC <- list(
  gdc_tpm  = list(url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.star_tpm.tsv.gz",
                  file = file.path(dl_dir, "TCGA-PRAD.star_tpm.tsv.gz")),
  gdc_pmap = list(url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
                  file = file.path(dl_dir, "gencode.v36.annotation.gtf.gene.probemap")),
  legacy   = list(url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PRAD.sampleMap%2FHiSeqV2.gz",
                  file = file.path(dl_dir, "TCGA.PRAD.sampleMap_HiSeqV2.gz"))
)
fetch <- function(src) {
  if (!file.exists(src$file)) {
    msg("downloading %s", src$url)
    ok <- tryCatch({ utils::download.file(src$url, src$file, mode = "wb", quiet = TRUE, method = "libcurl"); TRUE },
                   error = function(e) { msg("download failed: %s", conditionMessage(e)); FALSE })
    if (!ok) { unlink(src$file); return(FALSE) }
  }
  TRUE
}
head_last_modified <- function(url) {
  tryCatch({
    h <- base::curlGetHeaders(url)
    lm <- grep("^last-modified:", h, ignore.case = TRUE, value = TRUE)
    if (length(lm)) trimws(sub("^[^:]+:", "", lm[1])) else NA_character_
  }, error = function(e) NA_character_)
}
prov_entry <- function(src) {
  fi <- file.info(src$file)
  list(url = src$url, local_file = src$file, bytes = unname(fi$size), md5 = unname(tools::md5sum(src$file)),
       accessed_utc = format(fi$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),   # file mtime = download time
       server_last_modified = head_last_modified(src$url))
}
have_gdc    <- fetch(SRC$gdc_tpm) && fetch(SRC$gdc_pmap)
have_legacy <- fetch(SRC$legacy)
if (!have_gdc && !have_legacy) stop("Neither TCGA-PRAD source could be downloaded.")
primary_source <- if (have_gdc) "GDC_STAR_TPM" else "LEGACY_RSEM"
msg("primary source: %s", primary_source)

read_xena <- function(f) {
  dt <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE)
  m <- as.matrix(dt[, -1, drop = FALSE]); rownames(m) <- dt[[1]]; storage.mode(m) <- "double"; m
}
# barcode positions 14-15: 01 primary tumour, 06 metastasis, 11 solid tissue normal.
# One sample per patient and type (lexicographically first aliquot, i.e. vial A before B).
pick_samples <- function(ids, code) {
  d <- data.frame(id = ids, type = substr(ids, 14, 15), patient = substr(ids, 1, 12), stringsAsFactors = FALSE)
  d <- d[d$type == code, ]; d <- d[order(d$patient, d$id), ]
  list(keep = d$id[!duplicated(d$patient)], dropped = d$id[duplicated(d$patient)])
}

# ----------------------------------------------------------------------------
# 1. Gene lists (manuscript claims) and annotation
# ----------------------------------------------------------------------------
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
ann99 <- core$ann

CLAIMS <- list(
  activated = c("BCL11B", "TPGS1", "SSC4D", "USP3", "MDM4", "SRMS", "IGLL5", "SLX1A", "ZNF331", "LRIG1", "EHD1",
                "DIPK1B", "CAPN12", "CHMP4A", "ZNF548", "TRIM66", "PLEKHA3", "CCDC66", "FRMD4A"),
  repressed = c("NTNG1", "CSH2", "PIP4P2", "C3orf20", "SLX1B", "AC131160.1"),
  module_3D = c("EHD1", "SLX1A", "BCL11B", "SSC4D", "ZNF331", "NPHP4", "KIAA0895", "AGAP6", "TRPM4", "D2HGDH",
                "ZC3H12A", "KLHL3", "ULK4", "TPGS1", "FANCA", "TRPV3", "DIPK1B", "MCTS1", "IGLL5", "SRMS", "ARHGEF19"),
  module_3E = c("DIPK1B", "SSC4D", "FRMD4A", "IGLL5", "CNIH4", "LRIG1", "SRMS", "CCN1", "BCL11B", "TPGS1",
                "PPIP5K1", "HSPA6"))
ABSTRACT_GENES <- c("BCL11B", "TPGS1", "SSC4D", "USP3", "MDM4", "NTNG1", "CSH2", "PIP4P2")
KLRC3_NOTE <- "KLRC3"   # bottom row of the Fig 3E box (label truncated in the figure); sensitivity only
AR_GENES     <- c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "STEAP2", "PMEPA1")   # NKX3-1 deliberately excluded
PROL_GENES   <- c("MKI67", "TOP2A", "CCNB1", "CDK1")
STROMA_GENES <- c("MYH11", "CNN1", "DES", "ACTG2", "DCN", "LUM", "COL1A2", "COL3A1")  # smooth muscle + fibroblast
HALLMARKS  <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE",
                "HALLMARK_MYC_TARGETS_V1", "HALLMARK_APICAL_JUNCTION")
HM_LABEL   <- c(HALLMARK_E2F_TARGETS = "E2F targets", HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
                HALLMARK_MITOTIC_SPINDLE = "Mitotic spindle", HALLMARK_MYC_TARGETS_V1 = "MYC targets V1",
                HALLMARK_APICAL_JUNCTION = "Apical junction")
VARIANTS <- c("raw", "AR", "PROL", "AR_PROL", "AR_PROL_STROMA")      # partial-correlation covariate sets

claimed <- unique(unlist(CLAIMS))
sym2id <- function(s) {
  i <- match(s, ann99$symbol)
  if (anyNA(i)) stop("symbol(s) not in Ensembl 99 annotation: ", paste(s[is.na(i)], collapse = ", "))
  setNames(ann99$gene_id[i], s)
}
ID_NKX  <- sym2id("NKX3-1"); ID_CLM <- sym2id(claimed); ID_KLRC3 <- sym2id(KLRC3_NOTE)
ID_AR   <- sym2id(AR_GENES); ID_PROL <- sym2id(PROL_GENES); ID_STR <- sym2id(STROMA_GENES)

ginfo <- data.frame(symbol = claimed, gene_id = unname(ID_CLM), stringsAsFactors = FALSE)
for (k in names(CLAIMS)) ginfo[[paste0("in_", k)]] <- ginfo$symbol %in% CLAIMS[[k]]
ginfo$in_abstract <- ginfo$symbol %in% ABSTRACT_GENES
ginfo$expected_sign <- ifelse(ginfo$in_repressed, -1L, 1L)   # module genes cluster with NKX3-1 (up) -> +1
stopifnot(!any(ginfo$in_repressed & (ginfo$in_activated | ginfo$in_module_3D | ginfo$in_module_3E)))
ginfo$category <- ifelse(ginfo$in_activated, "Activated (text/Fig 1)",
                   ifelse(ginfo$in_repressed, "Repressed (text/Fig 1)", "Module only (Fig 3D/3E)"))

hm_all <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))
msig_db_version <- unique(hm_all$db_version)
HM_SETS <- lapply(split(hm_all$ensembl_gene, hm_all$gs_name), function(v) unique(v[!is.na(v) & v != ""]))
HM_SYM  <- lapply(split(hm_all$gene_symbol, hm_all$gs_name), unique)

# ----------------------------------------------------------------------------
# 2. Load TCGA-PRAD (GDC STAR TPM; fallback legacy RSEM). Parsed matrix cached in the temp dir.
# ----------------------------------------------------------------------------
cache <- file.path(work, "tcga_parsed.rds")
if (file.exists(cache)) {
  P <- readRDS(cache)
} else {
  P <- list()
  if (primary_source == "GDC_STAR_TPM") {
    m <- read_xena(SRC$gdc_tpm$file)
    pmap <- fread(SRC$gdc_pmap$file, data.table = FALSE)
    n_raw <- nrow(m)
    m <- m[!grepl("_PAR_Y$", rownames(m)), , drop = FALSE]
    full_id <- rownames(m)
    rownames(m) <- sub("\\.\\d+$", "", rownames(m))
    stopifnot(!anyDuplicated(rownames(m)))
    P$sym <- setNames(pmap$gene[match(full_id, pmap$id)], rownames(m))
    P$n_genes_raw <- n_raw; P$n_par_y_removed <- n_raw - nrow(m)
    P$unit <- "log2(TPM+1)"
  } else {
    m <- read_xena(SRC$legacy$file)
    eg <- suppressMessages(mapIds(org.Hs.eg.db, rownames(m), "ENSEMBL", "SYMBOL", multiVals = "first"))
    keep <- !is.na(eg) & !duplicated(eg); P$sym <- setNames(rownames(m)[keep], eg[keep])
    m <- m[keep, ]; rownames(m) <- eg[keep]
    P$n_genes_raw <- length(keep); P$n_par_y_removed <- 0L
    P$unit <- "log2(RSEM normalised count+1)"
  }
  ids <- colnames(m)
  P$n_samples_total <- length(ids)
  P$type_table <- as.list(table(substr(ids, 14, 15)))
  tum <- pick_samples(ids, "01"); nor <- pick_samples(ids, "11")
  P$tum_ids <- tum$keep; P$tum_dropped <- tum$dropped
  P$nor_ids <- nor$keep; P$nor_dropped <- nor$dropped
  P$X_tum <- m[, P$tum_ids]; P$X_nor <- m[, P$nor_ids]
  if (P$unit == "log2(TPM+1)") P$tpm_colsum_range <- range(colSums(2^m - 1))
  rm(m); gc()
  saveRDS(P, cache)
}
msg("tumours n=%d (dropped dup %d), normals n=%d", length(P$tum_ids), length(P$tum_dropped), length(P$nor_ids))
for (id in c(ID_NKX, ID_CLM, ID_AR, ID_PROL, ID_STR)) if (!id %in% rownames(P$X_tum)) stop("missing in TCGA: ", id)

# ----------------------------------------------------------------------------
# 3. Helpers
# ----------------------------------------------------------------------------
row_ranks <- function(X) { R <- rowRanks(X, ties.method = "average"); dimnames(R) <- dimnames(X); R }
unit_rows <- function(R) { Rc <- R - rowMeans(R); s <- sqrt(rowSums(Rc^2)); s[s == 0] <- NA; Rc / s }
resid_rows <- function(R, C) {       # residualise rows of R (genes x samples) on covariates C (samples x k)
  if (is.null(C)) return(R)
  q <- qr(cbind(1, C)); t(qr.resid(q, t(R)))
}
cor_vs <- function(U, u) as.vector(U %*% u)        # U: unit rows, u: unit vector -> Pearson r
p_from_r <- function(r, df) { t <- r * sqrt(df / pmax(1 - r^2, 1e-300)); 2 * pt(-abs(t), df) }
ci_bw <- function(r, n, k = 0) {    # Bonett & Wright (2000) Fisher-z interval for Spearman / partial Spearman
  se <- sqrt((1 + r^2 / 2) / (n - 3 - k)); z <- atanh(r)
  cbind(lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se))
}
pctl <- function(x, null) { s <- sort(null[!is.na(null)]); n <- length(s)   # mid-rank percentile
  vapply(x, function(v) if (is.na(v)) NA_real_ else 100 * (findInterval(v, s, left.open = TRUE) +
         0.5 * (findInterval(v, s) - findInterval(v, s, left.open = TRUE))) / n, numeric(1)) }
zrows <- function(X) { mu <- rowMeans(X); s <- rowSds(X); (X - mu) / s }
score_of <- function(Z, ids) colMeans(Z[intersect(ids, rownames(Z)), , drop = FALSE])
draw <- function(pool, n) pool[sample.int(length(pool), n)]
deciles_over <- function(x, ref_idx) {   # decile of x with breaks from x[ref_idx]
  br <- quantile(x[ref_idx], 0:10 / 10, na.rm = TRUE); br[1] <- -Inf; br[11] <- Inf
  cut(x, br, labels = FALSE)
}

# Analyse one cohort (tumour or normal): genome-wide Spearman / partial Spearman with NKX3-1
analyse_cohort <- function(X, label) {
  n <- ncol(X)
  is_tpm <- P$unit == "log2(TPM+1)"
  tpm_med  <- if (is_tpm) rowMedians(2^X - 1) else rowMedians(X)
  tpm_mean <- if (is_tpm) rowMeans(2^X - 1) else rowMeans(X)
  expressed <- rownames(X)[tpm_med >= if (is_tpm) 1 else log2(11)]
  nonzero <- rownames(X)[rowVars(X) > 0]
  detected <- intersect(rownames(X)[tpm_mean >= if (is_tpm) 0.1 else 1], nonzero)
  # correlations for all detected genes (mean TPM >= 0.1) so the "all genes" module sensitivity has matched
  # low-abundance pools; genome-wide null + genome-wide BH family = expressed genes (+ claimed genes)
  universe <- union(union(expressed, detected), intersect(c(ID_CLM, ID_KLRC3), nonzero))
  family <- union(expressed, intersect(c(ID_CLM, ID_KLRC3), nonzero))
  Z <- zrows(X[nonzero, ])
  ar <- score_of(Z, ID_AR); pr <- score_of(Z, ID_PROL); st <- score_of(Z, ID_STR)
  r_nkx <- rank(X[ID_NKX, ]); r_ar <- rank(ar); r_pr <- rank(pr); r_st <- rank(st)
  covs <- list(raw = NULL, AR = cbind(r_ar), PROL = cbind(r_pr), AR_PROL = cbind(r_ar, r_pr),
               AR_PROL_STROMA = cbind(r_ar, r_pr, r_st))
  R <- row_ranks(X[universe, ])
  res <- data.frame(gene_id = universe, symbol = unname(P$sym[universe]),
                    median_TPM = tpm_med[match(universe, rownames(X))], mean_TPM = tpm_mean[match(universe, rownames(X))],
                    mean_log2 = rowMeans(X[universe, ]), expressed = universe %in% expressed,
                    in_family = universe %in% family, stringsAsFactors = FALSE)
  U_by <- list()
  for (cv in names(covs)) {
    C <- covs[[cv]]; k <- if (is.null(C)) 0 else ncol(C)
    U <- unit_rows(resid_rows(R, C))
    y <- resid_rows(matrix(r_nkx, nrow = 1), C)[1, ]; y <- (y - mean(y)) / sqrt(sum((y - mean(y))^2))
    rho <- cor_vs(U, y)
    res[[paste0("rho_", cv)]] <- rho
    res[[paste0("p_", cv)]] <- p_from_r(rho, n - 2 - k)
    pa <- rep(NA_real_, length(rho)); pa[res$in_family] <- p.adjust(res[[paste0("p_", cv)]][res$in_family], "BH")
    res[[paste0("padj_genomewide_", cv)]] <- pa
    res[[paste0("pctl_", cv)]] <- pctl(rho, rho[res$expressed & res$gene_id != ID_NKX])
    U_by[[cv]] <- U
  }
  list(label = label, n = n, Z = Z, res = res, expressed = expressed, universe = universe,
       ar = ar, pr = pr, st = st, r_ar = r_ar, r_pr = r_pr, r_st = r_st, r_nkx = r_nkx, U = U_by,
       rho_nkx_ar = cor(r_nkx, r_ar), rho_nkx_pr = cor(r_nkx, r_pr), rho_nkx_st = cor(r_nkx, r_st),
       rho_ar_pr = cor(r_ar, r_pr))
}

msg("genome-wide correlations: tumours")
TUM <- analyse_cohort(P$X_tum, "Primary tumour")
msg("genome-wide correlations: normals")
NOR <- analyse_cohort(P$X_nor, "Solid tissue normal")
msg("tumour: %d expressed genes; rho(NKX3-1, AR) = %.3f; rho(NKX3-1, prolif) = %.3f; rho(NKX3-1, stroma) = %.3f",
    length(TUM$expressed), TUM$rho_nkx_ar, TUM$rho_nkx_pr, TUM$rho_nkx_st)

# cross-check the vectorised partial Spearman against ppcor
chk <- sapply(c(ID_CLM["BCL11B"], ID_CLM["EHD1"]), function(g) {
  pc <- ppcor::pcor.test(P$X_tum[ID_NKX, ], P$X_tum[g, ], cbind(TUM$ar, TUM$pr), method = "spearman")
  c(ppcor = pc$estimate, ours = TUM$res$rho_AR_PROL[TUM$res$gene_id == g])
})
stopifnot(max(abs(chk[1, ] - chk[2, ])) < 1e-8)

# ----------------------------------------------------------------------------
# 4. Analyses 1-2: claimed genes vs the genome-wide null
# ----------------------------------------------------------------------------
claimed_table <- function(C, tag) {
  r <- C$res[match(ginfo$gene_id, C$res$gene_id), ]
  out <- cbind(ginfo, r[, setdiff(colnames(r), c("gene_id", "symbol"))])
  for (cv in VARIANTS) {
    out[[paste0("padj_claimed_", cv)]] <- p.adjust(out[[paste0("p_", cv)]], "BH")
    out[[paste0("agree_", cv)]] <- sign(out[[paste0("rho_", cv)]]) == out$expected_sign
  }
  ci <- ci_bw(out$rho_raw, C$n, 0); out$rho_raw_lo <- ci[, 1]; out$rho_raw_hi <- ci[, 2]
  ci <- ci_bw(out$rho_AR_PROL, C$n, 2); out$rho_AR_PROL_lo <- ci[, 1]; out$rho_AR_PROL_hi <- ci[, 2]
  out$cohort <- tag; out$n_samples <- C$n
  out
}
CT <- claimed_table(TUM, "primary_tumour")
CN <- claimed_table(NOR, "solid_normal")

null_summary <- function(C) {
  do.call(rbind, lapply(VARIANTS, function(cv) {
    sel <- C$res$expressed & C$res$gene_id != ID_NKX
    v <- C$res[[paste0("rho_", cv)]][sel]
    q <- quantile(v, c(.025, .05, .25, .5, .75, .95, .975), na.rm = TRUE)
    data.frame(cohort = C$label, adjustment = cv, n_genes = sum(!is.na(v)), mean = mean(v, na.rm = TRUE),
               frac_positive = mean(v > 0, na.rm = TRUE),
               frac_BH05_genomewide = mean(C$res[[paste0("padj_genomewide_", cv)]][sel] < 0.05, na.rm = TRUE),
               q025 = q[1], q05 = q[2], q25 = q[3], q50 = q[4], q75 = q[5], q95 = q[6], q975 = q[7])
  }))
}
NULLS <- rbind(null_summary(TUM), null_summary(NOR))

# matched permutation (expression decile) for the signed mean rho and the agreement count
perm_signed <- function(C, genes, signs, cv, B = B_PERM, seed = SEED) {
  set.seed(seed)
  rho <- C$res[[paste0("rho_", cv)]]
  ref <- which(C$res$expressed & C$res$gene_id != ID_NKX)
  dec <- deciles_over(C$res$mean_log2, ref)
  gi <- match(genes, C$res$gene_id)
  pool_idx <- which(C$res$expressed & !(C$res$gene_id %in% c(ID_NKX, ID_CLM)))
  pools <- split(pool_idx, dec[pool_idx]); md <- split(seq_along(gi), dec[gi])
  obs_mean <- mean(signs * rho[gi]); obs_agree <- sum(sign(rho[gi]) == signs)
  nm <- numeric(B); na <- integer(B)
  for (b in seq_len(B)) {
    pick <- integer(length(gi))
    for (d in names(md)) pick[md[[d]]] <- draw(pools[[d]], length(md[[d]]))
    v <- rho[pick]; nm[b] <- mean(signs * v); na[b] <- sum(sign(v) == signs)
  }
  c(obs_signed_mean = obs_mean, null_signed_mean = mean(nm), p_perm_signed_mean = (1 + sum(nm >= obs_mean)) / (1 + B),
    obs_agree = obs_agree, null_agree_mean = mean(na), p_perm_agree = (1 + sum(na >= obs_agree)) / (1 + B))
}

sign_tests <- function(C, CTab) {
  sets <- list(activated = CTab$in_activated, repressed = CTab$in_repressed, module_3D = CTab$in_module_3D,
               module_3E = CTab$in_module_3E, all_claimed = rep(TRUE, nrow(CTab)), abstract_genes = CTab$in_abstract)
  rows <- list()
  for (cv in c("raw", "AR_PROL", "AR_PROL_STROMA")) for (s in names(sets)) for (scope in c("expressed", "all")) {
    sel <- sets[[s]] & !is.na(CTab[[paste0("rho_", cv)]])
    if (scope == "expressed") sel <- sel & CTab$expressed
    if (!any(sel)) next
    rho <- CTab[[paste0("rho_", cv)]][sel]; sg <- CTab$expected_sign[sel]
    agree <- sum(sign(rho) == sg); n <- sum(sel)
    nullv <- C$res[[paste0("rho_", cv)]][C$res$expressed & C$res$gene_id != ID_NKX]
    p0 <- mean(ifelse(sg > 0, mean(nullv > 0, na.rm = TRUE), mean(nullv < 0, na.rm = TRUE)))  # genome-wide base rate
    padj <- CTab[[paste0("padj_claimed_", cv)]][sel]
    row <- data.frame(cohort = C$label, adjustment = cv, gene_set = s, scope = scope, n_genes = n,
                      n_agree = agree, frac_agree = agree / n,
                      p_binom_vs_0.5 = binom.test(agree, n, 0.5, alternative = "greater")$p.value,
                      genomewide_base_rate = p0,
                      p_binom_vs_base_rate = binom.test(agree, n, p0, alternative = "greater")$p.value,
                      n_agree_BH05 = sum(sign(rho) == sg & padj < 0.05), n_opposite_BH05 = sum(sign(rho) != sg & padj < 0.05),
                      median_pctl_signed = median(ifelse(sg > 0, CTab[[paste0("pctl_", cv)]][sel],
                                                         100 - CTab[[paste0("pctl_", cv)]][sel])),
                      stringsAsFactors = FALSE)
    if (scope == "expressed" && n >= 3) row <- cbind(row, as.data.frame(as.list(perm_signed(C, CTab$gene_id[sel], sg, cv))))
    rows[[length(rows) + 1]] <- row
  }
  rbindlist(rows, fill = TRUE)
}
msg("sign tests + matched permutations")
ST <- rbind(sign_tests(TUM, CT), sign_tests(NOR, CN), fill = TRUE)

# ----------------------------------------------------------------------------
# 5. Analysis 3: module coherence vs expression-matched random gene sets
# ----------------------------------------------------------------------------
EXCL_IDS <- c(ID_NKX, ID_CLM, ID_AR, ID_PROL, ID_STR, ID_KLRC3)
coherence <- function(C, genes, set_name, cv = "raw", B = B_PERM, seed = SEED, scope = "expressed") {
  set.seed(seed)
  U <- C$U[[cv]]; res <- C$res                       # U rows aligned with res rows
  ok <- if (scope == "expressed") res$expressed else !is.na(U[, 1])   # "all" = detected + claimed genes
  gi <- match(intersect(genes, res$gene_id[ok]), res$gene_id)
  if (length(gi) < 3) return(NULL)
  okidx <- which(ok)
  dec <- deciles_over(res$mean_log2, okidx)
  pool_idx <- setdiff(okidx, which(res$gene_id %in% EXCL_IDS))
  pools <- split(pool_idx, dec[pool_idx]); md <- split(seq_along(gi), dec[gi])
  mp <- function(ix) { M <- tcrossprod(U[ix, , drop = FALSE]); mean(M[upper.tri(M)]) }
  obs <- mp(gi)
  nullv <- numeric(B)
  for (b in seq_len(B)) {
    pick <- integer(length(gi))
    for (d in names(md)) pick[md[[d]]] <- draw(pools[[d]], length(md[[d]]))
    nullv[b] <- mp(pick)
  }
  list(row = data.frame(cohort = C$label, gene_set = set_name, adjustment = cv, scope = scope,
                        n_genes = length(gi), genes = paste(res$symbol[gi], collapse = ","),
                        obs_mean_pairwise_rho = obs, null_mean = mean(nullv), null_sd = sd(nullv),
                        null_q95 = unname(quantile(nullv, .95)), null_q99 = unname(quantile(nullv, .99)),
                        z = (obs - mean(nullv)) / sd(nullv), percentile = 100 * mean(nullv < obs),
                        p_perm = (1 + sum(nullv >= obs)) / (1 + B), stringsAsFactors = FALSE),
       null = nullv)
}
msg("module coherence")
MOD_SETS <- list(module_3D = unname(ID_CLM[CLAIMS$module_3D]), module_3E = unname(ID_CLM[CLAIMS$module_3E]),
                 module_3D_plus_NKX3_1 = c(unname(ID_CLM[CLAIMS$module_3D]), ID_NKX),
                 module_3E_plus_NKX3_1 = c(unname(ID_CLM[CLAIMS$module_3E]), ID_NKX),
                 module_3E_plus_KLRC3 = c(unname(ID_CLM[CLAIMS$module_3E]), ID_KLRC3),
                 posctrl_AR_score_genes = unname(ID_AR), posctrl_proliferation_genes = unname(ID_PROL),
                 posctrl_stromal_genes = unname(ID_STR))
COH <- list(); COH_NULL <- list()
for (C in list(TUM, NOR)) for (s in names(MOD_SETS)) for (cv in c("raw", "AR_PROL", "AR_PROL_STROMA")) for (scope in c("expressed", "all")) {
  if (cv != "raw" && grepl("^posctrl", s)) next          # covariates are built from these genes
  if (scope == "all" && !grepl("^module_3[DE]$", s)) next
  cc <- coherence(C, MOD_SETS[[s]], s, cv, scope = scope)
  if (is.null(cc)) next
  COH[[length(COH) + 1]] <- cc$row
  if (cv == "raw" && scope == "expressed") COH_NULL[[paste(C$label, s)]] <- data.frame(cohort = C$label, gene_set = s, null = cc$null)
}
COH <- rbindlist(COH); COH_NULL <- rbindlist(COH_NULL)

pair_tab <- function(C, ids, set_name) {
  gi <- match(intersect(ids, C$res$gene_id[C$res$expressed]), C$res$gene_id); U <- C$U$raw
  M <- tcrossprod(U[gi, , drop = FALSE]); s <- C$res$symbol[gi]
  ij <- which(upper.tri(M), arr.ind = TRUE)
  data.frame(cohort = C$label, module = set_name, gene_a = s[ij[, 1]], gene_b = s[ij[, 2]], spearman_rho = M[ij])
}
PAIRS <- rbind(pair_tab(TUM, MOD_SETS$module_3D, "module_3D"), pair_tab(TUM, MOD_SETS$module_3E, "module_3E"))

# ----------------------------------------------------------------------------
# 6. Analysis 4: Hallmark programme scores vs NKX3-1
# ----------------------------------------------------------------------------
cov_mat <- function(C, cv) switch(cv, raw = NULL, AR = cbind(C$r_ar), AR_PROL = cbind(C$r_ar, C$r_pr),
                                  AR_STROMA = cbind(C$r_ar, C$r_st))
pathway_cor <- function(C) {
  rows <- list()
  for (h in HALLMARKS) {
    mem <- setdiff(intersect(HM_SETS[[h]], C$expressed), ID_NKX)
    r_sc <- rank(score_of(C$Z, mem))
    for (cv in c("raw", "AR", "AR_PROL", "AR_STROMA")) {
      Cm <- cov_mat(C, cv); k <- if (is.null(Cm)) 0 else ncol(Cm)
      a <- resid_rows(matrix(C$r_nkx, 1), Cm)[1, ]; b <- resid_rows(matrix(r_sc, 1), Cm)[1, ]
      r <- cor(a, b); ci <- ci_bw(r, C$n, k)
      rows[[length(rows) + 1]] <- data.frame(cohort = C$label, n_samples = C$n, pathway = h, label = HM_LABEL[[h]],
        n_genes_in_score = length(mem), n_genes_in_set = length(HM_SETS[[h]]), adjustment = cv,
        rho = r, ci_lo = ci[1], ci_hi = ci[2], p = p_from_r(r, C$n - 2 - k),
        overlap_with_prolif_genes = sum(ID_PROL %in% mem), stringsAsFactors = FALSE)
    }
  }
  out <- rbindlist(rows); out[, padj_BH := p.adjust(p, "BH"), by = .(cohort, adjustment)]; out
}
PW <- rbind(pathway_cor(TUM), pathway_cor(NOR))

# competitive version: the programme's correlation vs scores of random expression-matched gene sets of the
# same size (the TCGA genome-wide distribution of NKX3-1 correlations is shifted, see null summary)
pathway_competitive <- function(C, B = B_PERM, seed = SEED) {
  set.seed(seed)
  ex <- setdiff(C$expressed, ID_NKX); Zx <- C$Z[ex, ]
  dec <- deciles_over(C$res$mean_log2[match(ex, C$res$gene_id)], seq_along(ex))
  pools_all <- split(seq_along(ex), dec)
  q_ar <- qr(cbind(1, C$r_ar)); q_as <- qr(cbind(1, C$r_ar, C$r_st))
  y_raw <- C$r_nkx - mean(C$r_nkx); y_ar <- qr.resid(q_ar, C$r_nkx); y_as <- qr.resid(q_as, C$r_nkx)
  f_all <- function(sc) { r <- rank(sc)
    c(raw = cor(y_raw, r), AR = cor(y_ar, qr.resid(q_ar, r)), AR_STROMA = cor(y_as, qr.resid(q_as, r))) }
  rows <- list()
  for (h in HALLMARKS) {
    mem <- match(intersect(HM_SETS[[h]], ex), ex)
    o_all <- f_all(colMeans(Zx[mem, ]))
    md <- split(seq_along(mem), dec[mem]); pools <- lapply(pools_all, setdiff, mem)
    NV <- matrix(NA_real_, B, 3, dimnames = list(NULL, names(o_all)))
    for (b in seq_len(B)) {
      pick <- integer(length(mem))
      for (d in names(md)) pick[md[[d]]] <- draw(pools[[d]], length(md[[d]]))
      NV[b, ] <- f_all(colMeans(Zx[pick, , drop = FALSE]))
    }
    for (cv in names(o_all)) {
      o <- o_all[[cv]]; nv <- NV[, cv]
      rows[[length(rows) + 1]] <- data.frame(cohort = C$label, pathway = h, label = HM_LABEL[[h]], adjustment = cv,
        n_genes = length(mem), obs_rho = o, null_mean = mean(nv), null_q025 = unname(quantile(nv, .025)),
        null_q975 = unname(quantile(nv, .975)), percentile = 100 * mean(nv < o),
        p_lower_than_random = (1 + sum(nv <= o)) / (1 + B), p_higher_than_random = (1 + sum(nv >= o)) / (1 + B),
        stringsAsFactors = FALSE)
    }
  }
  rbindlist(rows)
}
msg("programme scores: competitive nulls")
PWC <- rbind(pathway_competitive(TUM), pathway_competitive(NOR))
PW_CTX <- rbindlist(lapply(HALLMARKS, function(h) {
  mem <- setdiff(intersect(HM_SETS[[h]], TUM$expressed), ID_NKX); sc <- rank(score_of(TUM$Z, mem))
  data.frame(pathway = h, rho_with_AR_score = cor(sc, TUM$r_ar), rho_with_prolif_score = cor(sc, TUM$r_pr),
             rho_with_stromal_score = cor(sc, TUM$r_st))
}))

# ----------------------------------------------------------------------------
# 7. Preranked GSEA: TCGA NKX3-1 correlation rankings and DU145 (corrected counts) Wald statistics.
#    DU145 Wald z show a monotone dependence on expression (low-expressed genes up, high-expressed down,
#    in every contrast, including CM vs control), so each DU145 ranking is also run after removing that
#    trend (robust loess of z on log10 baseMean, span 0.3).
# ----------------------------------------------------------------------------
run_fgsea <- function(stats, label) {
  set.seed(SEED)
  stats <- stats[!is.na(stats)]; stats <- stats[!duplicated(names(stats))]
  f <- suppressWarnings(fgsea(HM_SETS, stats, minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000, nproc = 1))
  data.frame(ranking = label, pathway = f$pathway, size = f$size, NES = f$NES, pval = f$pval, padj = f$padj,
             stringsAsFactors = FALSE)
}
msg("fgsea")
tr <- TUM$res[TUM$res$expressed & TUM$res$gene_id != ID_NKX, ]
GS <- rbind(
  run_fgsea(setNames(tr$rho_raw, tr$gene_id), "TCGA tumour: Spearman rho with NKX3-1"),
  run_fgsea(setNames(tr$rho_AR, tr$gene_id), "TCGA tumour: partial rho | AR score"),
  run_fgsea(setNames(tr$rho_AR_PROL, tr$gene_id), "TCGA tumour: partial rho | AR + proliferation"),
  run_fgsea(setNames(tr$rho_AR_PROL_STROMA, tr$gene_id), "TCGA tumour: partial rho | AR + proliferation + stroma"))
DU_LAB <- c(NKX_in_ctrl = "DU145 NKX3.1 vs control (basal)", NKX_in_CM = "DU145 NKX3.1 vs control (inflammation)",
            NKX_main = "DU145 NKX3.1 vs control (additive main effect)", CM_in_vector = "DU145 CM vs control (vector)")
DU_TREND <- list()
for (cn in names(DU_LAB)) {
  r <- core$fc$res[[cn]]; r <- r[r$gene_id != ID_NKX & !is.na(r$stat), ]
  r$lbm <- log10(r$baseMean)
  lo <- loess(stat ~ lbm, data = r, span = 0.3, family = "symmetric", degree = 1)
  r$stat_detrended <- r$stat - predict(lo, r)
  r$dec <- cut(r$lbm, quantile(r$lbm, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  DU_TREND[[cn]] <- data.frame(contrast = cn, baseMean_decile = 1:10,
    median_baseMean = tapply(r$baseMean, r$dec, median), median_wald_z = tapply(r$stat, r$dec, median),
    frac_negative = tapply(r$stat < 0, r$dec, mean), median_log2FC = tapply(r$log2FC_MLE, r$dec, median))
  GS <- rbind(GS, run_fgsea(setNames(r$stat, r$gene_id), DU_LAB[[cn]]),
              run_fgsea(setNames(r$stat_detrended, r$gene_id), paste0(DU_LAB[[cn]], ", expression trend removed")))
}
DU_TREND <- rbindlist(DU_TREND)
GSW <- dcast(as.data.table(GS), pathway ~ ranking, value.var = c("NES", "padj"))

# ----------------------------------------------------------------------------
# 8. Sensitivity: legacy RSEM (hg19, UCSC symbols): claimed-gene rho, global concordance, programme scores
# ----------------------------------------------------------------------------
RSEM <- NULL
if (have_legacy && primary_source == "GDC_STAR_TPM") {
  msg("legacy RSEM sensitivity")
  L <- read_xena(SRC$legacy$file)
  L <- L[, pick_samples(colnames(L), "01")$keep]; L <- L[rowVars(L) > 0, ]
  sym_map <- function(gid) {      # Ensembl -> Entrez (Ensembl 99 annotation) -> aliases -> first present in HiSeqV2
    eg <- ann99$entrez[match(gid, ann99$gene_id)]
    cand <- c(ann99$symbol[match(gid, ann99$gene_id)],
              if (!is.na(eg)) suppressMessages(AnnotationDbi::select(org.Hs.eg.db, eg, "ALIAS", "ENTREZID")$ALIAS))
    hit <- cand[cand %in% rownames(L)]; if (length(hit)) hit[1] else NA_character_
  }
  rl_nkx <- rank(L["NKX3-1", ])
  ULr <- unit_rows(row_ranks(L)); y <- (rl_nkx - mean(rl_nkx)) / sqrt(sum((rl_nkx - mean(rl_nkx))^2))
  rho_L <- setNames(cor_vs(ULr, y), rownames(L))
  s_star <- setNames(TUM$res$rho_raw, TUM$res$symbol)[TUM$res$expressed]
  common <- setdiff(intersect(names(s_star)[!duplicated(names(s_star))], names(rho_L)), "NKX3-1")
  RSEM <- list(n = ncol(L), n_common = length(common),
               concordance_pearson = cor(s_star[common], rho_L[common]),
               concordance_spearman = cor(s_star[common], rho_L[common], method = "spearman"),
               null_median_rho_rsem = median(rho_L[common]), null_median_rho_star = median(s_star[common]),
               claimed = data.frame(gene_id = ginfo$gene_id, symbol = ginfo$symbol,
                                    legacy_symbol = vapply(ginfo$gene_id, sym_map, character(1)), stringsAsFactors = FALSE))
  RSEM$claimed$rho_raw_legacyRSEM <- unname(rho_L[RSEM$claimed$legacy_symbol])
  RSEM$claimed$pctl_raw_legacyRSEM <- pctl(RSEM$claimed$rho_raw_legacyRSEM, rho_L[common])
  # programme scores in RSEM (genes with median log2 value >= log2(11), i.e. >= 10 normalised counts)
  ZL <- zrows(L); exL <- rownames(L)[rowMedians(L) >= log2(11)]
  ar_L <- rank(colMeans(ZL[intersect(AR_GENES, rownames(ZL)), ])); q_arL <- qr(cbind(1, ar_L))
  st_L <- rank(colMeans(ZL[intersect(STROMA_GENES, rownames(ZL)), ])); q_asL <- qr(cbind(1, ar_L, st_L))
  RSEM$programmes <- rbindlist(lapply(HALLMARKS, function(h) {
    mem <- setdiff(intersect(HM_SYM[[h]], exL), "NKX3-1"); sc <- rank(colMeans(ZL[mem, ]))
    data.frame(pathway = h, label = HM_LABEL[[h]], n_genes = length(mem), rho_raw = cor(rl_nkx, sc),
               rho_AR = cor(qr.resid(q_arL, rl_nkx), qr.resid(q_arL, sc)),
               rho_AR_STROMA = cor(qr.resid(q_asL, rl_nkx), qr.resid(q_asL, sc)))
  }))
  rm(L, ULr, ZL); gc()
}
if (!is.null(RSEM)) CT <- merge(CT, RSEM$claimed[, c("gene_id", "legacy_symbol", "rho_raw_legacyRSEM", "pctl_raw_legacyRSEM")],
                                by = "gene_id", all.x = TRUE, sort = FALSE)
KL <- TUM$res[TUM$res$gene_id == ID_KLRC3, ]

# ----------------------------------------------------------------------------
# 9. Tables
# ----------------------------------------------------------------------------
fwrite_tsv <- function(x, f) fwrite(x, file.path(tab_dir, f), sep = "\t", na = "NA", quote = FALSE)
ord_cols <- c("symbol", "gene_id", "category", "expected_sign", "in_activated", "in_repressed", "in_module_3D",
              "in_module_3E", "in_abstract", "median_TPM", "expressed",
              "rho_raw", "rho_raw_lo", "rho_raw_hi", "p_raw", "padj_claimed_raw", "padj_genomewide_raw", "pctl_raw", "agree_raw",
              "rho_AR", "p_AR", "padj_claimed_AR", "pctl_AR",
              "rho_PROL", "p_PROL", "padj_claimed_PROL", "pctl_PROL",
              "rho_AR_PROL", "rho_AR_PROL_lo", "rho_AR_PROL_hi", "p_AR_PROL", "padj_claimed_AR_PROL",
              "padj_genomewide_AR_PROL", "pctl_AR_PROL", "agree_AR_PROL",
              "rho_AR_PROL_STROMA", "p_AR_PROL_STROMA", "padj_claimed_AR_PROL_STROMA", "pctl_AR_PROL_STROMA", "agree_AR_PROL_STROMA")
cat_levels <- c("Activated (text/Fig 1)", "Repressed (text/Fig 1)", "Module only (Fig 3D/3E)")
ct_out <- as.data.frame(CT)[, c(ord_cols, intersect(c("legacy_symbol", "rho_raw_legacyRSEM", "pctl_raw_legacyRSEM"), colnames(CT)))]
ct_out <- ct_out[order(factor(ct_out$category, cat_levels), -ct_out$rho_raw), ]
cn_out <- as.data.frame(CN)[, ord_cols]; cn_out <- cn_out[match(ct_out$gene_id, cn_out$gene_id), ]
fwrite_tsv(ct_out, "tcga_claimed_genes_vs_NKX3-1_primary_tumours.tsv")
fwrite_tsv(cn_out, "tcga_claimed_genes_vs_NKX3-1_solid_normals.tsv")
fwrite_tsv(ST, "tcga_direction_sign_tests.tsv")
fwrite_tsv(NULLS, "tcga_genomewide_null_summary.tsv")
fwrite_tsv(COH, "tcga_module_coherence.tsv")
fwrite_tsv(PAIRS, "tcga_module_pairwise_spearman_primary_tumours.tsv")
fwrite_tsv(PW, "tcga_hallmark_scores_vs_NKX3-1.tsv")
fwrite_tsv(PWC, "tcga_hallmark_scores_vs_NKX3-1_competitive_null.tsv")
fwrite_tsv(PW_CTX, "tcga_hallmark_scores_vs_covariates_primary_tumours.tsv")
COV_CTX <- TUM$res[TUM$res$gene_id %in% c(ID_AR, ID_PROL, ID_STR), c("gene_id", "symbol", "median_TPM", "rho_raw", "pctl_raw")]
COV_CTX$score <- ifelse(COV_CTX$gene_id %in% ID_AR, "AR_score", ifelse(COV_CTX$gene_id %in% ID_PROL, "proliferation_score", "stromal_score"))
fwrite_tsv(COV_CTX, "tcga_covariate_genes_vs_NKX3-1_primary_tumours.tsv")
fwrite_tsv(GS, "tcga_vs_DU145_hallmark_fgsea_long.tsv")
fwrite_tsv(GSW, "tcga_vs_DU145_hallmark_fgsea_wide.tsv")
fwrite_tsv(DU_TREND, "DU145_wald_z_by_expression_decile.tsv")
if (!is.null(RSEM)) fwrite_tsv(RSEM$programmes, "tcga_legacyRSEM_hallmark_scores_vs_NKX3-1.tsv")
gw <- TUM$res[TUM$res$in_family, c("gene_id", "symbol", "median_TPM", "expressed", "rho_raw", "p_raw", "padj_genomewide_raw",
                  "pctl_raw", "rho_AR", "p_AR", "pctl_AR", "rho_AR_PROL", "p_AR_PROL", "padj_genomewide_AR_PROL", "pctl_AR_PROL",
                  "rho_AR_PROL_STROMA", "pctl_AR_PROL_STROMA")]
gw <- gw[order(-gw$rho_raw), ]
fwrite(gw, file.path(tab_dir, "tcga_genomewide_NKX3-1_spearman_primary_tumours.tsv.gz"), sep = "\t", na = "NA", quote = FALSE)

# ----------------------------------------------------------------------------
# 10. Provenance, results object, session info
# ----------------------------------------------------------------------------
prov <- list(
  task = "TCGA-PRAD independent-dataset check (Reviewer 1, major 5)",
  scripts = c("/mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga.R",
              "/mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_geo.py",
              "/mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_geo_rwpe1.R",
              "/mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_figures.R"),
  run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  primary_source = primary_source,
  sources = lapply(SRC[c(if (have_gdc) c("gdc_tpm", "gdc_pmap"), if (have_legacy) "legacy")], prov_entry),
  data_description = list(
    gdc_tpm = "UCSC Xena GDC hub, TCGA-PRAD gene expression RNAseq 'STAR - TPM' (GDC harmonised, GRCh38, GENCODE v36); values log2(TPM+1)",
    legacy = "UCSC Xena legacy TCGA hub, TCGA.PRAD.sampleMap/HiSeqV2 (RSEM normalised counts, hg19 UCSC gene symbols); values log2(norm_count+1)"),
  unit = P$unit, verification_tpm_colsum_range = P$tpm_colsum_range,
  n_genes_in_file = P$n_genes_raw, n_par_y_rows_removed = P$n_par_y_removed,
  n_samples_in_file = P$n_samples_total, sample_type_codes_in_file = P$type_table,
  n_primary_tumour_used = length(P$tum_ids), primary_duplicates_removed = P$tum_dropped,
  n_solid_normal_used = length(P$nor_ids), normal_duplicates_removed = P$nor_dropped,
  sample_selection = "Barcode positions 14-15: 01 = primary tumour (analysed), 11 = solid tissue normal (secondary), 06 = metastasis (excluded). One sample per patient and type (lexicographically first aliquot).",
  expression_filter = "expressed = median TPM >= 1 within cohort (tumours and normals separately); genome-wide null = expressed genes except NKX3-1; genome-wide BH family = expressed genes + claimed genes; correlations also computed for detected genes (mean TPM >= 0.1) for the 'all genes' module-coherence sensitivity",
  n_expressed_tumour = length(TUM$expressed), n_expressed_normal = length(NOR$expressed),
  gene_matching = "Ensembl gene IDs (version stripped) from the Ensembl 99 annotation used for DU145; all claimed/score genes present in GENCODE v36",
  msigdbr_version = as.character(packageVersion("msigdbr")), msigdb_db_version = msig_db_version,
  covariates = list(AR_score = AR_GENES, proliferation_score = PROL_GENES, stromal_score = STROMA_GENES,
                    definition = "mean of per-gene z-scores of log2(TPM+1) within cohort; partial Spearman = Pearson correlation of rank residuals after regressing ranks on covariate ranks"),
  seed = SEED, permutations = B_PERM,
  legacy_RSEM_sensitivity = if (!is.null(RSEM)) RSEM[c("n", "n_common", "concordance_pearson", "concordance_spearman",
                                                       "null_median_rho_rsem", "null_median_rho_star")] else NULL,
  context = list(rho_NKX3_1_vs_AR_score = TUM$rho_nkx_ar, rho_NKX3_1_vs_prolif_score = TUM$rho_nkx_pr,
                 rho_NKX3_1_vs_stromal_score = TUM$rho_nkx_st, rho_AR_vs_prolif = TUM$rho_ar_pr,
                 normals = list(rho_NKX3_1_vs_AR_score = NOR$rho_nkx_ar, rho_NKX3_1_vs_prolif_score = NOR$rho_nkx_pr,
                                rho_NKX3_1_vs_stromal_score = NOR$rho_nkx_st),
                 NKX3_1_median_TPM_tumour = median(2^P$X_tum[ID_NKX, ] - 1),
                 NKX3_1_median_TPM_normal = median(2^P$X_nor[ID_NKX, ] - 1),
                 KLRC3_tumour = as.list(KL[, c("median_TPM", "expressed", "rho_raw", "pctl_raw", "rho_AR_PROL", "pctl_AR_PROL")]))
)
writeLines(toJSON(prov, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"), file.path(out, "tcga_provenance.json"))
saveRDS(list(CT = as.data.frame(CT), CN = as.data.frame(CN), ST = ST, COH = COH, COH_NULL = COH_NULL, PW = PW, PWC = PWC,
             GS = GS, NULLS = NULLS, DU_TREND = DU_TREND, RSEM_programmes = if (!is.null(RSEM)) RSEM$programmes,
             prov = prov,
             scatter = list(X = P$X_tum[c(ID_NKX, ID_CLM[ABSTRACT_GENES]), ], sym = as.list(ID_CLM[ABSTRACT_GENES]),
                            id_nkx = unname(ID_NKX), ar = TUM$ar, pr = TUM$pr, st = TUM$st)),
        file.path(out, "tcga_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_15_tcga.txt"))
msg("analysis done")
