#!/usr/bin/env Rscript
# Reliability task, step 3: supplementary table for every gene named in the manuscript (text + figures),
# for both quantifications, + read-evidence classification + DESeq2 tests of paralogue/locus family totals.
#   corrected = core$fc  (STAR genome alignments, featureCounts: unique fragments, reverse-stranded, Ensembl 99 exons)
#   original  = core$sal (manuscript input: nf-core STAR->Salmon, tximport lengthScaledTPM, protein-coding tx2gene)
# Needs: 13_reliability_02_bam_audit.py outputs. Sign convention: treatment minus reference.
# Run: /home/biolab/miniconda3/envs/nkx3_r/bin/Rscript 13_reliability_03_table.R

suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(BiocParallel); library(dplyr); library(readr); library(jsonlite)
})
set.seed(42)
register(MulticoreParam(16))
proj <- "/mnt/ssd1/Projects/NKX3_wd"
out  <- file.path(proj, "results/revision_2026/reliability")
nf   <- "/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon"
core <- readRDS(file.path(proj, "results/revision_2026/core/core.rds"))
cd <- core$cd
grp_lab <- c(vector_ctrl = "Control-basal", vector_CM = "Control-inflammation",
             NKX3.1_ctrl = "NKX3.1-basal", NKX3.1_CM = "NKX3.1-inflammation")
contr <- c(CM_in_vector = "CM vs control (vector cells)", NKX_in_ctrl = "NKX3.1 vs control (basal)",
           NKX_in_CM = "NKX3.1 vs control (inflammation)", interaction = "Interaction (NKX3.1 x CM)",
           NKX_main = "NKX3.1 main effect (additive model)", CM_main = "CM main effect (additive model)")
ctx <- c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM", "interaction"); mains <- c("NKX_main", "CM_main")
stats_keep <- c(log2FC_MLE = "log2FC_MLE", log2FC_apeglm = "log2FC_shrunk", lfcSE_apeglm = "lfcSE_shrunk",
                Wald_stat = "stat", padj = "padj", padj_LFC1.5 = "padj_LFC1.5")

named <- read_tsv(file.path(out, "named_genes_curated.tsv"), show_col_types = FALSE)
ev    <- read_tsv(file.path(out, "read_evidence_per_sample.tsv"), show_col_types = FALSE)
fam_s <- read_tsv(file.path(out, "family_totals_per_sample.tsv"), show_col_types = FALSE)
fams  <- fromJSON(file.path(out, "families_definition.json"))

# ------------------------------------------------------------------ per-quantification gene block
quant_block <- function(obj, raw, tag, prefilter_label) {
  sf  <- sizeFactors(obj$dds_full)[cd$sample_id]
  ids <- named$gene_id
  inm <- ids %in% rownames(raw)
  res <- data.frame(gene_id = ids, check.names = FALSE)
  res[[paste0(tag, "_prefilter")]] <- ifelse(!inm, prefilter_label,
                                      ifelse(ids %in% names(obj$keep)[obj$keep], "tested",
                                             "filtered (<10 reads in >=3 samples)"))
  nrm <- matrix(NA_real_, length(ids), nrow(cd), dimnames = list(ids, cd$sample_id))
  rw  <- nrm
  nrm[inm, ] <- sweep(raw[ids[inm], cd$sample_id, drop = FALSE], 2, sf, "/")
  rw[inm, ]  <- raw[ids[inm], cd$sample_id, drop = FALSE]
  res[[paste0(tag, "_baseMean")]] <- rowMeans(nrm)
  for (g in names(grp_lab)) {
    s <- cd$sample_id[cd$group == g]
    res[[sprintf("%s_normcount_mean_%s", tag, grp_lab[g])]] <- rowMeans(nrm[, s, drop = FALSE])
    res[[sprintf("%s_normcount_sd_%s", tag, grp_lab[g])]]   <- apply(nrm[, s, drop = FALSE], 1, sd)
  }
  for (g in names(grp_lab)) {
    s <- cd$sample_id[cd$group == g]
    res[[sprintf("%s_rawcount_mean_%s", tag, grp_lab[g])]] <- rowMeans(rw[, s, drop = FALSE])
  }
  for (cn in c(ctx, mains)) {
    r <- obj$res[[cn]]
    for (st in names(stats_keep)) {
      res[[sprintf("%s_%s_%s", tag, cn, st)]] <- r[match(ids, r$gene_id), stats_keep[[st]]]
    }
  }
  res
}
fcb  <- quant_block(core$fc,  core$counts_fc,  "corrected", "not in count matrix")
salb <- quant_block(core$sal, core$counts_sal, "original",
                    "not in original count matrix (protein-coding tx2gene + CPM>=1 in >=3 samples filter)")

# ------------------------------------------------------------------ as published (original pairwise DESeq2, v15)
v15 <- file.path(proj, "results/rna_downstreams_out_v15")
pub_files <- c(CM_in_vector = "dgea_HM-vector-C_vs_HM-vector-CM_DESeq2_results.tsv",
               NKX_in_ctrl  = "dgea_HM-NKX3.1-C_vs_HM-vector-C_DESeq2_results.tsv",
               NKX_in_CM    = "dgea_HM-NKX3.1-CM_vs_HM-vector-CM_DESeq2_results.tsv")
pub_sign <- c(CM_in_vector = -1, NKX_in_ctrl = 1, NKX_in_CM = 1)   # v15 CM table was control minus CM
pub <- data.frame(gene_id = named$gene_id)
for (cn in names(pub_files)) {
  t <- read_tsv(file.path(v15, pub_files[[cn]]), show_col_types = FALSE)
  m <- match(named$gene_id, t$feature_id)
  pub[[sprintf("published_%s_log2FC", cn)]] <- pub_sign[[cn]] * t$log2FoldChange[m]
  pub[[sprintf("published_%s_fold", cn)]]   <- 2^(pub_sign[[cn]] * t$log2FoldChange[m])
  pub[[sprintf("published_%s_padj", cn)]]   <- t$padj[m]
}

# ------------------------------------------------------------------ family / locus totals: DESeq2 tests
fam_mat <- function(col) {
  w <- tidyr::pivot_wider(fam_s[, c("family", "sample_id", col)], names_from = sample_id, values_from = all_of(col))
  m <- as.matrix(w[, cd$sample_id]); rownames(m) <- paste0("FAMILY::", w$family); round(m)
}
test_rows <- function(base, add, sf, label) {
  m <- rbind(base, add)
  fits <- list()
  lv <- list(vc = list(NKX = c("vector", "NKX3.1"), CM = c("ctrl", "CM")),
             vcm = list(NKX = c("vector", "NKX3.1"), CM = c("CM", "ctrl")))
  for (k in names(lv)) {
    cdx <- cd; cdx$NKX <- factor(as.character(cdx$NKX), levels = lv[[k]]$NKX)
    cdx$CM <- factor(as.character(cdx$CM), levels = lv[[k]]$CM)
    d <- DESeqDataSetFromMatrix(m, cdx, ~ NKX + CM + NKX:CM)
    sizeFactors(d) <- sf
    fits[[k]] <- DESeq(d, quiet = TRUE, parallel = TRUE)
  }
  spec <- list(CM_in_vector = c("vc", "CM_CM_vs_ctrl"), NKX_in_ctrl = c("vc", "NKX_NKX3.1_vs_vector"),
               NKX_in_CM = c("vcm", "NKX_NKX3.1_vs_vector"), interaction = c("vc", "NKXNKX3.1.CMCM"))
  outl <- list()
  for (cn in names(spec)) {
    d <- fits[[spec[[cn]][1]]]; co <- spec[[cn]][2]
    r <- results(d, name = co, alpha = 0.05)
    s <- lfcShrink(d, coef = co, type = "apeglm", res = r, quiet = TRUE, parallel = TRUE)
    i <- rownames(add)
    outl[[cn]] <- data.frame(family = sub("^FAMILY::", "", i), quantification = label, contrast = cn,
                             log2FC_MLE = r[i, "log2FoldChange"], log2FC_apeglm = s[i, "log2FoldChange"],
                             lfcSE_apeglm = s[i, "lfcSE"], Wald_stat = r[i, "stat"], pvalue = r[i, "pvalue"],
                             padj = r[i, "padj"])
  }
  dl <- DESeqDataSetFromMatrix(m, cd, ~ group); sizeFactors(dl) <- sf
  dl <- DESeq(dl, test = "LRT", reduced = ~ 1, quiet = TRUE, parallel = TRUE)
  rl <- results(dl)
  lrt <- data.frame(family = sub("^FAMILY::", "", rownames(add)), quantification = label,
                    LRT_any_group_difference_pvalue = rl[rownames(add), "pvalue"],
                    LRT_any_group_difference_padj = rl[rownames(add), "padj"])
  list(wald = bind_rows(outl), lrt = lrt)
}
message("Family tests: genome fragments")
base_fc <- round(core$counts_fc[core$fc$keep, cd$sample_id])
ft_g <- test_rows(base_fc, fam_mat("genome_fragments_distinct"), sizeFactors(core$fc$dds_full)[cd$sample_id],
                  "genome (distinct fragments, unique + multi-mapped)")
message("Family tests: Salmon NumReads sums")
sal_nf <- read.delim(file.path(nf, "salmon.merged.gene_counts.tsv"), check.names = FALSE)
colnames(sal_nf) <- sub("^X", "", gsub(".", "-", colnames(sal_nf), fixed = TRUE))
smat <- round(as.matrix(sal_nf[, cd$sample_id])); rownames(smat) <- sal_nf$gene_id
smat <- smat[rowSums(smat >= 10) >= 3, ]
sf_sal <- estimateSizeFactorsForMatrix(smat)
ft_s <- test_rows(smat, fam_mat("salmon_numreads_sum"), sf_sal, "Salmon (nf-core NumReads, sum of members)")
fam_tests <- bind_rows(ft_g$wald, ft_s$wald)
fam_lrt   <- bind_rows(ft_g$lrt, ft_s$lrt)
fam_means <- fam_s %>% group_by(family, members, group) %>%
  summarise(genome = mean(genome_fragments_distinct), genome_unique = mean(genome_unique_sum),
            salmon = mean(salmon_numreads_sum), featurecounts = mean(featurecounts_unique_sum), .groups = "drop")
fam_wide <- fam_means %>% mutate(group = grp_lab[group]) %>%
  tidyr::pivot_wider(names_from = group, values_from = c(genome, genome_unique, salmon, featurecounts),
                     names_glue = "{.value}_mean_{group}")
fam_out <- fam_wide %>% left_join(fam_lrt %>% filter(grepl("^genome", quantification)) %>%
                                    select(family, genome_LRT_padj = LRT_any_group_difference_padj,
                                           genome_LRT_pvalue = LRT_any_group_difference_pvalue), by = "family") %>%
  left_join(fam_lrt %>% filter(grepl("^Salmon", quantification)) %>%
              select(family, salmon_LRT_padj = LRT_any_group_difference_padj,
                     salmon_LRT_pvalue = LRT_any_group_difference_pvalue), by = "family")
for (cn in ctx) for (q in c("genome", "Salmon")) {
  t <- fam_tests %>% filter(contrast == cn, grepl(paste0("^", q), quantification))
  pre <- sprintf("%s_%s", tolower(q), cn)
  fam_out[[paste0(pre, "_log2FC_apeglm")]] <- t$log2FC_apeglm[match(fam_out$family, t$family)]
  fam_out[[paste0(pre, "_log2FC_MLE")]]    <- t$log2FC_MLE[match(fam_out$family, t$family)]
  fam_out[[paste0(pre, "_padj")]]          <- t$padj[match(fam_out$family, t$family)]
}
write_tsv(fam_out, file.path(out, "family_totals_tests.tsv"))
write_tsv(fam_tests, file.path(out, "family_totals_tests_long.tsv"))

# ------------------------------------------------------------------ read-evidence classification
sf_fc <- sizeFactors(core$fc$dds_full)
ev$sf <- sf_fc[ev$sample_id]
gm <- function(col, norm = FALSE) {
  v <- ev[[col]]; if (norm) v <- v / ev$sf
  t <- tapply(v, list(ev$gene_id, ev$group), mean)
  t[, names(grp_lab), drop = FALSE]
}
S <- gm("salmon_numreads"); Go <- gm("sense_exonic_all"); Gc <- gm("sense_contained_all")
Gu <- gm("sense_exonic_unique"); Fc <- gm("featurecounts_unique")
Sn <- gm("salmon_numreads", TRUE); Gn <- gm("sense_contained_all", TRUE)
tot <- ev %>% group_by(gene_id) %>% summarise(across(c(sense_exonic_all, sense_exonic_multi, sense_exonic_shared,
                                                        sense_exonic_dup, sense_contained_all, antisense_exonic_all,
                                                        sense_intronic_all, salmon_numreads, featurecounts_unique), sum))
tot <- as.data.frame(tot); rownames(tot) <- tot$gene_id
ctr_def <- list(CM_in_vector = c("vector_CM", "vector_ctrl"), NKX_in_ctrl = c("NKX3.1_ctrl", "vector_ctrl"),
                NKX_in_CM = c("NKX3.1_CM", "vector_CM"))
fmt <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
fam_lookup <- function(gid) {
  f <- named$families[named$gene_id == gid]
  if (is.na(f) || f == "") return(character(0))
  strsplit(f, "; ")[[1]]
}
classify <- function(gid) {
  s <- S[gid, ]; go <- Go[gid, ]; gc <- Gc[gid, ]
  claims <- named$reported_claims[named$gene_id == gid]
  claims <- if (is.na(claims) || claims == "") character(0) else strsplit(claims, "; ")[[1]]
  mf <- tot[gid, "sense_exonic_multi"] / max(tot[gid, "sense_exonic_all"], 1)
  sh <- tot[gid, "sense_exonic_shared"] / max(tot[gid, "sense_exonic_all"], 1)
  # directions for reported contrasts (size-factor-normalised group means; pseudo-count 0.5)
  dirs <- sapply(claims, function(cl) {
    p <- strsplit(cl, ":")[[1]]; cdf <- ctr_def[[p[1]]]
    dS <- log2((Sn[gid, cdf[1]] + 0.5) / (Sn[gid, cdf[2]] + 0.5)); dG <- log2((Gn[gid, cdf[1]] + 0.5) / (Gn[gid, cdf[2]] + 0.5))
    c(dS = dS, dG = dG)
  })
  dir_ok <- if (length(claims)) all(apply(dirs, 2, function(v) (v["dS"] * v["dG"] > 0) || (abs(v["dS"]) < 0.2 && abs(v["dG"]) < 0.2))) else TRUE
  excess <- (s >= 10) & (s > 1.5 * go) & (s - go >= 10)
  big <- pmax(s, gc) >= 10
  ratio <- (s + 0.5) / (gc + 0.5)
  agree <- all(!big | (ratio >= 1 / 1.5 & ratio <= 1.5))
  cls <- NA; reason <- NA
  if (max(pmax(s, go)) < 10) {
    cls <- "low expression"
    reason <- sprintf("all sources < 10 fragments per sample on average (Salmon max %s, genome max %s)", fmt(max(s)), fmt(max(go)))
  } else if (any(excess)) {
    cls <- "not supported by read evidence"
    gx <- names(which(excess))
    reason <- sprintf("Salmon assigns %s fragments in %s but only %s sense fragments overlap the gene's exons in the genome alignments (unique + multi-mapped); Salmon-only fragments %s",
                      paste(fmt(s[gx]), collapse = "/"), paste(grp_lab[gx], collapse = "/"), paste(fmt(go[gx]), collapse = "/"),
                      if (all(gx %in% c("NKX3.1_ctrl", "NKX3.1_CM"))) "confined to NKX3.1-transfected samples"
                      else if (all(gx %in% c("vector_ctrl", "vector_CM"))) "confined to empty-vector samples" else "in several groups")
  } else if (max(go) >= 10 && (mf >= 0.5 || sh >= 0.5)) {
    cls <- "mappability-limited"
    reason <- sprintf("%.0f%% of locus fragments multi-mapped (NH>1), %.0f%% shared with exons of another same-strand gene; gene-level estimate not identifiable",
                      100 * mf, 100 * sh)
  } else if (agree && dir_ok) {
    cls <- "supported"
    reason <- sprintf("Salmon within 1.5-fold of exon-contained genome fragments in every group with >=10 fragments (ratios %s); directions agree",
                      paste(sprintf("%.2f", ratio[big]), collapse = "/"))
  } else {
    cls <- "not supported by read evidence"
    reason <- if (!agree) sprintf("Salmon/genome ratio outside 1.5-fold (%s) without multi-mapping or overlapping-gene explanation",
                                  paste(sprintf("%.2f", ratio[big]), collapse = "/"))
              else "reported direction differs between Salmon and genome alignments"
  }
  data.frame(gene_id = gid, read_evidence_class = cls, read_evidence_reason = reason,
             multi_mapped_fraction = mf, shared_fraction = sh,
             duplicate_fraction = tot[gid, "sense_exonic_dup"] / max(tot[gid, "sense_exonic_all"], 1),
             salmon_to_genome_ratio_max = if (any(big)) max(ratio[big]) else NA_real_,
             salmon_to_genome_ratio_min = if (any(big)) min(ratio[big]) else NA_real_,
             directions_agree = dir_ok,
             reported_dir_salmon = if (length(claims)) paste(sprintf("%s:%+.2f", sub(":.*", "", claims), dirs["dS", ]), collapse = "; ") else "",
             reported_dir_genome = if (length(claims)) paste(sprintf("%s:%+.2f", sub(":.*", "", claims), dirs["dG", ]), collapse = "; ") else "")
}
cls <- suppressWarnings(bind_rows(lapply(named$gene_id, classify)))

# family-level verdict for mappability-limited / shared genes
fam_verdict <- sapply(named$gene_id, function(gid) {
  fs <- fam_lookup(gid); if (!length(fs)) return("")
  claims <- named$reported_claims[named$gene_id == gid]
  claims <- if (is.na(claims) || claims == "") character(0) else sub(":.*", "", strsplit(claims, "; ")[[1]])
  if (!length(claims)) claims <- c("CM_in_vector", "NKX_in_ctrl", "NKX_in_CM")
  paste(sapply(fs, function(f) {
    t <- fam_tests %>% filter(family == f, grepl("^genome", quantification), contrast %in% claims)
    mx <- max(fam_means$genome[fam_means$family == f])
    if (!nrow(t) || all(is.na(t$padj)) || mx < 10)
      return(sprintf("%s total (genome): too few fragments to test (max group mean %.1f per sample)", f, mx))
    sprintf("%s total (genome): %s", f, paste(sprintf("%s log2FC(MLE) %+.2f, padj %s", t$contrast, t$log2FC_MLE,
                                                   formatC(t$padj, format = "g", digits = 2)), collapse = "; "))
  }), collapse = " | ")
})
cls$family_level_test <- fam_verdict[cls$gene_id]

# evidence group means (raw) for the summary sheet
evsum <- data.frame(gene_id = named$gene_id)
for (nm in c("salmon", "genome_overlap", "genome_contained", "genome_unique", "featurecounts")) {
  M <- switch(nm, salmon = S, genome_overlap = Go, genome_contained = Gc, genome_unique = Gu, featurecounts = Fc)
  for (g in names(grp_lab)) evsum[[sprintf("%s_mean_%s", nm, grp_lab[g])]] <- M[named$gene_id, g]
}
tt <- tot[named$gene_id, ]
evsum$antisense_exonic_total <- tt$antisense_exonic_all; evsum$sense_intronic_total <- tt$sense_intronic_all

# corrected-analysis verdict on the reported contrasts
verdict <- sapply(named$gene_id, function(gid) {
  claims <- named$reported_claims[named$gene_id == gid]
  if (is.na(claims) || claims == "") return("")
  paste(sapply(strsplit(claims, "; ")[[1]], function(cl) {
    p <- strsplit(cl, ":")[[1]]; r <- core$fc$res[[p[1]]]; i <- match(gid, r$gene_id)
    if (is.na(i)) return(sprintf("%s: not tested (pre-filter)", p[1]))
    sig <- !is.na(r$padj[i]) && r$padj[i] < 0.05
    same <- sign(r$log2FC_shrunk[i]) == ifelse(p[2] == "up", 1, -1)
    sprintf("%s (reported %s): log2FC MLE %+.2f / apeglm %+.2f, padj %s -> %s", p[1], p[2], r$log2FC_MLE[i],
            r$log2FC_shrunk[i], formatC(r$padj[i], format = "g", digits = 2),
            ifelse(sig & same, "reproduced", ifelse(sig, "significant, opposite direction", "not significant")))
  }), collapse = "; ")
})

final <- named %>% left_join(cls, by = "gene_id") %>%
  mutate(corrected_reported_contrast_result = verdict[gene_id]) %>%
  left_join(evsum, by = "gene_id") %>% left_join(fcb, by = "gene_id") %>% left_join(salb, by = "gene_id") %>%
  left_join(pub, by = "gene_id")
write_tsv(final, file.path(out, "named_genes_reliability_full.tsv"), na = "")

# compact class counts
cnt <- final %>% mutate(tier_short = substr(tier, 1, 1)) %>% count(tier_short, read_evidence_class) %>%
  tidyr::pivot_wider(names_from = read_evidence_class, values_from = n, values_fill = 0)
write_tsv(cnt, file.path(out, "class_counts_by_tier.tsv"))
print(as.data.frame(cnt))
print(final %>% filter(substr(tier, 1, 1) == "A") %>% select(symbol, read_evidence_class) %>% as.data.frame())
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_13_reliability_03.txt"))
