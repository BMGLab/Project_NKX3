#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tximport)
  library(DESeq2)
  library(edgeR)
  library(pheatmap)
})

# ---------------------- CLI args ----------------------
opt <- OptionParser(option_list = list(
  make_option("--data_dir", type="character", default="data/star_salmon",
              help="Root folder containing sample subfolders with Salmon quant.sf files [default %default]"),
  make_option("--sample_table", type="character", default=NA,
            help="Path to CSV with columns: sample_id,condition,file. If set, bypasses auto-discovery."),
  make_option("--out_dir", type="character", default="results/rna_downstreams_out",
              help="Output directory for figures and objects [default %default]"),
  make_option("--min_cpm", type="double", default=1.0,
              help="CPM threshold for filtering [default %default]"),
  make_option("--min_samps_frac", type="double", default=0.25,
              help="Keep features with CPM>=min_cpm in at least this fraction of samples [default %default]"),
  make_option("--tx2gene", type="character", default=NA,
              help="Path to tx2gene mapping (2 columns: transcript_id,gene_id). If not provided, will try to auto-locate *.tx2gene.{csv,tsv}."),
  make_option("--var_q_low",  type="double", default=0.05,
            help="Lower variance quantile to keep (0–1). Set 0 to disable lower cut. [default %default]"),
make_option("--var_q_high", type="double", default=0.995,
            help="Upper variance quantile to keep (0–1). Set 1 to disable upper cut. [default %default]"),
make_option("--var_mode",   type="character", default="both",
            help="Variance filter mode: 'both', 'lower', 'upper', 'none'. [default %default]"),
make_option("--var_min_keep", type="integer", default=5000,
            help="Minimum number of genes to retain after variance filtering; fallback widens if needed. [default %default]")
))
args <- parse_args(opt)

data_dir <- normalizePath(args$data_dir, mustWork = TRUE)
out_dir  <- args$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

message("Data dir: ", data_dir)
message("Out dir : ", normalizePath(out_dir, mustWork = FALSE))

# ---------------------- Samples: either user-supplied CSV or auto-discovery ----------------------
sample_table <- NULL

if (!is.na(args$sample_table)) {
  message("Using user-supplied sample table: ", args$sample_table)
  sample_table <- readr::read_csv(args$sample_table, show_col_types = FALSE)

  # Basic validation
  req_cols <- c("sample_id", "condition", "file")
  if (!all(req_cols %in% colnames(sample_table))) {
    stop("sample_table is missing required columns. Found: ",
         paste(colnames(sample_table), collapse = ", "),
         "\nRequired: ", paste(req_cols, collapse = ", "))
  }

  # Ensure character types and factorize condition
  sample_table <- sample_table %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      condition = as.factor(condition),
      file      = as.character(file)
    )

  # Check file existence
  missing_files <- sample_table$file[!file.exists(sample_table$file)]
  if (length(missing_files) > 0) {
    stop("These 'file' paths from sample_table do not exist:\n",
         paste(missing_files, collapse = "\n"))
  }

  # Enforce unique sample IDs
  if (anyDuplicated(sample_table$sample_id)) {
    dupes <- unique(sample_table$sample_id[duplicated(sample_table$sample_id)])
    stop("Duplicate sample_id(s) in sample_table: ", paste(dupes, collapse = ", "),
         "\nPlease make sample_id unique.")
  }

  # Use the given order
  quant_files <- sample_table$file
  names(quant_files) <- sample_table$sample_id
  sample_ids <- sample_table$sample_id

} else {
  # ---- Fallback: auto-discover quant.sf ----
  quant_files <- list.files(data_dir, pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
  if (length(quant_files) == 0) stop("No quant.sf files found under: ", data_dir)

  sample_ids <- basename(dirname(quant_files))
  infer_cond <- function(x) {
    m <- stringr::str_match(x, "^([^_\\-]+)")
    ifelse(is.na(m[,2]), x, m[,2])
  }
  sample_table <- tibble::tibble(
    sample_id = sample_ids,
    file = quant_files,
    condition = factor(infer_cond(sample_ids))
  )

  # Write template for manual edits
  readr::write_csv(sample_table %>% dplyr::select(sample_id, condition, file),
                   file.path(out_dir, "sample_table.csv"))
  message("Wrote sample table: ", file.path(out_dir, "sample_table.csv"))
}


# ---------------------- Locate tx2gene mapping (optional) ----------------------
tx2gene_path <- args$tx2gene
if (is.na(tx2gene_path)) {
  cand <- c(
    list.files(data_dir, pattern = "\\.tx2gene\\.(csv|tsv)$", recursive = TRUE, full.names = TRUE),
    list.files(dirname(data_dir), pattern = "\\.tx2gene\\.(csv|tsv)$", recursive = TRUE, full.names = TRUE)
  )
  tx2gene_path <- if (length(cand) > 0) cand[1] else NA
}
use_tx2gene <- FALSE
tx2gene <- NULL
if (!is.na(tx2gene_path) && file.exists(tx2gene_path)) {
  message("Using tx2gene mapping: ", tx2gene_path)
  if (grepl("\\.tsv$", tx2gene_path, ignore.case = TRUE)) {
    tx2gene <- readr::read_tsv(tx2gene_path, col_names = c("transcript","gene"), show_col_types = FALSE)
  } else {
    tx2gene <- readr::read_csv(tx2gene_path, col_names = c("transcript","gene"), show_col_types = FALSE)
  }
  stopifnot(all(c("transcript","gene") %in% colnames(tx2gene)))
  use_tx2gene <- TRUE
} else {
  message("No tx2gene provided/found. Proceeding at transcript level (DESeq2 can run, but gene-level is preferred).")
}

# ---------------------- Import with tximport ----------------------
names(quant_files) <- sample_ids

txi <- if (use_tx2gene) {
  tximport(
    files = quant_files,
    type = "salmon",
    tx2gene = as.data.frame(tx2gene),
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
  )
} else {
  # No tx2gene: work at transcript level
  tximport(
    files = quant_files,
    type = "salmon",
    txOut = TRUE,                 # <<< important
    countsFromAbundance = "no",   # keep raw estimated counts
    ignoreTxVersion = TRUE
  )
}


# ---------------------- Basic QC: library sizes, detected features ----------------------
counts <- round(txi$counts)
if (is.null(rownames(counts))) {
  stop("Counts matrix has no rownames; tximport likely failed. Check inputs.")
}
stopifnot(identical(colnames(counts), sample_ids))

lib_sizes <- colSums(counts)
detected  <- colSums(counts > 0)

qc_df <- tibble(sample_id = names(lib_sizes),
                library_size = as.numeric(lib_sizes),
                detected_features = as.numeric(detected)) %>%
  left_join(sample_table %>% select(sample_id, condition), by = "sample_id")

readr::write_csv(qc_df, file.path(out_dir, "qc_library_sizes.csv"))

# Library sizes plot
p_libsize <- ggplot(qc_df, aes(x = reorder(sample_id, library_size), y = library_size, fill = condition)) +
  geom_col() + coord_flip() +
  labs(x = "Sample", y = "Library size (sum of counts)", title = "Library sizes") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "qc_library_sizes.pdf"), p_libsize, width = 8, height = 6, device = "pdf")

# Detected features plot
p_detected <- ggplot(qc_df, aes(x = reorder(sample_id, detected_features), y = detected_features, fill = condition)) +
  geom_col() + coord_flip() +
  labs(x = "Sample", y = "Detected features (>0 counts)", title = "Detected features per sample") +
  theme_bw(base_size = 12)
ggsave(file.path(out_dir, "qc_detected_features.pdf"), p_detected, width = 8, height = 6, device = "pdf")

# ---------------------- Filtering low-count features (CPM-based) ----------------------
dge <- DGEList(counts = counts)
cpm_mat <- cpm(dge)
min_samps <- ceiling(ncol(counts) * args$min_samps_frac)
# Harsher filter parameters
args$min_cpm <- 1.0
args$min_samps_frac <- 0.25

min_samps <- ceiling(ncol(counts) * args$min_samps_frac)
keep <- rowSums(cpm_mat >= args$min_cpm) >= min_samps

message(sprintf(
  "Harsher filtering: keeping %d/%d features (CPM >= %.2f in >= %d samples).",
  sum(keep), nrow(counts), args$min_cpm, min_samps
))
counts_filt <- counts[keep, , drop = FALSE]


# ---------------------- Build and align colData ----------------------
# Keep the sample order from counts_filt
samples <- colnames(counts_filt)

# Derive coldata from the user-supplied sample_table
coldata <- sample_table |>
  dplyr::filter(sample_id %in% samples) |>
  dplyr::distinct(sample_id, .keep_all = TRUE) |>
  dplyr::select(sample_id, condition) |>
  as.data.frame()

# Row names must be sample IDs, and order must match counts_filt
row.names(coldata) <- coldata$sample_id
coldata <- coldata[samples, , drop = FALSE]

# Clean factor levels (drop unused, keep order stable)
coldata$condition <- droplevels(factor(coldata$condition))

# Final sanity check
stopifnot(identical(colnames(counts_filt), row.names(coldata)))

# Save filtered counts
saveRDS(counts_filt, file.path(out_dir, "counts_filtered.rds"))
write_tsv(as.data.frame(counts_filt) %>% tibble::rownames_to_column("feature_id"),
          file.path(out_dir, "counts_filtered.tsv"))

# ---------------------- Create DESeq2 dataset ----------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ condition
)

# Estimate size factors (library-size normalization)
dds <- DESeq2::estimateSizeFactors(dds)

# ---------------------- VST normalization (use local fit) ----------------------
# 'fitType="local"' avoids the parametric-trend message and is robust
vsd <- DESeq2::vst(dds, blind = TRUE, fitType = "local")

# Always materialize this matrix immediately and check it
vsd_mat <- SummarizedExperiment::assay(vsd)
if (is.null(dim(vsd_mat)) || ncol(vsd_mat) < 2) {
  stop("VST returned an empty or 1-column matrix; cannot proceed with correlation/PCA.")
}

# ---------------------- Variance-based filtering (loosened & configurable) ----------------------
gene_var <- apply(vsd_mat, 1, var, na.rm = TRUE)

# Clamp/validate inputs
qL <- max(0, min(1, args$var_q_low))
qH <- max(0, min(1, args$var_q_high))
mode <- match.arg(tolower(args$var_mode), c("both","lower","upper","none"))

# Compute quantiles (only once)
qs <- quantile(gene_var, probs = sort(unique(c(qL, qH))), na.rm = TRUE)
# Compute thresholds directly (no name indexing)
q_low  <- if (qL > 0)  as.numeric(quantile(gene_var, probs = qL,  na.rm = TRUE)) else -Inf
q_high <- if (qH < 1)  as.numeric(quantile(gene_var, probs = qH,  na.rm = TRUE)) else  Inf

# Build keep mask based on mode
keep_var <- switch(
  mode,
  "none"  = rep(TRUE, length(gene_var)),
  "lower" = gene_var >= q_low,
  "upper" = gene_var <= q_high,
  "both"  = (gene_var >= q_low) & (gene_var <= q_high)
)

# Ensure we keep at least var_min_keep genes; if not, relax to lower-only then none
n_keep <- sum(keep_var, na.rm = TRUE)
if (n_keep < args$var_min_keep && mode %in% c("both","upper")) {
  message(sprintf("Variance filter kept %d < %d genes; relaxing to lower-only cutoff.", n_keep, args$var_min_keep))
  keep_var <- gene_var >= q_low
  n_keep <- sum(keep_var, na.rm = TRUE)
}
if (n_keep < args$var_min_keep && mode != "none") {
  message(sprintf("Variance filter kept %d < %d genes; disabling variance filter.", n_keep, args$var_min_keep))
  keep_var <- rep(TRUE, length(gene_var))
  n_keep <- length(keep_var)
}

message(sprintf(
  "Variance filter (%s): keeping %d/%d genes with thresholds [%s q=%.3f -> %.4g, %s q=%.3f -> %.4g].",
  mode, n_keep, length(gene_var),
  ifelse(mode %in% c('both','lower'), "low", "off"), qL, ifelse(mode %in% c('both','lower'), q_low, NA_real_),
  ifelse(mode %in% c('both','upper'), "high", "off"), qH, ifelse(mode %in% c('both','upper'), q_high, NA_real_)
))

#
#vsd_mat <- vsd_mat[keep_var, , drop = FALSE]



# Save objects for downstream DE
saveRDS(dds, file.path(out_dir, "dds_prepped_for_deseq2.rds"))
saveRDS(vsd, file.path(out_dir, "vsd.rds"))

# ---------------------- Build Ensembl ID → Gene Symbol mapping ----------------------
library(rtracklayer)

gtf_path <- "data/Homo_sapiens.GRCh38.99.gtf"  # adjust path if needed
if (file.exists(gtf_path)) {
  message("Loading GTF for gene name mapping: ", gtf_path)
  gtf_gr <- rtracklayer::import(gtf_path)
  
  # Extract gene-level entries
  gene_meta <- as.data.frame(mcols(gtf_gr[gtf_gr$type == "gene"]))
  
  # Build mapping table
  if (all(c("gene_id", "gene_name") %in% colnames(gene_meta))) {
    gene_annot <- gene_meta[, c("gene_id", "gene_name")]
    gene_annot$gene_id <- sub("\\.\\d+$", "", gene_annot$gene_id)  # strip version numbers
    gene_annot <- gene_annot[!duplicated(gene_annot$gene_id), ]
  } else {
    warning("Could not find 'gene_id' or 'gene_name' in GTF; skipping symbol mapping.")
    gene_annot <- NULL
  }
} else {
  warning("GTF not found at ", gtf_path, "; skipping gene name mapping.")
  gene_annot <- NULL
}


# ---------------------- Replace Ensembl IDs with gene symbols ----------------------
if (exists("gene_annot") && !is.null(gene_annot)) {
  ids <- rownames(vsd_mat)
  ids_nover <- sub("\\.\\d+$", "", ids)
  symbols <- gene_annot$gene_name[match(ids_nover, gene_annot$gene_id)]
  
  # Keep gene symbol if available, else fallback to Ensembl ID
  clean_names <- ifelse(!is.na(symbols) & symbols != "", symbols, ids)
  
  rownames(vsd_mat) <- clean_names
}


# ===================== Two-group comparisons: scatter + z vs LFC =====================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

# ---- Parameters you can tweak
min_expr_any <- 2.0        # require at least one group's mean VST >= this
min_abs_lfc  <- 0.25      # ~≥2x on log2-ish scale
min_abs_z    <- 2.0      # |z| threshold
top_label    <- 20       # label this many most extreme genes by |z|

# ---- Helpers (handle Unicode dashes and safe filenames)
normalize_dash <- function(x) {
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x, perl = TRUE) # various dashes → '-'
  trimws(x)
}
safe_fname <- function(x) {
  x <- normalize_dash(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x
}

# Guess ID type (SYMBOL vs ENSEMBL), map → ENTREZ with Bioconductor
detect_id_type <- function(ids) {
  ids_nv <- sub("\\.\\d+$", "", ids)
  frac_ens <- mean(grepl("^ENSG\\d+$", ids_nv))
  if (is.na(frac_ens)) return("SYMBOL")
  if (frac_ens > 0.5) "ENSEMBL" else "SYMBOL"
}

map_to_entrez <- function(x_ids, fromType) {
  x_ids <- sub("\\.\\d+$", "", x_ids)
  suppressWarnings(
    bitr(x_ids, fromType = fromType, toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
}


# ---- Robust msigdbr → TERM2GENE (works across msigdbr versions) ----
msig_t2g <- function(collection, subcategory = NULL) {
  df <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
  if (nrow(df) == 0) {
    message("msigdbr: no rows for collection=", collection, ". Skipping.")
    return(NULL)
  }

  # Optional subcategory filter only if column exists and value is present
  if (!is.null(subcategory)) {
    if (!"gs_subcat" %in% colnames(df)) {
      message("msigdbr: 'gs_subcat' column not present for collection ", collection,
              ". Skipping requested subcategory '", subcategory, "'.")
      return(NULL)
    }
    if (!subcategory %in% unique(df$gs_subcat)) {
      message("msigdbr: subcategory '", subcategory, "' not found in collection ", collection,
              ". Available: ", paste(na.omit(unique(df$gs_subcat)), collapse=", "),
              ". Skipping this subcategory.")
      return(NULL)
    }
    df <- df[df$gs_subcat == subcategory, , drop = FALSE]
    if (nrow(df) == 0) {
      message("msigdbr: after filtering subcategory '", subcategory, "' got 0 rows. Skipping.")
      return(NULL)
    }
  }

  # Prefer ENTREZ if present
  if ("entrez_gene" %in% colnames(df)) {
    t2g <- df[, c("gs_name", "entrez_gene")]
    t2g <- t2g[!is.na(t2g$entrez_gene), , drop = FALSE]
    colnames(t2g) <- c("term", "gene")
    return(unique(t2g))
  }

  # Else map SYMBOL → ENTREZ
  if ("gene_symbol" %in% colnames(df)) {
    sym <- unique(df[, c("gs_name", "gene_symbol")])
    map <- suppressWarnings(
      clusterProfiler::bitr(sym$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    )
    m <- merge(sym, map, by.x="gene_symbol", by.y="SYMBOL", all.x=TRUE)
    t2g <- m[!is.na(m$ENTREZID), c("gs_name","ENTREZID")]
    colnames(t2g) <- c("term","gene")
    return(unique(t2g))
  }

  # Else map ENSEMBL → ENTREZ
  if ("ensembl_gene" %in% colnames(df)) {
    ens <- unique(df[, c("gs_name", "ensembl_gene")])
    map <- suppressWarnings(
      clusterProfiler::bitr(ens$ensembl_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    )
    m <- merge(ens, map, by.x="ensembl_gene", by.y="ENSEMBL", all.x=TRUE)
    t2g <- m[!is.na(m$ENTREZID), c("gs_name","ENTREZID")]
    colnames(t2g) <- c("term","gene")
    return(unique(t2g))
  }

  message("msigdbr: no usable ID columns (entrez_gene/gene_symbol/ensembl_gene). Skipping.")
  return(NULL)
}




# Normalize condition labels (robust to Unicode minus)
cond_vec <- factor(normalize_dash(as.character(coldata$condition)))
stopifnot(all(colnames(vsd_mat) == rownames(coldata)))
levels(cond_vec)  # for your info

# ---- Define your comparisons (use the labels as they *appear* in your coldata, dashes OK)
comparisons <- list(
  c("HM-vector-C",   "HM-vector-CM"),
  c("HM-NKX3.1-C",   "HM-vector-C"),
  c("HM-NKX3.1-CM",  "HM-vector-CM")
)

# ---- Core function: compute stats and plots for one comparison
run_comparison <- function(g1, g2) {
  g1n <- normalize_dash(g1); g2n <- normalize_dash(g2)
  i1 <- which(cond_vec == g1n); i2 <- which(cond_vec == g2n)
  if (length(i1) < 2 || length(i2) < 2) {
    stop(sprintf("Need ≥2 samples in each group. Found: %s=%d, %s=%d", g1n, length(i1), g2n, length(i2)))
  }

  # Per-gene means/vars
  mu1 <- rowMeans(vsd_mat[, i1, drop = FALSE])
  mu2 <- rowMeans(vsd_mat[, i2, drop = FALSE])
  v1  <- apply(vsd_mat[, i1, drop = FALSE], 1, var)
  v2  <- apply(vsd_mat[, i2, drop = FALSE], 1, var)
  n1  <- length(i1); n2 <- length(i2)

  # Effect sizes
  lfc <- mu1 - mu2
  se  <- sqrt(pmax(v1, 1e-6)/n1 + pmax(v2, 1e-6)/n2)  # Welch-SE on VST scale
  z   <- lfc / se
  p   <- 2 * pnorm(-abs(z))
  q   <- p.adjust(p, method = "BH")

  # Significance mask
  expressed_ok <- pmax(mu1, mu2) >= min_expr_any
  sig <- expressed_ok & (abs(lfc) >= min_abs_lfc) & (abs(z) >= min_abs_z)

  res <- data.frame(
    gene = rownames(vsd_mat),
    mean_g1 = mu1,
    mean_g2 = mu2,
    lfc_vst = lfc,
    se = se,
    z = z,
    pval = p,
    padj = q,
    expressed_ok = expressed_ok,
    significant = sig,
    stringsAsFactors = FALSE
  )
  # Order by |z|
  res <- res[order(-abs(res$z), -abs(res$lfc_vst)), ]

  # ---- Save results
  base <- sprintf("cmp_%s_vs_%s", safe_fname(g1n), safe_fname(g2n))
  readr::write_tsv(res, file.path(out_dir, paste0(base, "_results.tsv")))

  # ---- PLOT 1: mean(Cond1) vs mean(Cond2) scatter
  # label top |z|
  lab_idx <- head(seq_len(nrow(res)), top_label)
  pdf(file.path(out_dir, paste0(base, "_scatter_mean_vs_mean.pdf")), width = 7.5, height = 6, onefile = FALSE)
  p <- ggplot(res, aes(x = mean_g2, y = mean_g1)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(aes(alpha = significant), size = 1) +
    scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.2)) +
    labs(
      title = sprintf("Mean(%s) vs Mean(%s) — VST", g1n, g2n),
      x = sprintf("Mean VST: %s", g2n),
      y = sprintf("Mean VST: %s", g1n)
    ) +
    theme_bw(base_size = 12) +
    ggrepel::geom_text_repel(
      data = res[lab_idx, ],
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3
    )
  print(p)
  dev.off()

  # ---- PLOT 2: z vs LFC (volcano-like)
  pdf(file.path(out_dir, paste0(base, "_volcano_z_vs_lfc.pdf")), width = 7.5, height = 6, onefile = FALSE)
  q <- ggplot(res, aes(x = lfc_vst, y = z, color = significant)) +
    geom_hline(yintercept = c(-min_abs_z, min_abs_z), linetype = 3) +
    geom_vline(xintercept = c(-min_abs_lfc, min_abs_lfc), linetype = 3) +
    geom_point(size = 1) +
    scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey70")) +
    labs(
      title = sprintf("z vs LFC (VST): %s vs %s", g1n, g2n),
      x = sprintf("LFC (VST): %s − %s", g1n, g2n),
      y = "z (Welch approx)"
    ) +
    theme_bw(base_size = 12) +
    ggrepel::geom_text_repel(
      data = res[lab_idx, ],
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      show.legend = FALSE
    )
  print(q)
  dev.off()

    # ---- PLOT 3: Heatmap of significant genes (both directions) ----
  max_heatmap_genes <- 200     # maximum genes to display
  min_sig_for_heat  <- 10      # minimum to trigger plotting
  base_height_per_gene <- 0.15 # inches per gene for height adjustment

  sig_idx <- which(res$significant)
  if (length(sig_idx) >= min_sig_for_heat) {
    sig_order <- order(-abs(res$z[sig_idx]), -abs(res$lfc_vst[sig_idx]))
    keep_genes <- unique(res$gene[sig_idx][sig_order])
    if (length(keep_genes) > max_heatmap_genes)
      keep_genes <- keep_genes[seq_len(max_heatmap_genes)]

    # extract matrix (genes × samples in these two groups)
    i1 <- which(cond_vec == g1n)
    i2 <- which(cond_vec == g2n)
    mat <- vsd_mat[rownames(vsd_mat) %in% keep_genes, c(i1, i2), drop = FALSE]

    # z-score by gene
    mat <- t(scale(t(mat)))
    rownames(mat) <- make.unique(rownames(mat))

    # annotation for columns
    ann_df <- data.frame(
      condition = factor(c(rep(g1n, length(i1)), rep(g2n, length(i2))),
                         levels = c(g1n, g2n))
    )
    rownames(ann_df) <- colnames(mat)
    pal <- RColorBrewer::brewer.pal(3, "Set2")[1:2]
    ann_cols <- list(condition = setNames(pal, c(g1n, g2n)))

    # ---- adaptive height ----
    n_genes <- nrow(mat)
    heat_h  <- max(6, min(20, n_genes * base_height_per_gene))  # clamp between 6" and 20"

    # ---- blue-white-red palette ----
    heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

    pdf(file.path(out_dir, paste0(base, "_heatmap_sig_genes.pdf")),
        width = 6, height = heat_h, onefile = FALSE)
    pheatmap::pheatmap(
      mat,
      color = heat_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_method = "complete",
      main = sprintf("Significant genes (|z|≥%.1f, |LFC|≥%.1f): %s vs %s",
                     min_abs_z, min_abs_lfc, g1n, g2n),
      annotation_col = ann_df,
      annotation_colors = ann_cols,
      fontsize_row = ifelse(n_genes > 50, 6, 8),
      fontsize_col = 10,
      border_color = NA
    )
    dev.off()

    message(sprintf("[%s vs %s] Heatmap plotted for %d genes (height %.1f inches).",
                    g1n, g2n, n_genes, heat_h))
  } else {
    message(sprintf("[%s vs %s] Not enough significant genes (%d < %d) to draw heatmap.",
                    g1n, g2n, length(sig_idx), min_sig_for_heat))
  }

  # ===================== DESeq2 DGEA + heatmap for this comparison =====================
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(readr)
})

# ---- tunables
alpha_dgea   <- 0.05   # FDR (padj) cutoff
min_dgea_lfc <- 0.20   # absolute log2FC threshold for "significant" genes
max_hm_genes <- 200    # cap number of rows in heatmap (most sig first)
min_sig_to_plot <- 5  # don't draw heatmap if fewer sig genes

# sub-objects for this comparison
dds_sub <- dds[, c(i1, i2)]
dds_sub$condition <- droplevels(dds_sub$condition)

# run DESeq2 (re-level not needed if using 'contrast' below, but harmless)
dds_sub$condition <- relevel(dds_sub$condition, ref = g2n)
dds_sub <- DESeq(dds_sub, quiet = TRUE)

# results: g1 vs g2 (log2FC = g1 - g2)
res_dq <- results(dds_sub,
                  contrast = c("condition", g1n, g2n),
                  alpha = alpha_dgea,
                  tidy = FALSE)
# add gene symbols if available from GTF mapping (gene_annot built earlier)
res_df <- as.data.frame(res_dq)
res_df$feature_id <- rownames(res_df)

if (exists("gene_annot") && !is.null(gene_annot) && all(c("gene_id","gene_name") %in% colnames(gene_annot))) {
  # strip versions to match
  ga <- gene_annot
  ga$gene_id <- sub("\\.\\d+$", "", ga$gene_id)
  res_df$gene_id_nov <- sub("\\.\\d+$", "", res_df$feature_id)
  sym <- gene_annot$gene_name[match(res_df$gene_id_nov, ga$gene_id)]
  res_df$gene_symbol <- ifelse(!is.na(sym) & sym != "", sym, res_df$feature_id)
} else {
  res_df$gene_symbol <- res_df$feature_id
}

# order by padj, then abs(log2FC)
o <- order(res_df$padj, -abs(res_df$log2FoldChange), na.last = TRUE)
res_df <- res_df[o, ]

# write all results
base <- sprintf("dgea_%s_vs_%s", safe_fname(g1n), safe_fname(g2n))
readr::write_tsv(res_df,
                 file.path(out_dir, paste0(base, "_DESeq2_results.tsv")))

# pick significant set
sig_idx <- which(!is.na(res_df$padj) &
                 res_df$padj <= alpha_dgea &
                 abs(res_df$log2FoldChange) >= min_dgea_lfc)

if (length(sig_idx) >= min_sig_to_plot) {

  # keep up to max_hm_genes most significant
  keep_rows <- sig_idx[seq_len(min(length(sig_idx), max_hm_genes))]
  sig_ids   <- res_df$feature_id[keep_rows]  # these match rownames in vsd assay

  # pull VST matrix from vsd (use original IDs to subset safely)
  vst_all <- SummarizedExperiment::assay(vsd)  # rownames = original feature IDs
  # subset to sig genes and samples in the two groups
  mat <- vst_all[rownames(vst_all) %in% sig_ids, c(i1, i2), drop = FALSE]

  # row labels: prefer symbols if we have them
  lab_map <- res_df[, c("feature_id", "gene_symbol")]
  rowlabs <- lab_map$gene_symbol[match(rownames(mat), lab_map$feature_id)]
  rowlabs[is.na(rowlabs) | rowlabs == ""] <- rownames(mat)
  rownames(mat) <- make.unique(rowlabs)

  # z-score per gene (row) for heatmap
  mat_z <- t(scale(t(mat)))

  # annotation for columns
  ann_df <- data.frame(
    condition = factor(c(rep(g1n, length(i1)), rep(g2n, length(i2))),
                       levels = c(g1n, g2n))
  )
  rownames(ann_df) <- colnames(mat_z)

  # colors (same blue-white-red as in your other heatmaps)
  heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

  # adaptive height
  base_height_per_gene <- 0.15
  heat_h  <- max(6, min(20, nrow(mat_z) * base_height_per_gene))

  # draw and save
  pdf(file.path(out_dir, paste0(base, "_sig_heatmap.pdf")),
      width = 6, height = heat_h, onefile = FALSE)
  pheatmap::pheatmap(
    mat_z,
    color = heat_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "complete",
    main = sprintf("Significant DE genes (padj≤%.02f, |LFC|≥%.2f): %s vs %s",
                   alpha_dgea, min_dgea_lfc, g1n, g2n),
    annotation_col = ann_df,
    fontsize_row = ifelse(nrow(mat_z) > 50, 6, 8),
    fontsize_col = 10,
    border_color = NA
  )
  dev.off()

  message(sprintf("[DESeq2 %s vs %s] Heatmap plotted for %d genes (height %.1f inches).",
                  g1n, g2n, nrow(mat_z), heat_h))
} else {
  message(sprintf("[DESeq2 %s vs %s] Not enough significant genes (%d < %d) to draw heatmap.",
                  g1n, g2n, length(sig_idx), min_sig_to_plot))
}
# ===================== end DESeq2 DGEA + heatmap block =====================



    # ===================== Enrichment: GO ORA (UP/DOWN) + GO GSEA =====================
  suppressPackageStartupMessages({ library(enrichplot) })
  
  suppressPackageStartupMessages({
  library(msigdbr)          # CRAN
  library(clusterProfiler)  # Bioconductor
  library(org.Hs.eg.db)     # Bioconductor
  library(ggplot2)
  library(readr)
})

  # --- Helper: detect ID type (SYMBOL vs ENSEMBL), and map → ENTREZ
  detect_id_type <- function(ids) {
    ids_nv <- sub("\\.\\d+$", "", ids)
    frac_ens <- mean(grepl("^ENSG\\d+$", ids_nv))
    if (is.na(frac_ens)) return("SYMBOL")
    if (frac_ens > 0.5) "ENSEMBL" else "SYMBOL"
  }

  map_to_entrez <- function(x_ids, fromType) {
    x_ids <- sub("\\.\\d+$", "", x_ids)
    suppressWarnings(
      bitr(x_ids, fromType = fromType, toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
  }

  id_type <- detect_id_type(res$gene)

  # ---------- ORA: UP and DOWN (GO:BP) ----------
  # Select sig UP/DOWN sets using your thresholds (already in `res$significant` + sign of LFC)
  up_genes   <- res$gene[res$significant & res$lfc_vst > 0]
  down_genes <- res$gene[res$significant & res$lfc_vst < 0]

  # Map to ENTREZ
  up_map   <- if (length(up_genes))   map_to_entrez(up_genes,   id_type) else NULL
  down_map <- if (length(down_genes)) map_to_entrez(down_genes, id_type) else NULL

  do_enrichGO <- function(entrez_vec, label) {
    if (is.null(entrez_vec) || nrow(entrez_vec) == 0) return(NULL)
    ego <- enrichGO(
      gene          = unique(entrez_vec$ENTREZID),
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
    # write table
    write_tsv(as.data.frame(ego), file.path(out_dir, paste0(base, "_GO_BP_ORA_", label, ".tsv")))
    # dotplot PDF
    p <- dotplot(ego, showCategory = 20, font.size = 10) +
      ggtitle(sprintf("GO:BP ORA (%s) — %s vs %s", label, g1n, g2n))
    ggsave(file.path(out_dir, paste0(base, "_GO_BP_ORA_", label, ".pdf")),
           plot = p, device = "pdf", width = 7.5, height = 6)
    invisible(ego)
  }

  ego_up   <- do_enrichGO(up_map,   "UP")
  ego_down <- do_enrichGO(down_map, "DOWN")

  # ---------- GSEA: GO:BP (dotplot only) ----------
  all_ids <- sub("\\.\\d+$", "", res$gene)
  all_map <- map_to_entrez(all_ids, id_type)

  if (!is.null(all_map) && nrow(all_map) > 0) {
    ranks_df <- merge(all_map, data.frame(gene = all_ids, lfc = res$lfc_vst),
                      by.x = id_type, by.y = "gene", all.x = TRUE)
    ranks_df <- ranks_df[!is.na(ranks_df$ENTREZID) & !is.na(ranks_df$lfc), ]
    ranks_df <- ranks_df[order(abs(ranks_df$lfc), decreasing = TRUE), ]
    ranks_df <- ranks_df[!duplicated(ranks_df$ENTREZID), ]
    ranks <- ranks_df$lfc; names(ranks) <- ranks_df$ENTREZID
    ranks <- sort(ranks, decreasing = TRUE)

    if (length(ranks) >= 10) {
      gsea_bp <- clusterProfiler::gseGO(
        geneList     = ranks,
        OrgDb        = org.Hs.eg.db,
        ont          = "BP",
        keyType      = "ENTREZID",
        pvalueCutoff = 0.05,
        verbose      = FALSE
      )

      # --- GO:BP GSEA dotplot (no geneID/core_enrichment check) ---
      if (!is.null(gsea_bp) && nrow(as.data.frame(gsea_bp)) > 0) {
        readr::write_tsv(as.data.frame(gsea_bp),
                        file.path(out_dir, paste0(base, "_GO_BP_GSEA.tsv")))
        p_dot <- enrichplot::dotplot(gsea_bp, showCategory = 20) +
          ggplot2::ggtitle(sprintf("GO:BP GSEA — %s vs %s (ranked by LFC)", g1n, g2n))

        n_terms <- nrow(as.data.frame(gsea_bp))  # or as.data.frame(ego), etc.
        dot_h <- max(4, min(0.25 * n_terms, 12))  # 0.25" per term, capped between 4" and 12"

        ggplot2::ggsave(file.path(out_dir, paste0(base, "_GO_BP_GSEA_dotplot.pdf")),
                        plot = p_dot, device = "pdf", width = 7.5, height = dot_h)
      } else {
        message(sprintf("[%s vs %s] GO:BP GSEA returned no significant terms.", g1n, g2n))
      }
    } else {
      message(sprintf("[%s vs %s] Not enough ranked genes for GSEA.", g1n, g2n))
    }
  } else {
    message(sprintf("[%s vs %s] Could not map genes to ENTREZ; skipping GO:BP GSEA.", g1n, g2n))
  }



    # ===================== MSigDB Enrichment: ORA (UP/DOWN) and GSEA =====================
  id_type <- detect_id_type(res$gene)

  # --- define UP/DOWN gene sets (significant by your thresholds)
  up_genes   <- res$gene[res$significant & res$lfc_vst > 0]
  down_genes <- res$gene[res$significant & res$lfc_vst < 0]

  up_map   <- if (length(up_genes))   map_to_entrez(up_genes,   id_type) else NULL
  down_map <- if (length(down_genes)) map_to_entrez(down_genes, id_type) else NULL

  # --- ranked list for GSEA (ENTREZ ~ LFC)
  all_ids <- sub("\\.\\d+$", "", res$gene)
  all_map <- map_to_entrez(all_ids, id_type)
  ranks <- NULL
  if (!is.null(all_map) && nrow(all_map) > 0) {
    ranks_df <- merge(all_map, data.frame(gene = all_ids, lfc = res$lfc_vst),
                      by.x = id_type, by.y = "gene", all.x = TRUE)
    ranks_df <- ranks_df[!is.na(ranks_df$ENTREZID) & !is.na(ranks_df$lfc), ]
    # collapse duplicate ENTREZ by max |lfc|
    ranks_df <- ranks_df[order(abs(ranks_df$lfc), decreasing = TRUE), ]
    ranks_df <- ranks_df[!duplicated(ranks_df$ENTREZID), ]
    ranks <- ranks_df$lfc
    names(ranks) <- ranks_df$ENTREZID
    ranks <- sort(ranks, decreasing = TRUE)
  }

  # ---------- ORA helper using MSigDB TERM2GENE ----------
    do_msig_ora <- function(t2g, label) {
      if (is.null(t2g) || nrow(t2g) == 0) return(invisible(NULL))
      # UP
      if (!is.null(up_map) && nrow(up_map) > 0) {
        ego_up <- enricher(gene = unique(up_map$ENTREZID), TERM2GENE = t2g,
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
        if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
          write_tsv(as.data.frame(ego_up), file.path(out_dir, paste0(base, "_", label, "_ORA_UP.tsv")))
          p <- dotplot(ego_up, showCategory = 20, font.size = 10) +
            ggtitle(sprintf("%s ORA (UP) — %s vs %s", label, g1n, g2n))
          ggsave(file.path(out_dir, paste0(base, "_", label, "_ORA_UP.pdf")), plot = p,
                device = "pdf", width = 7.5, height = 6)
        }
      }
      # DOWN
      if (!is.null(down_map) && nrow(down_map) > 0) {
        ego_dn <- enricher(gene = unique(down_map$ENTREZID), TERM2GENE = t2g,
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
        if (!is.null(ego_dn) && nrow(as.data.frame(ego_dn)) > 0) {
          write_tsv(as.data.frame(ego_dn), file.path(out_dir, paste0(base, "_", label, "_ORA_DOWN.tsv")))
          p <- dotplot(ego_dn, showCategory = 20, font.size = 10) +
            ggtitle(sprintf("%s ORA (DOWN) — %s vs %s", label, g1n, g2n))
          ggsave(file.path(out_dir, paste0(base, "_", label, "_ORA_DOWN.pdf")), plot = p,
                device = "pdf", width = 7.5, height = 6)
        }
      }
    }

    # ---------- GSEA helper using MSigDB TERM2GENE (dotplot only) ----------
  do_msig_gsea <- function(t2g, label) {
    if (is.null(t2g) || nrow(t2g) == 0) return(invisible(NULL))
    if (!is.null(ranks) && length(ranks) >= 10) {
      gsea_res <- clusterProfiler::GSEA(geneList = ranks, TERM2GENE = t2g,
                                        pvalueCutoff = 0.05, verbose = FALSE)
      if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
        readr::write_tsv(as.data.frame(gsea_res),
                        file.path(out_dir, paste0(base, "_", label, "_GSEA.tsv")))
        p_dot <- enrichplot::dotplot(gsea_res, showCategory = 20) +
          ggplot2::ggtitle(sprintf("%s GSEA — %s vs %s (ranked by LFC)", label, g1n, g2n))
        
        n_terms <- nrow(as.data.frame(gsea_res))  # or as.data.frame(ego), etc.
        dot_h <- max(4, min(0.25 * n_terms, 12))  # 0.25" per term, capped between 4" and 12"

        ggplot2::ggsave(file.path(out_dir, paste0(base, "_", label, "_GSEA_dotplot.pdf")),
                        plot = p_dot, device = "pdf", width = 7.5, height = dot_h)
      }

    }
    invisible(NULL)
  }


  # ---------- Run MSigDB analyses ----------
  # Hallmarks (H)
  t2g_H <- msig_t2g("H")
  do_msig_ora(t2g_H, "MSIG_Hallmark")
  do_msig_gsea(t2g_H, "MSIG_Hallmark")

  # Canonical Pathways (C2) — if subcategories missing, fall back to full C2
  t2g_CP_KEGG  <- msig_t2g("C2", subcategory = "CP:KEGG")
  t2g_CP_REACT <- msig_t2g("C2", subcategory = "CP:REACTOME")
  t2g_CP_WP    <- msig_t2g("C2", subcategory = "CP:WIKIPATHWAYS")

  # If all three are NULL, use entire C2 as a fallback
  if (is.null(t2g_CP_KEGG) && is.null(t2g_CP_REACT) && is.null(t2g_CP_WP)) {
    message("MSigDB C2 subcategories unavailable in this msigdbr build; using whole C2 as fallback.")
    t2g_C2 <- msig_t2g("C2")
    if (!is.null(t2g_C2)) {
      do_msig_ora(t2g_C2, "MSIG_C2")
      do_msig_gsea(t2g_C2, "MSIG_C2")
    }
  } else {
    if (!is.null(t2g_CP_KEGG))  { do_msig_ora(t2g_CP_KEGG,  "MSIG_CP_KEGG");        do_msig_gsea(t2g_CP_KEGG,  "MSIG_CP_KEGG") }
    if (!is.null(t2g_CP_REACT)) { do_msig_ora(t2g_CP_REACT, "MSIG_CP_REACTOME");    do_msig_gsea(t2g_CP_REACT, "MSIG_CP_REACTOME") }
    if (!is.null(t2g_CP_WP))    { do_msig_ora(t2g_CP_WP,    "MSIG_CP_WIKIPATHWAYS"); do_msig_gsea(t2g_CP_WP,    "MSIG_CP_WIKIPATHWAYS") }
  }


  invisible(res)
}

# ---- Run all requested comparisons
for (cmp in comparisons) {
  run_comparison(cmp[1], cmp[2])
}

message("Saved per-comparison TSVs and PDFs (scatter + volcano) in: ", normalizePath(out_dir, mustWork = FALSE))




# ===================== Distinct-by-condition (UP & DOWN) → s2s correlation heatmap =====================
suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(readr)
  library(dplyr)
})

# ----- tunable parameters
n_top_each   <- 200   # top genes per direction (UP & DOWN) for each condition
min_expr_any <- 2     # require either group-mean or others-mean >= this (screen out very low expression)
min_lfc_vst  <- 0.25   # absolute VST mean difference (≈ ≥2x on log2-ish scale)
min_z        <- 2.0   # absolute z threshold

stopifnot(all(colnames(vsd_mat) == rownames(coldata)))
cond <- droplevels(coldata$condition)
cond_levels <- levels(cond)
message("Conditions detected: ", paste(cond_levels, collapse = ", "))

distinct_for_condition <- function(cond_label, vsd_mat, cond) {
  in_grp  <- which(cond == cond_label)
  others  <- setdiff(seq_len(ncol(vsd_mat)), in_grp)
  stopifnot(length(in_grp) >= 2, length(others) >= 2)

  mu_grp   <- rowMeans(vsd_mat[, in_grp, drop = FALSE])
  mu_other <- rowMeans(vsd_mat[, others, drop = FALSE])
  sd_other <- apply(vsd_mat[, others, drop = FALSE], 1, sd)

  lfc_vst <- mu_grp - mu_other
  z_rest  <- lfc_vst / pmax(sd_other, 1e-6)

  # keep genes that are expressed in at least one side and have substantial |effect|
  expressed_ok <- (pmax(mu_grp, mu_other) >= min_expr_any)
  up_keep      <- expressed_ok & (lfc_vst >=  min_lfc_vst) & ( z_rest >=  min_z)
  down_keep    <- expressed_ok & (lfc_vst <= -min_lfc_vst) & ( z_rest <= -min_z)

  res <- data.frame(
    gene         = rownames(vsd_mat),
    mean_group   = mu_grp,
    mean_others  = mu_other,
    sd_others    = sd_other,
    lfc_vst      = lfc_vst,
    z            = z_rest,
    condition    = cond_label,
    keep_up      = up_keep,
    keep_down    = down_keep,
    stringsAsFactors = FALSE
  )

  # rank within direction
  res_up   <- res[res$keep_up, , drop = FALSE][order(-res$z[res$keep_up], -res$lfc_vst[res$keep_up]), , drop = FALSE]
  res_down <- res[res$keep_down, , drop = FALSE][order( res$z[res$keep_down],  res$lfc_vst[res$keep_down]), , drop = FALSE] # more negative first

  list(all = res,
       top_up = head(res_up, n_top_each),
       top_down = head(res_down, n_top_each))
}

per_cond <- lapply(cond_levels, distinct_for_condition, vsd_mat = vsd_mat, cond = cond)
names(per_cond) <- cond_levels

# write per-condition tables; collect union
genes_union <- character(0)
for (cl in cond_levels) {
  res_all  <- per_cond[[cl]]$all
  res_up   <- per_cond[[cl]]$top_up
  res_down <- per_cond[[cl]]$top_down

  readr::write_tsv(res_all,  file.path(out_dir, paste0("distinct_by_condition_ALL_",     make.names(cl), ".tsv")))
  readr::write_tsv(res_up,   file.path(out_dir, paste0("distinct_by_condition_TOP_UP_",  make.names(cl), ".tsv")))
  readr::write_tsv(res_down, file.path(out_dir, paste0("distinct_by_condition_TOP_DOWN_",make.names(cl), ".tsv")))

  genes_union <- unique(c(genes_union, res_up$gene, res_down$gene))
}

if (length(genes_union) < 5) {
  stop("Too few distinct genes across conditions (", length(genes_union),
       "). Relax thresholds or increase n_top_each.")
}
message("Distinct genes union size (UP+DOWN, per condition): ", length(genes_union))

# subset and z-score per gene
mat_sub <- vsd_mat[rownames(vsd_mat) %in% genes_union, , drop = FALSE]
# ensure rownames are unique (symbols can collide); required by pheatmap
rownames(mat_sub) <- make.unique(rownames(mat_sub))
mat_sub <- t(scale(t(mat_sub)))

# sample↔sample correlation using only distinct genes
sample_cor_sel <- cor(mat_sub, method = "pearson", use = "pairwise.complete.obs")

# annotation
annotation_df <- data.frame(condition = cond, row.names = rownames(coldata))
condition_levels <- levels(annotation_df$condition)
pal <- brewer.pal(min(length(condition_levels), 8), "Set2")
if (length(condition_levels) > 8) pal <- colorRampPalette(pal)(length(condition_levels))
annotation_colors <- list(condition = setNames(pal, condition_levels))

# plot to PDF
pdf(file.path(out_dir, "sample_correlation_heatmap_distinct_by_condition_UP_AND_DOWN.pdf"),
    width = 9, height = 8, onefile = FALSE)
pheatmap::pheatmap(
  sample_cor_sel,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = paste0("Sample correlation using distinct genes per condition (UP & DOWN, top ", n_top_each, " each)"),
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  annotation_colors = annotation_colors,
  fontsize = 10,
  border_color = NA
)
dev.off()

# save union gene list
writeLines(genes_union, file.path(out_dir, paste0("distinct_genes_union_UP_AND_DOWN_top", n_top_each, "_per_condition.txt")))
message("Saved: sample_correlation_heatmap_distinct_by_condition_UP_AND_DOWN.pdf and per-condition tables.")




# ---------------------- PCA ----------------------
pca <- prcomp(t(mat_sub), scale. = FALSE) # using mat_sub instead of vsd_mat is also an option
pca_df <- as.data.frame(pca$x) |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::left_join(sample_table |>
                     dplyr::select(sample_id, condition),
                   by = "sample_id")

percentVar <- (pca$sdev^2) / sum(pca$sdev^2)

p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, label = sample_id, color = condition)) +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE) +
  ggplot2::labs(
    title = "PCA on VST-normalized counts",
    x = paste0("PC1 (", scales::percent(percentVar[1]), ")"),
    y = paste0("PC2 (", scales::percent(percentVar[2]), ")")
  ) +
  ggplot2::theme_bw(base_size = 12)

ggplot2::ggsave(file.path(out_dir, "pca_vst.pdf"), p_pca, width = 7.5, height = 6, device = "pdf")

# ---------------------- Export normalized counts (VST) ----------------------
vst_tbl <- as.data.frame(vsd_mat) |>
  tibble::rownames_to_column("feature_id")
readr::write_tsv(vst_tbl, file.path(out_dir, "vst_normalized_counts.tsv"))



# ---------------------- Session info & README ----------------------
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

readme <- c(
  "# RNA-seq downstreams (Salmon -> DESeq2)",
  "",
  "Artifacts:",
  "- sample_table.csv — edit `condition` if needed, then rerun for your design.",
  "- counts_filtered.rds / .tsv — CPM-filtered count matrix.",
  "- dds_prepped_for_deseq2.rds — DESeqDataSet ready for DESeq2::DESeq().",
  "- vsd.rds — VST-transformed DESeqDataSet for QC.",
  "- qc_library_sizes.pdf / qc_detected_features.pdf — basic QC plots.",
  "- sample_correlation_heatmap.pdf — Pearson correlation on VST.",
  "- pca_vst.pdf — PCA on VST.",
  "- vst_normalized_counts.tsv — VST expression matrix.",
  "",
  "Notes:",
  "- If you provide a tx2gene mapping (2 columns: transcript,gene), the pipeline runs at gene level using tximport with `lengthScaledTPM`.", # nolint: line_length_linter.
  "- Without tx2gene, analysis is at transcript level; gene-level DE is recommended.",
  "- Default design is `~ condition`; edit sample_table.csv (and/or modify script) for more complex designs.",
  sprintf("- Filtering keeps features with CPM >= %.2f in at least %d samples.", args$min_cpm, min_samps)
)
writeLines(readme, con = file.path(out_dir, "README.txt"))

message("Done. Outputs written to: ", normalizePath(out_dir, mustWork = FALSE))


# This is how to run this script from command line:

# '''
# Rscript scripts/rna_downstreams_salmon.R \
#   --sample_table results/sample_table.csv \
#   --tx2gene data/annotation.tx2gene_proteincoding.tsv \
#   --out_dir results/rna_downstreams_out
# '''