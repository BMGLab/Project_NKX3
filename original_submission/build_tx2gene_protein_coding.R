#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
gtf <- if (length(args) >= 1) args[1] else "data/Homo_sapiens.GRCh38.99.gtf"
out <- if (length(args) >= 2) args[2] else "results/annotation.tx2gene_proteincoding.tsv"

message("GTF: ", gtf)
message("OUT: ", out)
if (!file.exists(gtf)) stop("GTF file not found: ", gtf)

message("Importing GTF (this can take a minute)...")
gr <- rtracklayer::import(gtf)

# ----- Select only transcript-level rows (or mRNA)
tx_idx <- which(mcols(gr)$type %in% c("transcript", "mRNA"))
if (length(tx_idx) == 0L) stop("No 'transcript' or 'mRNA' rows found in GTF.")
meta <- as.data.frame(mcols(gr)[tx_idx, , drop = FALSE])

avail <- colnames(meta)
message("Available attribute columns in GTF: ", paste(avail, collapse = ", "))

# Column name selection (first match wins)
pick <- function(cands, avail) {
  hit <- intersect(cands, avail)
  if (length(hit) == 0L) NA_character_ else hit[1]
}
tx_col  <- pick(c("transcript_id","transcript_id_version","transcript_name"), avail)
gene_col<- pick(c("gene_id","gene_name"), avail)
bio_col <- pick(c("transcript_biotype","transcript_type","gene_biotype","gene_type"), avail)

message("Detected columns -> ",
        "transcript: ", ifelse(is.na(tx_col),"NONE",tx_col), ", ",
        "gene: ",       ifelse(is.na(gene_col),"NONE",gene_col), ", ",
        "biotype: ",    ifelse(is.na(bio_col),"NONE",bio_col))

if (is.na(tx_col) || is.na(gene_col)) {
  stop("Missing transcript or gene column in GTF. ",
       "Looked for transcript in: transcript_id/_version/_name; ",
       "gene in: gene_id/gene_name.")
}
if (is.na(bio_col)) {
  stop("Missing biotype column. Tried transcript_biotype/transcript_type/gene_biotype/gene_type.")
}

# Build raw vectors
transcript <- as.character(meta[[tx_col]])
gene       <- as.character(meta[[gene_col]])
biotype    <- as.character(meta[[bio_col]])

# Drop NAs (vectorized; no ifs)
ok <- !is.na(transcript) & !is.na(gene) & !is.na(biotype)
transcript <- transcript[ok]
gene       <- gene[ok]
biotype    <- biotype[ok]

# Protein-coding filter
pc <- biotype %in% c("protein_coding", "protein-coding")
transcript <- transcript[pc]
gene       <- gene[pc]

if (length(transcript) == 0L) {
  uniq_bio <- unique(biotype)
  stop("No protein-coding entries found. Biotype examples: ",
       paste(head(uniq_bio, 20), collapse = ", "))
}

# Strip Ensembl version suffixes
strip <- function(x) sub("\\.\\d+$","", x)
transcript <- strip(transcript)
gene       <- strip(gene)

# Unique transcript→gene pairs; write out
df <- unique(data.frame(transcript = transcript, gene = gene, stringsAsFactors = FALSE))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
# Option B: explicit logical flag
readr::write_tsv(df, out, col_names = TRUE)
message("Wrote: ", out, "  (", nrow(df), " transcript→gene rows; protein-coding only)")
