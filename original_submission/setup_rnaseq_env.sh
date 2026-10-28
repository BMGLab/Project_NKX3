#!/usr/bin/env bash
# =============================================================
# RNA-seq environment setup script for downstreams analysis
# Author: Yasin Kaymaz (Project NKX3)
# =============================================================

set -e  # Exit on error
set -o pipefail

ENV_NAME="rnaseq_env"
R_VERSION="4.4.1"

echo "==============================================="
echo " RNA-seq downstreams environment setup"
echo " Conda environment: $ENV_NAME"
echo " R version: $R_VERSION"
echo "==============================================="

# -------------------------------------------------------------
# 1. Check for conda or mamba
# -------------------------------------------------------------
if ! command -v conda &> /dev/null; then
  echo "❌ Conda not found. Please install Miniconda or Anaconda first."
  exit 1
fi

# Use mamba if available, else install it
if ! command -v mamba &> /dev/null; then
  echo "⚙️  Installing mamba into base environment..."
  conda install -n base -c conda-forge mamba -y
fi

# -------------------------------------------------------------
# 2. Create new environment with R and key tools
# -------------------------------------------------------------
echo "🚀 Creating environment '$ENV_NAME'..."
mamba create -y -n $ENV_NAME \
  -c conda-forge -c bioconda \
  r-base=$R_VERSION \
  r-essentials \
  r-tidyverse \
  r-data.table \
  r-optparse \
  r-stringr \
  r-tidyr \
  r-readr \
  r-dplyr \
  r-ggplot2 \
  r-pheatmap \
  r-tximport \
  bioconductor-deseq2 \
  bioconductor-edger \
  bioconductor-tximeta \
  bioconductor-biocparallel \
  bioconductor-genomicranges \
  bioconductor-summarizedexperiment \
  r-ggrepel \
  r-scales

echo "✅ Conda environment '$ENV_NAME' created."

# -------------------------------------------------------------
# 3. Activate environment
# -------------------------------------------------------------
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# -------------------------------------------------------------
# 4. Install any remaining R packages from CRAN/Bioconductor
# -------------------------------------------------------------
echo "📦 Installing any missing R packages via Rscript..."
Rscript - <<'EOF'
install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing ", p)
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

bioc_install_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing Bioconductor pkg: ", p)
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }
  }
}

cran_pkgs <- c("optparse", "readr", "dplyr", "tidyr", "stringr", 
               "ggplot2", "pheatmap", "ggrepel", "scales")

bioc_pkgs <- c("DESeq2", "edgeR", "tximport", "tximeta", 
               "SummarizedExperiment", "BiocParallel", "GenomicRanges")

install_if_missing(cran_pkgs)
bioc_install_if_missing(bioc_pkgs)

message("✅ All required R packages are installed.")
EOF

# -------------------------------------------------------------
# 5. Environment summary
# -------------------------------------------------------------
echo "==============================================="
echo "✅ Environment setup complete!"
echo "To activate it in future sessions, run:"
echo "   conda activate $ENV_NAME"
echo
echo "Then you can run:"
echo "   Rscript rna_downstreams_salmon.R --data_dir <path> --out_dir ./rna_downstreams_out"
echo "==============================================="
