#!/bin/bash
# Reliability task: per-gene read-evidence audit + reviewer-requested supplementary table + coverage figure.
# Runs end-to-end from a clean shell (tested with: env -i HOME=$HOME PATH=/usr/bin:/bin bash 13_reliability_00_run_all.sh).
# Inputs are read-only; outputs -> /mnt/ssd1/Projects/NKX3_wd/results/revision_2026/reliability/ ; temp -> work/reliability/
set -euo pipefail
S=/mnt/ssd1/Projects/NKX3_wd/scripts/revision
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
R=/home/biolab/miniconda3/envs/nkx3_r/bin/Rscript
TMP=/mnt/ssd1/Projects/NKX3_wd/work/reliability
mkdir -p "$TMP"; cd "$TMP"
$PY -I $S/13_reliability_01_gtf_models.py        # Ensembl 99 gene/exon models (temp pickle)
$PY -I $S/13_reliability_02_bam_audit.py         # BAM audit (<=4 concurrent readers), Salmon/featureCounts comparison, coverage arrays
$R      $S/13_reliability_03_table.R             # supplementary table columns, classification, family DESeq2 tests
$PY -I $S/13_reliability_04_excel.py             # Excel workbook + flat TSVs
$PY -I $S/13_reliability_05_coverage_figure.py   # coverage figure (composite + single panels)
