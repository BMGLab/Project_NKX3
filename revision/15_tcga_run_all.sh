#!/usr/bin/env bash
# 15_tcga_run_all.sh -- TCGA-PRAD + public NKX3.1 perturbation data check (Reviewer 1, major 5), end to end.
# Downloads (public data) go to work/; results to results/revision_2026/tcga/.
#   bash /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_run_all.sh
set -euo pipefail
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
RS=/home/biolab/miniconda3/envs/nkx3_r/bin/Rscript
S=/mnt/ssd1/Projects/NKX3_wd/scripts/revision
LOG=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/tcga/logs
mkdir -p "$LOG"
"$PY" -I "$S/15_tcga_geo.py"        > "$LOG/15_tcga_geo.log" 2>&1          # GEO search + RWPE1 matrices
"$RS" "$S/15_tcga.R"                > "$LOG/15_tcga.log" 2>&1              # TCGA-PRAD analyses 1-4
"$RS" "$S/15_tcga_geo_rwpe1.R"      > "$LOG/15_tcga_geo_rwpe1.log" 2>&1    # RWPE1 limma (needs TCGA table)
"$RS" "$S/15_tcga_figures.R"        > "$LOG/15_tcga_figures.log" 2>&1      # figures A-F + combined pages
"$PY" -I "$S/15_tcga_tables.py"     > "$LOG/15_tcga_tables.log" 2>&1       # workbook + provenance txt
echo "done: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
