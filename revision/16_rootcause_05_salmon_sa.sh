#!/usr/bin/env bash
# 16_rootcause_05_salmon_sa.sh
# Fix test (2): Salmon 1.10.1 selective-alignment (mapping) mode, no STAR.
# Index built on the nf-core-equivalent genome.transcripts.fa (RSEM extraction
# of the filtered Ensembl 99 GTF), k = 31, no decoys.
# Input: raw FASTQ trimmed as Trim Galore did (16_rootcause_01b_trim_raw.sh),
# i.e. all reads, including those STAR leaves unmapped.
#
# usage: 16_rootcause_05_salmon_sa.sh <sample_id>
set -euo pipefail
SAMPLE=${1:-24-004-7_S14}
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
SALMON=/home/biolab/miniconda3/envs/nkx3_salmon/bin/salmon
TXFA=$WORK/ref/rsem/genome.transcripts.fa
FGTF=$WORK/ref/Homo_sapiens.GRCh38.dna.primary_assembly.filtered.gtf
IDX=$WORK/ref/salmon_idx_k31
FQ=$WORK/fastq_trim/$SAMPLE
QD=$WORK/salmon_SA_rawtrim/$SAMPLE
RD=$OUT/runs/salmonSA_rawtrim/$SAMPLE
mkdir -p "$(dirname "$QD")" "$RD" "$OUT/logs"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

if [ ! -s "$IDX/versionInfo.json" ]; then
  echo "[$(date)] salmon index (k31, no decoys)"
  $SALMON index -t "$TXFA" -i "$IDX" -k 31 -p 16 > "$OUT/logs/05_salmon_index.log" 2>&1
fi
if [ ! -s "$QD/quant.sf" ]; then
  echo "[$(date)] salmon quant (selective alignment) $SAMPLE"
  $SALMON quant -i "$IDX" --libType=ISR --geneMap "$FGTF" --threads ${SALMON_THREADS:-16} \
     -1 "$FQ/${SAMPLE}_1_val_1.fq.gz" -2 "$FQ/${SAMPLE}_2_val_2.fq.gz" -o "$QD"
fi
cp "$QD/quant.sf" "$QD/cmd_info.json" "$QD/aux_info/meta_info.json" "$QD/logs/salmon_quant.log" "$RD/"
gzip -f "$RD/quant.sf"
echo "[$(date)] done"
