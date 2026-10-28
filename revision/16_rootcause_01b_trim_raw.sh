#!/usr/bin/env bash
# 16_rootcause_01b_trim_raw.sh
# Re-create the nf-core 3.14.0 TRIMGALORE output (*_val_1/2.fq.gz) from the raw
# FASTQ, so that STAR sees the same input as in the original run, INCLUDING
# the reads that STAR leaves unmapped (absent from the genome BAM).
# Trim Galore 0.6.7 (cutadapt 3.4) is not installed here; we reproduce what its
# trimming report shows it did, with cutadapt 5.2:
#   per read file, single-end: cutadapt -e 0.1 -q 20 -O 1 -a CTGTCTCTTATA
#   (Nextera adapter auto-detected by Trim Galore), then pair validation:
#   drop the pair if either read < 20 nt.
# usage: 16_rootcause_01b_trim_raw.sh <sample_id>
set -euo pipefail
SAMPLE=${1:-24-004-7_S14}
RAW=/mnt/lun2/rawdata/BilgeDebelec/BAP-RNAseq
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
CUTADAPT=/home/biolab/miniconda3/envs/mdsc-trim/bin/cutadapt
FQ=$WORK/fastq_trim/$SAMPLE
mkdir -p "$FQ" "$OUT/logs"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"
ADAPTER=CTGTCTCTTATA

for m in 1 2; do
  if [ ! -s "$FQ/${SAMPLE}_${m}_trimmed.fq.gz" ] && [ ! -s "$FQ/${SAMPLE}_${m}_val_${m}.fq.gz" ]; then
    echo "[$(date)] cutadapt R$m"
    $CUTADAPT -j ${CUTADAPT_CORES:-16} -e 0.1 -q 20 -O 1 -a $ADAPTER --compression-level 1 \
      -o "$FQ/${SAMPLE}_${m}_trimmed.fq.gz" "$RAW/${SAMPLE}_R${m}_001.fastq.gz" \
      > "$OUT/logs/01b_cutadapt_${SAMPLE}_R${m}.log"
  fi
done
if [ ! -s "$FQ/${SAMPLE}_2_val_2.fq.gz" ]; then
  echo "[$(date)] pair validation (length >= 20 for both mates)"
  $CUTADAPT -j ${CUTADAPT_CORES:-16} -m 20 --pair-filter=any --compression-level 1 \
    -o "$FQ/${SAMPLE}_1_val_1.fq.gz" -p "$FQ/${SAMPLE}_2_val_2.fq.gz" \
    "$FQ/${SAMPLE}_1_trimmed.fq.gz" "$FQ/${SAMPLE}_2_trimmed.fq.gz" \
    > "$OUT/logs/01b_validate_${SAMPLE}.log"
  rm -f "$FQ/${SAMPLE}_1_trimmed.fq.gz" "$FQ/${SAMPLE}_2_trimmed.fq.gz"
fi
grep -E "Total reads processed|Reads with adapters|Reads written|Total read pairs processed|Pairs that were too short|Pairs written" \
  "$OUT/logs/01b_cutadapt_${SAMPLE}_R1.log" "$OUT/logs/01b_cutadapt_${SAMPLE}_R2.log" "$OUT/logs/01b_validate_${SAMPLE}.log" || true
echo "[$(date)] done"
