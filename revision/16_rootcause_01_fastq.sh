#!/usr/bin/env bash
# 16_rootcause_01_fastq.sh
# Regenerate paired FASTQ for one sample from its nf-core STAR genome BAM
# (primary alignments only). Reads that STAR left unmapped (~8% "too short"
# + "other") are NOT in the genome BAM and are therefore absent here; they did
# not reach the nf-core toTranscriptome BAM either, so they cannot affect Salmon.
#
# usage: 16_rootcause_01_fastq.sh <sample_id>      (default 24-004-7_S14)
set -euo pipefail
SAMPLE=${1:-24-004-7_S14}
SRC=/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon/${SAMPLE}.markdup.sorted.bam
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
SAMTOOLS=/usr/local/bin/samtools
FQ=$WORK/fastq/$SAMPLE
mkdir -p "$FQ" "$OUT/logs"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

echo "[$(date)] collate + fastq from $SRC"
# Non-fast collate: groups mates by name in hash order, i.e. a pseudo-random
# order relative to genome coordinates (like the original FASTQ order), so
# Salmon's online/burn-in phase sees reads in a non-positional order.
$SAMTOOLS collate -@ 8 -l 1 -n 128 -O -u -T "$TMPDIR/${SAMPLE}.collate" "$SRC" \
 | $SAMTOOLS fastq -@ 8 -n -F 0x900 -c 1 \
     -1 "$FQ/${SAMPLE}_R1.fq.gz" -2 "$FQ/${SAMPLE}_R2.fq.gz" \
     -s "$FQ/${SAMPLE}_singletons.fq.gz" -0 "$FQ/${SAMPLE}_other.fq.gz" - \
     2> "$OUT/logs/01_fastq_${SAMPLE}.samtools_fastq.log"
cat "$OUT/logs/01_fastq_${SAMPLE}.samtools_fastq.log"

# record counts
{
  echo -e "file\treads"
  for f in "$FQ"/${SAMPLE}_*.fq.gz; do
    n=$(zcat "$f" | awk 'END{print NR/4}')
    echo -e "$(basename "$f")\t$n"
  done
} | tee "$OUT/logs/01_fastq_${SAMPLE}.counts.tsv"
echo "[$(date)] done"
