#!/usr/bin/env bash
# 16_rootcause_06_filter_fix.sh
# Diagnostic fix (3): remove from STAR's Aligned.toTranscriptome.out.bam every
# read that is absent from STAR's own genome BAM (i.e. reads STAR reported as
# unmapped, e.g. "too short"), then re-run Salmon with the nf-core options.
# If the artefact disappears, the inflated transcripts are fed exclusively by
# reads STAR itself rejected.
# usage: 16_rootcause_06_filter_fix.sh <sample_id> <run_tag>   (run_tag e.g. repro_rawtrim)
set -euo pipefail
SAMPLE=$1; TAG=$2
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
SAMTOOLS=/usr/local/bin/samtools
SALMON=/home/biolab/miniconda3/envs/nkx3_salmon/bin/salmon
FGTF=$WORK/ref/Homo_sapiens.GRCh38.dna.primary_assembly.filtered.gtf
TXFA=$WORK/ref/rsem/genome.transcripts.fa
SD=$WORK/star_${TAG}/$SAMPLE
FD=$WORK/filtered_${TAG}/$SAMPLE
QD=$WORK/salmon_${TAG}_genomeMappedOnly/$SAMPLE
RD=$OUT/runs/${TAG}_genomeMappedOnly/$SAMPLE
mkdir -p "$FD" "$(dirname "$QD")" "$RD"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

if [ ! -s "$QD/quant.sf" ]; then
  echo "[$(date)] filter toTranscriptome BAM to reads present in the genome BAM"
  $SAMTOOLS view -@ 8 -b -N "$SD/genome_names.txt" -o "$FD/tx_genome_mapped_only.bam" "$SD/${SAMPLE}.Aligned.toTranscriptome.out.bam"
  echo "[$(date)] Salmon alignment mode (nf-core options)"
  $SALMON quant --geneMap "$FGTF" --threads 6 --libType=ISR -t "$TXFA" \
     -a "$FD/tx_genome_mapped_only.bam" -o "$QD"
  rm -f "$FD/tx_genome_mapped_only.bam"
fi
cp "$QD/quant.sf" "$QD/cmd_info.json" "$QD/aux_info/meta_info.json" "$QD/logs/salmon_quant.log" "$RD/"
gzip -f "$RD/quant.sf"
echo "[$(date)] done"
