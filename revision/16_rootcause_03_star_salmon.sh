#!/usr/bin/env bash
# 16_rootcause_03_star_salmon.sh
# Re-run the nf-core/rnaseq 3.14.0 STAR_ALIGN + QUANTIFY_STAR_SALMON:SALMON_QUANT
# commands for one sample.
#
# usage: 16_rootcause_03_star_salmon.sh <sample_id> <mode> <fastq_source>
#   mode = repro     : STAR 2.7.9a, options identical to nf-core (--quantTranscriptomeBan Singleend)
#   mode = banISS    : STAR 2.7.9a, identical except --quantTranscriptomeBan IndelSoftclipSingleend
#                      (STAR's default; RSEM-compatible transcriptomic output)
#   mode = star2710b : STAR 2.7.10b, options identical to nf-core (Singleend)
#   fastq_source = rawtrim : raw FASTQ trimmed as Trim Galore did (16_rootcause_01b_trim_raw.sh)
#                  bamfq   : FASTQ regenerated from the nf-core genome BAM (16_rootcause_01_fastq.sh)
# The pre-built STAR index (STAR 2.7.4a, unfiltered Ensembl 99 GTF) is used
# read-only; on-the-fly sjdb insertion files go to the run's own _STARgenome dir.
set -euo pipefail
SAMPLE=${1:-24-004-7_S14}
MODE=${2:-repro}
FQSRC=${3:-rawtrim}
THREADS=${STAR_THREADS:-16}   # thread count does not change alignments (only output order)
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
SALMON=/home/biolab/miniconda3/envs/nkx3_salmon/bin/salmon
IDX=/home/sharedFolder/humanSTARindex
FGTF=$WORK/ref/Homo_sapiens.GRCh38.dna.primary_assembly.filtered.gtf
TXFA=$WORK/ref/rsem/genome.transcripts.fa
case $FQSRC in
  rawtrim) R1=$WORK/fastq_trim/$SAMPLE/${SAMPLE}_1_val_1.fq.gz; R2=$WORK/fastq_trim/$SAMPLE/${SAMPLE}_2_val_2.fq.gz ;;
  bamfq)   R1=$WORK/fastq/$SAMPLE/${SAMPLE}_R1.fq.gz;          R2=$WORK/fastq/$SAMPLE/${SAMPLE}_R2.fq.gz ;;
  *) echo "unknown fastq source $FQSRC"; exit 1 ;;
esac
case $MODE in
  repro)     BAN=Singleend;              STAR=/home/biolab/miniconda3/envs/nkx3_star/bin/STAR ;;
  banISS)    BAN=IndelSoftclipSingleend; STAR=/home/biolab/miniconda3/envs/nkx3_star/bin/STAR ;;
  star2710b) BAN=Singleend;              STAR=/home/biolab/miniconda3/envs/bench_star/bin/STAR ;;
  *) echo "unknown mode $MODE"; exit 1 ;;
esac
TAG=${MODE}_${FQSRC}
SD=$WORK/star_${TAG}/$SAMPLE
QD=$WORK/salmon_${TAG}/$SAMPLE
RD=$OUT/runs/${TAG}/$SAMPLE
mkdir -p "$SD" "$(dirname "$QD")" "$RD"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

# GENOME_BAM=keep (default) writes Aligned.out.bam; GENOME_BAM=names streams the
# identical genome BAM to stdout (--outStd BAM_Unsorted) and keeps only the sorted
# read names (saves ~6 GB of scratch per run; alignments are unaffected).
GENOME_BAM=${GENOME_BAM:-keep}
STD_OPT=(); [ "$GENOME_BAM" = names ] && STD_OPT=(--outStd BAM_Unsorted)
if [ ! -s "$SD/${SAMPLE}.Log.final.out" ]; then
  echo "[$(date)] STAR $($STAR --version) ($TAG, --quantTranscriptomeBan $BAN, genome BAM: $GENOME_BAM)"
  cd "$SD"
  run_star() {
  $STAR --genomeDir $IDX \
    --readFilesIn "$R1" "$R2" \
    --runThreadN $THREADS \
    --outFileNamePrefix "$SD/${SAMPLE}." \
    --sjdbGTFfile "$FGTF" \
    --outSAMattrRGline ID:$SAMPLE SM:$SAMPLE \
    --quantMode TranscriptomeSAM \
    --twopassMode Basic \
    --outSAMtype BAM Unsorted \
    --readFilesCommand zcat \
    --runRNGseed 0 \
    --outFilterMultimapNmax 20 \
    --alignSJDBoverhangMin 1 \
    --outSAMattributes NH HI AS NM MD \
    --quantTranscriptomeBan $BAN \
    --outSAMstrandField intronMotif "${STD_OPT[@]}"
  }
  if [ "$GENOME_BAM" = names ]; then
    run_star | samtools view -@4 - | cut -f1 | uniq | sort -u -S 8G --parallel=8 -T "$TMPDIR" > "$SD/genome_names.txt"
  else
    run_star
  fi
fi
cp "$SD/${SAMPLE}.Log.final.out" "$SD/${SAMPLE}.Log.out" "$RD/"

if [ ! -s "$QD/quant.sf" ]; then
  echo "[$(date)] Salmon alignment mode ($TAG)"
  # identical to nf-core 3.14.0 QUANTIFY_STAR_SALMON:SALMON_QUANT (cmd_info.json)
  $SALMON quant \
    --geneMap "$FGTF" \
    --threads 6 \
    --libType=ISR \
    -t "$TXFA" \
    -a "$SD/${SAMPLE}.Aligned.toTranscriptome.out.bam" \
    -o "$QD"
fi
cp "$QD/quant.sf" "$QD/cmd_info.json" "$QD/aux_info/meta_info.json" "$QD/logs/salmon_quant.log" "$RD/"
gzip -f "$RD/quant.sf"

# read-level bookkeeping: distinct read names in each BAM
if [ ! -s "$RD/read_name_counts.tsv" ]; then
  if [ ! -s "$SD/genome_names.txt" ]; then
    samtools view -@4 "$SD/${SAMPLE}.Aligned.out.bam" | cut -f1 | uniq | sort -u -S 8G --parallel=8 -T "$TMPDIR" > "$SD/genome_names.txt"
  fi
  G=$(wc -l < "$SD/genome_names.txt")
  T=$(samtools view -@4 "$SD/${SAMPLE}.Aligned.toTranscriptome.out.bam" | cut -f1 | uniq | sort -u -S 8G --parallel=8 -T "$TMPDIR" | tee "$SD/tx_names.txt" | wc -l)
  TO=$(comm -13 "$SD/genome_names.txt" "$SD/tx_names.txt" | tee "$SD/tx_only_names.txt" | wc -l)
  printf "metric\tvalue\nreads_in_genome_BAM\t%s\nreads_in_toTranscriptome_BAM\t%s\nreads_in_toTranscriptome_but_absent_from_genome_BAM\t%s\n" "$G" "$T" "$TO" > "$RD/read_name_counts.tsv"
  cat "$RD/read_name_counts.tsv"
fi
echo "[$(date)] done $SAMPLE $TAG"
