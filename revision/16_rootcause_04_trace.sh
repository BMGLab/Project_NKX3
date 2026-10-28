#!/usr/bin/env bash
# 16_rootcause_04_trace.sh
# Trace the reads that STAR's Aligned.toTranscriptome.out.bam assigns to
# flagged transcripts, and profile ALL toTranscriptome reads that are absent
# from the genome BAM (i.e. reads STAR reported as unmapped).
#   (1) every transcriptomic record of reads touching a target transcript
#   (2) the same reads in the genome BAM (present / absent)
#   (3) their FASTQ records (full mate sequences)
#   (4) python summary: alignment length, soft clips, NH, implied genomic locus,
#       plasmid-derived sequence motifs
#   (5) global per-gene profile: reads in toTranscriptome BAM vs reads that are
#       absent from the genome BAM
# usage: 16_rootcause_04_trace.sh <sample_id> <run_tag> <comma-separated transcript IDs>
set -euo pipefail
SAMPLE=$1; TAG=$2; TARGETS=$3
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
SAMTOOLS=/usr/local/bin/samtools
FGTF=$WORK/ref/Homo_sapiens.GRCh38.dna.primary_assembly.filtered.gtf
SD=$WORK/star_${TAG}/$SAMPLE
TD=$WORK/trace_${TAG}/$SAMPLE
RD=$OUT/trace/${SAMPLE}_${TAG}
FQ1=$WORK/fastq_trim/$SAMPLE/${SAMPLE}_1_val_1.fq.gz
FQ2=$WORK/fastq_trim/$SAMPLE/${SAMPLE}_2_val_2.fq.gz
TXBAM=$SD/${SAMPLE}.Aligned.toTranscriptome.out.bam
GBAM=$SD/${SAMPLE}.Aligned.out.bam
mkdir -p "$TD" "$RD"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

# transcript -> gene table (once)
TX2G=$WORK/ref/tx2gene_name.tsv
if [ ! -s "$TX2G" ]; then
  awk -F'\t' '$3=="transcript"{match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16);
     match($9,/gene_id "[^"]+"/); g=substr($9,RSTART+9,RLENGTH-10);
     n=g; if (match($9,/gene_name "[^"]+"/)) n=substr($9,RSTART+11,RLENGTH-12);
     print t"\t"g"\t"n}' "$FGTF" > "$TX2G"
fi

echo "[$(date)] (1) transcriptomic groups touching targets"
if [ ! -s "$TD/tx_groups.sam" ]; then
  $SAMTOOLS view -@ 8 "$TXBAM" | awk -F'\t' -v T="$TARGETS" '
    BEGIN{n=split(T,a,","); for(i=1;i<=n;i++) t[a[i]]=1}
    {if($1!=q){ if(hit) printf "%s", buf; buf=""; hit=0; q=$1} buf=buf $0 "\n"; if($3 in t) hit=1}
    END{if(hit) printf "%s", buf}' > "$TD/tx_groups.sam"
fi
cut -f1 "$TD/tx_groups.sam" | sort -u > "$TD/names.txt"
wc -l < "$TD/names.txt" | awk '{print "reads touching targets:",$1}'

echo "[$(date)] (2) genome BAM records of these reads"
# presence in the genome BAM is decided from tx_only_names.txt (reads present in the
# toTranscriptome BAM but absent from the genome BAM; written by 16_rootcause_03);
# genome records are extracted only if the genome BAM was kept.
comm -12 "$TD/names.txt" "$SD/tx_only_names.txt" > "$TD/names_absent_from_genome.txt"
if [ -s "$GBAM" ]; then
  $SAMTOOLS view -@ 8 -N "$TD/names.txt" "$GBAM" > "$TD/genome_records.sam"
else
  : > "$TD/genome_records.sam"
fi

echo "[$(date)] (3) FASTQ records"
for m in 1 2; do
  F=$FQ1; [ $m = 2 ] && F=$FQ2
  zcat "$F" | awk 'NR==FNR{n["@"$1]; next} FNR%4==1{split($1,a," "); p=(a[1] in n)} p' "$TD/names.txt" - > "$TD/reads_R$m.fq"
done
grep -F -f <(tr ',' '\n' <<< "$TARGETS" | sed 's/.*/transcript_id "&"/') "$FGTF" | awk -F'\t' '$3=="exon"' > "$TD/targets_exons.gtf"

echo "[$(date)] (4) per-read summary"
$PY -I - "$TD" "$RD" "$TARGETS" "$TX2G" <<'EOF'
import sys, re, collections
td, rd, targets_csv, tx2g_path = sys.argv[1:5]
targets = targets_csv.split(',')
tx2g = {}
for l in open(tx2g_path):
    t, g, n = l.rstrip('\n').split('\t'); tx2g[t] = (g, n)
# exon structures of targets (transcript order 5'->3')
ex = collections.defaultdict(list); txchr = {}; txstr = {}
for l in open(f"{td}/targets_exons.gtf"):
    p = l.split('\t'); t = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
    ex[t].append((int(p[3]), int(p[4]))); txchr[t] = p[0]; txstr[t] = p[6]
for t in ex:
    ex[t].sort(reverse=(txstr[t] == '-'))
def tx2genome(t, pos):
    """1-based transcript coordinate -> genome coordinate"""
    cum = 0
    for s, e in ex[t]:
        L = e - s + 1
        if pos <= cum + L:
            off = pos - cum - 1
            return s + off if txstr[t] == '+' else e - off
        cum += L
    return None
cig_re = re.compile(r'(\d+)([MIDNSHP=X])')
def cig_stats(c):
    m = s = ref = 0; ops = cig_re.findall(c)
    for n, op in ops:
        n = int(n)
        if op in 'M=X': m += n; ref += n
        elif op in 'DN': ref += n
        elif op == 'S': s += n
    return m, s, ref
# FASTQ lengths / sequences
seqs = {}
for m in (1, 2):
    with open(f"{td}/reads_R{m}.fq") as f:
        while True:
            h = f.readline()
            if not h: break
            s = f.readline().strip(); f.readline(); f.readline()
            seqs.setdefault(h[1:].split()[0], {})[m] = s
absent = set(l.strip() for l in open(f"{td}/names_absent_from_genome.txt"))
# genome records
gen = collections.defaultdict(list)
for l in open(f"{td}/genome_records.sam"):
    p = l.split('\t')
    nh = re.search(r'NH:i:(\d+)', l)
    gen[p[0]].append((int(p[1]), p[2], int(p[3]), p[5], int(nh.group(1)) if nh else -1))
# transcriptomic records grouped by read
txr = collections.defaultdict(list)
for l in open(f"{td}/tx_groups.sam"):
    p = l.rstrip('\n').split('\t')
    nh = re.search(r'NH:i:(\d+)', l)
    txr[p[0]].append(dict(flag=int(p[1]), tx=p[2], pos=int(p[3]), cigar=p[5], seq=p[9],
                          nh=int(nh.group(1)) if nh else -1))
# plasmid/vector motifs. The first four were observed in this dataset fused to
# the NKX3.1 CDS (soft-clipped flanks of NKX3-1 reads), i.e. they are
# demonstrably plasmid-derived in this experiment; the rest are standard
# pcDNA3.1-family elements used only as annotation.
def rc(s): return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
motifs = {
  'attB2_EGFP_junction(NKX3.1-3prime)': 'ACCCAGCTTTCTTGTACAAAGTGGTTGATGCTGTTAACATGGTGAGCAAGGGCGAGGAG',
  'EGFP_Nterm': 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGG',
  'attB1(NKX3.1-5prime)': 'CAAAAAAGCAGGCT',
  'TOPO_site': 'AAGGGCGAATTC',
  'bGH_polyA': 'CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCC',
  'SV40_polyA': 'AACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGC',
  'NeoR_Nterm': 'ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTGG',
  'T7_promoter': 'TAATACGACTCACTATAGGG',
}
K = 16
mk = {}
for name, s in motifs.items():
    for i in range(0, len(s) - K + 1):
        mk[s[i:i+K]] = name; mk[rc(s[i:i+K])] = name
def motif_hits(seq):
    hits = set()
    for i in range(0, len(seq) - K + 1):
        h = mk.get(seq[i:i+K])
        if h: hits.add(h)
    return hits
rows = []
for r, recs in txr.items():
    on_t = [x for x in recs if x['tx'] in targets]
    if not on_t: continue
    L1 = len(seqs.get(r, {}).get(1, '')); L2 = len(seqs.get(r, {}).get(2, ''))
    hits = motif_hits(seqs.get(r, {}).get(1, '')) | motif_hits(seqs.get(r, {}).get(2, ''))
    g = gen.get(r, [])
    gprim = [x for x in g if not (x[0] & 0x100)]
    for t in sorted(set(x['tx'] for x in on_t)):
        mates = [x for x in on_t if x['tx'] == t]
        c1 = ';'.join(x['cigar'] for x in mates if x['flag'] & 0x40)
        c2 = ';'.join(x['cigar'] for x in mates if x['flag'] & 0x80)
        M = sum(cig_stats(x['cigar'])[0] for x in mates)
        S = sum(cig_stats(x['cigar'])[1] for x in mates)
        starts = [x['pos'] for x in mates]; ends = [x['pos'] + cig_stats(x['cigar'])[2] - 1 for x in mates]
        gs = [tx2genome(t, min(starts)), tx2genome(t, max(ends))]
        locus = f"{txchr[t]}:{min(gs)}-{max(gs)}" if None not in gs else 'NA'
        rows.append(dict(read=r, transcript=t, gene=tx2g.get(t, ('', ''))[1], NH_tx=mates[0]['nh'],
            n_tx_alignments=len(set((x['tx'], x['pos']) for x in recs if x['flag'] & 0x40 or not x['flag'] & 0x1)),
            other_tx_genes=','.join(sorted(set(tx2g.get(x['tx'], ('', ''))[1] for x in recs) - {tx2g.get(t, ('', ''))[1]})),
            len_R1=L1, len_R2=L2, cigar_tx_R1=c1, cigar_tx_R2=c2, aligned_bases=M, softclipped_bases=S,
            aligned_frac_of_pair=round(M / (L1 + L2), 3) if L1 + L2 else 'NA',
            below_STAR_0p66_rule=(M < 0.66 * (L1 + L2 - 1)) if L1 + L2 else 'NA',
            implied_genomic_locus=locus, in_genome_BAM=(r not in absent),
            genome_primary=';'.join(f"{x[1]}:{x[2]}:{x[3]}:NH{x[4]}" for x in gprim[:2]),
            plasmid_motifs=','.join(sorted(hits)),
            seq_R1=seqs.get(r, {}).get(1, ''), seq_R2=seqs.get(r, {}).get(2, '')))
cols = list(rows[0].keys()) if rows else []
with open(f"{rd}/traced_reads.tsv", 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for x in rows: o.write('\t'.join(str(x[c]) for c in cols) + '\n')
# per transcript summary
summ = []
for t in targets:
    rr = [x for x in rows if x['transcript'] == t]
    n = len(rr)
    def frac(f): return round(sum(1 for x in rr if f(x)) / n, 3) if n else 'NA'
    loci = collections.Counter(x['implied_genomic_locus'].split(':')[0] + ':' + str(int(x['implied_genomic_locus'].split(':')[1].split('-')[0]) // 1000) + 'kb'
                               for x in rr if x['implied_genomic_locus'] != 'NA')
    mot = collections.Counter(m for x in rr for m in x['plasmid_motifs'].split(',') if m)
    al = sorted(x['aligned_bases'] for x in rr)
    summ.append(dict(transcript=t, gene=tx2g.get(t, ('', ''))[1], reads_in_toTranscriptome=n,
        absent_from_genome_BAM=sum(1 for x in rr if not x['in_genome_BAM']),
        frac_absent_from_genome_BAM=frac(lambda x: not x['in_genome_BAM']),
        median_aligned_bases=al[len(al)//2] if al else 'NA',
        frac_below_STAR_0p66_rule=frac(lambda x: x['below_STAR_0p66_rule'] is True),
        frac_with_plasmid_motif=frac(lambda x: bool(x['plasmid_motifs'])),
        top_motifs=';'.join(f"{k}:{v}" for k, v in mot.most_common(4)),
        frac_NH_tx_gt1=frac(lambda x: x['NH_tx'] > 1),
        top_implied_loci=';'.join(f"{k}:{v}" for k, v in loci.most_common(3))))
cols = list(summ[0].keys())
with open(f"{rd}/traced_per_transcript.tsv", 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for x in summ: o.write('\t'.join(str(x[c]) for c in cols) + '\n')
print(open(f"{rd}/traced_per_transcript.tsv").read())
EOF

echo "[$(date)] (5) global per-gene profile of reads absent from the genome BAM"
if [ -s "$SD/tx_only_names.txt" ] && [ ! -s "$RD/per_gene_txonly_profile.tsv.gz" ]; then
  $SAMTOOLS view -@ 8 "$TXBAM" | cut -f1,3 | awk -F'\t' -v OFS='\t' '
    FILENAME==ARGV[1]{g[$1]=$2; nm[$2]=$3; next}
    FILENAME==ARGV[2]{o[$1]; next}
    {gene=g[$2]; k=$1 SUBSEP gene; if (k!=last) {all[gene]++; if ($1 in o) only[gene]++; last=k}}
    END{print "gene_id","gene_name","reads_toTranscriptome","reads_absent_from_genome_BAM","frac_absent";
        for (x in all) print x, nm[x], all[x], (x in only ? only[x] : 0), sprintf("%.4f", (x in only ? only[x] : 0)/all[x])}' \
    "$TX2G" "$SD/tx_only_names.txt" - | sort -t$'\t' -k4,4nr | gzip > "$RD/per_gene_txonly_profile.tsv.gz"
  zcat "$RD/per_gene_txonly_profile.tsv.gz" | head -25 || true
fi
# keep small raw evidence
gzip -c "$TD/tx_groups.sam" > "$RD/tx_groups_touching_targets.sam.gz"
gzip -c "$TD/genome_records.sam" > "$RD/genome_records_of_traced_reads.sam.gz"
echo "[$(date)] done"
