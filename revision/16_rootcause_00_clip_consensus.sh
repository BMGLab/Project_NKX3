#!/usr/bin/env bash
# 16_rootcause_00_clip_consensus.sh
# Recover the plasmid sequences fused to the NKX3.1 CDS: consensus of the
# soft-clipped parts of primary alignments at the NKX3-1 locus (chr8, minus
# strand) in the nf-core genome BAM of one NKX3.1 library (indexed region query).
# usage: 16_rootcause_00_clip_consensus.sh [sample_id]
set -euo pipefail
SAMPLE=${1:-24-004-7_S14}
BAM=/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon/${SAMPLE}.markdup.sorted.bam
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
mkdir -p "$OUT/tables"
$PY -I - "$BAM" > "$OUT/tables/00_NKX3-1_junction_clip_consensus_${SAMPLE}.txt" <<'EOF'
import pysam, collections, sys
bam = pysam.AlignmentFile(sys.argv[1])
def rc(s): return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
left, right = collections.Counter(), collections.Counter()
L, R = collections.defaultdict(list), collections.defaultdict(list)
for a in bam.fetch('8', 23680900, 23683100):
    if a.is_secondary or a.is_supplementary: continue
    c = a.cigartuples
    if c[0][0] == 4 and c[0][1] >= 12:
        left[a.reference_start + 1] += 1; L[a.reference_start + 1].append(a.query_sequence[:c[0][1]][::-1])
    if c[-1][0] == 4 and c[-1][1] >= 12:
        right[a.reference_end] += 1; R[a.reference_end].append(a.query_sequence[-c[-1][1]:])
def consensus(seqs, minfrac=0.7, mindepth=5):
    out = []
    for i in range(max(len(s) for s in seqs)):
        col = collections.Counter(s[i] for s in seqs if len(s) > i); n = sum(col.values())
        if n < mindepth: break
        b, k = col.most_common(1)[0]; out.append(b if k / n >= minfrac else 'N')
    return ''.join(out)
print("NKX3-1 (ENSG00000167034) is on the minus strand: CDS start (ATG) at chr8:23,682,889; stop codon chr8:23,681,221-23,681,223.")
print("Top left-soft-clip anchors (alignment start, 1-based):", left.most_common(4))
print("Top right-soft-clip anchors (alignment end, 1-based):", right.most_common(4))
p, n = left.most_common(1)[0]
c = consensus(L[p])[::-1]
print(f"\n3' junction: {n} reads start at chr8:{p} (20 nt downstream of the NKX3.1 stop codon) with a left soft clip.")
print(f"  plasmid sequence downstream of NKX3.1, mRNA sense (junction on the left):\n  {rc(c)}")
p, n = right.most_common(1)[0]
c = consensus(R[p])
print(f"\n5' junction: {n} reads end at chr8:{p} (A of the NKX3.1 ATG) with a right soft clip.")
print(f"  plasmid sequence upstream of the NKX3.1 ATG, mRNA sense (junction on the right):\n  {rc(c)}")
EOF
cat "$OUT/tables/00_NKX3-1_junction_clip_consensus_${SAMPLE}.txt"
