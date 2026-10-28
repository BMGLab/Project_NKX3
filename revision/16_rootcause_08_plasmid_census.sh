#!/usr/bin/env bash
# 16_rootcause_08_plasmid_census.sh
# Count read pairs whose mate-1 or mate-2 sequence contains short plasmid/vector
# motifs (exact match, either orientation) in the trimmed FASTQ of a sample.
#   group "observed": sequences seen in THIS dataset inside reads that STAR left
#     unmapped but projected onto human transcripts (BCL11B, SSC4D, DPEP3, SRMS,
#     ...), or fused to the NKX3.1 CDS in soft-clipped NKX3-1 reads of
#     24-004-7_S14 (i.e. demonstrably construct-derived);
#   group "reference": standard elements of pcDNA3.1/HisMax-type vectors (as
#     described in the manuscript Methods), for annotation only;
#   group "NKX3-1": NKX3-1 CDS start/end (endogenous + plasmid NKX3.1).
# usage: 16_rootcause_08_plasmid_census.sh <sample_id>
set -euo pipefail
SAMPLE=$1
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
FQ=$WORK/fastq_trim/$SAMPLE
mkdir -p "$OUT/tables"
$PY -I - "$FQ/${SAMPLE}_1_val_1.fq.gz" "$FQ/${SAMPLE}_2_val_2.fq.gz" "$SAMPLE" "$OUT/tables/08_plasmid_motif_census_${SAMPLE}.tsv" <<'EOF'
import sys, subprocess
r1, r2, sample, out = sys.argv[1:5]
motifs = [  # (group, name, sequence)
 ('observed', 'CMV/T7_5UTR_to_attB1', 'GACTCACTATAGGGAGACCCAAGCTGGCTAGTTAAGCTGAGCATCAACAAGTTTGTACAAAAAAGCAGGCT'),
 ('observed', 'attB1', 'ACAAGTTTGTACAAAAAAGCAGGCT'),
 ('observed', 'junction_vector_to_NKX3-1_ATG', 'TTCGCCCTTGGATCCATGCTCAGGG'),
 ('observed', 'junction_NKX3-1_3UTR20_to_vector', 'GGTGACAACCGAATTCAAGGGCG'),
 ('observed', 'attB2', 'ACCCAGCTTTCTTGTACAAAGTGG'),
 ('observed', 'attB2_to_FP_ATG', 'GTGGTTGATGCTGTTAACATGGTGAGC'),
 ('observed', 'FP_Nterm(GFP-family)', 'ATGGTGAGCAAGGGCGAGGAG'),
 ('observed', 'FP_mid(GFP-family)', 'CTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCC'),
 ('observed', 'FP_aa203_Y(YFP-type)', 'CACTACCTGAGCTACCAGTCCGCCCTG'),
 ('reference', 'FP_aa203_T(EGFP-type)', 'CACTACCTGAGCACCCAGTCCGCCCTG'),
 ('observed', 'FP_Cterm_stop_PmeI', 'GGCATGGACGAGCTGTACAAGTAATGATAAGTTTAAAC'),
 ('observed', 'BlastR_bsd_ORF', 'CATCTCTGAAGACTACAGCGTCGCCAGCGCAGCTCTCTCTAGCGACGGCCGC'),
 ('reference', 'NeoR_Nterm', 'ATGATTGAACAAGATGGATTGCACGC'),
 ('reference', 'bGH_polyA_region', 'CTGTGCCTTCTAGTTGCCAGCCATC'),
 ('reference', 'SV40_polyA', 'TTGTTTATTGCAGCTTATAATGG'),
 ('reference', 'T7_promoter', 'TAATACGACTCACTATAGGG'),
 ('reference', 'Xpress_epitope(DLYDDDDK)', 'CTGTACGACGATGACGATAAG'),
 ('reference', 'His6_tag_CAT', 'CATCATCATCATCATCAT'),
 ('NKX3-1', 'NKX3-1_CDS_start', 'ATGCTCAGGGTTCCGGAGCCGC'),
 ('NKX3-1', 'NKX3-1_CDS_end', 'GGAGCCCAGCTTTTTGGTAATG'),
]
def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
pats = [(g, k, v, rc(v)) for g, k, v in motifs]
counts = {k: 0 for _, k, _ in motifs}; n = 0
p1 = subprocess.Popen(['zcat', r1], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
p2 = subprocess.Popen(['zcat', r2], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
f1, f2 = p1.stdout, p2.stdout
while True:
    h = f1.readline()
    if not h: break
    s1 = f1.readline(); f1.readline(); f1.readline()
    f2.readline(); s2 = f2.readline(); f2.readline(); f2.readline()
    n += 1
    s = s1 + '|' + s2
    for g, k, a, b in pats:
        if a in s or b in s: counts[k] += 1
with open(out, 'w') as o:
    o.write('sample\tgroup\tmotif\tsequence\tread_pairs_with_motif\tper_million_pairs\n')
    for g, k, v in motifs:
        o.write(f"{sample}\t{g}\t{k}\t{v}\t{counts[k]}\t{counts[k] / n * 1e6:.1f}\n")
    o.write(f"{sample}\t\tTOTAL_PAIRS\t\t{n}\t\n")
print(open(out).read())
EOF
