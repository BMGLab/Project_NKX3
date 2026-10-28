#!/usr/bin/env bash
# 16_rootcause_12_plasmid_screen_all.sh
# Quick per-sample screen of construct sequences in all 12 libraries: the first
# 2,000,000 reads of each raw R1 and R2 FASTQ (sequencer order) are streamed and
# searched for exact construct motifs (either orientation). Values are reads
# containing the motif per million reads screened. No temporary files.
set -euo pipefail
RAW=/mnt/lun2/rawdata/BilgeDebelec/BAP-RNAseq
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
mkdir -p "$OUT/tables"
$PY -I - "$RAW" "$OUT/tables/12_plasmid_screen_all_samples.tsv" 2000000 <<'EOF'
import sys, subprocess
raw, out, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
group = {'24-004-1_S56': 'Control-basal', '24-004-2_S52': 'Control-basal', '24-004-3_S41': 'Control-basal',
         '24-004-4_S27': 'Control-inflammation', '24-004-5_S19': 'Control-inflammation', '24-004-6_S17': 'Control-inflammation',
         '24-004-7_S14': 'NKX3.1-basal', '24-004-8_S53': 'NKX3.1-basal', '24-004-9_S42': 'NKX3.1-basal',
         '24-004-10_S20': 'NKX3.1-inflammation', '24-004-11_S18': 'NKX3.1-inflammation', '24004-12_S15': 'NKX3.1-inflammation'}
motifs = [('CMV/T7_5UTR_to_attB1', 'GACTCACTATAGGGAGACCCAAGCTGGCTAGTTAAGCTGAGCATCAACAAGTTTGTACAAAAAAGCAGGCT'),
          ('attB2_to_FP_ATG', 'GTGGTTGATGCTGTTAACATGGTGAGC'),
          ('FP_mid(GFP-family)', 'CTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCC'),
          ('FP_aa203_Y(YFP-type)', 'CACTACCTGAGCTACCAGTCCGCCCTG'),
          ('BlastR_bsd_ORF', 'CATCTCTGAAGACTACAGCGTCGCCAGCGCAGCTCTCTCTAGCGACGGCCGC'),
          ('Xpress_epitope(DLYDDDDK)', 'CTGTACGACGATGACGATAAG'),
          ('bGH_polyA_region', 'CTGTGCCTTCTAGTTGCCAGCCATC'),
          # observed in control-library reads that STAR left unmapped (assigned to C3orf20, LINC02626, NTNG1/SOX2-OT, CSH2)
          ('HisMax_5UTR_MCS_SP163', 'GCTAGCGTTTAAACTTAAGCTTAGCGCAGAGGCTTGGGGCAGCCGAGCGG'),
          ('lacZ_5prime(alpha_region)', 'CAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGG'),
          ('lacZ_internal', 'CTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGATGG'),
          ('MCS_to_bGH(XhoI-XbaI-ApaI-PmeI)', 'CTCGAGTCTAGAGGGCCCGTTTAAACCCGCTGATCAGCC'),
          ('NKX3-1_CDS_start', 'ATGCTCAGGGTTCCGGAGCCGC')]
def rc(s): return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
pats = [(m, s, rc(s)) for m, s in motifs]
with open(out, 'w') as o:
    o.write('sample\tgroup\treads_screened\t' + '\t'.join(m for m, _ in motifs) + '\n')
    for s in group:
        cnt = dict.fromkeys([m for m, _ in motifs], 0); n = 0
        for mate in (1, 2):
            p = subprocess.Popen(f"zcat {raw}/{s}_R{mate}_001.fastq.gz | head -n {4 * N}", shell=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
            for i, line in enumerate(p.stdout):
                if i % 4 != 1: continue
                n += 1
                for m, a, b in pats:
                    if a in line or b in line: cnt[m] += 1
            p.wait()
        o.write(f"{s}\t{group[s]}\t{n}\t" + '\t'.join(f"{cnt[m] / n * 1e6:.1f}" for m, _ in motifs) + '\n')
        o.flush()
print(open(out).read())
EOF
