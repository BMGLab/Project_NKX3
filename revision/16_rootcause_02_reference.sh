#!/usr/bin/env bash
# 16_rootcause_02_reference.sh
# Recreate the nf-core/rnaseq 3.14.0 reference intermediates:
#  (1) GTF_FILTER: keep GTF lines whose seqname is in the genome FASTA and that
#      carry a transcript_id (same rule as nf-core 3.14.0 bin/filter_gtf.py)
#  (2) MAKE_TRANSCRIPTS_FASTA: rsem-prepare-reference --gtf <filtered> <fasta> rsem/genome
#      -> genome.transcripts.fa  (nf-core used RSEM 1.3.1; we have RSEM 1.3.3,
#      whose reference-extraction code is unchanged for this purpose)
#  (3) compare target names / order / lengths with the nf-core quant.sf
set -euo pipefail
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
REF=/home/sharedFolder/humanSTARindex
GTF=$REF/Homo_sapiens.GRCh38.99.gtf
FA=$REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
SALBIN=/home/biolab/miniconda3/envs/nkx3_salmon/bin
NFQ=/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon/24-004-7_S14/quant.sf
R=$WORK/ref
mkdir -p "$R/rsem" "$OUT/logs" "$OUT/tables"
export TMPDIR=$WORK/tmp; mkdir -p "$TMPDIR"

FGTF=$R/Homo_sapiens.GRCh38.dna.primary_assembly.filtered.gtf
if [ ! -s "$FGTF" ]; then
  echo "[$(date)] GTF filter"
  $PY -I - "$FA.fai" "$GTF" "$FGTF" <<'EOF'
import re, sys
fai, gtf_in, gtf_out = sys.argv[1:4]
seqs = {l.split('\t')[0] for l in open(fai)}
n_in = n_out = 0
with open(gtf_in) as g, open(gtf_out, 'w') as o:
    for line in g:
        n_in += 1
        if line.split('\t')[0] in seqs and re.search(r'transcript_id "([^"]+)"', line):
            o.write(line); n_out += 1
print(f"GTF lines in: {n_in}; kept: {n_out}")
EOF
fi
awk -F'\t' '$3=="transcript"' "$FGTF" | wc -l | awk '{print "transcripts in filtered GTF:",$1}'

if [ ! -s "$R/rsem/genome.transcripts.fa" ]; then
  echo "[$(date)] rsem-prepare-reference"
  ( cd "$R" && PATH=$SALBIN:$PATH rsem-prepare-reference --gtf "$FGTF" "$FA" rsem/genome ) \
     > "$OUT/logs/02_rsem_prepare_reference.log" 2>&1
fi
$SALBIN/rsem-prepare-reference --version 2>/dev/null || true

echo "[$(date)] compare genome.transcripts.fa with nf-core quant.sf"
$PY -I - "$R/rsem/genome.transcripts.fa" "$NFQ" "$OUT/tables/02_transcript_fasta_vs_nfcore_quantsf.tsv" <<'EOF'
import sys
fa, q, out = sys.argv[1:4]
names, lens = [], {}
cur = None
for line in open(fa):
    if line.startswith('>'):
        cur = line[1:].split()[0]; names.append(cur); lens[cur] = 0
    else:
        lens[cur] += len(line.strip())
qn, ql = [], {}
with open(q) as f:
    next(f)
    for line in f:
        p = line.split('\t'); qn.append(p[0]); ql[p[0]] = int(p[1])
same_set = set(names) == set(qn)
same_order = names == qn
pos_fa = {n: i for i, n in enumerate(names)}
same_pos = sum(1 for i, n in enumerate(qn) if pos_fa.get(n) == i)
len_mismatch = [n for n in qn if n in lens and lens[n] != ql[n]]
with open(out, 'w') as o:
    o.write("metric\tvalue\n")
    o.write(f"n_targets_transcripts_fa\t{len(names)}\n")
    o.write(f"n_targets_nfcore_quant_sf\t{len(qn)}\n")
    o.write(f"same_name_set\t{same_set}\n")
    o.write(f"same_order\t{same_order}\n")
    o.write(f"targets_at_identical_index\t{same_pos}\n")
    o.write(f"length_mismatches_fa_vs_quantsf\t{len(len_mismatch)}\n")
    o.write(f"first_5_fa\t{','.join(names[:5])}\n")
    o.write(f"first_5_quantsf\t{','.join(qn[:5])}\n")
    o.write(f"examples_length_mismatch\t{','.join(len_mismatch[:10])}\n")
print(open(out).read())
EOF
echo "[$(date)] done"
