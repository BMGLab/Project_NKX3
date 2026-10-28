#!/usr/bin/env bash
# 16_rootcause_10_examples.sh
# Example reads (absent from STAR's genome BAM) that the nf-core-style STAR
# toTranscriptome BAM assigns to the flagged genes: per mate, the transcriptomic
# CIGAR, the aligned stretch (UPPER case) vs soft-clipped bases (lower case),
# mismatches in the aligned stretch vs the transcript, and the construct
# reference the read pair matches (exact 25-mers; 16_rootcause_09 references).
# usage: 16_rootcause_10_examples.sh <sample_id> <run_tag> [reads_per_gene]
set -euo pipefail
SAMPLE=$1; TAG=$2; NPER=${3:-3}
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
$PY -I - "$OUT/trace/${SAMPLE}_${TAG}" "$OUT/construct/construct_references.fa" "$WORK/ref/rsem/genome.transcripts.fa" \
   "$OUT/tables/10_example_reads_${SAMPLE}.tsv" "$SAMPLE" "$NPER" <<'EOF'
import sys, gzip, re, collections
tdir, ref_fa, tx_fa, out, sample, nper = sys.argv[1:7]; nper = int(nper)
genes = ['BCL11B', 'SSC4D', 'FAM86MP', 'DPEP3', 'IGLL5', 'SRMS', 'TPGS1', 'DIPK1B',
         'NTNG1', 'CSH2', 'C3orf20', 'SOX2-OT', 'ZNF595', 'NRF1', 'LINC02626']
def rc(s): return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
def fasta(path, keep=None):
    d, name = {}, None
    for l in open(path):
        if l.startswith('>'):
            name = l[1:].split()[0]
            if keep is None or name in keep: d[name] = []
            else: name = None
        elif name: d[name].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}
refs = fasta(ref_fa)
K = 25; kidx = {}
for n, s in refs.items():
    for i in range(len(s) - K + 1):
        k = s[i:i + K]
        if 'N' in k: continue
        kidx.setdefault(k, n); kidx.setdefault(rc(k), n)
rows = [l.rstrip('\n').split('\t') for l in open(f"{tdir}/traced_reads.tsv")]
hdr = rows[0]; rows = [dict(zip(hdr, r)) for r in rows[1:]]
chosen = collections.OrderedDict()
for g in genes:
    rr = [r for r in rows if r['gene'] == g and r['in_genome_BAM'] == 'False']
    seen = set()
    for r in rr:
        if r['read'] in seen: continue
        seen.add(r['read']); chosen[(g, r['read'])] = r
        if len(seen) >= nper: break
need_tx = set(r['transcript'] for r in chosen.values())
txseq = fasta(tx_fa, need_tx)
recs = collections.defaultdict(list)
names = set(k[1] for k in chosen)
with gzip.open(f"{tdir}/tx_groups_touching_targets.sam.gz", 'rt') as f:
    for l in f:
        p = l.split('\t')
        if p[0] in names: recs[p[0]].append(p)
cig_re = re.compile(r'(\d+)([MIDNSHP=X])')
out_rows = []
for (g, rname), r in chosen.items():
    hits = collections.Counter()
    for s in (r['seq_R1'], r['seq_R2']):
        for i in range(len(s) - K + 1):
            h = kidx.get(s[i:i + K])
            if h: hits[h] += 1
    best = hits.most_common(1)[0] if hits else ('none', 0)
    for p in recs[rname]:
        if p[2] != r['transcript']: continue
        flag, pos, cig, seq = int(p[1]), int(p[3]), p[5], p[9]
        mate = 'R1' if flag & 0x40 else 'R2'
        disp, ri, ti, aligned, mm = [], 0, pos - 1, 0, 0
        t = txseq.get(p[2], '')
        for n, op in cig_re.findall(cig):
            n = int(n)
            if op == 'S': disp.append(seq[ri:ri + n].lower()); ri += n
            elif op in 'M=X':
                seg = seq[ri:ri + n]; disp.append(seg.upper())
                mm += sum(1 for a, b in zip(seg, t[ti:ti + n]) if a != b)
                aligned += n; ri += n; ti += n
            elif op == 'I': disp.append(seq[ri:ri + n].lower()); ri += n
            elif op in 'DN': ti += n
        out_rows.append([sample, g, p[2], rname, mate, str(len(seq)), str(pos), cig, str(aligned), str(mm),
                         best[0], str(best[1]), ''.join(disp)])
with open(out, 'w') as o:
    o.write('\t'.join(['sample', 'gene', 'transcript', 'read', 'mate', 'read_len', 'tx_pos', 'tx_cigar',
                       'aligned_nt', 'mismatches_in_aligned', 'construct_reference_best', 'construct_25mer_hits_pair',
                       'sequence_transcript_orientation(ALIGNED_upper|softclip_lower)']) + '\n')
    for x in out_rows: o.write('\t'.join(x) + '\n')
print(f"{len(out_rows)} rows -> {out}")
for x in out_rows[:12]: print('\t'.join(x[:12]))
EOF
