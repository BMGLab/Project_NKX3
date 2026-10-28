#!/usr/bin/env bash
# 16_rootcause_11_bookkeeping.sh
# One table of read accounting for the nf-core original run and every re-run:
# STAR Log.final.out categories, distinct read names in the genome BAM and in
# Aligned.toTranscriptome.out.bam, reads present only in the latter, and the
# number of fragments Salmon processed.
set -euo pipefail
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
NF=/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
$PY -I - "$OUT" "$NF" <<'EOF'
import sys, os, re, json, glob
out, nf = sys.argv[1:3]
def star_log(p):
    d = {}
    for l in open(p):
        if '|' not in l: continue
        k, v = [x.strip() for x in l.split('|', 1)]
        d[k] = v
    g = lambda k: d.get(k, 'NA').replace('%', '')
    return {'STAR_input_pairs': g('Number of input reads'), 'uniquely_mapped': g('Uniquely mapped reads number'),
            'multi_mapped': g('Number of reads mapped to multiple loci'), 'too_many_loci': g('Number of reads mapped to too many loci'),
            'unmapped_too_short': g('Number of reads unmapped: too short'), 'unmapped_other': g('Number of reads unmapped: other')}
def star_version(p):
    for l in open(p):
        if l.startswith('STAR version='): return l.strip().split('=')[1]
    return 'NA'
def ban(p):
    for l in open(p):
        m = re.search(r'--quantTranscriptomeBan (\S+)', l)
        if m: return m.group(1)
    return 'NA'
rows = []
desc = {'repro_rawtrim': 'raw FASTQ, Trim Galore-equivalent trimming; nf-core options',
        'repro_bamFASTQ': 'FASTQ regenerated from the nf-core genome BAM (STAR-unmapped reads absent); nf-core options',
        'banISS_rawtrim': 'as repro_rawtrim but --quantTranscriptomeBan IndelSoftclipSingleend',
        'star2710b_rawtrim': 'as repro_rawtrim but STAR 2.7.10b',
        'repro_rawtrim_genomeMappedOnly': 'repro_rawtrim toTranscriptome BAM restricted to reads in the genome BAM',
        'salmonSA_rawtrim': 'Salmon selective alignment on the trimmed FASTQ (no STAR)'}
for s in sorted(os.listdir(nf)):
    if not os.path.exists(f"{nf}/log/{s}.Log.final.out") if os.path.isdir(f"{nf}/{s}") else True: continue
for s in ['24-004-7_S14', '24-004-1_S56']:
    r = dict(run='nfcore_original', sample=s, description='nf-core/rnaseq 3.14.0 output (March 2024)',
             STAR_version=star_version(f"{nf}/log/{s}.Log.out"), quantTranscriptomeBan=ban(f"{nf}/log/{s}.Log.out"))
    r.update(star_log(f"{nf}/log/{s}.Log.final.out"))
    r.update(reads_in_genome_BAM='NA (not re-counted)', reads_in_toTranscriptome_BAM='NA (BAM not kept by nf-core)',
             reads_in_toTranscriptome_but_absent_from_genome_BAM='NA',
             Salmon_fragments_processed=json.load(open(f"{nf}/{s}/aux_info/meta_info.json"))['num_processed'])
    rows.append(r)
    for d in sorted(glob.glob(f"{out}/runs/*/{s}")):
        tag = d.split('/')[-2]
        r = dict(run=tag, sample=s, description=desc.get(tag, ''))
        lo = f"{d}/{s}.Log.out"; lf = f"{d}/{s}.Log.final.out"
        r['STAR_version'] = star_version(lo) if os.path.exists(lo) else ('2.7.9a (reused)' if 'genomeMappedOnly' in tag else 'none')
        r['quantTranscriptomeBan'] = ban(lo) if os.path.exists(lo) else ('Singleend (reused)' if 'genomeMappedOnly' in tag else 'none')
        r.update(star_log(lf) if os.path.exists(lf) else {})
        rc = f"{d}/read_name_counts.tsv"
        if os.path.exists(rc):
            for l in list(open(rc))[1:]:
                k, v = l.rstrip('\n').split('\t'); r[k] = v
        mi = f"{d}/meta_info.json"
        if os.path.exists(mi):
            m = json.load(open(mi)); r['Salmon_fragments_processed'] = m['num_processed']
            if tag.startswith('salmonSA'): r['Salmon_fragments_mapped'] = m['num_mapped']
        rows.append(r)
cols = ['run', 'sample', 'description', 'STAR_version', 'quantTranscriptomeBan', 'STAR_input_pairs', 'uniquely_mapped', 'multi_mapped',
        'too_many_loci', 'unmapped_too_short', 'unmapped_other', 'reads_in_genome_BAM', 'reads_in_toTranscriptome_BAM',
        'reads_in_toTranscriptome_but_absent_from_genome_BAM', 'Salmon_fragments_processed', 'Salmon_fragments_mapped']
with open(f"{out}/tables/11_read_bookkeeping.tsv", 'w') as o:
    o.write('\t'.join(cols) + '\n')
    for r in rows: o.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
print(open(f"{out}/tables/11_read_bookkeeping.tsv").read())
EOF
