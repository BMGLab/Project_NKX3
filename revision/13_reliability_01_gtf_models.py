#!/usr/bin/env python
"""Parse the Ensembl 99 GTF (the annotation used by STAR, Salmon and featureCounts) into
gene / transcript / exon models and cache them for the reliability audit.
Run: /home/biolab/miniconda3/envs/nkx3_py/bin/python -I 13_reliability_01_gtf_models.py
Output (temp): /mnt/ssd1/Projects/NKX3_wd/work/reliability/gtf_models.pkl
"""
import re, pickle, collections, os, sys
GTF = '/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf'
TMP = '/mnt/ssd1/Projects/NKX3_wd/work/reliability'
os.makedirs(TMP, exist_ok=True)
attr_re = {k: re.compile(k + ' "([^"]+)"') for k in
           ['gene_id', 'gene_name', 'gene_biotype', 'transcript_id', 'transcript_name', 'transcript_biotype', 'tag']}

genes = {}                                    # gene_id -> dict
tx = collections.defaultdict(lambda: {'exons': [], 'cds': []})
with open(GTF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        feat = p[2]
        if feat not in ('gene', 'exon', 'CDS'):
            continue
        a = p[8]
        gid = attr_re['gene_id'].search(a).group(1)
        if feat == 'gene':
            m = attr_re['gene_name'].search(a)
            genes[gid] = dict(gene_id=gid, symbol=m.group(1) if m else gid,
                              biotype=attr_re['gene_biotype'].search(a).group(1),
                              chrom=p[0], start=int(p[3]), end=int(p[4]), strand=p[6])
            continue
        tid = attr_re['transcript_id'].search(a).group(1)
        t = tx[tid]
        if 'gene_id' not in t:
            tn = attr_re['transcript_name'].search(a); tb = attr_re['transcript_biotype'].search(a)
            t.update(gene_id=gid, chrom=p[0], strand=p[6], name=tn.group(1) if tn else tid,
                     biotype=tb.group(1) if tb else '', tags=set())
        for m in re.finditer('tag "([^"]+)"', a):
            t['tags'].add(m.group(1))
        iv = (int(p[3]) - 1, int(p[4]))      # 0-based half-open
        (t['exons'] if feat == 'exon' else t['cds']).append(iv)

def merge(iv):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]

for g in genes.values():
    g['transcripts'] = []
for tid, t in tx.items():
    t['exons'] = sorted(t['exons']); t['cds'] = sorted(t['cds']); t['tags'] = sorted(t['tags'])
    t['length'] = sum(e - s for s, e in t['exons'])
    genes[t['gene_id']]['transcripts'].append(tid)
for g in genes.values():
    g['exons_union'] = merge([iv for tid in g['transcripts'] for iv in tx[tid]['exons']])
    g['exonic_length'] = sum(e - s for s, e in g['exons_union'])

pickle.dump({'genes': genes, 'tx': dict(tx)}, open(os.path.join(TMP, 'gtf_models.pkl'), 'wb'))
print('genes', len(genes), 'transcripts', len(tx))
