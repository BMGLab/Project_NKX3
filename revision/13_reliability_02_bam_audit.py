#!/usr/bin/env python
"""Reliability task, step 2: per-gene read-evidence audit from the nf-core STAR genome BAMs.

For every gene named in the manuscript text / figures (plus loci that share reads with them),
count fragments (one count per read name, i.e. per read pair) per sample, strand-aware for the
reverse-stranded library (Illumina Stranded mRNA; Salmon libType ISR: read 1 antisense to the
transcript), split into:
  sense / antisense, exonic (any aligned block overlaps the union of the gene's Ensembl 99 exons)
  / intronic (inside the gene span but no exon overlap), unique (NH:i:1) / multi-mapped (NH>1,
  secondary alignments included), duplicate-flagged (0x400, Picard MarkDuplicates), and 'shared'
  (also overlaps exons of another gene on the same strand, i.e. not uniquely attributable).
Family / locus totals count distinct read names across all member genes (so a multi-mapped
fragment placed at SLX1A and SLX1B is counted once).
Also computes strand-specific per-base coverage for the coverage figure (all alignments and
NH==1 only, both mates, spliced blocks only), normalised later to reads per million uniquely
mapped read pairs (STAR Log.final.out).

Inputs are read-only. At most 4 BAM readers run concurrently (slow NFS).
Run: /home/biolab/miniconda3/envs/nkx3_py/bin/python -I 13_reliability_02_bam_audit.py
"""
import os, sys, re, json, pickle, bisect, collections, time
import multiprocessing as mp
import numpy as np
import pandas as pd
import pysam

PROJ = '/mnt/ssd1/Projects/NKX3_wd'
OUT = f'{PROJ}/results/revision_2026/reliability'
TMP = '/mnt/ssd1/Projects/NKX3_wd/work/reliability'
NF = '/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon'
os.makedirs(OUT, exist_ok=True); os.makedirs(TMP, exist_ok=True)

# ----------------------------------------------------------------------------- samples
SAMPLES = [  # (sample_id, group key, display label) in group order
    ('24-004-1_S56', 'vector_ctrl', 'Control-basal'), ('24-004-2_S52', 'vector_ctrl', 'Control-basal'),
    ('24-004-3_S41', 'vector_ctrl', 'Control-basal'),
    ('24-004-4_S27', 'vector_CM', 'Control-inflammation'), ('24-004-5_S19', 'vector_CM', 'Control-inflammation'),
    ('24-004-6_S17', 'vector_CM', 'Control-inflammation'),
    ('24-004-7_S14', 'NKX3.1_ctrl', 'NKX3.1-basal'), ('24-004-8_S53', 'NKX3.1_ctrl', 'NKX3.1-basal'),
    ('24-004-9_S42', 'NKX3.1_ctrl', 'NKX3.1-basal'),
    ('24-004-10_S20', 'NKX3.1_CM', 'NKX3.1-inflammation'), ('24-004-11_S18', 'NKX3.1_CM', 'NKX3.1-inflammation'),
    ('24004-12_S15', 'NKX3.1_CM', 'NKX3.1-inflammation')]
SAMPLE_IDS = [s for s, _, _ in SAMPLES]

# ----------------------------------------------------------------------------- curated gene list
# Where each gene is named in the submitted manuscript and what was claimed.
# Contrast keys (sign = treatment minus reference): CM_in_vector = CM vs control (vector cells);
# NKX_in_ctrl = NKX3.1 vs vector (basal); NKX_in_CM = NKX3.1 vs vector (inflammation).
T1 = 'A: reported result / marker (Abstract, Results, Discussion, Fig 1E-G, Supp Fig 2A-C)'
T2 = 'B: other gene named in manuscript text (context, not reported as DEG)'
T2b = 'B: named by Reviewer 1 (epithelial marker; not in manuscript)'
T3 = 'C: figure label only (Fig 1B-D / Supp Fig 2D-F top-20 |z| labels, Fig 3D-E module members)'
CUR = [
 # symbol, tier, where named, reported claims "contrast:direction" separated by ';'
 ('NTNG1', T1, 'Abstract; Results; Discussion; Fig 1C,F,D,G; Supp Fig 2B,C', 'NKX_in_ctrl:down; NKX_in_CM:down'),
 ('CSH2', T1, 'Abstract; Results; Discussion; Fig 1F,G,D; Supp Fig 2B,C', 'NKX_in_ctrl:down; NKX_in_CM:down'),
 ('PIP4P2', T1, 'Abstract; Results; Fig 1C,F; Supp Fig 2B', 'NKX_in_ctrl:down'),
 ('BCL11B', T1, 'Abstract; Results; Discussion; Fig 1C,D,F,G; Supp Fig 2B,C; Fig 3D,E', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('TPGS1', T1, 'Abstract; Results; Discussion; Fig 1D,F,G; Supp Fig 2B,C; Fig 3D,E', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('SSC4D', T1, 'Abstract; Results; Discussion; Fig 1C,F,G; Supp Fig 2B,C; Fig 3D,E', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('USP3', T1, 'Abstract; Discussion; Fig 1G; Supp Fig 2C', 'NKX_in_CM:up'),
 ('MDM4', T1, 'Abstract; Results; Discussion; Fig 1G; Supp Fig 2C', 'NKX_in_CM:up'),
 ('IGLL5', T1, 'Results; Fig 1D,F,G; Supp Fig 2B,C; Fig 3D,E', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('SLX1A', T1, 'Results; Discussion; Fig 1C,F; Supp Fig 2B; Fig 3D', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('SLX1B', T1, 'Fig 1F; Supp Fig 2B', 'NKX_in_ctrl:down'),
 ('RGPD5', T1, 'Results; Discussion; Fig 1E; Supp Fig 2A', 'CM_in_vector:up'),
 ('COMMD3-BMI1', T1, 'Fig 1E; Supp Fig 2A', 'CM_in_vector:down'),
 ('AC131160.1', T1, 'Fig 1F; Supp Fig 2B', 'NKX_in_ctrl:down'),
 ('NKX3-1', T1, 'Throughout; Fig 1A-G; Supp Fig 2; Fig 3', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('ASDURF', T1, 'Results; Fig 1B,E; Supp Fig 2A', 'CM_in_vector:up'),
 ('ARHGEF19', T1, 'Results; Discussion; Fig 1B,C,E,F; Supp Fig 2A,B; Fig 3D', 'CM_in_vector:up; NKX_in_ctrl:up'),
 ('FRMD4A', T1, 'Results (as down and as up); Discussion (up); Fig 1B,E,G; Supp Fig 2A,C; Fig 3E', 'CM_in_vector:down; NKX_in_CM:up'),
 ('SNX17', T1, 'Results; Discussion; Fig 1B,E; Supp Fig 2A', 'CM_in_vector:down'),
 ('NAA25', T1, 'Results; Fig 1E; Supp Fig 2A', 'CM_in_vector:down'),
 ('CCDC66', T1, 'Results; Discussion; Fig 1E,G; Supp Fig 2A,C', 'CM_in_vector:down; NKX_in_CM:up'),
 ('CAPN12', T1, 'Fig 1F; Supp Fig 2B', 'NKX_in_ctrl:up'),
 ('CHMP4A', T1, 'Fig 1F; Supp Fig 2B', 'NKX_in_ctrl:up'),
 ('EHD1', T1, 'Fig 1C,F; Supp Fig 2B; Fig 3D', 'NKX_in_ctrl:up'),
 ('ZNF548', T1, 'Fig 1F; Supp Fig 2B', 'NKX_in_ctrl:up'),
 ('DIPK1B', T1, 'Fig 1F,G; Supp Fig 2B,C; Fig 3D,E', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('SRMS', T1, 'Results; Discussion; Fig 1D,F,G; Supp Fig 2B,C; Fig 3D,E', 'NKX_in_ctrl:up; NKX_in_CM:up'),
 ('C3orf20', T1, 'Results; Fig 1G; Supp Fig 2C', 'NKX_in_CM:down'),
 ('TRIM66', T1, 'Results; Fig 1D,G; Supp Fig 2C', 'NKX_in_CM:up'),
 ('PLEKHA3', T1, 'Results; Fig 1G; Supp Fig 2C', 'NKX_in_CM:up'),
 ('LRIG1', T1, 'Results; Discussion; Fig 1D,G; Supp Fig 2C; Fig 3E', 'NKX_in_CM:up'),
 ('ZNF331', T1, 'Results; Discussion; Fig 3D', 'NKX_in_ctrl:up'),
 ('CDK4', T1, 'Discussion (NKX3.1 downregulates)', 'NKX_in_CM:down'),
 ('MCM2', T1, 'Discussion ("MCM complex genes")', 'NKX_in_CM:down'),
 ('MCM3', T1, 'Discussion ("MCM complex genes")', 'NKX_in_CM:down'),
 ('MCM4', T1, 'Discussion ("MCM complex genes")', 'NKX_in_CM:down'),
 ('MCM5', T1, 'Discussion ("MCM complex genes")', 'NKX_in_CM:down'),
 ('MCM6', T1, 'Discussion ("MCM complex genes")', 'NKX_in_CM:down'),
 ('MCM7', T1, 'Discussion ("MCM complex genes")', 'NKX_in_CM:down'),
 ('CD44', T1, 'Abstract; Results; Discussion; Fig 4D-E (flow cytometry, protein)', ''),
 ('PROM1', T1, 'Results; Discussion; Fig 4D-E (CD133, flow cytometry, protein)', ''),
 ('CDH2', T1, 'Results; Discussion; Fig 5A (N-cadherin, western blot, protein)', ''),
 ('AR', T2, 'Introduction; Results (DU145 AR-negative)', ''),
 ('MYC', T2, 'Results; Discussion (MYC targets gene set)', ''),
 ('EZH2', T2, 'Discussion (Lu EZH2 targets gene set)', ''),
 ('TP53', T2, 'Results; Discussion (p53 pathway; p53-mutant context)', ''),
 ('POU5F1', T2, 'Introduction (Oct4)', ''),
 ('TNF', T2, 'Introduction; Methods (TNF-alpha in CM)', ''),
 ('TGFB1', T2, 'Introduction; Methods (TGF-beta1 in CM)', ''),
 ('IL6', T2, 'Introduction', ''),
 ('VEGFA', T2, 'Introduction (VEGF)', ''),
 ('IL1B', T2, 'Introduction (IL-1beta)', ''),
 ('GAPDH', T2, 'Methods; Fig 1A, 5A (loading control)', ''),
 ('CDH1', T2b, 'Reviewer 1 major 3 (E-cadherin)', ''),
 ('TJP1', T2b, 'Reviewer 1 major 3 (ZO-1)', ''),
]
FIG3D = ['EHD1', 'SLX1A', 'NKX3-1', 'BCL11B', 'SSC4D', 'ZNF331', 'NPHP4', 'KIAA0895', 'AGAP6', 'TRPM4', 'D2HGDH',
         'ZC3H12A', 'KLHL3', 'ULK4', 'TPGS1', 'FANCA', 'TRPV3', 'DIPK1B', 'MCTS1', 'IGLL5', 'SRMS', 'ARHGEF19']
FIG3E = ['DIPK1B', 'SSC4D', 'FRMD4A', 'IGLL5', 'CNIH4', 'LRIG1', 'SRMS', 'CCN1', 'NKX3-1', 'BCL11B', 'TPGS1',
         'PPIP5K1', 'HSPA6', 'KLRC3']
V15 = f'{PROJ}/results/rna_downstreams_out_v15'
CMP = {  # original Welch-z tables behind Fig 1B-D / Supp Fig 2D-F; labels = top 20 by |z|
    'CM_in_vector': (f'{V15}/cmp_HM-vector-C_vs_HM-vector-CM_results.tsv', -1, 'Fig 1B; Supp Fig 2D'),
    'NKX_in_ctrl': (f'{V15}/cmp_HM-NKX3.1-C_vs_HM-vector-C_results.tsv', 1, 'Fig 1C; Supp Fig 2E'),
    'NKX_in_CM': (f'{V15}/cmp_HM-NKX3.1-CM_vs_HM-vector-CM_results.tsv', 1, 'Fig 1D; Supp Fig 2F'),
}

# Paralogue families named by the reviewer / task, plus locus families built below from
# same-strand exon overlaps (readthrough / nested genes).
PARALOG_FAMILIES = {
    'SLX1A+SLX1B': ['SLX1A', 'SLX1B'],
    'RGPD1-RGPD8 (RGPD7 not annotated in Ensembl 99)': ['RGPD1', 'RGPD2', 'RGPD3', 'RGPD4', 'RGPD5', 'RGPD6', 'RGPD8'],
    'GH1+GH2+CSH1+CSH2+CSHL1': ['GH1', 'GH2', 'CSH1', 'CSH2', 'CSHL1'],
    'COMMD3+BMI1+COMMD3-BMI1': ['COMMD3', 'BMI1', 'COMMD3-BMI1'],
}
# Coverage-figure windows: (panel title, chrom, start, end (1-based, inclusive), sense strand, genes to draw)
COV = [
    ('NKX3-1', '8', 23677697, 23683938, '-', ['NKX3-1']),
    ('BCL11B', '14', 99167287, 99274197, '-', ['BCL11B']),
    ('SSC4D', '7', 76388334, 76410697, '-', ['SSC4D', 'ZP3']),
    ('NTNG1', '1', 107139007, 107484458, '+', ['NTNG1']),
    ('GH/CSH cluster (CSH2)', '17', 63869000, 63921000, '-', ['CSH2', 'GH2', 'CSH1', 'CSHL1', 'GH1']),
    ('SLX1A', '16', 30192887, 30198561, '+', ['SLX1A', 'BOLA2B']),
    ('SLX1B', '16', 29453501, 29459219, '+', ['SLX1B', 'BOLA2']),
    ('TPGS1', '19', 506497, 520654, '+', ['TPGS1']),
    ('DIPK1B', '9', 136711572, 136725742, '+', ['DIPK1B', 'SNHG7']),
    ('SRMS', '20', 63539810, 63548504, '-', ['SRMS']),
    ('IGLL5', '22', 22886780, 22897111, '+', ['IGLL5', 'IGLJ1', 'IGLC1']),
]

# ----------------------------------------------------------------------------- annotation
M = pickle.load(open(f'{TMP}/gtf_models.pkl', 'rb'))
GENES, TX = M['genes'], M['tx']
BYSYM = collections.defaultdict(list)
for g in GENES.values():
    BYSYM[g['symbol']].append(g['gene_id'])
BYCHR = collections.defaultdict(list)
for g in GENES.values():
    BYCHR[g['chrom']].append(g)


def gid_of(sym):
    ids = BYSYM[sym]
    if len(ids) != 1:
        raise SystemExit(f'symbol {sym} maps to {ids}')
    return ids[0]


def overlaps_any(a, b):
    return any(s1 < e2 and s2 < e1 for s1, e1 in a for s2, e2 in b)


def merge(iv):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def same_strand_partners(gid):
    g = GENES[gid]
    out = []
    for h in BYCHR[g['chrom']]:
        if h['gene_id'] == gid or h['strand'] != g['strand'] or h['end'] < g['start'] or h['start'] > g['end']:
            continue
        if overlaps_any(g['exons_union'], h['exons_union']):
            out.append(h['gene_id'])
    return out


# build the curated table
rows = []
seen = set()
for sym, tier, where, claims in CUR:
    rows.append(dict(symbol=sym, tier=tier, where_named=where, reported_claims=claims)); seen.add(sym)
fig_lab = collections.defaultdict(list); fig_claim = collections.defaultdict(list)
for key, (path, sgn, panel) in CMP.items():
    t = pd.read_csv(path, sep='\t').head(20)
    for _, r in t.iterrows():
        d = 'up' if sgn * r['lfc_vst'] > 0 else 'down'
        fig_lab[r['gene']].append(panel); fig_claim[r['gene']].append(f'{key}:{d}')
for s in FIG3D:
    fig_lab[s].append('Fig 3D'); fig_claim[s].append('NKX_in_ctrl:up')
for s in FIG3E:
    fig_lab[s].append('Fig 3E'); fig_claim[s].append('NKX_in_CM:up')
for r in rows:   # add figure-label info to tier A/B genes that are also labelled
    if r['symbol'] in fig_lab:
        r['figure_labels'] = '; '.join(dict.fromkeys(fig_lab[r['symbol']]))
for s in fig_lab:
    if s in seen:
        continue
    cl = '; '.join(dict.fromkeys(fig_claim[s]))
    rows.append(dict(symbol=s, tier=T3, where_named='; '.join(dict.fromkeys(fig_lab[s])), reported_claims=cl,
                     figure_labels='; '.join(dict.fromkeys(fig_lab[s]))))
    seen.add(s)
named = pd.DataFrame(rows)
named['figure_labels'] = named['figure_labels'].fillna('')
named['gene_id'] = [gid_of(s) for s in named.symbol]
named['biotype'] = [GENES[g]['biotype'] for g in named.gene_id]
named['locus'] = [f"chr{GENES[g]['chrom']}:{GENES[g]['start']}-{GENES[g]['end']}({GENES[g]['strand']})" for g in named.gene_id]
named['same_strand_exon_overlap_genes'] = ['; '.join(GENES[h]['symbol'] for h in same_strand_partners(g)) for g in named.gene_id]

# families: paralogue families + locus families (named gene + same-strand exon-overlapping genes)
FAMILIES = {}
for fname, mem in PARALOG_FAMILIES.items():
    FAMILIES[fname] = [gid_of(s) for s in mem]
for g in named.gene_id:
    p = same_strand_partners(g)
    if p:
        mem = [g] + p
        fname = '+'.join(GENES[x]['symbol'] for x in mem)
        if not any(set(mem) == set(v) for v in FAMILIES.values()):
            FAMILIES['locus: ' + fname] = mem
fam_of = collections.defaultdict(list)
for f, mem in FAMILIES.items():
    for g in mem:
        fam_of[g].append(f)
named['families'] = ['; '.join(fam_of.get(g, [])) for g in named.gene_id]

# audit targets: named genes + all family members + genes drawn in coverage panels
targets = list(dict.fromkeys(list(named.gene_id) + [g for mem in FAMILIES.values() for g in mem] +
                             [gid_of(s) for c in COV for s in c[5]]))
named.to_csv(f'{OUT}/named_genes_curated.tsv', sep='\t', index=False)
json.dump({k: [GENES[g]['symbol'] for g in v] for k, v in FAMILIES.items()},
          open(f'{OUT}/families_definition.json', 'w'), indent=1)
print(f'{len(named)} named genes; {len(FAMILIES)} families; {len(targets)} audited loci', flush=True)

# per-target interval info
TINFO = {}
for g in targets:
    G = GENES[g]
    other = []
    for h in BYCHR[G['chrom']]:
        if h['gene_id'] == g or h['strand'] != G['strand'] or h['end'] < G['start'] or h['start'] > G['end']:
            continue
        other.extend(h['exons_union'])
    TINFO[g] = dict(chrom=G['chrom'], start0=G['start'] - 1, end=G['end'], strand=G['strand'],
                    exons=G['exons_union'], other=merge(other))


def iv_hit(blocks, ivs, starts):
    """any block overlaps any interval (ivs sorted, non-overlapping; starts = [s for s,e in ivs])"""
    for bs, be in blocks:
        i = bisect.bisect_right(starts, be - 1) - 1   # last interval starting before block end
        while i >= 0:
            s, e = ivs[i]
            if e <= bs:
                break
            if s < be and bs < e:
                return True
            i -= 1
    return False


TOL = 3   # bases an aligned block may extend past an exon edge and still count as exon-contained


def contained(blocks, ivs, starts):
    """every block lies inside one (merged) exon interval, +/- TOL"""
    for bs, be in blocks:
        i = bisect.bisect_right(starts, bs + TOL) - 1
        if i < 0 or not (ivs[i][0] - TOL <= bs and be <= ivs[i][1] + TOL):
            return False
    return True


def frag_strand(r):
    # reverse-stranded paired-end library: read 1 aligns antisense to the transcript
    if r.is_read2:
        return '-' if r.is_reverse else '+'
    return '+' if r.is_reverse else '-'


SENSE_EX, ANTI_EX, SENSE_IN, ANTI_IN, SHARED, DUP, PRIM = 1, 2, 4, 8, 16, 32, 64
CATS = ['sense_exonic_all', 'sense_exonic_unique', 'sense_exonic_multi', 'sense_exonic_dup', 'sense_exonic_unique_nodup',
        'sense_exonic_shared', 'sense_exonic_unique_specific', 'sense_exonic_multi_primary_here',
        'sense_contained_all', 'sense_contained_unique_specific',
        'antisense_exonic_all', 'antisense_exonic_unique', 'sense_intronic_all', 'sense_intronic_unique',
        'sense_exonic_mean_NM']


def audit_sample(sample_id):
    t0 = time.time()
    bam = pysam.AlignmentFile(f'{NF}/{sample_id}.markdup.sorted.bam')
    counts, names_sense = {}, {}
    need_names = set(g for mem in FAMILIES.values() for g in mem)
    nh_hist = {}
    for g in targets:
        I = TINFO[g]
        ex, exs = I['exons'], [s for s, _ in I['exons']]
        ot, ots = I['other'], [s for s, _ in I['other']]
        frags = {}
        for r in bam.fetch(I['chrom'], I['start0'], I['end']):
            if r.is_unmapped or r.is_supplementary or r.is_qcfail:
                continue
            blocks = [(s, e) for s, e in r.get_blocks() if s < I['end'] and e > I['start0']]
            if not blocks:
                continue          # spliced read spanning the gene without aligned bases in it
            sense = frag_strand(r) == I['strand']
            bits = 0
            if iv_hit(blocks, ex, exs):
                bits |= SENSE_EX if sense else ANTI_EX
                if sense and ot and iv_hit(blocks, ot, ots):
                    bits |= SHARED
            else:
                bits |= SENSE_IN if sense else ANTI_IN
            if r.is_duplicate:
                bits |= DUP
            if not r.is_secondary:
                bits |= PRIM
            nh = r.get_tag('NH'); nm = r.get_tag('NM') if r.has_tag('NM') else 0
            # transcript-compatible proxy: all aligned blocks inside the gene's exons and mate inside the gene span
            cont = bool(bits & SENSE_EX) and len(blocks) == len(r.get_blocks()) and contained(blocks, ex, exs) and \
                (not r.is_paired or r.mate_is_unmapped or (r.next_reference_name == I['chrom'] and
                                                           I['start0'] <= r.next_reference_start < I['end']))
            hi = r.get_tag('HI') if r.has_tag('HI') else 0
            f = frags.get(r.query_name)
            if f is None:
                frags[r.query_name] = [bits, nh, nm, {hi: cont}]
            else:
                f[0] |= bits; f[1] = max(f[1], nh); f[2] = max(f[2], nm)
                f[3][hi] = f[3].get(hi, True) and cont
        c = dict.fromkeys(CATS, 0)
        sn = set(); nms = []; hist = collections.Counter()
        for name, (bits, nh, nm, his) in frags.items():
            if bits & SENSE_EX:
                c['sense_exonic_all'] += 1; sn.add(name); nms.append(nm)
                if any(his.values()):
                    c['sense_contained_all'] += 1
                    if nh == 1 and not bits & SHARED:
                        c['sense_contained_unique_specific'] += 1
                hist[min(nh, 21)] += 1
                if nh == 1:
                    c['sense_exonic_unique'] += 1
                    if not bits & DUP:
                        c['sense_exonic_unique_nodup'] += 1
                    if not bits & SHARED:
                        c['sense_exonic_unique_specific'] += 1
                else:
                    c['sense_exonic_multi'] += 1
                    if bits & PRIM:
                        c['sense_exonic_multi_primary_here'] += 1
                if bits & DUP:
                    c['sense_exonic_dup'] += 1
                if bits & SHARED:
                    c['sense_exonic_shared'] += 1
            elif bits & SENSE_IN:
                c['sense_intronic_all'] += 1
                if nh == 1:
                    c['sense_intronic_unique'] += 1
            if bits & ANTI_EX:
                c['antisense_exonic_all'] += 1
                if nh == 1:
                    c['antisense_exonic_unique'] += 1
        c['sense_exonic_mean_NM'] = float(np.mean(nms)) if nms else np.nan
        counts[g] = c
        nh_hist[g] = dict(hist)
        if g in need_names:
            names_sense[g] = sn
    fam = {}
    for fname, mem in FAMILIES.items():
        u = set().union(*[names_sense[g] for g in mem])
        fam[fname] = len(u)
    # coverage
    cov = {}
    for title, chrom, s1, e1, strand, _ in COV:
        L = e1 - s1 + 1
        arr = {k: np.zeros(L + 1, dtype=np.int64) for k in ('sense_all', 'sense_unique', 'anti_all')}
        for r in bam.fetch(chrom, s1 - 1, e1):
            if r.is_unmapped or r.is_supplementary or r.is_qcfail:
                continue
            sense = frag_strand(r) == strand
            keys = ['sense_all'] + (['sense_unique'] if r.get_tag('NH') == 1 else []) if sense else ['anti_all']
            for bs, be in r.get_blocks():
                a, b = max(bs - (s1 - 1), 0), min(be - (s1 - 1), L)
                if a < b:
                    for k in keys:
                        arr[k][a] += 1; arr[k][b] -= 1
        cov[title] = {k: np.cumsum(v)[:L].astype(np.float32) for k, v in arr.items()}
    bam.close()
    print(f'  {sample_id} done in {time.time() - t0:.0f}s', flush=True)
    return sample_id, counts, fam, cov, nh_hist


def star_stats(sample_id):
    d = {}
    for line in open(f'{NF}/log/{sample_id}.Log.final.out'):
        if '|' in line:
            k, v = [x.strip() for x in line.split('|', 1)]
            d[k] = v
    return dict(input_pairs=int(d['Number of input reads']), unique_pairs=int(d['Uniquely mapped reads number']),
                multi_pairs=int(d['Number of reads mapped to multiple loci']),
                toomany_pairs=int(d['Number of reads mapped to too many loci']))


def _write_versions():
    import sys, importlib, platform
    with open(f'{OUT}/sessionInfo_13_reliability_02_bam_audit_python.txt', 'w') as fh:
        fh.write(f'python {sys.version}\nplatform {platform.platform()}\n')
        for m in ['pysam', 'numpy', 'pandas']:
            fh.write(f'{m} {importlib.import_module(m).__version__}\n')


if __name__ == '__main__':
    _write_versions()
    t0 = time.time()
    with mp.Pool(4) as pool:            # <= 4 concurrent BAM readers (slow NFS)
        results = pool.map(audit_sample, SAMPLE_IDS, chunksize=1)
    print(f'BAM audit finished in {time.time() - t0:.0f}s', flush=True)
    stats = pd.DataFrame({s: star_stats(s) for s in SAMPLE_IDS}).T
    stats.index.name = 'sample_id'
    stats.to_csv(f'{OUT}/star_mapping_stats.tsv', sep='\t')

    # long table gene x sample
    grp = {s: (g, l) for s, g, l in SAMPLES}
    recs = []
    for sid, counts, fam, cov, nhh in results:
        for g, c in counts.items():
            recs.append(dict(gene_id=g, symbol=GENES[g]['symbol'], sample_id=sid, group=grp[sid][0], label=grp[sid][1], **c))
    long = pd.DataFrame(recs)
    # Salmon (nf-core) gene counts and featureCounts (unique) for the same samples
    sal = pd.read_csv(f'{NF}/salmon.merged.gene_counts.tsv', sep='\t')
    colmap = {c: c[1:].replace('.', '-') for c in sal.columns if c.startswith('X')}
    sal = sal.rename(columns=colmap).set_index('gene_id')
    fc = pd.read_csv(f'{PROJ}/results/revision_checks_2026-10/featureCounts_unique_gene_counts.tsv', sep='\t').set_index('Geneid')
    long['salmon_numreads'] = [float(sal.at[g, s]) if g in sal.index else np.nan for g, s in zip(long.gene_id, long.sample_id)]
    long['featurecounts_unique'] = [int(fc.at[g, s]) for g, s in zip(long.gene_id, long.sample_id)]
    long.to_csv(f'{OUT}/read_evidence_per_sample.tsv', sep='\t', index=False)

    # family totals per sample
    frec = []
    for sid, counts, fam, cov, nhh in results:
        for f, n in fam.items():
            mem = FAMILIES[f]
            frec.append(dict(family=f, members='; '.join(GENES[g]['symbol'] for g in mem), sample_id=sid,
                             group=grp[sid][0], label=grp[sid][1], genome_fragments_distinct=n,
                             genome_unique_sum=sum(counts[g]['sense_exonic_unique'] for g in mem),
                             salmon_numreads_sum=float(sum(sal.at[g, sid] for g in mem if g in sal.index)),
                             featurecounts_unique_sum=int(sum(fc.at[g, sid] for g in mem))))
    pd.DataFrame(frec).to_csv(f'{OUT}/family_totals_per_sample.tsv', sep='\t', index=False)

    # NH distribution of sense-exonic fragments (pooled over samples)
    nhrec = collections.defaultdict(collections.Counter)
    for sid, counts, fam, cov, nhh in results:
        for g, h in nhh.items():
            nhrec[g].update(h)
    pd.DataFrame([dict(gene_id=g, symbol=GENES[g]['symbol'], **{f'NH_{k if k < 21 else ">20"}': v for k, v in sorted(h.items())})
                  for g, h in nhrec.items()]).fillna(0).to_csv(f'{TMP}/nh_distribution.tsv', sep='\t', index=False)

    # Salmon transcript-level counts for audited genes (to show which transcript carries the signal)
    txc = pd.read_csv(f'{NF}/salmon.merged.transcript_counts.tsv', sep='\t')
    txc = txc.rename(columns=colmap)
    txc = txc[txc.gene_id.isin(targets)].copy()
    txc.insert(1, 'symbol', [GENES[g]['symbol'] for g in txc.gene_id])
    txc.insert(2, 'transcript_name', [TX[t]['name'] for t in txc.tx])
    txc.insert(3, 'transcript_biotype', [TX[t]['biotype'] for t in txc.tx])
    txc.insert(4, 'transcript_length', [TX[t]['length'] for t in txc.tx])
    txc = txc[['tx', 'gene_id', 'symbol', 'transcript_name', 'transcript_biotype', 'transcript_length'] + SAMPLE_IDS]
    txc.to_csv(f'{OUT}/salmon_transcript_counts_audited_genes.tsv', sep='\t', index=False)

    # coverage arrays (temp) + window metadata
    covd = {sid: cov for sid, counts, fam, cov, nhh in results}
    pickle.dump({'cov': covd, 'windows': COV, 'stats': stats.to_dict('index'), 'samples': SAMPLES},
                open(f'{TMP}/coverage_arrays.pkl', 'wb'))
    print('done', flush=True)
