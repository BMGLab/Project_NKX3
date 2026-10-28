#!/usr/bin/env python
"""Reliability task, step 5: strand-specific coverage at the audited loci (supplementary figure candidate).

For each locus: per-group mean sense-strand coverage (both mates, spliced blocks only) from the nf-core STAR
genome BAMs, normalised per sample to reads per million uniquely mapped read pairs (STAR Log.final.out),
then averaged over the 3 replicates. Fill = all alignments (unique + multi-mapped); dark line = uniquely
mapped alignments only (drawn where >=20% of fragments are multi-mapped); grey row = antisense coverage
(mean of all 12 samples). Dashed line = depth implied by the original nf-core Salmon NumReads if those
fragments existed (sum over transcripts of NumReads x mean mapped length / transcript length, laid over
each transcript's exons). Long introns / intergenic gaps (>3% of the window, windows >12 kb) are shortened
to a fixed width and marked '//'.
Coverage arrays come from 13_reliability_02_bam_audit.py. The nf-core bigWigs are unscaled raw depth that
includes multi-mapped and duplicate alignments; reverse.bigWig carries minus-strand genes (verified with
NKX3-1) and equals the 'all alignments' arrays used here.
Run: /home/biolab/miniconda3/envs/nkx3_py/bin/python -I 13_reliability_05_coverage_figure.py
"""
import os, pickle, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors as mcolors
from matplotlib.patches import Rectangle

OUT = '/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/reliability'
FIG = f'{OUT}/coverage_figure'
TMP = '/mnt/ssd1/Projects/NKX3_wd/work/reliability'
NF = '/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon'
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'axes.linewidth': 0.5, 'svg.fonttype': 'none'})
GROUPS = [('vector_ctrl', 'Control-basal', '#909FC9'), ('vector_CM', 'Control-inflammation', '#E28DBB'),
          ('NKX3.1_ctrl', 'NKX3.1-basal', '#57C5D1'), ('NKX3.1_CM', 'NKX3.1-inflammation', '#F5906A')]
INK, MUTED, ANTI, OTHER = '#222222', '#666666', '#B5B5B5', '#8A8A8A'
NBINS = 600
DASH = (0, (2.2, 1.3))


def darker(hexc, f=0.62):
    r, g, b = mcolors.to_rgb(hexc)
    return (r * f, g * f, b * f)


C = pickle.load(open(f'{TMP}/coverage_arrays.pkl', 'rb'))
M = pickle.load(open(f'{TMP}/gtf_models.pkl', 'rb'))
GENES, TX = M['genes'], M['tx']
BYSYM = {g['symbol']: gid for gid, g in GENES.items()}
samples = C['samples']
uniq = {s: C['stats'][s]['unique_pairs'] for s, _, _ in samples}
mlen = {}
for s, _, _ in samples:
    for line in open(f'{NF}/log/{s}.Log.final.out'):
        if 'Average mapped length' in line:
            mlen[s] = float(line.split('|')[1])
ev = pd.read_csv(f'{OUT}/read_evidence_per_sample.tsv', sep='\t')
txc = pd.read_csv(f'{OUT}/salmon_transcript_counts_audited_genes.tsv', sep='\t')
full = pd.read_csv(f'{OUT}/named_genes_reliability_full.tsv', sep='\t')
cls_of = dict(zip(full.symbol, full.read_evidence_class))
mfrac_of = dict(zip(full.symbol, full.multi_mapped_fraction))
TARGET = {'GH/CSH cluster (CSH2)': 'CSH2'}


def gmean(sym, col, grp):
    return float(ev[(ev.symbol == sym) & (ev.group == grp)][col].mean())


def implied_profile(gene_syms, strand, s1, e1, grp):
    """group mean of sum_tx NumReads * mapped length / transcript length (per million unique pairs), per base"""
    L = e1 - s1 + 1
    prof = np.zeros(L)
    sids = [s for s, g, _ in samples if g == grp]
    for sym in gene_syms:
        gid = BYSYM[sym]
        if GENES[gid]['strand'] != strand:
            continue
        for _, r in txc[txc.gene_id == gid].iterrows():
            depth = np.mean([r[s] * mlen[s] / r['transcript_length'] / (uniq[s] / 1e6) for s in sids])
            if depth <= 0:
                continue
            for a, b in TX[r['tx']]['exons']:
                a0, b0 = max(a - (s1 - 1), 0), min(b - (s1 - 1), L)
                if a0 < b0:
                    prof[a0:b0] += depth
    return prof


def merge_iv(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def make_segments(win):
    """Display coordinates: long introns / intergenic gaps shortened to a fixed width (marked //)."""
    title, chrom, s1, e1, strand, gene_syms = win
    L = e1 - s1 + 1
    if L <= 6000:
        return [(0, L, 'keep', 0, L)], L
    pad = 40
    iv = []
    for sym in gene_syms:
        for a, b in GENES[BYSYM[sym]]['exons_union']:
            a0, b0 = max(a - (s1 - 1) - pad, 0), min(b - (s1 - 1) + pad, L)
            if a0 < b0:
                iv.append([a0, b0])
    iv = merge_iv(iv)
    T = int(0.03 * L)
    segs, pos = [], 0
    for a, b in iv:
        if a - pos > T:
            segs.append((pos, a, 'gap'))
        elif a > pos:
            segs.append((pos, a, 'keep'))
        segs.append((a, b, 'keep')); pos = b
    if L - pos > T:
        segs.append((pos, L, 'gap'))
    elif L > pos:
        segs.append((pos, L, 'keep'))
    kept = sum(b - a for a, b, k in segs if k == 'keep')
    D = max(50, int(0.035 * kept))
    out, d = [], 0
    for a, b, k in segs:
        w = (b - a) if k == 'keep' else D
        out.append((a, b, k, d, d + w)); d += w
    return out, d


def fwd(segs, g):
    for a, b, k, d0, d1 in segs:
        if a <= g <= b:
            return d0 + (g - a) * (d1 - d0) / max(b - a, 1)
    return segs[-1][4] if g > segs[-1][1] else 0.0


def inv(segs, x):
    for a, b, k, d0, d1 in segs:
        if d0 <= x <= d1:
            return a + (x - d0) * (b - a) / max(d1 - d0, 1e-9)
    return segs[-1][1]


def binned(arr, segs, dtot, nb=NBINS):
    xs = np.linspace(0, dtot, nb + 1)
    g = [int(round(inv(segs, x))) for x in xs]
    vals = np.array([arr[a:max(b, a + 1)].mean() for a, b in zip(g[:-1], g[1:])])
    return vals, (xs[:-1] + xs[1:]) / 2


def nice_scale(L):
    for v in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]:
        if v >= L / 6:
            return v
    return 100000


def fmt_n(x):
    return f'{x:,.0f}' if x >= 10 or x == 0 else f'{x:.1f}'


def fmt_r(v):
    return f'{v:,.0f}' if v >= 10 else (f'{v:.1f}' if v >= 1 else f'{v:.2g}')


def style_track(ax, xmax, ymax):
    ax.set_xlim(0, xmax); ax.set_ylim(0, ymax); ax.set_xticks([]); ax.set_yticks([])
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color('#9a9a9a'); ax.spines['bottom'].set_linewidth(0.4)


def plateaus(prof, segs):
    """runs of constant non-zero implied depth (genomic) -> (value, display start, display end)"""
    out = []
    nz = prof > 0
    if not nz.any():
        return out
    change = np.flatnonzero(np.diff(prof) != 0) + 1
    starts = np.concatenate([[0], change]); ends = np.concatenate([change, [len(prof)]])
    for a, b in zip(starts, ends):
        if prof[a] > 0:
            out.append((float(prof[a]), fwd(segs, a), fwd(segs, b)))
    return out


def assign_lanes(gene_syms, target, segs, dtot, s1, L):
    """target on the upper lane; other genes greedily on the first lane where gene span + label do not collide"""
    lanes_y = [0.86, 0.34]
    placed = {0: [], 1: []}
    out = {}
    order = [target] + [g for g in gene_syms if g != target]
    for sym in order:
        g = GENES[BYSYM[sym]]
        d0 = fwd(segs, max(g['start'] - s1, 0)); d1 = fwd(segs, min(g['end'] - s1 + 1, L))
        lw = 0.026 * len(sym) * dtot            # approximate label width in display units
        ext = (min(d0, (d0 + d1) / 2 - lw / 2), max(d1, (d0 + d1) / 2 + lw / 2))
        for li in ([0] if sym == target else [0, 1]):
            if all(ext[1] < e0 or ext[0] > e1 for e0, e1 in placed[li]):
                break
        else:
            li = 1
        placed[li].append(ext); out[sym] = (lanes_y[li], d0, d1)
    return out


def draw_panel(fig, spec, win):
    title, chrom, s1, e1, strand, gene_syms = win
    target = TARGET.get(title, title)
    L = e1 - s1 + 1
    segs, dtot = make_segments(win)
    gs = gridspec.GridSpecFromSubplotSpec(8, 3, subplot_spec=spec, width_ratios=[0.5, 1, 0.4],
                                          height_ratios=[1.55, 1, 1, 1, 1, 0.75, 1.9, 0.62], hspace=0.12, wspace=0.03)
    tr, ut, im = {}, {}, {}
    for gk, _, _ in GROUPS:
        sids = [s for s, g, _ in samples if g == gk]
        a = np.mean([C['cov'][s][title]['sense_all'] / (uniq[s] / 1e6) for s in sids], axis=0)
        u = np.mean([C['cov'][s][title]['sense_unique'] / (uniq[s] / 1e6) for s in sids], axis=0)
        tr[gk], x = binned(a, segs, dtot)
        ut[gk], _ = binned(u, segs, dtot)
        prof = implied_profile(gene_syms, strand, s1, e1, gk)
        im[gk] = plateaus(prof, segs)
    anti = np.mean([C['cov'][s][title]['anti_all'] / (uniq[s] / 1e6) for s, _, _ in samples], axis=0)
    anti_b, _ = binned(anti, segs, dtot)
    ymax = max(max(t.max() for t in tr.values()), max([v for p in im.values() for v, _, _ in p] + [0.05])) * 1.1
    mf = mfrac_of.get(target, 0) or 0
    # ---- group tracks (created first so the title can be aligned to them)
    axes_tr, axes_num = [], []
    for i, (gk, glab, col) in enumerate(GROUPS):
        axl = fig.add_subplot(gs[1 + i, 0]); axl.axis('off')
        axl.add_patch(Rectangle((0.0, 0.25), 0.065, 0.5, color=col, transform=axl.transAxes, lw=0))
        axl.text(0.095, 0.5, glab, fontsize=5.7, color=INK, va='center', ha='left', transform=axl.transAxes)
        ax = fig.add_subplot(gs[1 + i, 1]); axes_tr.append(ax)
        ax.fill_between(x, 0, tr[gk], step='mid', color=col, lw=0)
        ax.step(x, tr[gk], where='mid', color=darker(col), lw=0.35)
        if mf >= 0.2:
            ax.step(x, ut[gk], where='mid', color=INK, lw=0.55)
        for v, xa, xb in im[gk]:
            ax.hlines(v, xa, xb, color=INK, lw=0.8, linestyles=DASH)
        style_track(ax, dtot, ymax)
        axn = fig.add_subplot(gs[1 + i, 2]); axn.axis('off'); axes_num.append(axn)
        s_ = gmean(target, 'salmon_numreads', gk); g_ = gmean(target, 'sense_exonic_all', gk)
        txt = f'{fmt_n(s_)} | {fmt_n(g_)}'
        if mf >= 0.2:
            txt += f' (u {fmt_n(gmean(target, "sense_exonic_unique", gk))})'
        axn.text(1.0, 0.5, txt, fontsize=5.7, color=INK, va='center', ha='right', transform=axn.transAxes)
    # ---- title row
    axt = fig.add_subplot(gs[0, :]); axt.axis('off')
    sub = 'positive control (supported)' if target == 'NKX3-1' else cls_of.get(target, '')
    axt.text(0.0, 1.0, title, fontsize=8, fontweight='bold', fontstyle='italic' if title == target else 'normal',
             color=INK, va='top', ha='left', transform=axt.transAxes)
    axt.text(0.0, 0.02, sub, fontsize=6.2, color=INK, va='bottom', ha='left', transform=axt.transAxes)
    pt, pn, pa = axes_tr[0].get_position(), axes_num[0].get_position(), axt.get_position()
    fig.text(pt.x1, pa.y0 + 0.02 * pa.height, f'0-{fmt_r(ymax)} RPM', fontsize=5.6, color=MUTED, va='bottom', ha='right')
    fig.text(pn.x1, pa.y0 + 0.02 * pa.height, 'Salmon | genome', fontsize=5.6, color=MUTED, va='bottom', ha='right')
    # ---- antisense
    axl = fig.add_subplot(gs[5, 0]); axl.axis('off')
    axl.text(0.095, 0.5, 'antisense (all)', fontsize=5.5, color=MUTED, va='center', ha='left', transform=axl.transAxes)
    ax = fig.add_subplot(gs[5, 1])
    amax = max(anti_b.max(), 0.05) * 1.15
    ax.fill_between(x, 0, anti_b, step='mid', color=ANTI, lw=0)
    style_track(ax, dtot, amax)
    axn = fig.add_subplot(gs[5, 2]); axn.axis('off')
    axn.text(1.0, 0.5, f'0-{fmt_r(amax)} RPM', fontsize=5.2, color=MUTED, va='center', ha='right', transform=axn.transAxes)
    # ---- gene model (display coordinates)
    axg = fig.add_subplot(gs[6, 1]); axg.set_xlim(0, dtot); axg.set_ylim(0, 1); axg.axis('off')
    lanes = assign_lanes(gene_syms, target, segs, dtot, s1, L)
    for sym in gene_syms:
        g = GENES[BYSYM[sym]]
        lane, d0, d1 = lanes[sym]
        colg = INK if sym == target else OTHER
        axg.plot([d0, d1], [lane, lane], color=colg, lw=0.6, solid_capstyle='butt')
        for a, b in g['exons_union']:
            a0, b0 = max(a - (s1 - 1), 0), min(b - (s1 - 1), L)
            if a0 < b0:
                xa, xb = fwd(segs, a0), fwd(segs, b0)
                axg.add_patch(Rectangle((xa, lane - 0.1), max(xb - xa, dtot * 0.003), 0.2, color=colg, lw=0))
        n_ch = max(2, int(7 * (d1 - d0) / dtot))
        for xx in np.linspace(d0 + (d1 - d0) * 0.08, d1 - (d1 - d0) * 0.08, n_ch):
            dx = dtot * 0.007 * (1 if g['strand'] == '+' else -1)
            axg.plot([xx - dx, xx, xx - dx], [lane + 0.06, lane, lane - 0.06], color=colg, lw=0.4)
        axg.text((d0 + d1) / 2, lane - 0.12, sym, fontsize=5.6, fontstyle='italic', color=colg, ha='center', va='top')
    for a, b, k, d0, d1 in segs:   # shortened gaps
        if k == 'gap':
            xm = (d0 + d1) / 2
            for off in (-0.005, 0.005):
                axg.plot([xm + (off - 0.004) * dtot, xm + (off + 0.004) * dtot], [0.9, 1.0], color=MUTED, lw=0.5,
                         clip_on=False)
    # ---- footer: compact locus (left), shortened-intron note and scale bar (track column)
    compressed = any(k == 'gap' for a, b, k, d0, d1 in segs)
    axl = fig.add_subplot(gs[7, 0]); axl.axis('off')
    axl.text(0.0, 0.5, f'chr{chrom}:{s1 / 1e6:.3f}-{e1 / 1e6:.3f} Mb ({strand})', fontsize=5.2, color=MUTED,
             va='center', ha='left', transform=axl.transAxes)
    axf = fig.add_subplot(gs[7, 1]); axf.axis('off'); axf.set_xlim(0, dtot); axf.set_ylim(0, 1)
    kept = sum(b - a for a, b, k, _, _ in segs if k == 'keep')
    sb = nice_scale(kept)
    axf.plot([dtot - sb, dtot], [0.5, 0.5], color=INK, lw=0.8, solid_capstyle='butt')
    axf.text(dtot - sb - dtot * 0.02, 0.5, f'{sb // 1000} kb' if sb >= 1000 else f'{sb} bp', fontsize=5.2, color=INK,
             ha='right', va='center')


# figure_qa (2026-10-21): key text names the Salmon count shown (nf-core NumReads, not the length-scaled
# DESeq2 input plotted in Fig 1E-G diagnostic / verify_reads); same number of key lines; presentation only
def legend_panel(fig, spec):
    ax = fig.add_subplot(spec); ax.axis('off'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    y = 0.98
    ax.text(0.0, y, 'Key', fontsize=7, fontweight='bold', color=INK, va='top'); y -= 0.12
    for j, (_, glab, col) in enumerate(GROUPS):
        xx = 0.0 if j % 2 == 0 else 0.5
        ax.add_patch(Rectangle((xx, y - 0.04), 0.04, 0.04, color=col, lw=0))
        ax.text(xx + 0.055, y - 0.02, glab, fontsize=6, color=INK, va='center')
        if j % 2 == 1:
            y -= 0.075
    y -= 0.03
    items = [('fill', 'Coverage, sense strand, all alignments (unique + multi-mapped), reads per million uniquely mapped pairs (RPM); mean of 3 replicates; one scale per panel'),
             ('line', 'Uniquely mapped alignments only (where >= 20% are multi-mapped)'),
             ('dash', 'Depth implied by the original Salmon NumReads'),
             ('anti', 'Antisense coverage, mean of 12 samples (own scale)'),
             ('model', 'Ensembl 99 exons: target gene black, others grey; // shortened intron'),
             ('num', 'Mean fragments per sample: Salmon NumReads (nf-core, not length-scaled) | genome (sense, exon-overlapping, unique + multi-mapped; u = uniquely mapped)')]
    lh = 0.066
    for kind, txt in items:
        lines = textwrap.wrap(txt, 76)
        yc = y - lh * 0.5
        if kind == 'fill':
            ax.add_patch(Rectangle((0.0, yc - 0.02), 0.045, 0.04, color='#c4c4c4', lw=0))
        elif kind == 'line':
            ax.plot([0.0, 0.045], [yc, yc], color=INK, lw=0.6)
        elif kind == 'dash':
            ax.plot([0.0, 0.045], [yc, yc], color=INK, lw=0.8, ls=DASH)
        elif kind == 'anti':
            ax.add_patch(Rectangle((0.0, yc - 0.02), 0.045, 0.04, color=ANTI, lw=0))
        elif kind == 'model':
            ax.plot([0.0, 0.045], [yc, yc], color=INK, lw=0.6)
            ax.add_patch(Rectangle((0.01, yc - 0.016), 0.02, 0.032, color=INK, lw=0))
        elif kind == 'num':
            ax.text(0.0, yc, '12 | 3', fontsize=5.4, color=INK, va='center')
        ax.text(0.07, y, '\n'.join(lines), fontsize=5.8, color=INK, va='top', linespacing=1.12)
        y -= lh * len(lines) + 0.02


W = {w[0]: w for w in C['windows']}
ORDER = ['NKX3-1', 'BCL11B', 'SSC4D', 'TPGS1', 'DIPK1B', 'SRMS', 'IGLL5', 'NTNG1', 'GH/CSH cluster (CSH2)', 'SLX1A', 'SLX1B']

def _write_versions():
    import sys, importlib, platform
    with open(f'{OUT}/sessionInfo_13_reliability_05_coverage_figure_python.txt', 'w') as fh:
        fh.write(f'python {sys.version}\nplatform {platform.platform()}\n')
        for m in ['matplotlib', 'numpy', 'pandas']:
            fh.write(f'{m} {importlib.import_module(m).__version__}\n')


if __name__ == '__main__':
    _write_versions()
    for t in ORDER:                                    # single panels for the PowerPoint deck
        fig = plt.figure(figsize=(3.45, 1.9))
        gsp = gridspec.GridSpec(1, 1, left=0.015, right=0.985, top=0.97, bottom=0.02)
        draw_panel(fig, gsp[0], W[t])
        fn = t.replace('/', '_').replace(' ', '_').replace('(', '').replace(')', '')
        fig.savefig(f'{FIG}/coverage_{fn}.pdf'); fig.savefig(f'{FIG}/coverage_{fn}.png', dpi=300)
        plt.close(fig)
    fig = plt.figure(figsize=(7.0, 9.6))               # composite: 2 columns x 6 rows, key in the last slot
    gsp = gridspec.GridSpec(6, 2, left=0.012, right=0.988, top=0.992, bottom=0.006, hspace=0.16, wspace=0.07)
    for k, t in enumerate(ORDER):
        draw_panel(fig, gsp[k // 2, k % 2], W[t])
    legend_panel(fig, gsp[5, 1])
    fig.savefig(f'{OUT}/SuppFig_coverage_named_loci.pdf')
    fig.savefig(f'{OUT}/SuppFig_coverage_named_loci.png', dpi=300)
    plt.close(fig)
    print('saved')
