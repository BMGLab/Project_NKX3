#!/usr/bin/env bash
# 16_rootcause_14_figure.sh
# Supplementary-style figure: gene-level read counts of the manuscript's
# headline genes in one NKX3.1-basal library (24-004-7_S14) and one
# Control-basal library (24-004-1_S56) under
#   (i)  the nf-core STAR 2.7.9a -> Salmon output used in the manuscript,
#   (ii) the same pipeline with the STAR-unmapped reads removed from the transcriptome
#        BAM (identical read set to STAR 2.7.10b output; verified on 24-004-7_S14), and
#   (iii) featureCounts on uniquely mapped genome alignments.
# Input: tables/07_gene_counts_by_method.tsv. Output: figures/14_counts_by_method.{pdf,png}
set -euo pipefail
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
mkdir -p "$OUT/figures"
$PY -I - "$OUT" <<'EOF'
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
out = sys.argv[1]
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans'],
                     'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
T = pd.read_csv(f"{out}/tables/07_gene_counts_by_method.tsv", sep='\t')
panels = [('24-004-7_S14', 'NKX3.1-basal library (24-004-7_S14)',
           ['BCL11B', 'SSC4D', 'IGLL5', 'SRMS', 'DPEP3', 'FAM86MP', 'TPGS1', 'DIPK1B', 'NKX3-1']),
          ('24-004-1_S56', 'Control-basal library (24-004-1_S56)',
           ['NTNG1', 'CSH2', 'C3orf20', 'SOX2-OT', 'ZNF595', 'NRF1', 'LINC02626', 'BCL11B', 'NKX3-1'])]
series = [('nfcore_original', 'nf-core STAR 2.7.9a → Salmon (manuscript)', '#eb6834', 'o', 'full'),
          ('fix_genomeMappedOnly', 'Same, STAR-unmapped reads removed (= STAR 2.7.10b)', '#2a78d6', 'D', 'full'),
          ('featureCounts_unique', 'featureCounts, uniquely mapped pairs', '#1baf7a', 's', 'none')]
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
for ax, (s, title, genes) in zip(axes, panels):
    d = T[T['sample'] == s].set_index('gene')
    genes = [g for g in genes if g in d.index]
    y = np.arange(len(genes))[::-1]
    for k, (col, lab, c, m, fill) in enumerate(series):
        if col not in d.columns: continue
        v = np.log10(d.loc[genes, col].astype(float).values + 1)
        off = (k - 1) * 0.22
        ax.scatter(v, y + off, s=22, marker=m, facecolors=c if fill == 'full' else 'white', edgecolors=c,
                   linewidths=1.0, zorder=3, label=lab)
    for yy in y: ax.axhline(yy, color='#e6e6e6', lw=0.5, zorder=0)
    ax.set_yticks(y); ax.set_yticklabels(genes, fontstyle='italic')
    ticks = [0, 1, 10, 100, 1000, 10000]
    ax.set_xticks(np.log10(np.array(ticks) + 1)); ax.set_xticklabels([f"{t:,}" for t in ticks])
    ax.set_xlim(-0.15, np.log10(30001))
    ax.set_xlabel('Gene-level read count in this library (log scale; 0 = none)')
    ax.set_title(title, fontsize=7.5, loc='left', fontweight='bold')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    ax.tick_params(axis='y', length=0)
handles = [Line2D([], [], marker=m, linestyle='none', markersize=5, markerfacecolor=c if f == 'full' else 'white',
                  markeredgecolor=c, markeredgewidth=1.0, label=lab) for _, lab, c, m, f in series]
fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 1.0),
           handletextpad=0.3, columnspacing=1.2)
fig.tight_layout(rect=(0, 0, 1, 0.93), w_pad=2.0)
fig.savefig(f"{out}/figures/14_counts_by_method.pdf")
fig.savefig(f"{out}/figures/14_counts_by_method.png", dpi=300, facecolor='white')
print('saved')
EOF
