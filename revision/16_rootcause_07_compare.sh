#!/usr/bin/env bash
# 16_rootcause_07_compare.sh
# Gene-level comparison (sum of Salmon NumReads over transcripts) of every
# quantification run for the requested samples vs nf-core and featureCounts.
# Writes:
#   tables/07_gene_counts_by_method.tsv      (listed genes x method x sample)
#   tables/07_global_agreement.tsv           (all genes; per method/sample)
#   tables/07_artefact_genes_<sample>.tsv    (genes losing >= 10 reads when the
#                                             STAR-unmapped reads are removed)
# usage: 16_rootcause_07_compare.sh <sample_id> [<sample_id> ...]
set -euo pipefail
WORK=/mnt/ssd1/Projects/NKX3_wd/work/rootcause
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
mkdir -p "$OUT/tables"
$PY -I - "$OUT" "$WORK/ref/tx2gene_name.tsv" "$@" <<'EOF'
import sys, os, gzip, numpy as np, pandas as pd
out, tx2g_path, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
NF = '/mnt/lun3/processed/BilgeDebelec/BAP-RNAseq/results.first/star_salmon'
FC = '/mnt/ssd1/Projects/NKX3_wd/results/revision_checks_2026-10/featureCounts_unique_gene_counts.tsv'
tx2g = pd.read_csv(tx2g_path, sep='\t', header=None, names=['tx', 'gene_id', 'gene_name'])
gname = tx2g.drop_duplicates('gene_id').set_index('gene_id')['gene_name']
fc = pd.read_csv(FC, sep='\t').set_index('Geneid')
methods = [  # (label, path template, description)
  ('nfcore_original', NF + '/{s}/quant.sf', 'nf-core 3.14.0 output used in the manuscript (STAR 2.7.9a Singleend -> Salmon alignment mode)'),
  ('repro_rawtrim', out + '/runs/repro_rawtrim/{s}/quant.sf.gz', 'our re-run: raw FASTQ, Trim Galore-equivalent trimming, STAR 2.7.9a Singleend -> Salmon (nf-core options)'),
  ('repro_bamFASTQ', out + '/runs/repro_bamFASTQ/{s}/quant.sf.gz', 'our re-run from FASTQ regenerated from the genome BAM (lacks STAR-unmapped reads)'),
  ('fix_genomeMappedOnly', out + '/runs/repro_rawtrim_genomeMappedOnly/{s}/quant.sf.gz', 'repro_rawtrim toTranscriptome BAM restricted to reads present in the genome BAM -> Salmon'),
  ('fix_STAR2.7.10b', out + '/runs/star2710b_rawtrim/{s}/quant.sf.gz', 'STAR 2.7.10b (bug fixed upstream), otherwise identical options -> Salmon'),
  ('fix_banIndelSoftclipSingleend', out + '/runs/banISS_rawtrim/{s}/quant.sf.gz', 'STAR 2.7.9a with --quantTranscriptomeBan IndelSoftclipSingleend (STAR default) -> Salmon'),
  ('fix_SalmonSA', out + '/runs/salmonSA_rawtrim/{s}/quant.sf.gz', 'Salmon 1.10.1 selective alignment (k31, no decoys), no STAR'),
]
genes = ['BCL11B', 'NTNG1', 'CSH2', 'SSC4D', 'SRMS', 'IGLL5', 'TPGS1', 'DIPK1B', 'NKX3-1',
         'SLX1A', 'FAM86MP', 'EHD1', 'MSH3', 'DPEP3', 'C3orf20', 'SOX2-OT', 'ZNF595', 'NRF1', 'LINC02626']
name2id = {}
for gid, n in gname.items():
    if n in genes and n not in name2id: name2id[n] = gid
def gene_counts(path):
    q = pd.read_csv(path, sep='\t')
    q = q.merge(tx2g[['tx', 'gene_id']], left_on='Name', right_on='tx', how='left')
    return q.groupby('gene_id')['NumReads'].sum()
long, glob = [], []
for s in samples:
    gc = {}
    for lab, tpl, desc in methods:
        p = tpl.format(s=s)
        if os.path.exists(p): gc[lab] = gene_counts(p)
    df = pd.DataFrame(gc).fillna(0)
    df['featureCounts_unique'] = fc[s].reindex(df.index).fillna(0) if s in fc.columns else np.nan
    for n in genes:
        gid = name2id.get(n)
        if gid is None: continue
        for lab in df.columns:
            long.append(dict(sample=s, gene=n, gene_id=gid, method=lab, count=round(float(df.loc[gid, lab]) if gid in df.index else 0.0, 1)))
    # global agreement vs featureCounts and vs nf-core
    for lab in df.columns:
        x = np.log2(df[lab] + 1); y = np.log2(df['featureCounts_unique'] + 1); z = np.log2(df['nfcore_original'] + 1)
        expressed = (df[lab] >= 10) | (df['featureCounts_unique'] >= 10)
        glob.append(dict(sample=s, method=lab, total_reads=round(df[lab].sum()),
            pearson_log2_vs_featureCounts=round(np.corrcoef(x[expressed], y[expressed])[0, 1], 4),
            pearson_log2_vs_nfcore=round(np.corrcoef(x[expressed], z[expressed])[0, 1], 4),
            genes_ge20_with_featureCounts_le2=int(((df[lab] >= 20) & (df['featureCounts_unique'] <= 2)).sum())))
    if 'repro_rawtrim' in df and 'fix_genomeMappedOnly' in df:
        d = df.copy(); d['gene_name'] = gname.reindex(d.index).values
        d['reads_from_STAR_unmapped'] = d['repro_rawtrim'] - d['fix_genomeMappedOnly']
        d['frac_from_STAR_unmapped'] = (d['reads_from_STAR_unmapped'] / d['repro_rawtrim'].replace(0, np.nan)).round(3)
        a = d[d['reads_from_STAR_unmapped'] >= 10].sort_values('reads_from_STAR_unmapped', ascending=False)
        a.round(1).to_csv(f"{out}/tables/07_artefact_genes_{s}.tsv", sep='\t', index_label='gene_id')
        print(f"{s}: genes losing >=10 reads when STAR-unmapped reads are removed: {len(a)}; "
              f">=50% of count: {int((a['frac_from_STAR_unmapped'] >= 0.5).sum())}")
L = pd.DataFrame(long)
W = L.pivot_table(index=['sample', 'gene', 'gene_id'], columns='method', values='count').reset_index()
order = ['sample', 'gene', 'gene_id', 'featureCounts_unique'] + [m[0] for m in methods if m[0] in W.columns]
W[order].to_csv(f"{out}/tables/07_gene_counts_by_method.tsv", sep='\t', index=False)
pd.DataFrame(glob).to_csv(f"{out}/tables/07_global_agreement.tsv", sep='\t', index=False)
pd.DataFrame([dict(method=m[0], description=m[2]) for m in methods] +
             [dict(method='featureCounts_unique', description='featureCounts 2.1.1 on the nf-core genome BAMs: uniquely mapped read pairs, reverse-stranded, Ensembl 99 exons')]
            ).to_csv(f"{out}/tables/07_method_legend.tsv", sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20)
print(W[order].to_string(index=False))
print(pd.DataFrame(glob).to_string(index=False))
EOF
