#!/usr/bin/env bash
# 16_rootcause_13_crossref_degs.sh
# Cross-reference the ORIGINAL pipeline DEGs (rna_downstreams_out_v15 dgea_* DESeq2
# tables; padj <= 0.05 & |LFC| >= 0.20 as in the manuscript) with the genes whose
# nf-core Salmon count, in the two re-processed libraries, is driven by reads that
# STAR itself left unmapped (07_artefact_genes_<sample>.tsv: count lost when the
# STAR-unmapped reads are removed), annotated with the construct class of those
# reads (09_per_gene_STARunmapped_reads_by_class.tsv; FP / bsd / NKX3.1-junction
# references only - control-vector sequences such as lacZ/bGH/HisMax leader are
# not in that reference set and appear as 'none'). A gene is flagged when >= 10 reads and
# >= 50% of its nf-core count disappear once the STAR-unmapped reads are removed AND the
# STAR-unmapped reads reaching it are >= 50% of that count.
set -euo pipefail
OUT=/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause
V15=/mnt/ssd1/Projects/NKX3_wd/results/rna_downstreams_out_v15
PY=/home/biolab/miniconda3/envs/nkx3_py/bin/python
$PY -I - "$OUT" "$V15" <<'EOF'
import sys, os, numpy as np, pandas as pd
out, v15 = sys.argv[1:3]
libs = [('24-004-7_S14', 'S14_NKX3.1-basal'), ('24-004-1_S56', 'S56_Control-basal')]
art, cls = {}, {}
for s, lab in libs:
    a = pd.read_csv(f"{out}/tables/07_artefact_genes_{s}.tsv", sep='\t').set_index('gene_id')
    art[lab] = a
    p = f"{out}/construct/{s}_repro_rawtrim/09_per_gene_STARunmapped_reads_by_class.tsv"
    cls[lab] = pd.read_csv(p, sep='\t').set_index('gene_name') if os.path.exists(p) else None
contrasts = {'NKX3.1 vs control (basal)': 'dgea_HM-NKX3.1-C_vs_HM-vector-C_DESeq2_results.tsv',
             'NKX3.1 vs control (inflammation)': 'dgea_HM-NKX3.1-CM_vs_HM-vector-CM_DESeq2_results.tsv',
             'control-basal minus control-inflammation (original sign)': 'dgea_HM-vector-C_vs_HM-vector-CM_DESeq2_results.tsv'}
rows = []
for cname, f in contrasts.items():
    d = pd.read_csv(f"{v15}/{f}", sep='\t')
    sig = d[(d.padj <= 0.05) & (d.log2FoldChange.abs() >= 0.20)]
    for _, r in sig.iterrows():
        row = dict(contrast=cname, gene=r.gene_symbol, gene_id=r.feature_id, log2FC_original=round(r.log2FoldChange, 2),
                   padj_original=f"{r.padj:.2e}", direction='up' if r.log2FoldChange > 0 else 'down')
        flags = []
        for s, lab in libs:
            a = art[lab]
            if r.feature_id in a.index:
                lost = float(a.loc[r.feature_id, 'reads_from_STAR_unmapped']); fr = float(a.loc[r.feature_id, 'frac_from_STAR_unmapped'])
            else:
                lost, fr = 0.0, 0.0
            row[f'nfcore_count_{lab}'] = None
            row[f'reads_from_STAR_unmapped_{lab}'] = round(lost, 1)
            row[f'frac_of_count_{lab}'] = round(fr, 3) if not np.isnan(fr) else 0.0
            c = cls[lab]
            row[f'construct_class_{lab}'] = ''
            if c is not None and r.gene_symbol in c.index:
                cc = c.loc[r.gene_symbol]
                parts = [f"{k}:{int(cc[k])}" for k in ['FP', 'bsd', 'neo', 'NKX3.1_construct_junctions', 'none'] if int(cc[k]) > 0]
                row[f'construct_class_{lab}'] = ';'.join(parts)
            # reads absent from the genome BAM that reach this gene in the toTranscriptome BAM
            tx_only = int(c.loc[r.gene_symbol, 'txonly_reads']) if (c is not None and r.gene_symbol in c.index) else 0
            nfc = float(a.loc[r.feature_id, 'nfcore_original']) if r.feature_id in a.index else 0.0
            row[f'STARunmapped_reads_reaching_gene_{lab}'] = tx_only
            # flag: >=10 reads and >=50% of the nf-core count are lost when STAR-unmapped reads are removed, AND
            # the STAR-unmapped reads reaching the gene amount to >=50% of its nf-core count (excludes pure EM
            # re-allocation between paralogues, e.g. SLX1A/SLX1B)
            if lost >= 10 and fr >= 0.5 and tx_only >= 0.5 * nfc: flags.append(lab)
        row['count_driven_by_STAR_unmapped_reads_in'] = ';'.join(flags)
        rows.append(row)
R = pd.DataFrame(rows).drop(columns=[c for c in pd.DataFrame(rows).columns if c.startswith('nfcore_count_')])
R.to_csv(f"{out}/tables/13_original_DEGs_vs_STARunmapped_artefact.tsv", sep='\t', index=False)
old = f"{out}/tables/13_original_DEGs_vs_construct_artefact.tsv"
if os.path.exists(old): os.remove(old)
pd.set_option('display.width', 250)
for cname in contrasts:
    x = R[R.contrast == cname]
    k = x[x.count_driven_by_STAR_unmapped_reads_in != '']
    print(f"{cname}: {len(x)} original DEGs; count driven (>=50%, >=10 reads) by STAR-unmapped reads in a re-processed library: {len(k)}")
    print(k[['gene', 'log2FC_original', 'padj_original', 'count_driven_by_STAR_unmapped_reads_in',
             'frac_of_count_S14_NKX3.1-basal', 'frac_of_count_S56_Control-basal']].to_string(index=False))
EOF
