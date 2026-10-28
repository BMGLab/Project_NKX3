"""60_plasmid_flags.py - genes whose Salmon counts (original nf-core 3.14.0 run) come mainly from plasmid reads.

Rule, applied to the two libraries traced at read level in the rootcause task (one per plasmid group):
  a gene is flagged when >= 10 plasmid-derived reads were assigned to it AND excluding the reads that STAR
  did not align to the genome removes >= 50% of its Salmon count.
  - 24-004-7_S14 (NKX3.1 plasmid): genome-wide k-mer classification of STAR-unmapped reads
    (fluorescent protein, bsd, NKX3.1-construct junctions; construct/09_per_gene_STARunmapped_reads_by_class.tsv);
    count change = nf-core Salmon vs STAR 2.7.10b re-run (tables/07_artefact_genes_*.tsv).
  - 24-004-1_S56 (control plasmid): targeted trace of the genes reported in the original analysis plus
    neighbours (trace/*/traced_reads.tsv; plasmid motifs or lacZ/bGH/HisMax origin by sequence);
    count change = nf-core Salmon vs transcriptome BAM restricted to genome-mapped reads.
NKX3-1 is not flagged: it is the transfected gene (about 20% of its Salmon count comes from construct reads).
Output: results/revision_2026/salmon_v3/plasmid_flagged_genes.tsv
"""
import pandas as pd
RC = '/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/rootcause'
OUT = '/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/salmon_v3/plasmid_flagged_genes.tsv'

# NKX3.1 plasmid library
c = pd.read_csv(f'{RC}/construct/24-004-7_S14_repro_rawtrim/09_per_gene_STARunmapped_reads_by_class.tsv', sep='\t')
c['plasmid_reads'] = c[['FP', 'bsd', 'neo', 'NKX3.1_construct_junctions']].sum(axis=1)
a = pd.read_csv(f'{RC}/tables/07_artefact_genes_24-004-7_S14.tsv', sep='\t')
s14 = c.merge(a[['gene_id', 'gene_name', 'nfcore_original', 'fix_STAR2.7.10b']], on='gene_name', how='left')
s14 = s14.rename(columns={'fix_STAR2.7.10b': 'salmon_without_unmapped'})
s14['library'] = '24-004-7_S14 (NKX3.1 plasmid)'

# control plasmid library (targeted trace)
t = pd.read_csv(f'{RC}/trace/24-004-1_S56_repro_rawtrim/traced_reads.tsv', sep='\t', low_memory=False)
t = t.drop_duplicates(['read', 'gene'])
t = t[~t['in_genome_BAM'].astype(str).eq('True')]
s56 = t.groupby('gene').size().rename('plasmid_reads').reset_index().rename(columns={'gene': 'gene_name'})
b = pd.read_csv(f'{RC}/tables/07_artefact_genes_24-004-1_S56.tsv', sep='\t')
fixcol = [x for x in b.columns if x.startswith('fix_')][0]
s56 = s56.merge(b[['gene_id', 'gene_name', 'nfcore_original', fixcol]], on='gene_name', how='left').rename(columns={fixcol: 'salmon_without_unmapped'})
s56['library'] = '24-004-1_S56 (control plasmid)'

d = pd.concat([s14[['library', 'gene_id', 'gene_name', 'nfcore_original', 'salmon_without_unmapped', 'plasmid_reads']],
               s56[['library', 'gene_id', 'gene_name', 'nfcore_original', 'salmon_without_unmapped', 'plasmid_reads']]])
d['frac_count_removed'] = 1 - d['salmon_without_unmapped'] / d['nfcore_original']
f = d[(d.plasmid_reads >= 10) & (d.frac_count_removed >= 0.5) & (d.gene_name != 'NKX3-1')]
f = f.sort_values(['library', 'nfcore_original'], ascending=[False, False]).round(3)
f.to_csv(OUT, sep='\t', index=False)
print(f.to_string(index=False))
print('fix column used for S56:', fixcol)
