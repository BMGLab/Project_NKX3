# Project_NKX3

Analysis code for the RNA-seq part of the NKX3.1 study in DU145 prostate cancer cells (Isel et al., Scientific Reports, in revision). DU145 cells were transfected with NKX3.1 or an empty vector and cultured with or without conditioned medium (CM). There are three libraries per group, 12 in total.

Raw reads are in the European Nucleotide Archive under [PRJEB111856](https://www.ebi.ac.uk/ena/browser/view/PRJEB111856). Processed results are in the paper's Supplementary Data 1–2 and Supplementary Table S1. This repository contains code only.

## Layout

```
original_submission/   downstream analysis of the original submission
revision/              analyses added during the revision, numbered in run order
envs/                  conda environment exports (R, Python, Salmon, STAR)
```

## Processing

Reads were processed with nf-core/rnaseq v3.14.0. In that pipeline, STAR (v2.7.9a) projects the alignments onto transcripts and Salmon quantifies these projections. Gene-level counts were imported with tximport (`original_submission/rna_downstreams_salmon.R`, with the protein-coding transcript-to-gene table from `build_tx2gene_protein_coding.R`).

STAR v2.7.9a also writes alignments of reads that fail its mapping filters to the transcriptome BAM; the defect was fixed in v2.7.10a. As a result, short partial alignments of reads from the transfected plasmids were counted against a few human genes. The `16_rootcause_*` scripts trace these reads. `60_plasmid_flags.py` marks the affected genes (†), which are then left out of the gene-set and principal component analyses.

## Revision scripts

| Script | What it does | Used for |
|---|---|---|
| `00_core_deseq2.R` | One DESeq2 model on all 12 libraries (~ NKX3.1 + CM + NKX3.1:CM), apeglm shrinkage, all contrasts, test against a 1.5-fold threshold | Supplementary Data 1 |
| `16_rootcause_*.sh` | Read-level trace of the transcriptome-BAM defect, plasmid read census, plasmid screen of all 12 libraries | Methods; Response Figure R1A |
| `60_plasmid_flags.py` | Flags genes whose counts come mainly from plasmid reads | Fig. 1, Table S1, Data 1 |
| `61_fig1BCD_original_style.R` | Signed z-statistic plots (treatment minus reference) | Fig. 1B–D |
| `13_reliability_*` | Builds the list of genes named in the original submission, with read-level checks | input to `62` |
| `62_salmon_named_genes.R` | baseMean, normalised counts, unshrunken and shrunken fold changes; paralogue and readthrough group totals | Supplementary Fig. 2, Table S1 |
| `63_gsea_salmon.R` | Preranked GSEA (Hallmark, C2, GO BP) on the DESeq2 Wald statistic; expression-trend and NKX3-1-removal sensitivity runs | Supplementary Data 2 |
| `69_pub_figures_salmon_fig2.R`, `65_pub_figures_salmon_compact_GO.R`, `64_pub_figures_salmon.R` | NES dot plots | Fig. 2, Fig. 3A–C, Supplementary Fig. 3 |
| `66_suppfig1_salmon.R` | PCA and sample correlation | Supplementary Fig. 1C–D |
| `67_camera_fry_salmon.R` | Gene-set tests that account for inter-gene correlation (camera, fry) | Limitations paragraph |
| `68_package_v3_tables.py` | Supplementary Table S1 and Data 1–2 workbooks | Supplementary tables |
| `70_response_fig_plasmid.R` | Plasmid reads and plasmid-derived genes per library | Response Figure R1 |
| `14_modules*.R`, `43_fig3_alternative.R`, `50_verify_fig3_*` | Tests of the former co-expression modules, including the original display procedure on mislabelled groupings | Response to Reviewer 1, point 5; Response Figure R2 |
| `15_tcga*` | TCGA-PRAD co-expression and public NKX3.1-overexpression data in RWPE1 cells (GSE81439, GSE115335) | Response Figures R3–R4 |

## Re-running

The scripts contain absolute paths of the analysis server (`/mnt/...`). Re-running them requires the raw reads from ENA, the nf-core/rnaseq output and path edits. The R scripts were run in `envs/nkx3_r.yml` and the Python scripts in `envs/nkx3_py.yml`. Most R scripts also write their `sessionInfo()` next to their outputs.

## License

MIT (see `LICENSE`).
