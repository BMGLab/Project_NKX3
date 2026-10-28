#!/usr/bin/env python
"""15_tcga_tables.py -- assemble the supplementary workbook and a human-readable provenance file for the
TCGA-PRAD / RWPE1 independent-dataset check. Reads only this task's own outputs.

Run (last step of 15_tcga_run_all.sh):
  /home/biolab/miniconda3/envs/nkx3_py/bin/python -I /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_tables.py
"""
import json, os, sys, platform
import pandas as pd

OUT = "/mnt/ssd1/Projects/NKX3_wd/results/revision_2026/tcga"
TAB = os.path.join(OUT, "tables")


def t(name):
    return pd.read_csv(os.path.join(TAB, name), sep="\t")


SHEETS = [
    ("TCGA_claimed_tumour", "tcga_claimed_genes_vs_NKX3-1_primary_tumours.tsv",
     "41 genes reported in the manuscript (activated/repressed in Abstract, Results, Fig 1; Fig 3D/3E modules) vs NKX3-1 in "
     "TCGA-PRAD primary tumours. rho_* = Spearman (raw) or partial Spearman adjusted for AR score (AR), proliferation score "
     "(PROL), both (AR_PROL) or both + stromal score (AR_PROL_STROMA). pctl_* = percentile among the NKX3-1 correlations of "
     "all expressed genes (median TPM >= 1). padj_claimed_* = BH within the 41 genes; padj_genomewide_* = BH over expressed "
     "+ claimed genes. agree_* = sign matches the reported direction. *_legacyRSEM = same on legacy RSEM data."),
    ("TCGA_claimed_normal", "tcga_claimed_genes_vs_NKX3-1_solid_normals.tsv", "As above, 52 solid tissue normal samples."),
    ("TCGA_sign_tests", "tcga_direction_sign_tests.tsv",
     "Direction agreement with the reported sign. p_binom_vs_0.5 and p_binom_vs_base_rate (genome-wide fraction of genes "
     "correlating in that direction) are one-sided binomial tests; p_perm_* compare the signed mean rho and the agreement "
     "count with 10,000 random gene sets matched on expression decile (expressed genes only)."),
    ("TCGA_null_summary", "tcga_genomewide_null_summary.tsv",
     "Distribution of NKX3-1 correlations over all expressed genes (the null against which genes are ranked)."),
    ("TCGA_module_coherence", "tcga_module_coherence.tsv",
     "Mean pairwise Spearman rho among module genes vs 10,000 random equal-size gene sets matched on expression decile. "
     "scope 'expressed' = genes with median TPM >= 1; 'all' = all detected module genes (mean TPM >= 0.1). posctrl_* = "
     "positive-control sets (AR-score genes, proliferation genes, stromal genes)."),
    ("TCGA_module_pairs", "tcga_module_pairwise_spearman_primary_tumours.tsv", "Pairwise Spearman rho among expressed module genes."),
    ("TCGA_programme_scores", "tcga_hallmark_scores_vs_NKX3-1.tsv",
     "Single-sample programme scores (mean z of expressed Hallmark genes, NKX3-1 excluded) vs NKX3-1; 95% CI Bonett-Wright."),
    ("TCGA_programme_competitive", "tcga_hallmark_scores_vs_NKX3-1_competitive_null.tsv",
     "Programme-score correlation vs scores of 10,000 random gene sets of equal size matched on expression decile."),
    ("TCGA_programme_covariates", "tcga_hallmark_scores_vs_covariates_primary_tumours.tsv",
     "Correlation of each programme score with the AR, proliferation and stromal scores (tumours)."),
    ("TCGA_RSEM_programmes", "tcga_legacyRSEM_hallmark_scores_vs_NKX3-1.tsv",
     "Programme-score correlations recomputed on the legacy RSEM (upper-quartile normalised) data."),
    ("TCGA_covariate_genes", "tcga_covariate_genes_vs_NKX3-1_primary_tumours.tsv", "Genes used for the covariate scores."),
    ("Hallmark_fgsea_TCGA_DU145", "tcga_vs_DU145_hallmark_fgsea_long.tsv",
     "Preranked fgsea (50 Hallmarks). TCGA rankings = NKX3-1 correlations; DU145 rankings = Wald z of the corrected "
     "(featureCounts) DESeq2 model, raw and with the expression-dependent trend removed (robust loess on log10 baseMean)."),
    ("DU145_trend", "DU145_wald_z_by_expression_decile.tsv",
     "Median DU145 Wald z by baseMean decile (corrected counts): low-expressed genes trend up, high-expressed genes down, "
     "in every contrast, including CM vs control."),
    ("RWPE1_claimed", "rwpe1_claimed_genes_limma.tsv",
     "Same 41 genes in public RWPE1 data: limma-trend on DESeq VST values, NKX3.1 WT minus vector (GSE81439, 6 vs 6; "
     "GSE115335, paraquat-treated, 6 vs 6). passes_filter = mean VST >= 4; pctl_t = percentile among tested genes."),
    ("RWPE1_sign_tests", "rwpe1_direction_sign_tests.tsv", "Direction agreement in RWPE1 (one-sided binomial tests)."),
    ("RWPE1_camera", "rwpe1_camera_claimed_sets_and_hallmarks.tsv", "limma camera tests for the claimed gene sets and five Hallmarks."),
    ("RWPE1_fgsea", "rwpe1_hallmark_fgsea.tsv", "Preranked fgsea on limma t, raw and expression-detrended."),
    ("RWPE1_trend", "rwpe1_limma_t_by_expression_decile.tsv", "Median limma t by expression decile in RWPE1."),
    ("RWPE1_concordance", "rwpe1_global_concordance.tsv",
     "Genome-wide Spearman correlation of gene statistics: RWPE1 NKX3.1 response vs DU145 NKX3.1 response and vs TCGA NKX3-1 correlations."),
    ("GEO_curated", "geo_nkx3-1_perturbation_datasets_curated.tsv", "Manually curated NKX3.1 perturbation datasets from the GEO search."),
    ("GEO_all_hits", "geo_nkx3-1_search_all_hits.tsv", "All GEO series returned by the E-utilities queries (see geo_esearch_esummary_raw.json)."),
]


def main():
    readme = pd.DataFrame([(s, f, d) for s, f, d in SHEETS], columns=["sheet", "source_file", "description"])
    head = pd.DataFrame({"sheet": ["(all)"], "source_file": ["-"], "description": [
        "Independent-dataset check for Reviewer 1, major point 5 (TCGA-PRAD) and public NKX3.1 perturbation data (RWPE1). "
        "Sign convention: treatment minus reference (NKX3.1 minus vector control); TCGA = correlation with NKX3-1 across tumours. "
        "Correlation across tumours does not establish regulation; TCGA-PRAD tumours are AR-positive with endogenous NKX3-1, "
        "unlike AR-negative DU145 cells with ectopic NKX3.1."]})
    readme = pd.concat([head, readme], ignore_index=True)
    xlsx = os.path.join(TAB, "Supplementary_Table_TCGA_RWPE1_NKX3-1_check.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        readme.to_excel(w, sheet_name="README", index=False)
        for sheet, f, _ in SHEETS:
            t(f).to_excel(w, sheet_name=sheet[:31], index=False)
    # human-readable provenance
    P = json.load(open(os.path.join(OUT, "tcga_provenance.json")))
    G = json.load(open(os.path.join(OUT, "geo_rwpe1_download_provenance.json")))
    R = json.load(open(os.path.join(OUT, "rwpe1_summary.json")))
    raw = json.load(open(os.path.join(OUT, "geo_esearch_esummary_raw.json")))
    L = ["PROVENANCE -- TCGA-PRAD / public NKX3.1 data check (Reviewer 1, major 5)", "=" * 72,
         f"Analysis run (UTC): {P['run_utc']}", f"Primary source: {P['primary_source']}", ""]
    for k, s in P["sources"].items():
        L += [f"[{k}] {s['url']}", f"    local: {s['local_file']}", f"    bytes: {s['bytes']}  md5: {s['md5']}",
              f"    downloaded (UTC): {s['accessed_utc']}   server Last-Modified: {s['server_last_modified']}"]
    L += ["", f"Values: {P['unit']}; column sums of TPM verified: {P.get('verification_tpm_colsum_range')}",
          f"Genes in file: {P['n_genes_in_file']} (PAR_Y rows removed: {P['n_par_y_rows_removed']})",
          f"Samples in file: {P['n_samples_in_file']}; sample-type codes: {P['sample_type_codes_in_file']}",
          f"Primary tumours used: {P['n_primary_tumour_used']} (duplicate aliquots removed: {', '.join(P['primary_duplicates_removed'])})",
          f"Solid tissue normals used: {P['n_solid_normal_used']} (duplicates removed: {', '.join(P['normal_duplicates_removed']) or 'none'})",
          f"Sample selection: {P['sample_selection']}", f"Expression filter: {P['expression_filter']}",
          f"Expressed genes: tumours {P['n_expressed_tumour']}, normals {P['n_expressed_normal']}",
          f"Gene matching: {P['gene_matching']}",
          f"MSigDB via msigdbr {P['msigdbr_version']} (db {P['msigdb_db_version']})",
          f"Covariates: {json.dumps(P['covariates'])}", f"Seed {P['seed']}; permutations {P['permutations']}",
          f"Legacy RSEM sensitivity: {json.dumps(P.get('legacy_RSEM_sensitivity'))}",
          f"Context: {json.dumps(P['context'])}", "",
          f"GEO search (NCBI E-utilities, db=gds) accessed {raw['accessed_utc']}:"]
    for q, v in raw["queries"].items():
        L.append(f"    {q}: {v['count']} hits  term = {v['term']}")
    L += ["", "GEO supplementary matrices analysed:"]
    for acc, s in G.items():
        L += [f"[{acc}] {s['url']}", f"    bytes: {s['bytes']}  md5: {s['md5']}  downloaded (UTC): {s['accessed_utc']}"]
    L += [f"RWPE1 analysis: {json.dumps(R)}", "",
          f"Python {platform.python_version()}, pandas {pd.__version__}", "Scripts:"] + [f"    {s}" for s in P["scripts"]] + \
         ["    /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_tables.py",
          "    /mnt/ssd1/Projects/NKX3_wd/scripts/revision/15_tcga_run_all.sh"]
    with open(os.path.join(OUT, "tcga_provenance.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote", xlsx)


if __name__ == "__main__":
    main()
