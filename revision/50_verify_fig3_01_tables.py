#!/usr/bin/env python
"""
50_verify_fig3_01_tables.py -- task verify_round2/fig3 (adversarial check of fig3_alternative)

Independent re-derivation (Python; no code shared with 43_fig3_alternative.R) of
  1. the rows of the ORIGINAL Fig 3 heatmaps as printed in the v15 PDFs (pdftotext -bbox), compared
     with fig3alt_original_fig3_rows_corrected_status.tsv (row order) and with the corrected DESeq2
     TSVs (tested / padj / direction counts);
  2. the candidate leading-edge heatmap rows: leading edge (core_enrichment) of the corrected
     Hallmark GSEA tables, ranked by |Wald statistic| from core/deseq2_fc_<contrast>.tsv, top-k per
     set (k raised until the union reaches 30), and the log2FC (MLE) / shrunken log2FC / padj that the
     figure annotates; NES printed in the panel headers;
  3. replicate-consistency numbers of the displayed genes from core/vst_fc.tsv.

Inputs are read-only. Outputs: results/revision_2026/verify_round2/fig3/v01_*.tsv
Run: /home/biolab/miniconda3/envs/nkx3_py/bin/python -I /mnt/ssd1/Projects/NKX3_wd/scripts/revision/50_verify_fig3_01_tables.py
"""
import math
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

PROJ = "/mnt/ssd1/Projects/NKX3_wd"
REV = f"{PROJ}/results/revision_2026"
F3A = f"{REV}/fig3_alternative"
OUT = f"{REV}/verify_round2/fig3"
WORK = ("/mnt/ssd1/Projects/NKX3_wd/work/verify_round2_fig3")
V15 = f"{PROJ}/results/rna_downstreams_out_v15"
os.makedirs(OUT, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def read_tsv(p, **kw):
    return pd.read_csv(p, sep="\t", **kw)


# ----------------------------------------------------------------------------------------------
# 1. original Fig 3 rows: PDF labels vs fig3alt row-status table; corrected status recomputed
# ----------------------------------------------------------------------------------------------
def pdf_row_labels(pdf):
    """Row labels (top to bottom) of a pheatmap PDF: words in the most common x column."""
    with tempfile.NamedTemporaryFile(suffix=".html", dir=WORK, delete=False) as tmp:
        tmpf = tmp.name
    subprocess.run(["pdftotext", "-bbox", pdf, tmpf], check=True)
    rows = []
    pat = re.compile(r'<word xMin="([0-9.]+)" yMin="([0-9.]+)" xMax="([0-9.]+)" yMax="([0-9.]+)">(.*?)</word>')
    with open(tmpf, encoding="utf-8") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                rows.append((float(m.group(1)), float(m.group(2)), m.group(5)))
    os.unlink(tmpf)
    xs = pd.Series([round(r[0], 2) for r in rows])
    col = xs.value_counts().idxmax()
    sel = [(y, w) for (x, y, w) in rows if round(x, 2) == col and w != "condition"]
    sel.sort()
    lab = [w.replace("−", "-") for (_, w) in sel]
    return lab


ORIG = {
    "A": dict(pdf="cmp_HM-vector-C_vs_HM-vector-CM_heatmap_sig_genes.pdf", contrast="CM_in_vector", sign=-1),
    "B": dict(pdf="cmp_HM-NKX3.1-C_vs_HM-vector-C_heatmap_sig_genes.pdf", contrast="NKX_in_ctrl", sign=1),
    "C": dict(pdf="cmp_HM-NKX3.1-CM_vs_HM-vector-CM_heatmap_sig_genes.pdf", contrast="NKX_in_CM", sign=1),
}
rs = read_tsv(f"{F3A}/fig3alt_original_fig3_rows_corrected_status.tsv")
orig_rows = []
for k, o in ORIG.items():
    lab = pdf_row_labels(f"{V15}/{o['pdf']}")
    sub = rs[rs.panel == k].sort_values("heatmap_row")
    # compare the printed labels with the table symbols as they are (both carry the same make.unique suffixes;
    # stripping ".1" would wrongly alter clone names such as AC011043.1)
    tab_syms = sub.symbol.tolist()
    same_len = len(lab) == len(tab_syms)
    exact = same_len and all(a == b for a, b in zip(lab, tab_syms))
    n_mismatch = sum(a != b for a, b in zip(lab, tab_syms)) + abs(len(lab) - len(tab_syms))
    # corrected status from the DESeq2 TSV (independent of fig3alt)
    d = read_tsv(f"{REV}/core/deseq2_fc_{o['contrast']}.tsv").set_index("gene_id")
    ids = sub.gene_id.tolist()
    tested = [g in d.index for g in ids]
    lfc = np.array([d.loc[g, "log2FC_MLE"] if g in d.index else np.nan for g in ids])
    padj = np.array([d.loc[g, "padj"] if g in d.index else np.nan for g in ids])
    orig_diff = sub.orig_welch_diff_vst_treat_minus_ref.to_numpy()
    same = np.sign(lfc) == np.sign(orig_diff)
    tested_arr = np.array(tested)
    sig_syms = [s for s, p in zip(sub.symbol, padj) if p == p and p < 0.05]
    # rows whose direction call in the fig3alt table differs from the unrounded TSV value
    tab_same = sub.same_direction.to_numpy()
    flips = [f"{s} (table LFC {t}, TSV LFC {v:.2e})" for s, t, v, a, b, tt in
             zip(sub.symbol, sub.corrected_log2FC_MLE, lfc, same, tab_same, tested) if tt and str(b) in ("True", "False")
             and bool(a) != (str(b) == "True")]
    # NA padj (independent filtering / outlier) among tested rows
    n_padj_na = int(np.sum(tested_arr & np.isnan(padj)))
    orig_rows.append(dict(
        panel=k, pdf=o["pdf"], pdf_rows=len(lab), table_rows=len(tab_syms), order_identical=exact,
        n_label_mismatches=n_mismatch, rows_not_in_deseq2_table=int((~tested_arr).sum()),
        rows_tested_padj_NA=n_padj_na,
        rows_padj_lt_0_05=len(sig_syms), which=",".join(sig_syms),
        tested_same_direction=int(np.sum(same & tested_arr)),
        tested_opposite_direction=int(np.sum(~same & tested_arr)),
        median_abs_orig_diff=round(float(np.median(np.abs(orig_diff))), 3),
        median_abs_corrected_lfc=round(float(np.nanmedian(np.abs(lfc))), 3),
        max_unique_reads_le5=int((sub.unique_reads_max_per_library <= 5).sum()),
        boxed_rows=int(sub.in_boxed_module.sum()),
        boxed_rows_not_tested=int((sub.in_boxed_module.to_numpy() & ~tested_arr).sum()),
        direction_calls_differing_from_fig3alt_table="; ".join(flips),
    ))
orig_rows = pd.DataFrame(orig_rows)
orig_rows.to_csv(f"{OUT}/v01_original_fig3_rows_recheck.tsv", sep="\t", index=False)
log("== original Fig 3 rows: PDF vs fig3alt table, corrected status recomputed from deseq2_fc TSVs")
log(orig_rows.to_string(index=False))

# ----------------------------------------------------------------------------------------------
# 2. candidate leading-edge heatmaps
# ----------------------------------------------------------------------------------------------
PANELS = {
    "A": dict(contrast="CM_in_vector", sets=["HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
              ref="vector_ctrl", trt="vector_CM"),
    "B": dict(contrast="NKX_in_ctrl", sets=["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"],
              ref="vector_ctrl", trt="NKX3.1_ctrl"),
    "C": dict(contrast="NKX_in_CM", sets=["HALLMARK_APICAL_JUNCTION"], ref="vector_CM", trt="NKX3.1_CM"),
}
# what the PNG shows (read off fig3alt_candidate_leading_edge_heatmaps_6samples.png by eye, top to bottom)
PNG_ROWS = {
    "A": "CDC20 ERH APEX1 KPNA2 HDGF CCT5 PTGES3 CCT3 EIF2S1 RAN GSPT1 ILF2 CCT2 CUL1 SNRPG "
         "ATP5MG ATP5MC3 NDUFB4 SUCLA2 NDUFA8 HCCS NDUFB6 MRPS15 AIFM1 ATP6V1G1 ATP1B1 SLC25A20 SLC25A4 COX7C ATP5PF",
    "B": "SLBP DCK CDKN2A RPA1 TRIP13 CDCA8 DCLRE1B TUBG1 HMGB2 NCAPD2 NUP205 MELK SSRP1 "
         "CTCF STMN1 PLK1 CDC25B GSPT1 "
         "CUL3 HOXC10 RASAL2 FOXN3 KIF5B TOP1 CBX1 HIF1A TRAIP SQLE MAP3K20 RAD54L NEK2",
    "C": "YWHAH MVD ALOX15B MYL9 ACTB JUP TUBG1 SHROOM2 CNN2 ACTG1 STX4 PFN1 ACTN4 GNAI2 MSN GRB7 VAV2 EXOC4 "
         "CRAT ACTC1 ACTN1 IRS1 ICAM5 BAIAP2 CTNNA1 ITGB4 TAOK2 RHOF GTF2F1 MYH9",
}
PNG_NES = {"HALLMARK_MYC_TARGETS_V1": -2.53, "HALLMARK_OXIDATIVE_PHOSPHORYLATION": -2.16,
           "HALLMARK_E2F_TARGETS": -1.78, "HALLMARK_G2M_CHECKPOINT": -1.82, "HALLMARK_APICAL_JUNCTION": -2.38}
N_ROWS = 30

le_tab = read_tsv(f"{F3A}/fig3alt_leading_edge_genes.tsv")
vst = read_tsv(f"{REV}/core/vst_fc.tsv")
vst = vst.set_index(vst.columns[0])
cd_groups = {  # sample -> group (core cd order)
    "24-004-1_S56": "vector_ctrl", "24-004-2_S52": "vector_ctrl", "24-004-3_S41": "vector_ctrl",
    "24-004-4_S27": "vector_CM", "24-004-5_S19": "vector_CM", "24-004-6_S17": "vector_CM",
    "24-004-7_S14": "NKX3.1_ctrl", "24-004-8_S53": "NKX3.1_ctrl", "24-004-9_S42": "NKX3.1_ctrl",
    "24-004-10_S20": "NKX3.1_CM", "24-004-11_S18": "NKX3.1_CM", "24004-12_S15": "NKX3.1_CM",
}
le_rows, panel_rows, gene_rows = [], [], []
for pn, p in PANELS.items():
    g = read_tsv(f"{REV}/gsea/tables/gsea_fc_hallmark_{p['contrast']}.tsv", dtype={"core_enrichment": str})
    d = read_tsv(f"{REV}/core/deseq2_fc_{p['contrast']}.tsv", dtype={"entrez": str})
    r = d[d.entrez.notna() & np.isfinite(d.stat)].copy()
    r["abs_stat"] = r.stat.abs()
    r = r.sort_values("abs_stat", ascending=False, kind="mergesort").drop_duplicates("entrez", keep="first")
    by_entrez = r.set_index("entrez")
    le_by_set = {}
    for s in p["sets"]:
        row = g[g.ID == s].iloc[0]
        le = row.core_enrichment.split("/")
        missing = [e for e in le if e not in by_entrez.index]
        m = by_entrez.loc[[e for e in le if e in by_entrez.index]].sort_values("abs_stat", ascending=False, kind="mergesort")
        le_by_set[s] = m.gene_id.tolist()
        le_rows.append(dict(panel=pn, contrast=p["contrast"], set=s, NES_table=round(row.NES, 4),
                            NES_png=PNG_NES[s], NES_matches=round(row.NES, 2) == PNG_NES[s],
                            FDR_table=row["p.adjust"], le_size=len(le), le_missing_from_ranking=len(missing),
                            le_genes_with_positive_stat=int((m.stat > 0).sum())))
    # top-k per set, k raised until the union reaches N_ROWS
    if len(p["sets"]) == 1:
        sel = le_by_set[p["sets"][0]][:N_ROWS]; k = N_ROWS
    else:
        k = math.ceil(N_ROWS / len(p["sets"]))
        while True:
            sel = list(dict.fromkeys(sum([v[:k] for v in le_by_set.values()], [])))
            if len(sel) >= N_ROWS or k >= max(len(v) for v in le_by_set.values()):
                break
            k += 1
    sym = d.set_index("gene_id").symbol
    dd = d.set_index("gene_id").loc[sel]
    # display order: set1-only, both, set2-only; within block by Wald stat ascending
    if len(p["sets"]) == 2:
        s1, s2 = p["sets"]
        def mem(gid):
            a, b = gid in le_by_set[s1], gid in le_by_set[s2]
            return 0 if (a and not b) else (1 if (a and b) else 2)
        order = sorted(sel, key=lambda gid: (mem(gid), dd.loc[gid, "stat"]))
    else:
        order = sorted(sel, key=lambda gid: dd.loc[gid, "stat"])
    my_syms = [sym[gid] for gid in order]
    png_syms = PNG_ROWS[pn].split()
    f3 = le_tab[(le_tab.panel == pn) & (le_tab.displayed)].sort_values("display_row")
    f3_syms = f3.symbol.tolist()
    # LFC annotation: fig3alt table vs DESeq2 TSV
    f3l = f3.set_index("gene_id").log2FC_MLE
    lfc_diff = float(np.max(np.abs(f3l.to_numpy() - d.set_index("gene_id").loc[f3l.index, "log2FC_MLE"].to_numpy())))
    # replicate consistency on vst_fc.tsv
    ref_s = [s for s, gg in cd_groups.items() if gg == p["ref"]]
    trt_s = [s for s, gg in cd_groups.items() if gg == p["trt"]]
    X = vst.loc[sel]
    mr, mt = X[ref_s].mean(axis=1), X[trt_s].mean(axis=1)
    sgn = np.sign(mt - mr)
    trt_ok = (np.sign(X[trt_s].sub(mr, axis=0)).eq(sgn, axis=0)).all(axis=1)
    ref_ok = (np.sign(X[ref_s].sub(mt, axis=0)).eq(-sgn, axis=0)).all(axis=1)
    sep = np.where(sgn < 0, X[trt_s].max(axis=1) < X[ref_s].min(axis=1), X[trt_s].min(axis=1) > X[ref_s].max(axis=1))
    lfc = dd.log2FC_MLE
    panel_rows.append(dict(
        panel=pn, contrast=p["contrast"], k_per_set=k, n_displayed=len(sel),
        my_rows_equal_fig3alt_table=my_syms == f3_syms, my_rows_equal_png=my_syms == png_syms,
        max_abs_lfc_diff_table_vs_tsv=lfc_diff,
        n_positive_lfc_displayed=int((lfc > 0).sum()),
        n_lfc_below_minus1_truncated=int((lfc < -1).sum()),
        median_log2FC_MLE=round(float(lfc.median()), 3),
        median_fold_change_pct=round(100 * (2 ** float(lfc.median()) - 1), 1),
        min_log2FC_MLE=round(float(lfc.min()), 3), max_log2FC_MLE=round(float(lfc.max()), 3),
        median_log2FC_shrunk=float(f"{dd.log2FC_shrunk.median():.2g}"),
        max_abs_log2FC_shrunk=float(f"{dd.log2FC_shrunk.abs().max():.2g}"),
        min_padj=float(f"{dd.padj.min():.2g}"), n_padj_lt_0_05=int((dd.padj < 0.05).sum()),
        n_padj_lt_0_1=int((dd.padj < 0.1).sum()),
        median_baseMean=round(float(dd.baseMean.median()), 0),
        pct_both_groups_all_replicates_agree=round(100 * float((trt_ok & ref_ok).mean()), 1),
        pct_complete_separation=round(100 * float(np.mean(sep)), 1),
    ))
    for gid in order:
        gene_rows.append(dict(panel=pn, symbol=sym[gid], gene_id=gid, stat=round(dd.loc[gid, "stat"], 3),
                              log2FC_MLE=round(dd.loc[gid, "log2FC_MLE"], 4),
                              log2FC_shrunk=float(f"{dd.loc[gid, 'log2FC_shrunk']:.3g}"),
                              padj=float(f"{dd.loc[gid, 'padj']:.3g}"), baseMean=round(dd.loc[gid, "baseMean"], 1)))
    if my_syms != png_syms:
        log(f"  panel {pn}: my order vs PNG differ: {[(a, b) for a, b in zip(my_syms, png_syms) if a != b]}")
le_rows = pd.DataFrame(le_rows)
panel_rows = pd.DataFrame(panel_rows)
gene_rows = pd.DataFrame(gene_rows)
le_rows.to_csv(f"{OUT}/v01_leading_edge_sets_recheck.tsv", sep="\t", index=False)
panel_rows.to_csv(f"{OUT}/v01_leading_edge_panels_recheck.tsv", sep="\t", index=False)
gene_rows.to_csv(f"{OUT}/v01_leading_edge_displayed_genes_recheck.tsv", sep="\t", index=False)
log("== leading-edge sets (NES in table vs PNG header)")
log(le_rows.to_string(index=False))
log("== panels (independent re-derivation vs fig3alt table and PNG)")
log(panel_rows.T.to_string())

# compare with fig3alt_replicate_consistency.tsv
rc = read_tsv(f"{F3A}/fig3alt_replicate_consistency.tsv")
log("== fig3alt_replicate_consistency.tsv (for comparison)")
log(rc[["panel", "median_log2FC_MLE", "median_fold_change_pct", "median_log2FC_shrunk_apeglm", "min_gene_padj",
        "pct_both_groups_all_replicates_agree", "pct_complete_separation"]].to_string(index=False))

with open(f"{OUT}/v01_log.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
    fh.write(f"\npython {sys.version}\npandas {pd.__version__} numpy {np.__version__}\n")
