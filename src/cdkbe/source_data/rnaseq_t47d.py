"""Figure 5f-5i: the T47D CDK4 RNA-seq panels.

From `cdk4_program_analysis/pub_fig_square.py`, which draws all four as one
2x2 figure (`pub_fig_square_4panel_mid.pdf`):

    5f  CDK4 transcript level per genotype
    5g  Hallmark E2F target score per sample
    5h  volcano, CDK4 W198* vs CBE NT
    5i  volcano, CDK4 W198* vs D196N
"""
from __future__ import annotations
import json

import numpy as np
import pandas as pd
import scipy.stats as sps

from .. import io, paths

BASE = paths.PROCESSED / "rnaseq" / "T47D_RNAseq"
ANALYSIS = BASE / "cdk4_program_analysis"
RESULTS = ANALYSIS / "results"

CTRL = "CBE_NT"
VARIANTS = ["CDK4_E56K", "CDK4_D196N", "CDK4_W198"]
GENOTYPES = [CTRL] + VARIANTS

BASEMEAN_MIN = 50          # a gene below this has too few reads to call
LFC_T, PADJ_T = 1.0, 0.05


def _expression():
    """-> (tpm, log2 tpm, sample->genotype, symbol->gene_id)."""
    tpm = io.read(BASE / "salmon.merged.gene_tpm.tsv", sep="\t", index_col=0)
    names = tpm.pop("gene_name")
    meta = pd.Series(tpm.columns.str.extract(r"^T47D_(.+)_veh_\d+$")[0].values,
                     index=tpm.columns)
    samples = [s for g in GENOTYPES for s in meta.index[meta == g]]
    tpm, meta = tpm[samples], meta[samples]
    log_tpm = np.log2(tpm + 1)
    order = tpm.mean(axis=1).sort_values(ascending=False).index
    s2i = pd.Series(order, index=names.reindex(order).values)
    return tpm, log_tpm, meta, s2i[~s2i.index.duplicated()]


def _holm(p):
    p = np.asarray(p, float)
    out = np.empty(len(p))
    running = 0.0
    for k, i in enumerate(np.argsort(p)):
        running = max(running, min(1.0, (len(p) - k) * p[i]))
        out[i] = running
    return out


def figure_5f() -> pd.DataFrame:
    """CDK4 transcript level (TPM) per genotype, with padj against the control."""
    tpm, _, meta, s2i = _expression()
    gene = s2i["CDK4"]
    res = {v: io.read(RESULTS / f"{v}_vs_{CTRL}_veh.csv").set_index("gene_id")
           for v in VARIANTS}
    rows = [dict(genotype=g, sample=s, CDK4_TPM=tpm.loc[gene, s],
                 padj_vs_CBE_NT=res[g].loc[gene, "padj"] if g in res else np.nan)
            for g in GENOTYPES for s in meta.index[meta == g]]
    return pd.DataFrame(rows)


def figure_5g() -> pd.DataFrame:
    """Hallmark E2F target score per sample, with Holm-adjusted p vs the control.

    The score is the mean z-score across the expressed E2F target genes; the
    test is one-sided (variants are expected to lower it) and Holm-corrected
    over the three variants.
    """
    _, log_tpm, meta, s2i = _expression()
    hall = json.loads((BASE / "T47D_CDK4_20260927_results" /
                       "hallmark_gene_sets_cache.json").read_text())
    ids = [s2i[s] for s in hall["E2F Targets"]
           if s in s2i.index and log_tpm.loc[s2i[s]].mean() > 1]
    m = log_tpm.loc[ids]
    score = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).mean()

    grp = lambda g: list(meta.index[meta == g])
    raw = []
    for v in VARIANTS:
        t, p2 = sps.ttest_ind(score[grp(v)], score[grp(CTRL)], equal_var=False)
        raw.append(p2 / 2 if t < 0 else 1 - p2 / 2)
    adj = dict(zip(VARIANTS, _holm(raw)))
    return pd.DataFrame([
        dict(genotype=g, sample=s, E2F_mean_zscore=score[s],
             holm_p_vs_CBE_NT=adj.get(g, np.nan), n_E2F_genes=len(ids))
        for g in GENOTYPES for s in grp(g)])


def _volcano(filename: str) -> pd.DataFrame:
    d = io.read(RESULTS / filename).set_index("gene_id")
    d = d.dropna(subset=["padj", "log2FoldChange"])
    d = d[d["baseMean"] >= BASEMEAN_MIN].copy()
    d["neg_log10_padj"] = -np.log10(d["padj"].clip(lower=1e-300))
    d["significant_up"] = (d.padj < PADJ_T) & (d.log2FoldChange > LFC_T)
    d["significant_down"] = (d.padj < PADJ_T) & (d.log2FoldChange < -LFC_T)
    cols = ["gene_name", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
            "neg_log10_padj", "significant_up", "significant_down"]
    return d.reset_index()[["gene_id"] + [c for c in cols if c in d.columns]]


def figure_5h() -> pd.DataFrame:
    """Volcano: CDK4 W198* vs CBE NT, vehicle."""
    return _volcano("CDK4_W198_vs_CBE_NT_veh_nofilt.csv")


def figure_5i() -> pd.DataFrame:
    """Volcano: CDK4 W198* vs D196N, vehicle."""
    return _volcano("CDK4_W198_vs_CDK4_D196N_veh_nofilt.csv")
