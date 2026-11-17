"""CDK9 L156F RNA-seq volcanoes (supplementary 5e).

From `05_RNA_seq/RNAseq_analysis.ipynb` cells 2-4. Five
DESeq2 contrasts side by side: the L156F and NT lines at 24 and 48 hours of
KI-CDK9d-32, plus the untreated L156F-vs-NT baseline.

Genes are filtered to `baseMean >= 20` before anything is called, so a gene too
lowly expressed to have a trustworthy fold change cannot become a point. A gene
is differentially expressed at |LFC| >= 1 and padj <= 0.01, and the counts in
each panel's legend are what the figure reports.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import io, paths

RNASEQ = paths.PROCESSED / "rnaseq"
BASEMEAN_CUTOFF = 20
LFC_CUTOFF, PADJ_CUTOFF = 1.0, 0.01

#: contrast -> the title its panel carries
CONTRASTS = {
    "L156F_24hours_vs_L156F_no_drug": "L156F KI-CDK9d-32\n24 hrs vs. No Drug",
    "L156F_48hours_vs_L156F_no_drug": "L156F KI-CDK9d-32\n48 hrs vs. No Drug",
    "NT_24hours_vs_NT_no_drug": "NT gRNA KI-CDK9d-32\n24 hrs vs. No Drug",
    "NT_48hours_vs_NT_no_drug": "NT gRNA KI-CDK9d-32\n48 hrs vs. No Drug",
    "L156F_no_drug_vs_NT_no_drug": "L156F vs. NT gRNA\n(No Drug)",
}


def gene_names() -> dict:
    tx = io.read(RNASEQ / "tx2gene.tsv", sep="\t")
    return dict(zip(tx["gene_id"], tx["gene_name"]))


def contrast(name: str, names: dict | None = None) -> pd.DataFrame:
    """One DESeq2 contrast, expression-filtered and labelled with gene names."""
    names = names if names is not None else gene_names()
    d = io.read(RNASEQ / f"{name}.csv").rename(columns={"Unnamed: 0": "gene_id"})
    d["gene_name"] = d["gene_id"].map(names)
    return d[d["baseMean"] >= BASEMEAN_CUTOFF].reset_index(drop=True)


def source() -> pd.DataFrame:
    """Every contrast's filtered genes with its DE call -- the panel's data."""
    names = gene_names()
    out = []
    for key in CONTRASTS:
        d = contrast(key, names).copy()
        d["contrast"] = key
        d["direction"] = np.where(
            (d["log2FoldChange"] >= LFC_CUTOFF) & (d["padj"] <= PADJ_CUTOFF), "up",
            np.where((d["log2FoldChange"] <= -LFC_CUTOFF) & (d["padj"] <= PADJ_CUTOFF),
                     "down", "ns"))
        out.append(d[["contrast", "gene_id", "gene_name", "baseMean",
                      "log2FoldChange", "padj", "direction"]])
    return pd.concat(out, ignore_index=True)


def counts() -> pd.DataFrame:
    """DE up/down per contrast -- the numbers each panel's legend prints."""
    d = source()
    g = (d[d.direction != "ns"].groupby(["contrast", "direction"]).size()
         .unstack(fill_value=0).reindex(list(CONTRASTS)).fillna(0).astype(int))
    return g.reset_index()


def supp_5e(outdir=None, **kw):
    """-> (figure, the per-contrast DE counts)."""
    names = gene_names()
    fig, ax = plt.subplots(1, 5, figsize=(16, 5), sharey=True, sharex=True)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(axis="x", which="major", labelsize=15)
        a.tick_params(axis="y", which="major", labelsize=15)
        a.set_xlabel("LFC", fontsize=15)
    ax[0].set_ylabel("$-log_{10}$(FDR)", fontsize=15)

    for i, (key, title) in enumerate(CONTRASTS.items()):
        d = contrast(key, names)
        up = d[(d["log2FoldChange"] >= LFC_CUTOFF) & (d["padj"] <= PADJ_CUTOFF)]
        down = d[(d["log2FoldChange"] <= -LFC_CUTOFF) & (d["padj"] <= PADJ_CUTOFF)]
        rest = d.drop(list(up.index) + list(down.index))
        # the data marks are rasterized; the axes and text stay vector
        ax[i].scatter(rest["log2FoldChange"], -np.log10(rest["padj"]),
                      color="tab:grey", alpha=.2, rasterized=True)
        ax[i].scatter(up["log2FoldChange"], -np.log10(up["padj"]), color="tab:red",
                      label=f"{len(up)} DE UP", alpha=.5, rasterized=True)
        ax[i].scatter(down["log2FoldChange"], -np.log10(down["padj"]), color="tab:blue",
                      label=f"{len(down)} DE DOWN", alpha=.5, rasterized=True)
        ax[i].legend(fontsize=12)
        ax[i].set_title(title, fontsize=14)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "Supp_5e_RNA_seq_comp.pdf"), dpi=300, transparent=True)
    return fig, counts()


FIGURES = {"S5e": supp_5e}
PANELS = {"S5e": counts}
