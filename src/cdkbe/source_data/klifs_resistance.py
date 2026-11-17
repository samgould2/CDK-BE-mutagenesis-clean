"""Clinically observed resistance mutations against the KLIFS pocket (supp 15b).

From `07_KLIFS/KLIFS_analysis.ipynb` cells 17-21.

Every OncoKB resistance mutation in the ten comparison kinases is asked one
question: does its codon fall on one of the 85 numbered KLIFS pocket positions
for that kinase? The answer is a stacked bar per gene -- inside the pocket on
the bottom, outside on top.

The KLIFS index comes from each kinase's own `KLIFS2UniProtIdx` map, so a
position is matched by residue number in that protein rather than by alignment
column; the `pocket_seq` check asserts the residue the mutation reports is the
one KLIFS has at that index, which is what makes the match trustworthy.

Note this uses `oncokb_klifs/`, the OncoKB comparison kinases, not the nine
CDK json files beside it.
"""
from __future__ import annotations

import json
from functools import reduce

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .. import io, paths

KLIFS = paths.REFERENCE / "KLIFS"
ONCOKB = KLIFS / "OncoKB_resistance.csv"


def klifs_index(genes, subdir="oncokb_klifs") -> pd.DataFrame:
    """KLIFS index -> residue number and pocket residue, one column pair per gene."""
    frames = []
    for gene in genes:
        with open(KLIFS / subdir / f"{gene}.json") as fh:
            data = json.load(fh)
        k = pd.DataFrame(list(data["KLIFS2UniProtIdx"].items()),
                         columns=["KLIFS_idx", gene])
        k[f"{gene}_seq"] = list(data["klifs"]["pocket_seq"])
        frames.append(k)
    return reduce(lambda a, b: pd.merge(a, b, on="KLIFS_idx"), frames)


def annotated() -> pd.DataFrame:
    """OncoKB resistance mutations, each marked inside or outside the pocket."""
    oncokb = io.read(ONCOKB)
    genes = sorted(oncokb["Gene"].unique())
    comb = klifs_index(genes)
    # declare the columns up front: assigning into a fresh all-NaN float column
    # leaves `KLIFS Catalytic Site` as object dtype, where `~` is bitwise
    # negation on Python bools rather than a logical not
    oncokb["KLIFS_Index"] = pd.Series([pd.NA] * len(oncokb), dtype="object")
    oncokb["KLIFS Catalytic Site"] = False
    for i, val in oncokb.iterrows():
        g, cod, mutation = val["Gene"], val["Codon"], val["Mutation"]
        hit = comb[comb[g] == cod]
        if len(hit):
            oncokb.loc[i, "KLIFS_Index"] = hit["KLIFS_idx"].values[0]
            oncokb.loc[i, "KLIFS Catalytic Site"] = True
            # the residue KLIFS has at that index must be the one the mutation
            # reports; otherwise the numbering does not line up
            assert hit[f"{g}_seq"].values[0] == mutation[0], f"{g} {mutation}"
        else:
            oncokb.loc[i, "KLIFS_Index"] = None
    return oncokb.astype({"KLIFS Catalytic Site": bool})


def counts() -> pd.DataFrame:
    """Per gene, how many unique resistance mutations fall in or out of the pocket."""
    uniq = annotated().drop_duplicates(subset=["Gene", "Mutation"])
    cat = (pd.DataFrame(uniq[["Gene", "KLIFS Catalytic Site"]]
                        .groupby("Gene").value_counts()).reset_index())
    return cat.pivot_table(index="Gene", columns="KLIFS Catalytic Site",
                           values="count").fillna(0)


def supp_15b(outdir=None, **kw):
    """-> (figure, the per-gene in/out counts)."""
    pivot = counts()
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.barplot(x=pivot.index, y=pivot[True], color="skyblue", ax=ax,
                edgecolor="black", label="Within KLIFS Catalytic Site")
    sns.barplot(x=pivot.index, y=pivot[False], bottom=pivot[True], color="lightcoral",
                ax=ax, edgecolor="black", label="Outside KLIFS Catalytic Site")
    ax.set_title("Presence of Resistance Mutations in KLIFS Catalytic Site", fontsize=15)
    ax.legend(fontsize=15, loc="upper right")
    ax.set_ylabel("No. of Resistance Mutations", fontsize=15)
    ax.set_xlabel("")
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=14)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "Supp_15b_oncokb_KLIFS_inclusion.pdf"), transparent=True)
    return fig, pivot.reset_index()


FIGURES = {"S15b": supp_15b}
PANELS = {"S15b": lambda: counts().reset_index()}
