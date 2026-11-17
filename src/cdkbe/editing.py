"""Sensor-editing tables: per-read calls -> one variant call per gRNA.

This logic was duplicated in `lollipop_plots_canonical`, `fig1_analysis`,
`editing_analysis_v2` and both `logo_plots_*_canonical` notebooks. The copies
agreed, but only because nobody had edited one of them; it lives here now.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd

from . import io, paths

#: the edit each editor is allowed to make; reads carrying anything else are dropped
EDITS = {"ABE": ("A", "G"), "CBE": ("C", "T")}

_DNA_CHANGE = re.compile(r"\+(\d+)([ACGT])>([ACGT])")


def first_number(text: str) -> int:
    """Codon number out of an HGVSp string ('L156F' -> 156)."""
    return [int(i) for i in re.findall(r"\d+", text)][0]


def hgvsp_from_full_canonical(sample: str, editing_dir=None) -> pd.DataFrame:
    """Per-gRNA HGVSp table from the canonical full editing file.

    Keeps only reads whose every edit is the editor's own (an ABE read carrying
    a C>T came from something other than the ABE, so it is not evidence about
    this guide's edit).
    """
    editing_dir = editing_dir or paths.EDITING
    ref, alt = EDITS[sample[:3]]
    full = io.read(f"{editing_dir}/{sample}_full_canonical.csv")
    edits = full["DNA Change"].str.findall(_DNA_CHANGE)
    keep = edits.apply(lambda e: all((r == ref) and (a == alt) for _, r, a in e))
    hgvsp = (full[keep].groupby(["gRNA_id", "HGVSp"], as_index=False)["#Reads"].sum()
             .sort_values(["gRNA_id", "#Reads"], ascending=[True, False])
             .reset_index(drop=True))
    hgvsp["%Reads"] = hgvsp["#Reads"] / hgvsp.groupby("gRNA_id")["#Reads"].transform("sum") * 100
    return hgvsp


def savs_from_full_canonical(sample: str, editing_dir=None,
                             own_edits_only: bool = False) -> pd.DataFrame:
    """Single amino-acid variants, splitting compound outcomes into their parts.

    A guide whose sensor shows `R46W_K48E` is evidence for both variants, so
    the compound outcome's reads are counted toward each one. The denominator
    stays that guide's **total** reads, not the sum of the split rows, so a
    compound edit contributes its full share to every variant it contains and
    the percentages do not sum to 100 for such guides -- which is correct: they
    are per-variant prevalence, not a partition.

    `own_edits_only` controls whether reads carrying a substitution the editor
    cannot make are dropped first. **The two shipped conventions differ**, so
    this is not a detail to guess at:

    * the subpool `SAVs_canonical/*.csv` tables were built *without* the
      filter -- `own_edits_only=False` reproduces them exactly (16,026 rows for
      CBE_subpool1), and that is what the main-figure panels used;
    * `logo_plots_KB_compound_mut.ipynb` applies the filter, so supplementary
      9b needs `own_edits_only=True`.

    See docs/DISCREPANCIES.md item 14.
    """
    if own_edits_only:
        hgvsp = hgvsp_from_full_canonical(sample, editing_dir)
    else:
        editing_dir = editing_dir or paths.EDITING
        full = io.read(f"{editing_dir}/{sample}_full_canonical.csv")
        hgvsp = (full.groupby(["gRNA_id", "HGVSp"], as_index=False)["#Reads"].sum()
                 .sort_values(["gRNA_id", "#Reads"], ascending=[True, False])
                 .reset_index(drop=True))
        hgvsp["%Reads"] = (hgvsp["#Reads"]
                           / hgvsp.groupby("gRNA_id")["#Reads"].transform("sum") * 100)
    totals = hgvsp.groupby("gRNA_id")["#Reads"].sum()
    sav = hgvsp.copy()
    sav["HGVSp"] = sav["HGVSp"].str.split("_")
    sav = sav.explode("HGVSp")
    sav = sav.groupby(["gRNA_id", "HGVSp"], as_index=False)["#Reads"].sum()
    sav["%Reads"] = 100 * sav["#Reads"] / sav["gRNA_id"].map(totals)
    sav["Codon"] = [h if h == "WT" else first_number(h) for h in sav["HGVSp"]]
    return (sav.sort_values(["gRNA_id", "#Reads"], ascending=[True, False])
            .reset_index(drop=True))


def top_hgvsp(edit: pd.DataFrame) -> pd.DataFrame:
    """Most abundant non-WT variant per gRNA, with its codon and the WT fraction.

    A gRNA whose reads are entirely WT keeps HGVSp 'WT' and codon 0, so it stays
    in the table and can be filtered on `WT_perc` rather than disappearing.
    """
    non_wt = edit[edit["HGVSp"] != "WT"]
    wt = (edit[edit["HGVSp"] == "WT"].rename(columns={"%Reads": "WT_perc"})
          [["gRNA_id", "WT_perc"]])
    # `kind="stable"` matters: where two variants tie on %Reads the original
    # `groupby.apply(nlargest)` kept whichever came first in the incoming order,
    # and an unstable sort silently picks a different variant for those guides.
    top = (non_wt.sort_values("%Reads", ascending=False, kind="stable")
           .drop_duplicates("gRNA_id").reset_index(drop=True))
    out = pd.merge(top, wt, on="gRNA_id", how="outer")
    out["HGVSp"] = out["HGVSp"].fillna("WT")
    out = out.fillna(0)
    out["Codon"] = [0 if h == "WT" else first_number(h) for h in out["HGVSp"]]
    return out


def codon_rates(savs: pd.DataFrame, library: pd.DataFrame) -> pd.DataFrame:
    """Per-gRNA mutation rate summed by codon (WT excluded, percentages intact)."""
    n = (savs[savs["Codon"] != "WT"]
         .groupby(["gRNA_id", "Codon"], as_index=False)["%Reads"].sum())
    n = pd.merge(n, library[["gRNA_id", "Gene"]], on="gRNA_id")
    n["Codon"] = n["Codon"].astype(int)
    return n


def edit_dict(samples, editing_dir=None) -> dict:
    """{sample: top-variant table} for the samples a figure needs."""
    return {s: top_hgvsp(hgvsp_from_full_canonical(s, editing_dir)) for s in samples}
