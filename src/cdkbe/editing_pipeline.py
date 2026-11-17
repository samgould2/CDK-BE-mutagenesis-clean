"""Per-sample sensor editing calls -> the per-screen editing tables.

From sensor_HGVSp_pipeline.ipynb, which is
a working log rather than a pipeline: it also contains the library design, the
inputs for the cluster CRISPResso jobs and several exploratory sections, it
defines `extract_numbers` four times and `HGVSp_combiner` twice, and the lines
that actually save the output tables are commented out. Only the reduction from
per-sample calls to the shipped tables is ported here.

Four tables come out of each screen:

    <screen>_full                every distinct edit, with its DNA change
    <screen>_full_canonical      the same, keeping only canonical edits
    <screen>_HGVSp               reads per (gRNA, protein consequence)
    <screen>_HGVSp_canonical     the same, canonical only

"Canonical" means `Canonical_edit` is true -- the edit is the one the editor
makes. `Canonical_window` records whether it also fell in the expected window,
but is NOT part of the filter.
"""
from __future__ import annotations
import os
from pathlib import Path

import pandas as pd

from . import counts, io, paths

#: the columns that identify one distinct edit
EDIT_KEYS = ["Edit", "HGVSp", "Num_edits", "DNA Change",
             "Canonical_edit", "Canonical_window", "gRNA_id"]

#: screen -> (raw editing directory, subpool) for the shipped editing tables
SCREEN_POOL = {
    "CBE_subpool1": ("CBE_subpool_1", "F1-R1"),
    "ABE_subpool1": ("ABE_subpool_1", "F1-R1"),
    "CBE_CDK12_13": ("CDK12_13", "F2-R2"),
    "ABE_CDK12_13": ("CDK12_13", "F2-R2"),
    "CBE_CDK2_4_6": ("CDK2_4_6", "F3-R3"),
    "ABE_CDK2_4_6": ("CDK2_4_6", "F3-R3"),
    "CBE_KB_compound_mut": ("KB_compound_mut", "F1-R1"),
}


def sample_files(screen: str, key: pd.DataFrame | None = None,
                 include_plasmid: bool = False) -> list[Path]:
    """The per-sample editing files that belong to one screen's table."""
    key = counts.sample_key() if key is None else key
    location, _ = SCREEN_POOL[screen]
    editor = screen.split("_")[0]
    rows = key[(key["Screen_data_location"] == location)]
    rows = rows[rows["Editor"].isin([editor] + (["PLASMID"] if include_plasmid else []))]
    d = paths.RAW / "screen_editing" / location
    present = sorted(os.listdir(d))
    out = []
    for fn in rows["File_Name"]:
        out += [d / f for f in present if fn in f]
    return out


def combine(files) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum read counts across samples -> (full edit table, per-HGVSp table)."""
    frames = [io.read(f) for f in files]
    full = (pd.concat(frames).groupby(EDIT_KEYS).sum().reset_index()
            .sort_values(["gRNA_id", "#Reads"], ascending=[True, False])
            .reset_index(drop=True))
    hgvsp = (full[["gRNA_id", "HGVSp", "#Reads"]]
             .groupby(["gRNA_id", "HGVSp"], as_index=False).sum()
             .sort_values(["gRNA_id", "#Reads"], ascending=[True, False])
             .reset_index(drop=True))
    return full, hgvsp


def add_percentages(df: pd.DataFrame, guides) -> pd.DataFrame:
    """%Reads within each gRNA, restricted to the given (targeting) guides."""
    sub = df[df["gRNA_id"].isin(set(guides))].copy()
    sub["%Reads"] = 100 * sub["#Reads"] / sub.groupby("gRNA_id")["#Reads"].transform("sum")
    return sub.reset_index(drop=True)


def targeting_guides(screen: str, library: pd.DataFrame | None = None) -> list[str]:
    library = counts.library() if library is None else library
    _, pool = SCREEN_POOL[screen]
    m = (library["Pool"] == pool) & (library["classification"] == "targeting")
    return library.loc[m, "gRNA_id"].tolist()


def build(screen: str, key: pd.DataFrame | None = None,
          library: pd.DataFrame | None = None) -> dict:
    """-> {'full', 'full_canonical', 'HGVSp', 'HGVSp_canonical'} for one screen."""
    files = sample_files(screen, key)
    if not files:
        raise FileNotFoundError(f"{screen}: no per-sample editing files found")
    full, hgvsp = combine(files)
    guides = targeting_guides(screen, library)
    full = add_percentages(full, guides)
    hgvsp = add_percentages(hgvsp, guides)

    canon_rows = full[full["Canonical_edit"]].reset_index(drop=True)
    canon_full = add_percentages(canon_rows.drop(columns="%Reads"), guides)
    canon_hgvsp = (canon_full[["gRNA_id", "HGVSp", "#Reads"]]
                   .groupby(["gRNA_id", "HGVSp"], as_index=False).sum()
                   .sort_values(["gRNA_id", "#Reads"], ascending=[True, False])
                   .reset_index(drop=True))
    canon_hgvsp = add_percentages(canon_hgvsp, guides)
    return {"full": full, "full_canonical": canon_full,
            "HGVSp": hgvsp, "HGVSp_canonical": canon_hgvsp}
