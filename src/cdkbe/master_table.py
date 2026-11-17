"""One long-format table per editor, stacking every screen and baseline.

From master_table_gen.ipynb, where loading
the nine screens x three baselines took 169 lines of near-identical
`pd.read_csv` calls followed by 80 lines assigning `Screen_ID` and `Editor` one
variable at a time. Both are `SOURCES` here.

The result is what `fig1_analysis` and the conservation modules read as
`{ABE,CBE}_master_LFC_table.csv`: one row per guide per baseline, with each
baseline's columns suffixed so they can live side by side.
"""
from __future__ import annotations
import pandas as pd

from . import io, paths

#: (Screen_ID, Editor, table stem, count flavour, strip the editor prefix?)
#:
#: Subpool 1's tables already name their columns without an editor prefix; the
#: others carry `ABE_`/`CBE_` and it is stripped so every screen agrees.
#: The compound-mutant screen uses proto counts, because its guides share
#: barcodes.
SOURCES = [
    ("Subpool1",  "ABE", "ABE_subpool1",    "barcode", False),
    ("Subpool1",  "CBE", "CBE_subpool1",    "barcode", False),
    ("Subpool2",  "ABE", "CDK12_13_ABE",    "barcode", True),
    ("Subpool2",  "CBE", "CDK12_13_CBE",    "barcode", True),
    ("Subpool3",  "ABE", "CDK2_4_6_ABE",    "barcode", True),
    ("Subpool3",  "CBE", "CDK2_4_6_CBE",    "barcode", True),
    ("SY-5609",   "ABE", "SY_5609_ABE",     "barcode", True),
    ("SY-5609",   "CBE", "SY_5609_CBE",     "barcode", True),
    ("Iterative", "CBE", "Compound_mutant", "proto",   False),
]

BASELINES = ("Plasmid", "DMSO", "T0")
KEEP = ["gRNA_id", "classification", "Screen_ID", "Editor"]


def load(stem: str, flavour: str, baseline: str, strip_prefix: bool) -> pd.DataFrame:
    d = io.read(paths.LFC / f"{flavour}_counts" / f"{baseline}_base" /
                f"{stem}_{flavour}_{baseline}_base.csv")
    if strip_prefix:
        d.columns = d.columns.str.replace(r"(ABE_|CBE_)", "", regex=True)
    return d


def tag_columns(d: pd.DataFrame, baseline: str) -> pd.DataFrame:
    """Suffix every per-baseline column so three baselines can sit side by side."""
    cols = ([c for c in d.columns if "_REP" in c]
            + [c for c in d.columns if "_median" in c]
            + [c for c in d.columns if "_avg" in c]
            + [c for c in d.columns if "FDR" in c]
            + ["base_RAW"])
    renamed = d.rename(columns={c: f"{c}_{baseline}" for c in cols})
    out = renamed[KEEP + [f"{c}_{baseline}" for c in cols]].copy()
    out["Base_Value"] = baseline
    return out


def build(editor: str) -> pd.DataFrame:
    """The master table for one editor: every screen, every baseline, stacked."""
    parts = []
    for baseline in BASELINES:
        for screen_id, ed, stem, flavour, strip in SOURCES:
            if ed != editor:
                continue
            d = load(stem, flavour, baseline, strip)
            d["Screen_ID"] = screen_id
            d["Editor"] = ed
            parts.append(tag_columns(d, baseline))
    return pd.concat(parts)


def path(editor: str):
    return paths.LFC / f"{editor}_master_LFC_table.csv"
