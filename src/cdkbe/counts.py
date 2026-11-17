"""Per-sample sequencing counts -> per-screen count and RPM tables.

From count_tables.ipynb, where the
same four lines were repeated for each of the nine screen/editor combinations.
The combinations are a table here (`SCREENS`) so adding one is a row.

Three count flavours come out of each screen, and they differ in how strictly a
read had to match:

    barcode   reads assigned by the sensor barcode
    proto     reads assigned by the protospacer
    matched   reads where barcode and protospacer agree

Most analyses use `barcode`; the compound-mutant screen uses `proto` (its
guides share barcodes).
"""
from __future__ import annotations
import os
from typing import Iterable

import numpy as np
import pandas as pd

from . import io, paths

#: (output name, Screen_data_location, editors to include, subpool, extra samples)
#: `extra` pulls in a sample that physically sat in another screen's run --
#: the plasmid library, and the CBE replicate that was sequenced with the ABE arm.
SCREENS = [
    # no editor filter: this screen's own plasmid sample belongs in its table
    dict(name="CBE_subpool1",    location="CBE_subpool_1", editors=None,     pool="F1-R1",
         extra=("ABE_subpool_1", ["CBE"])),
    dict(name="ABE_subpool1",    location="ABE_subpool_1", editors=["ABE"],  pool="F1-R1",
         extra=("CBE_subpool_1", ["PLASMID"])),
    dict(name="SY_5609_CBE",     location="SY_5609",       editors=["CBE"],  pool="F1-R1",
         extra=("CBE_subpool_1", ["PLASMID"])),
    dict(name="SY_5609_ABE",     location="SY_5609",       editors=["ABE"],  pool="F1-R1",
         extra=("CBE_subpool_1", ["PLASMID"])),
    dict(name="Compound_mutant", location="KB_compound_mut", editors=None,   pool="F1-R1",
         extra=("CBE_subpool_1", ["PLASMID"])),
    dict(name="CDK12_13_CBE",    location="CDK12_13",      editors=["CBE", "PLASMID"], pool="F2-R2", extra=None),
    dict(name="CDK12_13_ABE",    location="CDK12_13",      editors=["ABE", "PLASMID"], pool="F2-R2", extra=None),
    dict(name="CDK2_4_6_CBE",    location="CDK2_4_6",      editors=["CBE", "PLASMID"], pool="F3-R3", extra=None),
    dict(name="CDK2_4_6_ABE",    location="CDK2_4_6",      editors=["ABE", "PLASMID"], pool="F3-R3", extra=None),
]

FLAVOURS = {"barcode": "bc_count", "proto": "total_guide_count", "matched": "matched_guide_count"}


def sample_key() -> pd.DataFrame:
    return io.read_excel(paths.RAW / "sample_key.xlsx")


def load_raw_counts(key: pd.DataFrame | None = None) -> dict:
    """{File_Name: counts}, summing the lanes a sample was split across.

    A sample sequenced on two lanes has two files; their counts add. The
    original used a bare `except` to cover the one-lane case, which would also
    have swallowed a genuinely missing or malformed file -- this checks the
    number of files instead.
    """
    key = sample_key() if key is None else key
    out = {}
    for _, val in key.iterrows():
        loc, file_name = val["Screen_data_location"], val["File_Name"]
        d = paths.RAW / "screen_counts" / loc
        files = sorted(f for f in os.listdir(d) if file_name in f)
        if not files:
            raise FileNotFoundError(f"{file_name}: no count file under {d}")
        a = io.read(d / files[0])
        for extra in files[1:]:                      # same sample, further lanes
            b = io.read(d / extra)
            for c in FLAVOURS.values():
                a[c] = a[c] + b[c]
        out[file_name] = a
    return out


def samples_for(screen: dict, key: pd.DataFrame) -> pd.DataFrame:
    """The rows of the sample key that belong in one screen's count table."""
    m = key["Screen_data_location"] == screen["location"]
    if screen["editors"] is not None:
        m &= key["Editor"].isin(screen["editors"])
    rows = [key[m]]
    if screen["extra"]:
        loc, editors = screen["extra"]
        rows.append(key[(key["Screen_data_location"] == loc) & (key["Editor"].isin(editors))])
    return pd.concat(rows).reset_index(drop=True)


def count_tables(sample_table: pd.DataFrame, raw: dict, library: pd.DataFrame,
                 pool: str = "F1-R1") -> dict:
    """{flavour: table} for one screen, restricted to that subpool's guides.

    Every sample's counts are in the same guide order, which the assert below
    checks -- the tables are assembled column-wise, so a reordered sample would
    silently misalign every guide.
    """
    guides = set(library.loc[library["Pool"] == pool, "gRNA_id"])
    order, cols = None, {f: [] for f in FLAVOURS}
    names = []
    for _, val in sample_table.iterrows():
        sub = raw[val["File_Name"]]
        sub = sub[sub["Guide_ID"].isin(guides)]
        ids = list(sub["Guide_ID"])
        if order is None:
            order = ids
        elif ids != order:
            raise ValueError(f"{val['Sample']}: guide order differs from the first sample")
        for f, col in FLAVOURS.items():
            cols[f].append(list(sub[col]))
        names.append(val["Sample"])
    return {f: pd.DataFrame(dict(zip(["gRNA_id"] + names, [order] + v)))
            for f, v in cols.items()}


def rpm(table: pd.DataFrame, pseudocount: int = 0) -> pd.DataFrame:
    """Reads per million, per sample column.

    No pseudocount by default: adding one here would make it impossible to
    filter on genuinely-zero guides downstream, so it is added at the LFC step.
    """
    out = table.copy()
    for col in table.columns[1:]:
        v = np.asarray(table[col], dtype=float) + pseudocount
        out[col] = 1e6 * v / v.sum()
        assert np.isclose(out[col].sum(), 1e6), f"{col}: RPM does not sum to 1e6"
    return out


def library() -> pd.DataFrame:
    return io.read(paths.REFERENCE / "library" / "CDK_library_final.csv")
