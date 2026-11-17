"""The shipped source-data workbooks must be complete and self-describing.

`source_data_tables/*.xlsx` is what a reader gets instead of this repository:
one workbook per figure, one sheet per panel, holding the numbers that panel
plots. That only works if every panel with source data has a sheet, if each
sheet carries the labels as well as the values, and if the INDEX says what the
sheet is. None of that is true by construction -- the workbooks are built by a
script that can silently fall behind the registry, and did: `Figure_4.xlsx`
shipped with seven sheets while thirteen panels had source data, and sheets
built from a matrix lost their row labels to `to_excel(index=False)`.

So this checks the files on disk, not the builder.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from cdkbe.source_data import PANELS, SUPP_PANELS

TABLES = Path(__file__).resolve().parents[1] / "source_data_tables"

def workbook(panel: str) -> Path:
    """The workbook a panel belongs in, from its own name."""
    supp = re.match(r"S(\d+)", panel)
    if supp:
        return TABLES / f"Supplementary_Figure_{supp.group(1)}.xlsx"
    return TABLES / f"Figure_{panel[0]}.xlsx"


@pytest.fixture(scope="module")
def sheets():
    """{workbook stem: {sheet name: frame}} for every shipped workbook."""
    out = {}
    for path in sorted(TABLES.glob("*.xlsx")):
        book = pd.ExcelFile(path)
        out[path.stem] = {name: book.parse(name) for name in book.sheet_names}
    return out


def panels_in(sheets, stem) -> set:
    """Which panels a workbook covers, from its INDEX sheet."""
    index = sheets[stem]["INDEX"]
    return set(index["panel"].astype(str))


def test_every_panel_with_source_data_has_a_sheet(sheets):
    missing = []
    for panel in list(PANELS) + list(SUPP_PANELS):
        stem = workbook(panel).stem
        if stem not in sheets or panel not in panels_in(sheets, stem):
            missing.append(panel)
    assert missing == []


def test_every_sheet_belongs_to_a_registered_panel(sheets):
    known = set(PANELS) | set(SUPP_PANELS)
    stray = []
    for stem, book in sheets.items():
        if not stem.startswith(("Figure_", "Supplementary_Figure_")):
            continue          # the Sanger supplement is not a per-panel book
        stray += [p for p in panels_in(sheets, stem) if p not in known]
    assert stray == []


def test_figure_4_carries_every_panel_this_repository_can_supply(sheets):
    """The case that went stale, pinned by name."""
    covered = panels_in(sheets, "Figure_4")
    assert {"4a", "4b", "4c", "4d", "4e", "4f", "4g", "4k", "4l",
            "4n", "4o", "4p"} <= covered


def test_no_sheet_is_a_block_of_numbers_with_nothing_naming_the_rows(sheets):
    """A panel cannot be redrawn from values whose rows are unlabelled."""
    bare = []
    for stem, book in sheets.items():
        for name, frame in book.items():
            if name == "INDEX" or frame.empty:
                continue
            labelled = [c for c in frame.columns if frame[c].dtype == object]
            if not labelled:
                bare.append(f"{stem}:{name}")
    assert bare == []


def test_no_sheet_carries_a_column_with_no_header(sheets):
    """`Unnamed: 7` is either a dropped index or spreadsheet scratch."""
    unnamed = []
    for stem, book in sheets.items():
        for name, frame in book.items():
            unnamed += [f"{stem}:{name}:{c}" for c in frame.columns
                        if str(c).startswith("Unnamed:")]
    assert unnamed == []


def test_every_index_row_says_what_the_sheet_is(sheets):
    blank = []
    for stem, book in sheets.items():
        if "INDEX" not in book:
            continue
        index = book["INDEX"]
        for _, row in index.iterrows():
            if not str(row.get("description", "")).strip() or \
                    str(row.get("description")) == "nan":
                blank.append(f"{stem}:{row['sheet']}")
    assert blank == []


def test_the_index_row_counts_match_the_sheets_they_describe(sheets):
    for stem, book in sheets.items():
        if "INDEX" not in book:
            continue
        for _, row in book["INDEX"].iterrows():
            frame = book[str(row["sheet"])]
            assert len(frame) == row["rows"], f"{stem}:{row['sheet']}"
            assert len(frame.columns) == row["columns"], f"{stem}:{row['sheet']}"


def test_the_matrix_panels_kept_their_row_labels(sheets):
    """The correlation clustermaps are the case `index=False` was losing."""
    for name in ("S4a_ABE", "S4a_CBE", "S4b_ABE", "S4c_ABE"):
        frame = sheets["Supplementary_Figure_4"][name]
        assert frame.columns[0] == "label"
        assert frame["label"].dtype == object
        # a correlation matrix is square once its labels are counted
        assert len(frame) == len(frame.columns) - 1
