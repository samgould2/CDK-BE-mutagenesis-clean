#!/usr/bin/env python3
"""Check that each panel's source sheet carries what the panel prints.

    python3 scripts/audit_source_data.py             # main figures
    python3 scripts/audit_source_data.py supp        # supplementary
    python3 scripts/audit_source_data.py all

A workbook can be complete in the sense that every panel has a sheet and still
be useless: a sheet of p-values reproduces nothing if the bar heights are
somewhere else, and a column of `CDK6_G53S` does not explain a panel labelled
`CDK6 G53S_M54I`. So this draws each panel, pulls out every word it prints,
and reports the ones that appear nowhere in its sheets.

**Reading the output.** Most of what it reports is fine and always will be:
axis titles (`Median LFC (vs. DMSO)`), legend headings (`Palbociclib:`),
section words (`Controls`, `Intron`), and counts computed from the sheet and
printed on the figure (`343 up`). Those are wording, not data. What matters is
a **category or series label** in the list -- a cell line, a variant, a
condition, a gene -- because that is a row the sheet does not have.

The check is deliberately crude: it looks for each word as a substring of some
cell or column name, case-insensitively, and only sees text columns, so a
label spelled `CDK2_S261G` in the sheet and `CDK2 S261G` on the panel reads as
missing. It is a prompt to go and look, not a test, which is why it is a
script rather than part of the suite.

It found two real gaps when it was written: figure 4k and 4l shipped only
their significance table, not the gate frequencies they plot, and figure 1f
printed each gene's length without that length being in its sheet. Both are
fixed; see `docs/DISCREPANCIES.md` item 34.
"""
from __future__ import annotations

import io
import re
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

warnings.filterwarnings("ignore")
import matplotlib                                              # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
import pandas as pd                                            # noqa: E402

from cdkbe.source_data import DRAWN, SUPP_DRAWN                # noqa: E402

TABLES = ROOT / "source_data_tables"
_BOOKS: dict = {}


def workbook(panel: str):
    supp = re.match(r"S(\d+)", panel)
    stem = (f"Supplementary_Figure_{supp.group(1)}" if supp
            else f"Figure_{panel[0]}")
    if stem not in _BOOKS:
        _BOOKS[stem] = pd.ExcelFile(TABLES / f"{stem}.xlsx")
    return _BOOKS[stem]


def sheet_words(panel: str) -> set:
    """Every column name and text value in the sheets belonging to a panel."""
    book = workbook(panel)
    index = book.parse("INDEX")
    names = [str(s) for s, p in zip(index["sheet"], index["panel"].astype(str))
             if p == panel]
    words = set()
    for name in names:
        frame = book.parse(name)
        words |= {str(c).lower() for c in frame.columns}
        for column in frame.columns:
            if frame[column].dtype == object:
                words |= {str(v).lower() for v in frame[column].dropna().unique()}
    return {w for w in words if len(w) > 1}


def panel_words(drawer) -> list:
    """Every string the panel prints, in drawing order."""
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        result = drawer()
        figure = result[0] if isinstance(result, tuple) else result
        if not hasattr(figure, "savefig"):
            return None                 # a render, not a matplotlib figure
        buf = io.StringIO()
        figure.savefig(buf, format="svg")
        plt.close("all")
    return [t.strip() for t in
            re.findall(r"<text[^>]*>([^<]*)</text>", buf.getvalue()) if t.strip()]


def audit(drawn: dict) -> int:
    worst = 0
    for panel in sorted(drawn):
        try:
            printed = panel_words(drawn[panel])
        except Exception as exc:
            print(f"{panel:6s} SKIP   {type(exc).__name__}: {str(exc)[:70]}")
            continue
        if printed is None:
            print(f"{panel:6s} -      not a matplotlib figure")
            continue
        try:
            words = sheet_words(panel)
        except Exception as exc:
            print(f"{panel:6s} SHEET  {type(exc).__name__}: {str(exc)[:60]}")
            continue
        missing = sorted({t for t in printed
                          if re.search(r"[A-Za-z]{2}", t)
                          and not any(t.lower() in w or w in t.lower()
                                      for w in words)})
        worst = max(worst, len(missing))
        print(f"{panel:6s} {len(printed):4d} printed, {len(missing):3d} not in sheet"
              + (f"  -> {missing[:8]}" if missing else ""))
    return worst


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "main"
    if which in ("main", "all"):
        print("== main figures")
        audit(DRAWN)
    if which in ("supp", "all"):
        print("== supplementary figures")
        audit(SUPP_DRAWN)
