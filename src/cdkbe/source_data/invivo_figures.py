"""Draw the in vivo panels (figures 2m, 2n, 4n, 4o).

Thin wrappers over `cdkbe.invivo_panel.combined_panel`, holding the arguments
each published panel was drawn with so the figure notebooks stay short.
"""
from __future__ import annotations
import pandas as pd

from .. import invivo_panel, paths, plotting
from . import invivo

#: the two layout variants the manuscript used
SQUARE = dict(figsize=(6.0, 6.0), box_w=0.38, wrap_labels=True, force_square=True)
WIDE = dict(figsize=(8.4, 6.0), box_w=0.30)

COHORTS = {
    "CDK9 I157V": dict(drug="Atuveciclib", endpoint="Day 5"),
    "CDK7 L18F": dict(drug="SY-5609", endpoint="Day 7"),
}


def _panel(d, cohort, drug, y, ylabel, title, outfile=None, **kw):
    fig, stats = invivo_panel.combined_panel(
        d[d["Cohort"] == cohort], cohort, drug, y=y, ylabel=ylabel, title=title,
        outfile=str(outfile) if outfile else None, **kw)
    if "Cohort" not in stats.columns:
        stats.insert(0, "Cohort", cohort)
    return fig, stats


def figure_2m(outdir=None, **kw):
    """Tumour volume fold change, both cohorts. Returns [(fig, stats), ...]."""
    d = invivo.figure_2m()
    out = []
    for cohort, info in COHORTS.items():
        day = info["endpoint"]
        f = (outdir / f"Figure_2m_{cohort.split()[0]}.pdf") if outdir else None
        out.append(_panel(d, cohort, info["drug"], "FoldChange",
                          f"Tumor volume fold change\n({day} / Day 0)",
                          f"{day} endpoint", outfile=f, ref_line=1, **(kw or SQUARE)))
    return out


def figure_2n(outdir=None, **kw):
    """Ki67, both cohorts."""
    d = invivo.figure_2n()
    out = []
    for cohort, info in COHORTS.items():
        f = (outdir / f"Figure_2n_{cohort.split()[0]}.pdf") if outdir else None
        out.append(_panel(d, cohort, info["drug"], "Ki67 %", "% Ki67+ cells",
                          "Ki67", outfile=f, bottom=0, **(kw or SQUARE)))
    return out


def figure_4n(outdir=None, **kw):
    """CDK6 R46W tumour volume fold change, palbociclib."""
    d = invivo.figure_4n()
    f = (outdir / "Figure_4n_CDK6.pdf") if outdir else None
    return [_panel(d, "CDK6 R46W", "Palbociclib", "FoldChange",
                   "Tumor volume fold change\n(Day 7 / Day 1)", "Day 7 endpoint",
                   outfile=f, wt_label="NT", ref_line=1, **(kw or WIDE))]


def figure_4o(outdir=None, **kw):
    """CDK6 R46W Ki67, palbociclib."""
    d = invivo.figure_4o()
    f = (outdir / "Figure_4o_CDK6.pdf") if outdir else None
    return [_panel(d, "CDK6 R46W", "Palbociclib", "Ki67 %", "% Ki67+ cells", "Ki67",
                   outfile=f, wt_label="NT", bottom=0, **(kw or WIDE))]


FIGURES = {"2m": figure_2m, "2n": figure_2n, "4n": figure_4n, "4o": figure_4o}

#: These four return a LIST of (figure, stats), one entry per cohort, because
#: each published panel is drawn once per tumour model. The list is positional,
#: so `PARTS` names the entries -- in the same order the functions build them,
#: which is `COHORTS` order for 2m/2n and the single CDK6 cohort for 4n/4o.
#: `scripts/render_figures.py` reads this to name the files it writes; without
#: it a reader would be left to tell two cohorts apart by list index.
PARTS = {
    "2m": [c.split()[0] for c in COHORTS],
    "2n": [c.split()[0] for c in COHORTS],
    "4n": ["CDK6"],
    "4o": ["CDK6"],
}
