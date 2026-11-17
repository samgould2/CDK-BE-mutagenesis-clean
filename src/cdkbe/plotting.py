"""Figure style and palettes, taken verbatim from the original notebooks.

These values are NOT design choices made here -- they are the published
figures' own colours and fonts, centralised so a change applies everywhere
instead of to whichever notebook was edited last.
"""
from __future__ import annotations
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

#: per-CDK colours (editing_analysis_v2, fig1_analysis)
CDK_COLORS = {
    "CDK7": "cyan", "CDK8": "tab:blue", "CDK19": "cornflowerblue", "CDK9": "darkviolet",
    "CDK12": "tab:green", "CDK13": "darkolivegreen",
    "CDK2": "goldenrod", "CDK4": "tab:red", "CDK6": "tab:orange",
}
#: subpool-3 gene colours (logo_plots_CDK2_4_6_canonical) -- avoids purple/blue,
#: which already mean ABE/CBE, so gene and editor never collide
GENE_COLORS_S3 = {"CDK2": "#ccaa00", "CDK4": "#9e1b1b", "CDK6": "#f2941a"}
EDITOR_COLORS = {"ABE": "tab:purple", "CBE": "tab:blue"}
EDITOR_MARKERS = {"ABE": "o", "CBE": "s"}

CONDITION_PALETTES = {
    "subpool1": ["indigo", "darkviolet", "mediumorchid", "thistle", "magenta"],
    "subpool2": ["darkolivegreen", "tab:green", "lightgreen"],
    "subpool3": ["lightcoral", "indianred", "firebrick", "maroon", "darkorange", "goldenrod"],
}

GENES = {
    "subpool1": ["CDK7", "CDK8", "CDK19", "CDK9"],
    "subpool2": ["CDK12", "CDK13"],
    "subpool3": ["CDK2", "CDK4", "CDK6"],
}


def use_style() -> None:
    """Helvetica with editable PDF text, as every original notebook set."""
    plt.rc("font", family="Helvetica")
    mpl.rcParams["pdf.fonttype"] = 42      # keep text selectable in Illustrator
    mpl.rcParams["svg.fonttype"] = "none"


def save(fig, path, formats=("pdf", "png"), dpi=300, transparent=True, **kw):
    """Write a figure in several formats; returns the paths written."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    stem = p.with_suffix("")
    out = []
    for ext in formats:
        f = Path(f"{stem}.{ext}")
        fig.savefig(f, dpi=dpi, transparent=transparent, **kw)
        out.append(f)
    return out
