"""Figure style matched to Gould et al., *A functional map of CDK-drug interactions*.

The parent manuscript's figures were made with the conventions in
`data/CDK-BE-mutagenesis/notebooks/**` — verified by reading those notebooks rather
than guessed:

    plt.rc('font', family='Helvetica')          every notebook's first cell
    ABE  -> 'tab:purple'   (43 uses)            CBE -> 'tab:blue' (59 uses)
    highlight / hits -> 'tab:red'               ink -> 'black' (230 uses)
    controls -> 'tab:grey'                      diverging -> 'seismic'
    savefig(..., transparent=True)              rasters at dpi=300-400

Panel letters are lowercase, bold, flush left above the axes (a, b, c ...), which is
how the manuscript's Figures 1-7 are lettered. Titles are sentence case, small, and
sit next to the letter. There is no figure-level banner title — the parent figures
carry their explanation in the legend, not on the canvas.

This module deliberately does NOT reproduce the conventions of the earlier
BEAT-seq figures (large header, italic method line under every panel). Those were
useful working annotations; they are not the manuscript's style.

Usage:

    from analysis import style
    style.apply()
    fig, axes = plt.subplots(...)
    style.panel(ax, "a", "CDK4 loss arrests G1 only where CDK4 is essential")
    style.save(fig, "KEY1_variant_function")
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import paths

# --- palette ----------------------------------------------------------------
ABE = "tab:purple"
CBE = "tab:blue"
HIT = "tab:red"
INK = "black"
CTRL = "tab:grey"
MUTED = "#BFBFBF"
BAND = "#D6E4F0"      # shaded reference band, matches the manuscript's ATP-site shading
GOLD = "gold"
DIVERGING = "seismic"

EDITOR_COLOR = {"ABE": ABE, "CBE": CBE}
LINE_MARKER = {"A549": "o", "T47D": "s"}

# Gene order used throughout the parent manuscript: subpool 1, then 2, then 3.
CDK_ORDER = ["CDK7", "CDK8", "CDK9", "CDK19", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]

# The manuscript's own per-CDK colours, read out of
# data/CDK-BE-mutagenesis/notebooks/ (the dict in the structure notebooks, and the hex
# twins in 10_PyMOL/align_catalytic_v2.py). Use these rather than a generic colormap so
# the single-cell figures key to the same colours as Figures 1-7.
GENE_COLOR = {
    "CDK7": "cyan",              # #00ffff
    "CDK8": "tab:blue",          # #1f77b4
    "CDK9": "darkviolet",        # #9400d3
    "CDK19": "cornflowerblue",   # #6495ed
    "CDK12": "tab:green",        # #2ca02c
    "CDK13": "darkolivegreen",   # #556b2f
    "CDK2": "goldenrod",
    "CDK4": "tab:red",
    "CDK6": "tab:orange",
}


def apply(base: float = 9.0) -> None:
    """Set rcParams. `base` is the body font size in points."""
    mpl.rcParams.update({
        "font.family": "Helvetica",
        "font.size": base,
        "axes.titlesize": base + 1,
        "axes.labelsize": base,
        "xtick.labelsize": base - 1,
        "ytick.labelsize": base - 1,
        "legend.fontsize": base - 1,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.grid": False,
        "legend.frameon": False,
        "figure.facecolor": "none",
        "axes.facecolor": "none",
        "savefig.transparent": True,
        # live text in Illustrator, as the parent figures were delivered
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel(ax, letter: str, title: str = "", pad: float = 6.0,
          size: float | None = None) -> None:
    """Lowercase bold panel letter flush left, optional sentence-case title beside it."""
    size = size or mpl.rcParams["axes.titlesize"]
    ax.set_title(f"$\\bf{{{letter}}}$   {title}" if title else f"$\\bf{{{letter}}}$",
                 loc="left", pad=pad, fontsize=size)


def despine(ax, top: bool = True, right: bool = True) -> None:
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)


def save(fig, stem: str, dpi: int = 400) -> None:
    """Write PDF (vector, live text) + PNG, both transparent, as the parent repo did."""
    paths.ensure_dirs()
    for ext in ("pdf", "png"):
        fig.savefig(paths.FIGURES / f"{stem}.{ext}", dpi=dpi,
                    transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {paths.FIGURES / stem}.pdf / .png")
