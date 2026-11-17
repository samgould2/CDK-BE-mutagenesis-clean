"""Figure 4g and supplementary 13a: allele frequency across palbociclib selection.

Eight edited T47D (4g) or MCF7 (S13a) lines, each sequenced before and after
palbociclib selection. The same measurement as figure 4f, in a different cell
line and against a different drug, so `allele_frequency.py` is the neighbouring
module and the tables share its schema.

**Where the numbers come from.** The tables here were read back out of the
published figures, which are vector files written by matplotlib against a
labelled linear axis and so invert exactly. `scripts/recover_selection_tables.py`
does the reading, and will not write anything until the same inversion has
reproduced the T47D cell-cycle panels' bar heights, which this repository knows
independently from their FlowJo export. Fifteen of the sixteen values come back
as whole percentages, the same form figure 4f's table takes.

**What this module does not claim.** The style is read off the published SVG --
figure size, colours, marker and line widths, font sizes, legend placement --
so the redraw sits beside the original rather than replacing it, and the
published image stays under `data/raw/collaborator_panels/T47D_MCF7/`. These
are recovered values: a primary quantification of the Sanger traces would
replace them.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from .. import io, paths

#: panel -> (the recovered table, the cell line, the file the drawer writes)
SPEC = {
    "4g": (paths.PROCESSED / "validation" / "allele_frequency_4g_T47D.csv",
           "T47D", "Figure_4g_T47D_allele_frequency.pdf"),
    "S13a": (paths.PROCESSED / "validation" / "allele_frequency_S13a_MCF7.csv",
             "MCF7", "Supp_13a_MCF7_allele_frequency.pdf"),
}

#: the published panel's own colours, in its own order, read from the artwork
PALETTE = ["#534ab7", "#0f6e56", "#378add", "#d85a30",
           "#ba7517", "#993556", "#639922", "#d4537e"]

#: likewise its geometry and type sizes, in points and inches
FIGSIZE = (3.37, 3.30)
LINEWIDTH = 0.4
MARKERSIZE = 6.0          # the artwork's marker paths have radius 3pt
MARKEREDGEWIDTH = 0.5
TICK_FONTSIZE = 9.5
LABEL_FONTSIZE = 10.0
LEGEND_FONTSIZE = 7.5
LEGEND_TITLE = "Palbociclib:"
XTICKLABELS = ["Pre-\nSelection", "Post-\nSelection"]
YLABEL = "Allele Frequency (%)"
#: measured off both published panels, which agree to the last decimal
XLIM = (-0.15, 1.15)
YLIM = (0.0, 102.0)


def source(panel: str) -> pd.DataFrame:
    """One row per variant: its allele frequency before and after selection."""
    return io.read(SPEC[panel][0])


def panel(name: str, outdir=None, **_):
    """-> (figure, the table it plots)."""
    table, _, stem = SPEC[name]
    d = source(name).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    for i, row in d.iterrows():
        ax.plot([0, 1], [row["Pre-Selection"], row["Post-Selection"]],
                color=PALETTE[i % len(PALETTE)], linewidth=LINEWIDTH,
                marker="o", markersize=MARKERSIZE,
                markeredgecolor="black", markeredgewidth=MARKEREDGEWIDTH,
                label=row["Allele"])

    ax.set_xticks([0, 1])
    ax.set_xticklabels(XTICKLABELS, fontsize=TICK_FONTSIZE)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_yticks(range(0, 101, 10))
    ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
    ax.set_ylabel(YLABEL, fontsize=LABEL_FONTSIZE)
    ax.spines[["right", "top"]].set_visible(False)
    # the legend sits outside the axes on the right, as in the published panel
    ax.legend(title=LEGEND_TITLE, fontsize=LEGEND_FONTSIZE,
              title_fontsize=LEGEND_FONTSIZE, loc="upper left",
              bbox_to_anchor=(1.02, 1.0), frameon=False,
              handlelength=1.6, labelspacing=0.35, borderpad=0.0)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / stem), transparent=True, bbox_inches="tight")
    return fig, d


def _wrap(name: str):
    def drawer(outdir=None, **kw):
        return panel(name, outdir=outdir, **kw)
    drawer.__name__ = drawer.__qualname__ = f"draw_{name}"
    line = SPEC[name][1]
    drawer.__doc__ = (f"Allele frequency before and after palbociclib selection "
                      f"in {line}, one line per edited variant.")
    return drawer


MAIN_PANELS = {"4g": lambda: source("4g")}
MAIN_FIGURES = {"4g": _wrap("4g")}

#: the supplementary registry rejects a main-figure key, so `PANELS` and
#: `FIGURES` here carry only supplementary 13a and 4g stays in its own pair
PANELS = {"S13a": lambda: source("S13a")}
FIGURES = {"S13a": _wrap("S13a")}
PARTS = {k: [SPEC[k][2].rsplit(".", 1)[0]] for k in SPEC}
