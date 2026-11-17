"""Figure 4f: allele frequency before and after selection, per validated variant.

From `09_analysis_and_viz/analysis.ipynb` cell 137, whose
`fig.savefig` names `figures/CDK4_validation_v3.pdf` -- the file
`docs/source_data_mapping_main.csv` lists for figure 4f.

**The numbers were literals in that cell.** They are Sanger/allele-frequency
quantifications that were typed into the notebook as two Python lists; no
upstream table produces them, in this repository or in the original tree, and
nothing recomputes them. They are therefore reproduced here verbatim, as data
rather than as code, in `data/processed/validation/allele_frequency_4f.csv`.

That makes this module a transcription, not a derivation. If the quantification
is ever redone from the `.ab1` traces under
`data/reference/validation_guides/2025_1_15_CDK_validation_guide_seq/`, this
table is what the result has to be checked against.
See `docs/DISCREPANCIES.md` item 20.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from .. import io, paths

TABLE = paths.PROCESSED / "validation" / "allele_frequency_4f.csv"


def source() -> pd.DataFrame:
    """One row per variant: its allele frequency before and after selection."""
    return io.read(TABLE)


def draw(outdir=None, **kw):
    """Allele frequency before and after selection, one line per variant.

    -> (figure, the table it plots).
    """
    d = source()
    fig, ax = plt.subplots(figsize=(3.3, 5))
    cmap = plt.colormaps["tab20"].colors if hasattr(plt, "colormaps") \
        else plt.cm.tab20.colors
    for i, r in d.reset_index(drop=True).iterrows():
        pre, post = r["Pre-Selection"], r["Post-Selection"]
        ax.plot([0, 1], [pre, post], label=r["Allele"], linewidth=2,
                color=cmap[i % len(cmap)])
        ax.scatter([0, 1], [pre, post], edgecolor="black", linewidth=1, s=100,
                   color=cmap[i % len(cmap)])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Pre-\nSelection", "Post-\nSelection"])
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.set_xlim(-0.15, 1.25)
    ax.set_ylabel("Allele Frequency (%)", fontsize=15)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "CDK4_validation_v3.pdf"), transparent=True)
    return fig, d


draw.__name__ = draw.__qualname__ = "draw_4f"
draw.__doc__ = ("Allele frequency before and after selection, one line per "
                "validated CDK2/4/6 variant.")

PANELS = {"4f": source}
FIGURES = {"4f": draw}
