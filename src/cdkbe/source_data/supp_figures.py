"""Supplementary panels that reuse the main-figure machinery.

The twelve supplementary lollipops are the same plot as figures 2a/2h/3a with a
different drug: `09_analysis_and_viz/lollipop_plots_canonical.ipynb`
builds them all from one `plot_def`, so they are specified here as data rather
than written out again. Condition names, display titles and LFC cutoffs are
copied from that notebook's `PLOTS` list verbatim -- in particular

* Senexin B's condition keeps its space (`Senexin B_15000nM`); only the output
  filename replaced it with an underscore, and
* CDK12-IN-2 uses an LFC cutoff of 2.5 where the other subpool-2 drugs use 2,
  and the subpool-3 drugs all use 1.
"""
from __future__ import annotations

from .. import lfc_correlation as _lc
from .screen_figures import PANEL_SPEC, draw

#: supplementary panel -> the same 6-tuple PANEL_SPEC uses
#: (subpool, LFC table key, CBE condition, ABE condition, title, LFC cutoff)
SUPP_SPEC = {
    # --- supplementary 5: the other two subpool-1 CDK9 compounds -------------
    "S5A": ("subpool1", "subpool1_DMSO", "KI-CDK9d-32N_5000nM",
            "KI-CDK9d-32N_5000nM", "KI-CDK9d-32N 5,000nM", 2),
    "S5B": ("subpool1", "subpool1_DMSO", "KB-0742_1500nM",
            "KB-0742_1500nM", "KB-0742 1,500nM", 2),
    # --- supplementary 8: CDK8/19 ------------------------------------------
    "S8a": ("subpool1", "subpool1_DMSO", "SEL120_4000nM",
            "SEL120_4000nM", "SEL120 4,000nM", 2),
    "S8b": ("subpool1", "subpool1_DMSO", "Senexin B_15000nM",
            "Senexin B_15000nM", "Senexin B 15,000nM", 2),
    # --- supplementary 10: the other two CDK12/13 compounds ------------------
    "S10a": ("subpool2", "subpool2_DMSO", "CBE_BSJ-4-116", "ABE_BSJ-4-116",
             "BSJ-4-116", 2),
    "S10b": ("subpool2", "subpool2_DMSO", "CBE_CDK12-IN-2", "ABE_CDK12-IN-2",
             "CDK12-IN-2", 2.5),
    # --- supplementary 12: the six CDK4/6 inhibitors -------------------------
    "S12a": ("subpool3", "subpool3_DMSO", "CBE_RIBO", "ABE_RIBO", "RIBO", 1),
    "S12b": ("subpool3", "subpool3_DMSO", "CBE_PALBO", "ABE_PALBO", "PALBO", 1),
    "S12c": ("subpool3", "subpool3_DMSO", "CBE_ABEMA", "ABE_ABEMA", "ABEMA", 1),
    "S12d": ("subpool3", "subpool3_DMSO", "CBE_ATIRMO", "ABE_ATIRMO", "ATIRMO", 1),
    "S12e": ("subpool3", "subpool3_DMSO", "CBE_INX", "ABE_INX", "INX", 1),
    "S12f": ("subpool3", "subpool3_DMSO", "CBE_TAGTO", "ABE_TAGTO", "TAGTO", 1),
}

# `draw` reads PANEL_SPEC by key, so register the supplementary panels there
# rather than duplicating the loader.
PANEL_SPEC.update(SUPP_SPEC)


def _panel(key):
    def drawer(outdir=None, **kw):
        return draw(key, outdir=outdir, **kw)
    drawer.__name__ = f"supp_{key}"
    drawer.__doc__ = f"Supplementary {key[1:]}: lollipop, {SUPP_SPEC[key][4]}."
    return drawer


FIGURES = {k: _panel(k) for k in SUPP_SPEC}


def source(key):
    """-> {editor: the frame the panel plots}, without keeping the figure."""
    import matplotlib.pyplot as plt

    fig, kept = draw(key)
    plt.close(fig)
    return kept


PANELS = {k: (lambda key=k: source(key)) for k in SUPP_SPEC}


# ------------------------------------------------------------ supplementary 4
#: panel -> the subpool whose two editor clustermaps it shows. Each panel is a
#: pair (ABE and CBE); the mapping lists both PDFs against one letter.
S4_SUBPOOL = {"S4a": "subpool1", "S4b": "subpool2", "S4c": "subpool3"}


def _s4_figures(subpool):
    def drawer(outdir=None, **kw):
        """-> {editor: (figure, correlation matrix)}."""
        return {ed: _lc.panel(ed, subpool, outdir=outdir)
                for ed in _lc.EDITORS}
    drawer.__name__ = f"supp_S4_{subpool}"
    drawer.__doc__ = (f"Supplementary 4: LFC correlation clustermaps, "
                      f"{subpool}, one per editor.")
    return drawer


def _s4_source(subpool):
    def src():
        """-> {editor: correlation matrix}."""
        return {ed: _lc.correlation(ed, subpool)[0] for ed in _lc.EDITORS}
    return src


FIGURES.update({k: _s4_figures(sp) for k, sp in S4_SUBPOOL.items()})
PANELS.update({k: _s4_source(sp) for k, sp in S4_SUBPOOL.items()})
