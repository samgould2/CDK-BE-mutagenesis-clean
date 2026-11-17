"""Allele-specific MTT dose response (supplementary 9d-9g).

From `source_data/04_MTT_data/FACS_dose_response_2024_8_22.ipynb` cells
131-134, all four over the plate loaded in cell 124: `2025_8_10.csv`, with the
`Error` rows dropped. The live fit for those cells is `cdkbe.dose_response`
(the notebook's cell-125 definition); see that module for the normalisation.

Four lines on each curve panel -- NT, CDK9 L156F, CDK7 L18F and the double --
so the CDK9 and CDK7 resistance alleles can be read against each other and
against the double mutant in one experiment. Each point is the mean of three
wells with an SD error bar; the curve is the fit, drawn across a decade wider
than the data so the plateaus are visible.

9g is the same data as 9f shown as a bar chart over the top four KB-0742 doses,
where the curves are flat and the question is whether any line separates.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .. import dose_response as dr

PLATE = "2025_8_10.csv"
CELL_TYPES = ["NT", "L156F", "L18F", "L156F/L18F"]
LABELS = {"NT": "NT gRNA", "L156F": "CDK9 L156F",
          "L18F": "CDK7 L18F", "L156F/L18F": "L156F/L18F"}
COLOURS = {"NT": "black", "L156F": "tab:blue",
           "L18F": "tab:red", "L156F/L18F": "tab:purple"}

#: panel -> (drug, x-limits of the finished axis)
CURVES = {
    "S9d": ("SY-5609", (1e-10, 2e-6)),
    "S9e": ("KI-CDK9d-32", (1e-9, 4e-6)),
    "S9f": ("KB-0742", (9e-8, 10e-6)),
}


#: Figures 2d and 2k, from the same plate and the same cell-125 fit. Each entry
#: is (drug, x-limits, the two lines drawn, y-limits, the notebook's savefig
#: stem), read off cells 126, 127 and 130 -- whose `fig.savefig` names are the
#: files `docs/source_data_mapping_main.csv` lists for those panels, which is
#: what ties these arguments to the published figure.
MAIN_CURVES = {
    "2d": [
        ("KI-CDK9d-32", (1e-9, 1e-5), ["NT", "L156F"], (-5, 110),
         "L156F_KI-CDK9d-32_wide"),
        ("KB-0742", (9e-8, 10e-6), ["NT", "L156F"], (-5, 110),
         "L156F_KB-0742_wide"),
    ],
    "2k": [
        ("SY-5609", (1e-10, 1e-6), ["NT", "L18F"], (40, 105), "L18F_SY_wide"),
    ],
}

#: the slug each part of a main panel is named by, for the rendered file name
MAIN_PARTS = {k: [spec[0] for spec in v] for k, v in MAIN_CURVES.items()}

#: `scripts/render_figures.py` reads PARTS to name the files a panel produces.
#: 2d is one lettered panel built from two drug curves, so it writes two.
PARTS = MAIN_PARTS


def plate() -> pd.DataFrame:
    """The allele-specific plate, with the flagged wells dropped (cell 124)."""
    d = dr.plate(PLATE)
    return d[d["Error"] != True] if "Error" in d.columns else d  # noqa: E712


def curve_panel(drug: str, xlim, outdir=None, cell_types=None,
                ylim=(0, 105), figsize=(5, 6), stem=None):
    """-> (figure, per-line fit parameters).

    `cell_types`, `ylim` and `figsize` default to the four-line supplementary
    panels. Figures 2d and 2k are the SAME plate and the same fits with two of
    the four lines and a shorter canvas (notebook cells 126, 127 and 130), so
    they are this function with different arguments rather than a second
    implementation that could drift from it.
    """
    d = plate()
    cell_types = list(cell_types or CELL_TYPES)
    fig, ax = plt.subplots(figsize=figsize)
    x_extrap = np.linspace(0, 1e-4, 1_000_000)
    rows = []
    for ct in cell_types:
        a, b, c, dd, r2, subset = dr.fit(d, ct, drug)
        ax.plot(x_extrap, dr.func2(x_extrap, a, b, c, dd),
                label=f"{LABELS[ct]}, $r^{{2}}=$ {np.round(r2, 2)}",
                c=COLOURS[ct], linewidth=3)
        x = np.asarray(subset["Concentration"]) * 1e-9
        y = 100 * np.asarray(subset["Viability"])
        # wells come in replicate triples, in order; average each triple
        ax.errorbar(x[::3],
                    [np.average(y[i:i + 3]) for i in range(0, len(y), 3)],
                    yerr=[np.std(y[i:i + 3]) for i in range(0, len(y), 3)],
                    ls="none", c=COLOURS[ct], fmt="o", markersize=5, capsize=7)
        rows.append({"Cell Line": ct, "Drug": drug, "a": a, "b": b, "c_molar": c,
                     "d": dd, "LC50_nM": c * 1e9, "r_squared": r2})
    ax.set_xscale("log")
    ax.set_ylabel("Viability %", fontsize=15)
    ax.set_xlabel(f"{drug} dose (nM)", fontsize=15)
    ax.set_ylim(*ylim)
    ax.set_title(drug, fontsize=15)
    ax.grid()
    ax.set_xticks([1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
    ax.set_xticklabels([1, 10, 100, 1000, 10000])
    ax.set_xlim(*xlim)
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=15)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / f"{stem or f'all_{drug}'}.pdf"), transparent=True)
    return fig, pd.DataFrame(rows)


def bar_panel(outdir=None):
    """9g: KB-0742 viability at the top four doses, as bars."""
    d = plate()
    drug = "KB-0742"
    comb = pd.concat([dr.fit(d, ct, drug)[-1] for ct in CELL_TYPES])
    comb = comb.assign(Viability=100 * comb["Viability"])
    top = comb[comb["Concentration"] >= 1250]
    palette = [COLOURS[c] for c in CELL_TYPES]
    fig, ax = plt.subplots(figsize=(5, 5))
    kw = dict(data=top, x="Concentration", y="Viability", hue="Cell Line",
              palette=palette, legend=False)
    sns.barplot(ax=ax, errorbar="sd", edgecolor="black", linewidth=1, **kw)
    sns.stripplot(ax=ax, dodge=True, edgecolor="white", linewidth=1, s=10, **kw)
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.set_xticklabels([1.25, 1.5, 1.75, 2])     # nM -> micromolar, as published
    ax.set_ylabel("Viability %", fontsize=15)
    ax.set_xlabel(f"{drug} dose " + r"($\mu$M)", fontsize=15)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "all_KB_bar.pdf"), transparent=True)
    return fig, top


def _wrap(panel):
    drug, _ = CURVES[panel]
    f = lambda outdir=None, **kw: curve_panel(*CURVES[panel], outdir)   # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    f.__doc__ = (f"Allele-specific MTT dose response to {drug}: "
                 f"{', '.join(LABELS.values())}.")
    return f


FIGURES = {k: _wrap(k) for k in CURVES}


def _bar(outdir=None, **kw):
    """KB-0742 viability at the top four doses, as bars, one group per cell line."""
    return bar_panel(outdir)


FIGURES["S9g"] = _bar

def _main_wrap(panel):
    def drawer(outdir=None, **kw):
        out = {}
        for drug, xlim, cts, ylim, stem in MAIN_CURVES[panel]:
            out[drug] = curve_panel(drug, xlim, outdir, cell_types=cts,
                                    ylim=ylim, figsize=(5, 3), stem=stem)
        return out
    lines = ", ".join(LABELS[c] for c in MAIN_CURVES[panel][0][2])
    drugs = ", ".join(spec[0] for spec in MAIN_CURVES[panel])
    drawer.__name__ = drawer.__qualname__ = f"draw_{panel}"
    drawer.__doc__ = (f"Allele-specific MTT dose response: {lines}, "
                      f"against {drugs}.")
    return drawer


#: figures 2d and 2k. Kept apart from `FIGURES`/`PANELS`, which
#: `source_data._supp()` merges wholesale into the supplementary registry --
#: a main-figure key in there would be counted as a supplementary panel.
MAIN_FIGURES = {k: _main_wrap(k) for k in MAIN_CURVES}


def main_source(panel: str) -> pd.DataFrame:
    """The points figures 2d / 2k plot, one row per well."""
    frames = []
    for drug, _, cts, _, _ in MAIN_CURVES[panel]:
        frames.append(dr.points_table(plate(), cts, drug))
    return pd.concat(frames, ignore_index=True)


PANELS = {k: (lambda _k=k: dr.points_table(plate(), CELL_TYPES, CURVES[_k][0]))
          for k in CURVES}
MAIN_PANELS = {k: (lambda _p=k: main_source(_p)) for k in MAIN_CURVES}
PANELS["S9g"] = lambda: dr.points_table(plate(), CELL_TYPES, "KB-0742")
