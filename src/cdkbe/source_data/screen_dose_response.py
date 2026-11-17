"""Dose-response curves for the compounds used in each screen (supp 2a-2c).

From `FACS_dose_response_2024_8_22.ipynb` cells 157-168. One grid of 4PL fits
per subpool, each panel a drug with its LC50 marked.

**Each of the three uses a different background convention**, because the
notebook redefines `log_logistic_4_param` between them. This is not cosmetic --
it moves every viability value -- so the convention is named per panel rather
than shared:

* 2a (cell 158) subtracts **no background at all** (`no_cell_avg = 0`);
* 2b and 2c (cell 162) take it from wells whose **Cell Line** says there are no
  cells, where cell 125 -- used by figure 3d and supp 9d-g -- takes it from
  wells whose **Drug** is "No Cells".

2a and 2b reproduce every published LC50 exactly. 2c reproduces three of six;
see docs/DISCREPANCIES.md item 15.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import dose_response as dr

SPEC = {
    "S2a": dict(plate="2024_09_27.csv", cell_type="A549", background="none",
                a549_only=True, grid=(2, 3), figsize=(12, 6),
                drugs=["KI-CDK9d-32", "KI-CDK9d-32N", "KB-0742",
                       "SY-5609", "SEL120", "Senexin B"]),
    "S2b": dict(plate="2025_4_9_csv.csv", cell_type="A549",
                background="cell_line_no_cells", a549_only=False,
                grid=(1, 3), figsize=(12, 3.2),
                drugs=["BSJ-4-116", "CDK12-IN-2", "HQ461"]),
    "S2c": dict(plate="2025_9_9_CDK246.csv", cell_type="NT CBE",
                background="cell_line_no_cells", a549_only=False,
                grid=(2, 3), figsize=(12, 6),
                drugs=["Palbociclib", "Ribociclib", "Abemaciclib",
                       "Atirmociclib", "INX-315", "Tagtociclib"]),
}


def plate(key: str) -> pd.DataFrame:
    spec = SPEC[key]
    d = dr.plate(spec["plate"])
    if spec["a549_only"]:
        # cell 157: the plate carries A2780 too, and blanks are read as zero
        d = d[d["Cell Line"] == "A549"].fillna(0)
    return d


def parameters(key: str) -> pd.DataFrame:
    """Fitted 4PL parameters and LC50 per drug -- the panel's source data."""
    spec = SPEC[key]
    d = plate(key)
    rows = []
    for drug in spec["drugs"]:
        a, b, c, dd, r2, subset = dr.fit_with(d, spec["cell_type"], drug,
                                              spec["background"])
        rows.append(dict(Drug=drug, a=a, b=b, c_molar=c, d=dd,
                         LC50_nM=c * 1e9, r_squared=r2, n_wells=len(subset)))
    return pd.DataFrame(rows)


def draw(key: str, outdir=None, **kw):
    """-> (figure, the fitted parameters)."""
    spec = SPEC[key]
    d = plate(key)
    nrow, ncol = spec["grid"]
    fig, ax = plt.subplots(nrow, ncol, figsize=spec["figsize"],
                           sharey=True, sharex=True, squeeze=False)
    x_extrap = np.linspace(0, 1e-4, 1_000_000)
    rows = []
    for i, drug in enumerate(spec["drugs"]):
        a, b, c, dd, r2, subset = dr.fit_with(d, spec["cell_type"], drug,
                                              spec["background"])
        rows.append(dict(Drug=drug, a=a, b=b, c_molar=c, d=dd,
                         LC50_nM=c * 1e9, r_squared=r2, n_wells=len(subset)))
        axis = ax[i // ncol][i % ncol]
        axis.plot(x_extrap, dr.func2(x_extrap, a, b, c, dd))
        axis.set_title(f"{drug}, $r^2=${np.round(r2, 2)}", fontsize=15)
        x = np.asarray(subset["Concentration"]) * 1e-9
        y = 100 * np.asarray(subset["Viability"])
        axis.errorbar(x[::3],
                      [np.average(y[j:j + 3]) for j in range(0, len(y), 3)],
                      yerr=[np.std(y[j:j + 3]) for j in range(0, len(y), 3)],
                      ls="none", c="black", fmt="o", markersize=5, capsize=7)
        axis.plot([c, c], [0, 115], linestyle="dashed", color="black",
                  label=f"LC50:\n{np.round(c * 1e9, 0)} nM")
        axis.legend(fontsize=13, loc="lower left")
        axis.set_xscale("log")
        axis.set_xlim(1e-10, 1e-4)
        axis.set_ylim(0, 115)
        axis.grid()
        axis.set_xticks([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4])
        axis.set_xticklabels([.1, 1, 10, 100, 1000, "1e4", "1e5"])
        axis.spines[["right", "top"]].set_visible(False)
        axis.tick_params(axis="both", which="major", labelsize=15)
    for r in range(nrow):
        ax[r][0].set_ylabel("Viability %", fontsize=15)
    for cidx in range(ncol):
        ax[nrow - 1][cidx].set_xlabel("Dose (nM)", fontsize=15)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / f"Supp_{key[1:]}_DR_curves.pdf"), transparent=True)
    return fig, pd.DataFrame(rows)


def _wrap(panel):
    sp = SPEC[panel]
    f = lambda **kw: draw(panel, **kw)          # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    f.__doc__ = (f"4PL dose-response grid in {sp['cell_type']} from plate "
                 f"{sp['plate']}, LC50 marked: {', '.join(sp['drugs'])}. "
                 f"Background convention: {sp['background']}.")
    return f


FIGURES = {k: _wrap(k) for k in SPEC}
PANELS = {k: (lambda _k=k: parameters(_k)) for k in SPEC}
