"""CDK12/13 variant validation MTT (figure 3d).

From `source_data/04_MTT_data/FACS_dose_response_2024_8_22.ipynb` cells 136,
142 and 143: the plate `2025_8_31_CDK12_13.csv`, the extra-sum-of-squares LC50
test, and the 3 x 6 grid of dose-response panels.

Each panel is one variant (red) against the NT line made with the **same
editor** (black) -- ABE variants against `NT ABE`, CBE variants against
`NT CBE` -- with the LC50 fold shift and its Holm-adjusted significance in the
corner. Statistics live in `cdkbe.dose_response.lc50_table`.

The BSJ-4-116 row uses a tenfold lower x-range than the other two drugs; it is
the most potent of the three and its curves would otherwise sit against the
left spine.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .. import dose_response as dr

PLATE = "2025_8_31_CDK12_13.csv"
DRUGS = ["BSJ-4-116", "CDK12-IN-2", "HQ461"]

#: the six variants, in panel order, each with its matched NT control
CONTROL_OF = {
    "CDK12 I733V (3)": "NT ABE",
    "CDK12 G731K (4)": "NT CBE",
    "CDK12 M816I_D817N (5)": "NT CBE",
    "CDK12 G822K (6)": "NT CBE",
    "CDK12 G822E (7)": "NT CBE",
    "CDK13 I711V (8)": "NT ABE",
}
CELL_TYPES = list(CONTROL_OF) + ["NT ABE", "NT CBE"]
TITLES = {k: k.rsplit(" (", 1)[0].replace(" ", "\n") for k in CONTROL_OF}


def plate():
    return dr.plate(PLATE)


def stats():
    """LC50 fold shift and Holm-adjusted F test, one row per variant x drug."""
    return dr.lc50_table(plate(), DRUGS, CELL_TYPES, CONTROL_OF)


def _points(subset):
    """Replicate triples -> (x, mean, sd), in the order the wells were read."""
    x = np.asarray(subset["Concentration"]) * 1e-9
    y = 100 * np.asarray(subset["Viability"])
    return (x[::3],
            [np.average(y[i:i + 3]) for i in range(0, len(y), 3)],
            [np.std(y[i:i + 3]) for i in range(0, len(y), 3)])


def figure_3d(outdir=None, **kw):
    """-> (figure, statistics table)."""
    d = plate()
    st = stats()
    fits = {(drug, ct): dr.fit(d, ct, drug) for drug in DRUGS for ct in CELL_TYPES}
    x_extrap = np.linspace(0, 1e-4, 1_000_000)

    fig, ax = plt.subplots(3, 6, figsize=(15, 7), sharex=False, sharey=True)
    for k, drug in enumerate(DRUGS):
        for i, variant in enumerate(CONTROL_OF):
            a = ax[k][i]
            for who, colour in ((variant, "tab:red"), (CONTROL_OF[variant], "black")):
                *p, r2, subset = fits[(drug, who)]
                a.plot(x_extrap, dr.func2(x_extrap, *p), c=colour, linewidth=3,
                       label=f"$r^{{2}}=$ {np.round(r2, 2)}")
                xs, mu, sd = _points(subset)
                a.errorbar(xs, mu, yerr=sd, ls="none", c=colour, fmt="o",
                           markersize=5, capsize=7)
            a.set_xscale("log")
            a.set_ylim(0, 105)
            a.set_xticks([1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
            a.set_xticklabels(["1", "10", "100", "1000", "10K"])
            # BSJ-4-116 is the most potent of the three; its curves need a
            # decade lower range than the other two rows
            a.set_xlim(*((5e-9, 2e-6) if k == 0 else (5e-8, 2e-5)))
            a.spines[["right", "top"]].set_visible(False)
            a.tick_params(axis="both", which="major", labelsize=15)
            a.set_xlabel("")
            if k == 0:
                a.set_title(TITLES[variant], fontsize=18)
            row = st[(st["Drug"] == drug) & (st["Variant"] == variant)].iloc[0]
            a.text(0.97, 0.97, f"{row['fold_shift']:.1f}×\n{row['sig']}",
                   transform=a.transAxes, ha="right", va="top",
                   fontsize=14, linespacing=1.1)
        ax[k][0].set_ylabel(f"{drug}\nViability %", fontsize=18)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_3d_CDK12_13_validation_MTT.pdf"),
                    transparent=True)
    return fig, st


FIGURES = {"3d": figure_3d}
