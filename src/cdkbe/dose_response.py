"""Four-parameter log-logistic dose-response fits for the MTT plates.

From `source_data/04_MTT_data/FACS_dose_response_2024_8_22.ipynb`, which is a
working log rather than a pipeline: `log_logistic_4_param` is redefined eleven
times with different bounds and normalisations, and `allele_specific` is
reassigned nine times. The definition ported here is **cell 125**, which is the
one live for every panel this module draws (the next redefinition is cell 150,
after all of them).

Normalisation, in the original's order:

1. `no_cell_avg` -- the mean of every `No Cells` well *across the whole plate*,
   subtracted as background;
2. `no_drug_avg` -- that cell line's own DMSO wells, minus the same background;
3. viability = (measurement - background) / no_drug_avg.

So each line is normalised to its own vehicle, and the background is shared.

The fit is `curve_fit` on (a-d)/(1+(x/c)^b)+d with concentrations in molar,
starting from [90, 1, .005, 20] and with `a` bounded to [85, 100] -- the upper
asymptote is pinned near full viability so a line that never drops cannot be
fitted with an implausibly low plateau.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from . import paths

MTT = paths.RAW / "mtt_dose_response"

#: starting values and bounds, exactly as the live definition had them
P0 = np.array([90.0, 1.0, 0.005, 20.0])
BOUNDS = ([85, 0, 0, 0], [100, np.inf, np.inf, 100])


def func2(x, a, b, c, d):
    """4PL logistic: upper asymptote a, slope b, EC50 c (molar), lower d."""
    return ((a - d) / (1.0 + ((x / c) ** b))) + d


def plate(name: str) -> pd.DataFrame:
    """One MTT plate, with the long concentration column renamed."""
    return pd.read_csv(MTT / name).rename(
        columns={"Concentration (nM or % for DMSO)": "Concentration"})


def viability(df: pd.DataFrame, cell_type: str) -> pd.DataFrame:
    """Add a `Viability` column for one cell line, normalised to its own DMSO."""
    sub = df[df["Cell Line"] == cell_type].copy()
    no_cell = np.average(df[df["Drug"] == "No Cells"]["Measurement-Ref"])
    no_drug = np.average(sub[sub["Drug"] == "DMSO"]["Measurement-Ref"]) - no_cell
    sub["Viability"] = (sub["Measurement-Ref"] - no_cell) / no_drug
    return sub


def fit(df: pd.DataFrame, cell_type: str, drug: str):
    """-> (a, b, c, d, r_squared, subset) for one cell line and drug."""
    sub = viability(df, cell_type)
    subset = sub[(sub["Drug"] == drug) & (sub["Cell Line"] != "No cells")]
    x = np.asarray(subset["Concentration"]) * 1e-9        # nM -> M
    y = 100 * np.asarray(subset["Viability"])
    params, _ = curve_fit(func2, x, y, P0, bounds=BOUNDS)
    resid = func2(x, *params) - y
    r2 = 1.0 - (np.var(resid) / np.var(y))
    return (*params, r2, subset)


def curve_table(df: pd.DataFrame, cell_types, drug: str) -> pd.DataFrame:
    """Fitted parameters and LC50 per cell line -- the source data for a panel."""
    rows = []
    for ct in cell_types:
        a, b, c, d, r2, subset = fit(df, ct, drug)
        rows.append({"Cell Line": ct, "Drug": drug, "a": a, "b": b, "c_molar": c,
                     "d": d, "LC50_nM": c * 1e9, "r_squared": r2,
                     "n_wells": len(subset)})
    return pd.DataFrame(rows)


def points_table(df: pd.DataFrame, cell_types, drug: str) -> pd.DataFrame:
    """Every well behind a panel: concentration, viability, cell line."""
    out = []
    for ct in cell_types:
        a, b, c, d, r2, subset = fit(df, ct, drug)
        s = subset[["Cell Line", "Drug", "Concentration", "Measurement",
                    "Measurement-Ref", "Viability"]].copy()
        s["Viability_pct"] = 100 * s["Viability"]
        out.append(s)
    return pd.concat(out, ignore_index=True)


# --------------------------------------------------------- LC50 significance

def lc50_ftest(sub_var, sub_ctrl, params_var, params_ctrl):
    """Extra-sum-of-squares F test: do two lines share one LC50?

    Full model: the two separate 4PL fits already plotted (8 parameters).
    Reduced model: a joint fit where log10(LC50) is shared and a/b/d stay
    separate (7 parameters). The difference in residual sum of squares, scaled
    by the full model's, is an F statistic on (1, n-8) degrees of freedom.

    Sharing *log10* LC50 rather than LC50 keeps the shared parameter on the
    scale the concentrations are spaced on, which is what makes the reduced
    fit converge across four orders of magnitude.
    """
    import scipy.stats

    x1 = np.asarray(sub_var["Concentration"]) * 1e-9
    y1 = 100 * np.asarray(sub_var["Viability"])
    x2 = np.asarray(sub_ctrl["Concentration"]) * 1e-9
    y2 = 100 * np.asarray(sub_ctrl["Viability"])

    ss_full = (np.sum((func2(x1, *params_var) - y1) ** 2)
               + np.sum((func2(x2, *params_ctrl) - y2) ** 2))

    n1 = len(x1)
    x_all, y_all = np.concatenate([x1, x2]), np.concatenate([y1, y2])

    def joint(x, a1, b1, d1, a2, b2, d2, logc):
        c = 10 ** logc
        return np.concatenate([func2(x[:n1], a1, b1, c, d1),
                               func2(x[n1:], a2, b2, c, d2)])

    a1, b1, c1, d1 = params_var
    a2, b2, c2, d2 = params_ctrl
    p0 = [a1, b1, d1, a2, b2, d2, np.mean([np.log10(c1), np.log10(c2)])]
    lb = [85, 0, 0, 85, 0, 0, -12]
    ub = [100, np.inf, 100, 100, np.inf, 100, -2]
    p0 = np.clip(p0, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
    p_red, _ = curve_fit(joint, x_all, y_all, p0, bounds=[lb, ub], maxfev=20000)
    ss_red = np.sum((joint(x_all, *p_red) - y_all) ** 2)
    # guard: the separate fits can land in a worse local minimum than the joint one
    ss_red = max(ss_red, ss_full)

    df_full = len(x_all) - 8
    f = ((ss_red - ss_full) / 1) / (ss_full / df_full)
    return f, scipy.stats.f.sf(f, 1, df_full)


def holm(pvals):
    """Holm step-down adjusted p-values."""
    pvals = np.asarray(pvals)
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvals[idx])
        adj[idx] = min(1, running)
    return adj


def stars(p):
    for cut, mark in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")):
        if p < cut:
            return mark
    return "ns"


def lc50_table(df, drugs, cell_types, control_of):
    """Every variant against its own NT control, with Holm-adjusted p-values.

    `control_of` maps each variant to the NT line made with the *same editor*,
    so an ABE variant is compared against NT ABE and a CBE one against NT CBE.
    """
    fits = {(drug, ct): fit(df, ct, drug) for drug in drugs for ct in cell_types}
    rows = []
    for drug in drugs:
        for ct in control_of:
            ctrl = control_of[ct]
            *pv, _, sv = fits[(drug, ct)]
            *pc, _, sc = fits[(drug, ctrl)]
            f, p = lc50_ftest(sv, sc, pv, pc)
            rows.append(dict(Drug=drug, Variant=ct, Control=ctrl,
                             LC50_variant_nM=pv[2] * 1e9, LC50_control_nM=pc[2] * 1e9,
                             fold_shift=pv[2] / pc[2], F=f, p=p))
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out["p"])
    out["sig"] = out["p_holm"].apply(stars)
    return out


# ------------------------------------------------- per-plate normalisations
#
# `FACS_dose_response_2024_8_22.ipynb` redefines `log_logistic_4_param` eleven
# times, and the differences are not cosmetic: they change where the background
# comes from. Three conventions appear among the panels ported here.

#: name -> how that plate's no-cell background is found
BACKGROUND = {
    #: cell 125 (figure 3d, supp 9d-g): wells whose *Drug* is "No Cells"
    "drug_no_cells": lambda d: np.average(d[d["Drug"] == "No Cells"]["Measurement-Ref"]),
    #: cell 158 (supp 2a): no background subtraction at all
    "none": lambda d: 0.0,
    #: cell 162 (supp 2b, 2c): wells whose *Cell Line* says there are none.
    #: The notebook matched "NO CELLS" exactly; 2025_9_9_CDK246.csv spells it
    #: "No cells", so the match is made case-insensitively here -- see
    #: docs/DISCREPANCIES.md item 15.
    "cell_line_no_cells": lambda d: np.average(
        d[d["Cell Line"].astype(str).str.upper() == "NO CELLS"]["Measurement-Ref"]),
}


def viability_with(df: pd.DataFrame, cell_type: str, background: str) -> pd.DataFrame:
    """`viability`, with the background convention named explicitly."""
    sub = df[df["Cell Line"] == cell_type].copy()
    no_cell = BACKGROUND[background](df)
    no_drug = np.average(sub[sub["Drug"] == "DMSO"]["Measurement-Ref"]) - no_cell
    sub["Viability"] = (sub["Measurement-Ref"] - no_cell) / no_drug
    return sub


def fit_with(df: pd.DataFrame, cell_type: str, drug: str, background: str):
    """`fit`, with the background convention named explicitly."""
    sub = viability_with(df, cell_type, background)
    subset = sub[(sub["Drug"] == drug)
                 & (~sub["Cell Line"].astype(str).str.upper().eq("NO CELLS"))]
    x = np.asarray(subset["Concentration"]) * 1e-9
    y = 100 * np.asarray(subset["Viability"])
    params, _ = curve_fit(func2, x, y, P0, bounds=BOUNDS)
    resid = func2(x, *params) - y
    return (*params, 1.0 - (np.var(resid) / np.var(y)), subset)
