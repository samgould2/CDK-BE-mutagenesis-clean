"""B-ALL transplantation leukemia burden (supplementary 7f-7i).

From `11_in_vivo_analysis/invivo.ipynb` cells 2-10. Two
transplantation experiments, each read straight from its flow export:

    7f, 7g  CDK9 L156F   vs KI-CDK9d-32
    7h, 7i  CDK7 L18F    vs SY-5609

The first panel of each pair is the raw burden by genotype and treatment; the
second normalises each treated mouse to the **mean of its own genotype's
vehicle arm**, so the two genotypes are compared on the same 100% baseline
rather than against each other's absolute burden.

The notebook ran `ttest_ind` on the normalised values but printed the result
rather than drawing it, so the p-values are returned alongside the frame
instead of being written onto the panel.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns

from .. import paths

FLOW = paths.RAW / "flow"

#: experiment -> (file, sheet, mutant genotype, treatment arm, colour)
EXPERIMENTS = {
    "CDK9": ("Cdk9 L156F  B-ALL Transplantation Experiment 19 .xlsx.xls",
             "Sheet2", "CDK9 L156F", "KI-CDK9d-32", "tab:purple"),
    "CDK7": ("B-ALL mCdk7 L18F SY-5609 Competition D10.xls",
             "Sheet1", "CDK7 L18F", "SY-5609", "cyan"),
}

BURDEN = "Leukemia Burden"
NORMED = "% of Vehicle Leukemia Burden"


def raw(experiment: str) -> pd.DataFrame:
    """Per-mouse leukemia burden, as exported."""
    fn, sheet, *_ = EXPERIMENTS[experiment]
    return pd.read_excel(FLOW / fn, sheet)


def normalised(experiment: str) -> pd.DataFrame:
    """Treated mice as a percentage of their own genotype's vehicle mean."""
    d = raw(experiment)
    _, _, mutant, treatment, _ = EXPERIMENTS[experiment]
    out = []
    for geno in ("WT", mutant):
        vehicle = d[(d["Genotype"] == geno) & (d["Treatment"] == "Vehicle")][BURDEN]
        treated = d[(d["Genotype"] == geno) & (d["Treatment"] == treatment)][BURDEN]
        out.append(pd.DataFrame({"Genotype": geno,
                                 NORMED: 100 * np.asarray(treated) / vehicle.mean()}))
    return pd.concat(out, ignore_index=True)


def ttest(experiment: str):
    """WT vs mutant on the normalised burden -> (t, p); the notebook printed this."""
    d = normalised(experiment)
    _, _, mutant, _, _ = EXPERIMENTS[experiment]
    return scipy.stats.ttest_ind(d[d["Genotype"] == "WT"][NORMED],
                                 d[d["Genotype"] == mutant][NORMED])


def _style(ax, ylabel):
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.set_xlabel("Genotype", fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)


def panel_raw(experiment: str, outdir=None):
    """-> (figure, frame): burden by genotype and treatment."""
    d = raw(experiment)
    _, _, _, treatment, colour = EXPERIMENTS[experiment]
    fig, ax = plt.subplots(figsize=(5, 5))
    kw = dict(data=d, x="Genotype", y=BURDEN, hue="Treatment",
              edgecolor="black", linewidth=1, palette=["grey", colour])
    sns.barplot(ax=ax, errorbar="sd", **kw)
    sns.stripplot(ax=ax, dodge=True, size=13, legend=False, **kw)
    ax.legend(loc="upper right", fontsize=15, frameon=False)
    _style(ax, "% Leukemia Burden")
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / f"invivo_{experiment}.pdf"), transparent=True)
    return fig, d


def panel_normed(experiment: str, outdir=None):
    """-> (figure, frame): treated mice as a percentage of their vehicle arm."""
    d = normalised(experiment)
    _, _, _, treatment, colour = EXPERIMENTS[experiment]
    fig, ax = plt.subplots(figsize=(3, 5))
    kw = dict(data=d, x="Genotype", y=NORMED, edgecolor="black",
              linewidth=1, palette=[colour])
    sns.barplot(ax=ax, errorbar="sd", **kw)
    sns.stripplot(ax=ax, dodge=False, size=13, legend=False, **kw)
    _style(ax, NORMED)
    ax.set_title(f"{treatment} Treated Mice\n\n", fontsize=14.5)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / f"invivo_{experiment}_normed.pdf"), transparent=True)
    return fig, d


def _wrap(panel, draw_fn, experiment, what):
    _, _, variant, drug, _ = EXPERIMENTS[experiment]
    f = lambda outdir=None, **kw: draw_fn(experiment, outdir)     # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    f.__doc__ = f"{what} for the {variant} B-ALL transplant, vehicle vs {drug}."
    return f


FIGURES = {
    "S7f": _wrap("S7f", panel_raw, "CDK9", "Leukemia burden per mouse"),
    "S7g": _wrap("S7g", panel_normed, "CDK9", "Leukemia burden as % of vehicle"),
    "S7h": _wrap("S7h", panel_raw, "CDK7", "Leukemia burden per mouse"),
    "S7i": _wrap("S7i", panel_normed, "CDK7", "Leukemia burden as % of vehicle"),
}

PANELS = {
    "S7f": lambda: raw("CDK9"),
    "S7g": lambda: normalised("CDK9"),
    "S7h": lambda: raw("CDK7"),
    "S7i": lambda: normalised("CDK7"),
}
