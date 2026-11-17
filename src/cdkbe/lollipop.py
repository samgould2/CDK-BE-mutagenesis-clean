"""Residue-axis lollipop plots: median LFC per guide against codon position.

From lollipop_plots_canonical.ipynb, which
held three ~150-line plotters (`subpool1_plotter`, `subpool2_plotter`,
`subpool3_plotter`) that were identical apart from the gene list, the figure
size, the panel width ratios and the x ticks. Those four things are `SUBPOOLS`
here and there is one plotter.

Each panel is one gene. Marker area scales with -log10(FDR); guides that clear
both the FDR and LFC cutoffs are drawn opaque and labelled, the rest faint.
Shaded bands are the KLIFS pocket positions. The last panel holds the control
guide classes.
"""
from __future__ import annotations

import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from . import io, lfc, paths

try:                                     # optional: labels still draw without it
    from adjustText import adjust_text
except ImportError:                      # pragma: no cover
    adjust_text = None

EDITOR_COLOR = {"CBE": "tab:blue", "ABE": "tab:purple"}

#: gene order, figure size, panel width ratios (last = controls) and x ticks
SUBPOOLS = {
    "subpool1": dict(
        genes=["CDK7", "CDK8", "CDK19", "CDK9"], figsize=(15, 5),
        width_ratios=[352, 472, 511, 379, 300],
        xticks=[[0, 100, 200, 346], [0, 100, 200, 300, 400, 464],
                [0, 100, 200, 300, 400, 502], [0, 100, 200, 300, 372]]),
    "subpool2": dict(
        genes=["CDK12", "CDK13"], figsize=(15, 6),
        width_ratios=[1490, 1512, 500],
        xticks=[[0, 250, 500, 750, 1000, 1250, 1490],
                [0, 250, 500, 750, 1000, 1250, 1512]]),
    "subpool3": dict(
        genes=["CDK2", "CDK4", "CDK6"], figsize=(15, 6),
        width_ratios=[298, 303, 326, 150],
        xticks=[[0, 50, 100, 150, 200, 250, 298], [0, 50, 100, 150, 200, 250, 303],
                [0, 50, 100, 150, 200, 250, 326]]),
}


def find_ranges(numbers):
    """Consecutive runs in a list of positions -> [(start, end), ...]."""
    s = sorted(set(int(n) for n in numbers))
    if not s:
        return []
    out, start = [], s[0]
    for prev, cur in zip(s, s[1:]):
        if cur != prev + 1:
            out.append((start, prev))
            start = cur
    out.append((start, s[-1]))
    return out


def klifs_ranges(gene: str) -> list:
    k = io.read(paths.REFERENCE / "KLIFS" / "KLIFS_compiled.csv")
    return find_ranges(k[gene].dropna())


def control_classes(m: pd.DataFrame, lfc_col: str, wt_cutoff: float, editor: str) -> pd.DataFrame:
    """The three control guide classes drawn in the last panel.

    The order below is the x order of that panel -- the original got it from
    the order these frames happened to be concatenated in.
    """
    classes = (
        ("Intron", m[m["classification"] == "intron"]),
        ("NT", m[m["Gene"] == "NT"]),
        (f"WT\u2265\n{wt_cutoff}%", m[m["WT_perc"] >= wt_cutoff]),
    )
    parts = []
    for label, sub in classes:
        sub = sub.copy()
        sub["Guide Class"] = label
        parts.append(sub)
    out = pd.concat(parts)
    out["Editor"] = editor
    return out.rename(columns={lfc_col: "LFC"})


def lollipop(subpool: str, targ: dict, merged: dict, conditions: dict, title: str,
             lfc_cutoff: float, fdr_cutoff: float = lfc.FDR_CUTOFF,
             min_base_counts: int = lfc.MIN_BASE_COUNTS,
             wt_cutoff_main: float = lfc.WT_CUTOFF_MAIN, wt_cutoff_ctrl: float = 98,
             editors=("CBE", "ABE"), baseline: str = "DMSO", n_library: int | None = None):
    """Draw one lollipop figure.

    `targ`/`merged` are {editor: frame} as returned by `cdkbe.lfc.load_merged`,
    `conditions` is {editor: condition name} -- the two arms can be named
    differently in their respective tables.
    """
    cfg = SUBPOOLS[subpool]
    genes = cfg["genes"]
    fig, ax = plt.subplots(1, len(genes) + 1, figsize=cfg["figsize"],
                           sharey=True, width_ratios=cfg["width_ratios"])

    kept = {}
    for editor in editors:
        cond = conditions[editor]
        name, fdr = f"LFC_median_{cond}", f"FDR_{cond}"
        sub = lfc.guide_filter(targ[editor], genes, min_base_counts, wt_cutoff_main)
        kept[editor] = sub
        for i, gene in enumerate(genes):
            a = ax[i]
            g = sub[sub["Gene"] == gene]
            faint = g[(g[fdr] > fdr_cutoff) | (g[name] < lfc_cutoff)]
            hit = g[(g[fdr] <= fdr_cutoff) & (g[name] >= lfc_cutoff)]
            for part, alpha in ((faint, 0.2), (hit, 1.0)):
                a.scatter(part["Codon"], part[name],
                          s=20 * -np.log10(part[fdr]) + 1, alpha=alpha,
                          color=EDITOR_COLOR[editor], rasterized=True)
                # the stem of each lollipop
                a.vlines(part["Codon"], 0, part[name], linewidth=0.4,
                         color="black", alpha=0.2, rasterized=True)
            if len(hit):
                texts = [a.text(r["Codon"], r[name], r["HGVSp"], fontsize=14,
                                color=EDITOR_COLOR[editor],
                                path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")])
                         for _, r in hit.iterrows()]
                if adjust_text is not None:
                    adjust_text(texts, only_move={"points": "xy", "texts": "xy"}, ax=a,
                                arrowprops=dict(arrowstyle="-", linewidth=1, alpha=1,
                                                color=EDITOR_COLOR[editor]))

    for i, gene in enumerate(genes):
        a = ax[i]
        a.axhline(0, color="black", linestyle="dashed")
        a.set_title(gene, fontsize=15)
        a.set_xticks(cfg["xticks"][i])
        a.spines[["right", "top"]].set_visible(False)
        a.tick_params(axis="both", which="major", labelsize=13)
        for lo, hi in klifs_ranges(gene):
            a.axvspan(lo, hi, alpha=0.2)
    ax[0].set_ylabel(f"Median LFC (vs. {baseline})", fontsize=15)

    # ---- controls
    ctrl = pd.concat([
        control_classes(merged[e][merged[e]["base_RAW"] >= min_base_counts],
                        f"LFC_median_{conditions[e]}", wt_cutoff_ctrl, e)
        for e in editors])
    c = ax[-1]
    order = ["Intron", "NT", f"WT≥\n{wt_cutoff_ctrl}%"]
    pal = [EDITOR_COLOR[e] for e in editors]
    sns.stripplot(data=ctrl, x="Guide Class", y="LFC", ax=c, hue="Editor", order=order,
                  palette=pal, hue_order=list(editors), alpha=0.2, s=6, legend=False,
                  dodge=len(editors) > 1, rasterized=True)
    sns.boxplot(data=ctrl, x="Guide Class", y="LFC", ax=c, hue="Editor", order=order,
                palette=pal, hue_order=list(editors))
    c.axhline(0, color="black", linestyle="dashed")
    c.spines[["right", "top"]].set_visible(False)
    c.tick_params(axis="both", which="major", labelsize=13)
    c.set_xlabel("")
    c.set_ylabel("")
    c.set_title("Controls", fontsize=15)
    c.legend(fontsize=10)

    counts = ", ".join(f"{e}: {len(kept[e]):,}" + (f"/{n_library:,}" if n_library else "")
                       + " gRNAs" for e in editors)
    fig.suptitle(
        f"{title} | hits with FDR < {fdr_cutoff} & LFC ≥ {lfc_cutoff} & "
        f"base mean ≥ {min_base_counts} & ≥ {100 - wt_cutoff_main}% sensor editing\n{counts}",
        fontsize=14)
    fig.tight_layout()
    return fig, kept
