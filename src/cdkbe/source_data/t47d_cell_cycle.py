"""T47D cell-cycle panels: figure 4k, 4l and supplementary 13c, 13d.

Four panels, one experiment. T47D lines carrying a base-edited CDK variant
were treated for seven days with DMSO, palbociclib or serum starvation, then
stained and gated into G0, G1, S and G2/M. The FlowJo export of those four
gate frequencies is `data/raw/collaborator_panels/T47D_cell_cycle/`, 54 rows
-- six cell lines by three treatments by three replicates -- and everything
below is computed from it.

    4k, 4l    % S phase by treatment, with variant-vs-NT brackets
    S13c, S13d  the full G0/G1/S/G2M stack, DMSO and palbociclib only

The delivered package carries the notebook that drew them
(`cell_cycle_stacked_barplot.ipynb`) and its `figure_style.py`, both left in
place beside the data, and the published SVGs are vendored under
`tests/reference/t47d_cell_cycle/`. This module is that notebook's plotting
code moved into the package, so the panels come out of a test run rather than
out of somebody's working directory.

Two details are worth stating because they change the numbers:

* **The four gates are siblings, so they do not sum to 100%.** They are
  frequencies of the same parent population, and small gaps at the gate
  boundaries are normal. The stacked panels rescale each condition's four
  values to 100; the %S-phase panels use the gate frequency as measured.
* **The significance brackets are Welch's t-test, Bonferroni-corrected across
  the whole panel** -- every variant-by-treatment comparison drawn on one
  figure, so four tests per panel and `p_adj = 4p`, not two.

Several of the other T47D and MCF7 panels are included as their published
images rather than redrawn; `docs/COVERAGE.md` marks which. This experiment
came with its measurements, so it redraws.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

from .. import paths

FOLDER = paths.RAW / "collaborator_panels" / "T47D_cell_cycle"
DATA_CSV = FOLDER / "20260908_T47D_cell_cycle.csv"

PHASE_ORDER = ["G0", "G1", "S", "G2M"]
PHASE_LABELS = {"G0": "Sub-G1", "G1": "G1", "S": "S", "G2M": "G2/M"}
#: the published panels' fills, from the notebook's own override of
#: `figure_style.PHASE_COLORS`
PHASE_COLORS = {"G0": "#DE918F", "G1": "#B2CEE7", "S": "#668FCA",
                "G2M": "#C0A4C7"}
CELL_LINE_COLORS = {
    "ABE_NT": "#8C8C86", "CBE_NT": "#8C8C86",
    "CDK2_S261G": "#5B84B1", "CDK4_I51V": "#9A82B0",
    "CDK6_G53S": "#D9702E", "CDK6_R46W": "#7FB49B",
}

ABE_LINES = ["ABE_NT", "CDK2_S261G", "CDK4_I51V"]
CBE_LINES = ["CBE_NT", "CDK6_G53S", "CDK6_R46W"]
#: the labels the published panels print, which are not the sample names: each
#: CDK6 line carries a second edit that the sample name does not record
CELL_LINE_LABELS = {"CDK6_G53S": "CDK6 G53S_M54I",
                    "CDK6_R46W": "CDK6 R44C_R46W"}

TREATMENT_LABELS = {"Palbo": "Palbociclib", "ss": "serum free"}
#: serum starvation is measured but not drawn on these four panels
TREATMENT_ORDER = ["DMSO", "Palbo"]
SIG_TREATMENTS = ["DMSO", "Palbo"]

GROUP_WIDTH = 0.8
GENE_RE = re.compile(r"(CDK\d+)")
STACKED_FIGSIZE = (3.4, 2.6)
STACKED_LABEL_ROTATION = 45
SPHASE_FIG_HEIGHT = 2.1


# ----------------------------------------------------------------- reading --

def _parse_sample(name: str):
    parts = name.split("_")
    return "_".join(parts[:-2]), parts[-2], parts[-1]


def replicates() -> pd.DataFrame:
    """The FlowJo export, one row per well.

    -> Sample, CellLine, Treatment, Replicate, G0, G1, S, G2M
    """
    d = pd.read_csv(DATA_CSV)
    d.columns = [c.strip().rstrip(":").strip() for c in d.columns]
    d = d.dropna(axis=1, how="all")

    def col(needle):
        return next(c for c in d.columns if needle.lower() in c.lower())

    names = {col("Sample name"): "Sample", col("G1 |"): "G1", col("S |"): "S",
             col("G2M |"): "G2M", col("G0"): "G0"}
    d = d.drop(columns=[c for c in d.columns if c == "Sample"]).rename(columns=names)
    d["Sample"] = d["Sample"].astype(str).str.strip()
    d[PHASE_ORDER] = d[PHASE_ORDER].apply(pd.to_numeric, errors="coerce")
    parsed = d["Sample"].apply(_parse_sample)
    d["CellLine"] = parsed.apply(lambda p: p[0])
    d["Treatment"] = parsed.apply(lambda p: p[1])
    d["Replicate"] = parsed.apply(lambda p: p[2])
    return (d[["Sample", "CellLine", "Treatment", "Replicate"] + PHASE_ORDER]
            .dropna(subset=PHASE_ORDER).reset_index(drop=True))


def means() -> pd.DataFrame:
    """Mean gate frequency per cell line and treatment, over the replicates."""
    return (replicates().groupby(["CellLine", "Treatment"])[PHASE_ORDER]
            .mean().reset_index())


def _label(name: str, two_line: bool = False) -> str:
    label = CELL_LINE_LABELS.get(name, name.replace("_", " "))
    return label.replace(" ", "\n", 1) if two_line else label


# -------------------------------------------------------------- statistics --

def _stars(p: float) -> str:
    if p != p:
        return "n/a"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def significance(lines, control: str, treatments=None) -> pd.DataFrame:
    """Welch's t-test of %S phase, each variant against its own NT arm.

    Bonferroni over every comparison the panel draws -- `BONFERRONI_FAMILY =
    "figure"` in the delivered notebook -- so with two variants and two
    treatments the multiplier is four, not two.
    """
    d = replicates()
    treatments = list(treatments or SIG_TREATMENTS)
    rows = []
    for treatment in treatments:
        ref = d[(d["CellLine"] == control)
                & (d["Treatment"] == treatment)]["S"].to_numpy()
        for line in lines:
            if line == control:
                continue
            x = d[(d["CellLine"] == line)
                  & (d["Treatment"] == treatment)]["S"].to_numpy()
            p = (float("nan") if min(len(ref), len(x)) < 2 else
                 stats.ttest_ind(x, ref, equal_var=False,
                                 alternative="two-sided").pvalue)
            rows.append({"CellLine": line, "Treatment": treatment,
                         "vs_control": control, "n_control": len(ref),
                         "n_treatment": len(x), "p_value": p})
    out = pd.DataFrame(rows)
    m = int(out["p_value"].notna().sum())
    out["m"] = m
    out["p_adj"] = (out["p_value"] * m).clip(upper=1.0)
    out["sig"] = out["p_adj"].apply(_stars)
    return out


# ------------------------------------------------------------------ panels --

def _positions(groups, bars):
    width = GROUP_WIDTH / len(bars)
    return ({(g, b): gi + (bi - (len(bars) - 1) / 2) * width
             for gi, g in enumerate(groups) for bi, b in enumerate(bars)},
            width)


def stacked(lines, outdir=None, stem=None, treatments=None):
    """S13c/S13d: the whole stack, grouped by treatment.

    Each condition's four gate frequencies are rescaled to sum to 100 -- see
    this module's docstring for why they do not already.
    """
    treatments = list(treatments or TREATMENT_ORDER)
    grouped = means()
    x_pos, width = _positions(treatments, lines)

    fig, ax = plt.subplots(figsize=STACKED_FIGSIZE)
    labels = []
    for (treatment, line), x in x_pos.items():
        row = grouped[(grouped["CellLine"] == line)
                      & (grouped["Treatment"] == treatment)]
        if row.empty:
            continue
        vals = row.iloc[0][PHASE_ORDER].astype(float)
        vals = vals / vals.sum() * 100
        bottom = 0.0
        for phase in PHASE_ORDER:
            ax.bar(x, float(vals[phase]), bottom=bottom, width=width,
                   color=PHASE_COLORS[phase], edgecolor="black", linewidth=0.3)
            bottom += float(vals[phase])
        labels.append(ax.annotate(
            _label(line), xy=(x + width / 3, 0), xycoords="data",
            xytext=(0, -4), textcoords="offset points", ha="right", va="top",
            rotation=STACKED_LABEL_ROTATION, rotation_mode="anchor",
            fontsize=10))

    ax.set_xticks(range(len(treatments)))
    ax.set_xticklabels([TREATMENT_LABELS.get(t, t) for t in treatments],
                       fontsize=10)
    # put the treatment labels just under the longest rotated cell-line label,
    # measured rather than guessed, so the gap holds at any figure height
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest = min(a.get_window_extent(renderer).y0 for a in labels)
    axis = ax.get_window_extent(renderer).y0
    ax.tick_params(axis="x", pad=(axis - lowest) * 72 / fig.dpi + 4, length=0)
    ax.tick_params(axis="y", length=4, width=0.3)
    ax.set_xlim(-0.6, len(treatments) - 1 + 0.6)
    ax.set_ylabel("% of cells", fontsize=10)
    ax.set_ylim(0, 110)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(0.3)
    ax.spines["bottom"].set_linewidth(0.3)

    order = list(reversed(PHASE_ORDER))
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=PHASE_COLORS[p],
                             edgecolor="black", linewidth=0.3) for p in order]
    ax.legend(handles, [PHASE_LABELS[p] for p in order], loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=10)
    _save(fig, outdir, stem)
    return fig, grouped[grouped["Treatment"].isin(treatments)
                        & grouped["CellLine"].isin(lines)]


def s_phase(lines, control: str, outdir=None, stem=None, treatments=None):
    """4k/4l: %S phase by treatment, bars with their replicates and brackets."""
    treatments = list(treatments or TREATMENT_ORDER)
    d = replicates()
    x_pos, width = _positions(treatments, lines)
    summary = (d.groupby(["CellLine", "Treatment"])["S"]
               .agg(["mean", "std"]).reset_index())

    fig, ax = plt.subplots(figsize=(max(2.8, 1.1 * len(treatments)
                                        * (len(lines) / 3)), SPHASE_FIG_HEIGHT))
    for (treatment, line), x in x_pos.items():
        row = summary[(summary["CellLine"] == line)
                      & (summary["Treatment"] == treatment)]
        if row.empty:
            continue
        mean, sd = float(row["mean"].iloc[0]), float(row["std"].iloc[0])
        ax.bar(x, mean, width=width, facecolor="white", edgecolor="black",
               linewidth=0.3, zorder=2)
        ax.errorbar(x, mean, yerr=sd, ecolor="black", elinewidth=0.4,
                    capsize=3, capthick=0.6, fmt="none", zorder=2)
        vals = d[(d["CellLine"] == line)
                 & (d["Treatment"] == treatment)]["S"].to_numpy()
        for j, v in enumerate(vals):
            ax.scatter(x + (j - (len(vals) - 1) / 2) * width * 0.28, v, s=30,
                       color=CELL_LINE_COLORS.get(line, "#888888"),
                       edgecolor="black", linewidth=0.35, zorder=3)

    sig = significance(lines, control,
                       [t for t in SIG_TREATMENTS if t in treatments])
    tops = _brackets(ax, d, x_pos, lines, control, sig,
                     [t for t in SIG_TREATMENTS if t in treatments])

    ax.set_xticks(range(len(treatments)))
    ax.set_xticklabels([TREATMENT_LABELS.get(t, t) for t in treatments],
                       fontsize=8)
    ax.tick_params(axis="x", pad=8, length=0)
    ax.tick_params(axis="y", labelsize=7, width=0.3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.set_xlim(-0.6, len(treatments) - 1 + 0.6)
    ax.set_ylabel("% S phase cells", fontsize=9)
    data_max = float(d[d["CellLine"].isin(lines)
                       & d["Treatment"].isin(treatments)]["S"].max())
    ax.set_ylim(0, max((max(tops) + 6) if tops else 0.0, data_max + 4) * 1.1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(0.3)
    ax.spines["bottom"].set_linewidth(0.3)

    handles = [plt.Line2D([0], [0], marker="o", linestyle="",
                          markerfacecolor=CELL_LINE_COLORS.get(cl, "#888888"),
                          markeredgecolor="black", markersize=5,
                          markeredgewidth=0.35) for cl in lines]
    ax.legend(handles, [_label(cl) for cl in lines], loc="upper right",
              frameon=False, fontsize=7)
    _save(fig, outdir, stem)
    return fig, sig


def _brackets(ax, d, x_pos, lines, control, sig, treatments):
    """Significance brackets, narrowest span lowest, stacked per treatment."""
    step, gap, height, tops = 3, 3, 1, []
    for treatment in treatments:
        spans = sorted(
            ((line, x_pos[(treatment, control)], x_pos[(treatment, line)])
             for line in lines if line != control),
            key=lambda c: abs(c[2] - c[1]))
        previous = None
        for line, x1, x2 in spans:
            lo, hi = sorted([x1, x2])
            covered = [d[(d["CellLine"] == cl) & (d["Treatment"] == tr)]["S"].max()
                       for (tr, cl), x in x_pos.items() if lo - 1e-9 <= x <= hi + 1e-9]
            floor = (max(covered) if covered else 0.0) + gap
            y = floor if previous is None else max(floor, previous + step)
            previous = y
            row = sig[(sig["CellLine"] == line)
                      & (sig["Treatment"] == treatment)].iloc[0]
            ax.plot([lo, lo, hi, hi], [y, y + height, y + height, y],
                    color="black", linewidth=0.3)
            ax.text((lo + hi) / 2, y + height, row["sig"], ha="center",
                    va="bottom", fontsize=7)
            tops.append(y + height)
    return tops


# --------------------------------------------------------------- registry ---

_STEM = {
    "4k": "T47D_cell_cycle_stacked_Sphase_ABE",
    "4l": "T47D_cell_cycle_stacked_Sphase_CBE",
    "S13c": "T47D_cell_cycle_stacked_DMSO_Palbo_ABE",
    "S13d": "T47D_cell_cycle_stacked_DMSO_Palbo_CBE",
}
#: panel -> (the cell-line group, its non-targeting control)
GROUPS = {"4k": (ABE_LINES, "ABE_NT"), "4l": (CBE_LINES, "CBE_NT"),
          "S13c": (ABE_LINES, "ABE_NT"), "S13d": (CBE_LINES, "CBE_NT")}


def _save(fig, outdir, stem):
    if outdir is not None and stem:
        fig.savefig(str(Path(outdir) / f"{stem}.pdf"), bbox_inches="tight",
                    transparent=True)


def panel(name: str, outdir=None, **_):
    lines, control = GROUPS[name]
    if name in ("4k", "4l"):
        return s_phase(lines, control, outdir=outdir, stem=_STEM[name])
    return stacked(lines, outdir=outdir, stem=_STEM[name])


def source(name: str) -> dict:
    """-> the frames a reader needs to redraw this panel, keyed by sheet name.

    A significance table alone does not reproduce a panel: the bars come from
    the gate frequencies, and the brackets from the test. Both are returned,
    plus the published cell-line label beside the export's own spelling, so
    `CDK6_G53S` and the `CDK6 G53S_M54I` printed on the panel sit in one row.
    """
    lines, control = GROUPS[name]
    keep = lambda d: d[d["CellLine"].isin(lines)
                       & d["Treatment"].isin(TREATMENT_ORDER)].reset_index(drop=True)

    def labelled(d):
        d = d.copy()
        d.insert(1, "Label", [_label(c) for c in d["CellLine"]])
        return d

    out = {"replicates": labelled(keep(replicates())),
           "means": labelled(keep(means()))}
    if name in ("4k", "4l"):
        # the %S-phase panels plot the gate frequency as measured
        out["significance"] = labelled(significance(lines, control))
    else:
        # the stacked panels rescale each condition's four gates to 100
        g = out["means"].copy()
        total = g[PHASE_ORDER].sum(axis=1)
        for phase in PHASE_ORDER:
            g[phase] = g[phase] / total * 100
        out["means_rescaled_to_100"] = g
    return out


def _wrap(name: str):
    def drawer(outdir=None, **kw):
        return panel(name, outdir=outdir, **kw)
    drawer.__name__ = drawer.__qualname__ = f"panel_{name}"
    drawer.__doc__ = f"T47D cell-cycle panel {name}."
    return drawer


MAIN_PANELS = {k: (lambda _p=k: source(_p)) for k in ("4k", "4l")}
MAIN_FIGURES = {k: _wrap(k) for k in ("4k", "4l")}
PANELS = {k: (lambda _p=k: source(_p)) for k in ("S13c", "S13d")}
FIGURES = {k: _wrap(k) for k in ("S13c", "S13d")}

#: `scripts/render_figures.py` reads this to name the file a panel produces
PARTS = {k: [v] for k, v in _STEM.items()}
