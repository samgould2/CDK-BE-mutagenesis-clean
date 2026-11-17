"""BEAT-seq QC and differential-expression panels (supplementary 14a-14k).

These eleven panels were cut from five composite PDFs in the delivered bundle's
`plotting/figures/` tree. **No plotting code was delivered with that tree** --
the bundle shipped outputs only -- so each panel here is drawn from the tables
the bundle saved beside its PDF, and every number it puts on the page is
checked against the number printed on the published one.

Which letter is which sub-panel **is** established, contrary to what this module
said before: the composites carry their own sub-panel letters, printed in the
titles. `experiment_summary.pdf` is lettered a-d and those four are
supplementary 14a-14d in the same order; `de_summary.pdf` is also lettered a-d,
and those are supplementary 14g-14j. So:

    14a  Cells recovered            14g  DE calls
    14b  Guide assignment           14h  Not power-limited
    14c  Library coverage           14i  Same guide, both lines
    14d  Representation depth       14j  Active guide-arms
    14e  Cells per guide, per arm
    14f  Designed consequence       14k  Pooled correlation clustermap

What pins the port is that the delivered tables regenerate the annotations the
published panels print, to the digit:

* 93 / 95 / 94 / 93 % passing QC and 138,363 post-QC cells (14a);
* 73 / 70 / 86 / 88 % singlets and 101,805 singlets (14b);
* 91/92, 71/72, 91/92, 71/72 guides recovered (14c);
* medians 509 / 537 / 138 / 71 and means 500 / 549 / 136 / 70 (14d, 14e);
* 4/232 missense, 2/34 nonsense, 2/22 start_loss (14f);
* 2/82, 0/62, 3/82, 3/62 active guide-arms (14g);
* Spearman rho +0.06, -0.03, -0.05, +0.05 of cells against DE genes (14h);
* Spearman rho 0.26 between the two cell lines on the same guide (14i);
* the eight active guide-arms with their cell counts and edit rates (14j).

14e's histogram bins are the one thing not printed anywhere, and they are
recovered rather than chosen: 29 logarithmic bins across the pooled range
reproduce all four arms' bar heights exactly (`BINS`), which no neighbouring
bin count does.

**14k is not redrawn.** The pseudobulk matrix is here and so are its annotation
tracks, but the published clustermap is over "2,607 shared genes" while the
matrix holds 4,779 with no missing values, so the filter that chose those 2,607
is in nothing the bundle delivered and cannot be inferred. Guessing a gene set
would change every correlation on the page. See docs/DISCREPANCIES.md item 13.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import io, paths

FIGS = paths.PROCESSED / "beatseq" / "outputs" / "plotting" / "figures"
TABLES = (paths.PROCESSED / "beatseq" / "outputs" / "05_primary_analysis"
          / "outputs" / "tables")
PANEL_PDFS = paths.REFERENCE / "beatseq_panels"

#: the four arms, in the order every published panel puts them
ARMS = ["A549/ABE", "A549/CBE", "T47D/ABE", "T47D/CBE"]

#: arm colours, read off the pooled-correlation legend swatches
ARM_COLOUR = {"A549/ABE": "#17706B", "A549/CBE": "#77C4BA",
              "T47D/ABE": "#5B3A8C", "T47D/CBE": "#A98FD0"}

#: designed-consequence colours, read off the same legend
CLASS_COLOUR = {"missense": "#4C78A8", "nonsense": "#B2182B",
                "start_loss": "#8C564B", "control": "#4C9F70"}

BLUE = "#2980B9"        # the QC panels' one accent colour
GREY_MID = "#90A4AE"    # multiplets, and guide-arms with 1-9 DE genes
GREY_PALE = "#CFD8DC"   # unassigned cells, and guide-arms with none
GREY_DARK = "#37474F"   # the dashed "designed" outline, and its label
RED = "#C0392B"         # the DE >= 10 activity threshold
ANNOT = "0.35"          # the grey of the per-bar annotations

#: the DE-gene count at which a guide-arm is called active
ACTIVE_AT = 10

#: 14d/14e's axis: the decades the published panels tick
DEPTH_TICKS = [30, 100, 300, 1000]

#: 14e's bins. 29 logarithmic bins across the pooled range of `cells_per_guide`
#: reproduce every bar height in the published panel; 28 and 30 do not.
BINS = 29

#: supplementary panel -> (composite PDF, {sheet name: csv beside it})
COMPOSITES = {
    "experiment_summary": ("qc/experiment_summary.pdf", {
        "guide_assignment": "qc/guide_assignment_summary_stats.csv",
        "library_coverage": "qc/library_coverage.csv"}),
    "cells_per_guide": ("qc/cells_per_guide_per_arm.pdf", {
        "cells_per_guide": "qc/cells_per_guide.csv"}),
    "de_by_variant_class": ("de/de_by_variant_class.pdf", {
        "de_hits": "de/de_hits.csv"}),
    "de_summary": ("de/de_summary.pdf", {
        "de_summary_stats": "de/de_summary_stats.csv",
        "de_target_enrichment": "de/de_target_enrichment.csv"}),
    "pooled_corr": ("clustermaps/pooled_corr.pdf", {
        "clustermap_null": "clustermaps/clustermap_null.csv"}),
}

#: which composite each supplementary letter is a sub-panel of, and which
#: sub-panel letter it carries inside that composite
PANEL_COMPOSITE = {
    "S14a": ("experiment_summary", "a"), "S14b": ("experiment_summary", "b"),
    "S14c": ("experiment_summary", "c"), "S14d": ("experiment_summary", "d"),
    "S14e": ("cells_per_guide", None),
    "S14f": ("de_by_variant_class", None),
    "S14g": ("de_summary", "a"), "S14h": ("de_summary", "b"),
    "S14i": ("de_summary", "c"), "S14j": ("de_summary", "d"),
    "S14k": ("pooled_corr", None),
}


# ------------------------------------------------------------- the tables

def guide_assignment() -> pd.DataFrame:
    """Cells and guide calls per arm, in the published arm order."""
    d = io.read(FIGS / "qc" / "guide_assignment_summary_stats.csv")
    arms = d[d["level"] == "arm"].set_index("group")
    return arms.reindex(ARMS).reset_index().rename(columns={"group": "arm"})


def overall() -> pd.Series:
    """The experiment-wide row the first two titles quote."""
    d = io.read(FIGS / "qc" / "guide_assignment_summary_stats.csv")
    return d[d["level"] == "overall"].iloc[0]


def library_coverage() -> pd.DataFrame:
    """Guides recovered against guides designed, per arm.

    The delivered table splits each arm into targeting guides and non-targeting
    controls; the published panel plots their sum, which is where 91/92 comes
    from -- 82 of 82 targeting plus 9 of 10 controls.
    """
    d = io.read(FIGS / "qc" / "library_coverage.csv")
    g = d.groupby("arm")[["designed", "recovered"]].sum()
    return g.reindex(ARMS).reset_index()


def cells_per_guide() -> pd.DataFrame:
    """One row per guide per arm: how many cells carried it."""
    return io.read(FIGS / "qc" / "cells_per_guide.csv")


def depth_stats() -> pd.DataFrame:
    """The median and mean 14e prints, rounded the way the panel prints them."""
    c = cells_per_guide()
    g = c.groupby("arm")["n_cells"]
    out = pd.DataFrame({"n_guides": g.size(),
                        "median": g.median().map(round),
                        "mean": g.mean().map(round),
                        "max": g.max()})
    return out.reindex(ARMS)


def guide_annotation() -> pd.DataFrame:
    """Every guide-arm with its DE call, variant class, cells and edit rate.

    324 rows: the 288 targeting guide-arms plus 36 non-targeting controls.
    `guide` is the protospacer alone, which is what pairs a guide across the
    two cell lines in 14i.
    """
    d = io.read(TABLES / "de" / "pseudobulk_guide_annotation.csv").copy()
    d["guide"] = d["arm_guide"].str.split("|", n=1).str[1]
    return d


def targeting() -> pd.DataFrame:
    """The 288 targeting guide-arms: controls are not DE-called."""
    d = guide_annotation()
    return d[~d["is_control"]].reset_index(drop=True)


def de_summary_stats() -> pd.DataFrame:
    """Active / subthreshold / inactive counts per arm, plus 14h's rho."""
    d = io.read(FIGS / "de" / "de_summary_stats.csv")
    arms = d[d["level"] == "arm"].set_index("group")
    return arms.reindex(ARMS).reset_index().rename(columns={"group": "arm"})


def de_hits() -> pd.DataFrame:
    """The eight active guide-arms, longest bar first."""
    d = io.read(FIGS / "de" / "de_hits.csv")
    return d.sort_values("num_de_genes", ascending=False).reset_index(drop=True)


def composite(name: str) -> dict:
    """-> {sheet: frame} for one composite figure, as the bundle saved it."""
    _, sheets = COMPOSITES[name]
    return {sheet: io.read(FIGS / rel) for sheet, rel in sheets.items()}


def published_pdf(key: str):
    """The delivered composite PDF this panel was cut from."""
    return PANEL_PDFS / COMPOSITES[PANEL_COMPOSITE[key][0]][0]


# ------------------------------------------------------------- the drawing

def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _title(ax, letter, text):
    """The composites letter their sub-panels; a standalone panel does not.

    In the assembled supplementary figure each of these carries its own
    manuscript letter, so the composite-local one would be a second, wrong
    label sitting beside it.
    """
    ax.set_title(f"{letter}  {text}" if letter else text, loc="left")


def _arm_labels(ax):
    ax.set_xticks(range(len(ARMS)))
    ax.set_xticklabels(ARMS, rotation=30, ha="right")


def _depth_axis(ax, which="y"):
    """The shared log axis of 14d, 14e and 14h, ticked at the published decades."""
    getattr(ax, f"set_{which}scale")("log")
    a = ax.yaxis if which == "y" else ax.xaxis
    a.set_ticks(DEPTH_TICKS)
    a.set_ticklabels([f"{t:,}" for t in DEPTH_TICKS])
    a.set_minor_locator(plt.NullLocator())


# ------------------------------------------------------------------- 14a

def source_14a() -> pd.DataFrame:
    """Cells before and after QC, per arm."""
    d = guide_assignment()
    return d[["arm", "total_pre_qc", "total_post_qc", "pct_passing_qc"]]


def _draw_14a(ax, letter=None):
    d = source_14a()
    x = np.arange(len(d))
    w = 0.40                     # the published pair touches, with no gap
    ax.bar(x - w / 2, d["total_pre_qc"], w, color="#E6E6E6", label="Pre-QC")
    ax.bar(x + w / 2, d["total_post_qc"], w, color=BLUE, label="Post-QC")
    for i, row in d.iterrows():
        ax.text(i + w / 2, row["total_post_qc"], f"{row['pct_passing_qc']:.0f}%",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Cells")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:,.0f}")
    ax.margins(y=0.20)           # the headroom that puts a tick at 80,000
    ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    _arm_labels(ax)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    _title(ax, letter, f"Cells recovered ({overall()['total_post_qc']:,} post-QC)")
    _style(ax)


def draw_14a(outdir=None, **kw):
    """Cells recovered per arm, before and after quality control."""
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    _draw_14a(ax)
    fig.tight_layout()
    return fig, source_14a()


# ------------------------------------------------------------------- 14b

def source_14b() -> pd.DataFrame:
    """Singlet / multiplet / unassigned shares of post-QC cells, per arm."""
    d = guide_assignment()
    return d[["arm", "singlet_pct", "multiplet_pct", "unassigned_pct"]]


def _draw_14b(ax, letter=None):
    d = source_14b()
    x = np.arange(len(d))
    bottom = np.zeros(len(d))
    for col, colour, label in [("singlet_pct", BLUE, "singlet"),
                               ("multiplet_pct", GREY_MID, "multiplet"),
                               ("unassigned_pct", GREY_PALE, "unassigned")]:
        ax.bar(x, d[col], 0.65, bottom=bottom, color=colour, label=label)
        bottom += d[col].values
    for i, v in enumerate(d["singlet_pct"]):
        ax.text(i, 103, f"{v:.0f}%", ha="center", va="bottom", fontsize=9,
                color=BLUE, fontweight="bold")
    ax.set_ylabel("% of post-QC cells")
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylim(0, 140)          # room for the labels, and the legend above them
    _arm_labels(ax)
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left",
              handlelength=1.2, columnspacing=1.0)
    _title(ax, letter, f"Guide assignment ({overall()['singlet']:,} singlets)")
    _style(ax)


def draw_14b(outdir=None, **kw):
    """How post-QC cells split into singlets, multiplets and unassigned."""
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    _draw_14b(ax)
    fig.tight_layout()
    return fig, source_14b()


# ------------------------------------------------------------------- 14c

def source_14c() -> pd.DataFrame:
    """Guides recovered against guides designed, per arm."""
    return library_coverage()


def _draw_14c(ax, letter=None):
    d = source_14c()
    x = np.arange(len(d))
    ax.bar(x, d["designed"], 0.65, facecolor="none", edgecolor=GREY_DARK,
           linestyle="--", linewidth=1.2)
    ax.bar(x, d["recovered"], 0.65, color=BLUE)
    for i, row in d.iterrows():
        ax.text(i, row["designed"],
                f"{row['recovered']:.0f}/{row['designed']:.0f}",
                ha="center", va="bottom", fontsize=9, color=GREY_DARK,
                fontweight="bold")
    ax.set_ylabel("Guides")
    _arm_labels(ax)
    _title(ax, letter, "Library coverage (recovered / designed)")
    _style(ax)


def draw_14c(outdir=None, **kw):
    """Library coverage: guides seen against guides in the design."""
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    _draw_14c(ax)
    fig.tight_layout()
    return fig, source_14c()


# ------------------------------------------------------------------- 14d

def source_14d() -> pd.DataFrame:
    """Cells per guide, one row per guide per arm."""
    return cells_per_guide()


def _draw_14d(ax, letter=None, seed=0):
    d = source_14d()
    rng = np.random.default_rng(seed)
    groups = [d.loc[d["arm"] == a, "n_cells"].values for a in ARMS]
    bp = ax.boxplot(groups, positions=range(len(ARMS)), widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color=GREY_DARK, linewidth=1.2),
                    whiskerprops=dict(color=GREY_DARK),
                    capprops=dict(color=GREY_DARK),
                    boxprops=dict(edgecolor=GREY_DARK))
    for patch, arm in zip(bp["boxes"], ARMS):
        patch.set_facecolor(ARM_COLOUR[arm])
    for i, vals in enumerate(groups):
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(vals)), vals, s=4,
                   color="0.4", alpha=0.6, linewidths=0, zorder=3)
    _depth_axis(ax, "y")
    ax.set_ylabel("Cells per guide")
    _arm_labels(ax)
    _title(ax, letter, "Representation depth")
    _style(ax)


def draw_14d(outdir=None, **kw):
    """The spread of cells per guide within each arm."""
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    _draw_14d(ax)
    fig.tight_layout()
    return fig, depth_stats().reset_index()


# ------------------------------------------------------------------- 14e

def bin_edges() -> np.ndarray:
    """The 29 logarithmic bins shared by all four of 14e's histograms."""
    v = cells_per_guide()["n_cells"]
    return np.logspace(np.log10(v.min()), np.log10(v.max()), BINS + 1)


def source_14e() -> pd.DataFrame:
    """Per-arm guide counts, medians and means."""
    return depth_stats().reset_index()


def draw_14e(outdir=None, **kw):
    """Cells per guide as four histograms, one per arm.

    The median (dashed) and mean (solid) lines are drawn over each. The text
    block sits on whichever side the distribution leaves empty, which is the
    right only for T47D/CBE, whose deepest guide falls short of the middle of
    the shared axis.
    """
    c = cells_per_guide()
    stats = depth_stats()
    edges = bin_edges()
    mid = np.sqrt(edges[0] * edges[-1])
    fig, axes = plt.subplots(1, 4, figsize=(12, 2.8), sharex=True, sharey=True)
    for ax, arm in zip(axes, ARMS):
        vals = c.loc[c["arm"] == arm, "n_cells"]
        ax.hist(vals, bins=edges, color=ARM_COLOUR[arm])
        ax.axvline(stats.loc[arm, "median"], color="black", linestyle="--",
                   linewidth=1.2)
        ax.axvline(stats.loc[arm, "mean"], color="black", linewidth=1.2)
        right = stats.loc[arm, "max"] < mid
        ax.text(0.97 if right else 0.03, 0.95,
                f"median {stats.loc[arm, 'median']:,}\n"
                f"mean {stats.loc[arm, 'mean']:,}",
                transform=ax.transAxes, ha="right" if right else "left",
                va="top", fontsize=9)
        ax.set_title(f"{arm} (n={stats.loc[arm, 'n_guides']})")
        ax.set_xlabel("Cells per guide")
        _depth_axis(ax, "x")
        _style(ax)
    axes[0].set_ylabel("Guides")
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "cells_per_guide_per_arm.pdf"),
                    bbox_inches="tight")
    return fig, source_14e()


# ------------------------------------------------------------------- 14f

def source_14f() -> pd.DataFrame:
    """Active guide-arms per designed consequence."""
    t = targeting()
    g = t.groupby("variant_class").agg(active=("is_active", "sum"),
                                       n_guide_arms=("is_active", "size"))
    g["pct_active"] = 100 * g["active"] / g["n_guide_arms"]
    return g.reset_index()


def draw_14f(outdir=None, **kw):
    """What fraction of guide-arms were active, by designed consequence."""
    d = source_14f()
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    colours = [CLASS_COLOUR.get(c, GREY_MID) for c in d["variant_class"]]
    ax.bar(range(len(d)), d["pct_active"], 0.62, color=colours)
    for i, row in d.iterrows():
        ax.text(i, row["pct_active"],
                f"{row['active']:.0f}/{row['n_guide_arms']:.0f}",
                ha="center", va="bottom", fontsize=10)
    ax.margins(y=0.18)           # the headroom the published panel leaves
    ax.set_ylabel("% of guide-arms active")
    ax.set_xticks(range(len(d)))
    ax.set_xticklabels(d["variant_class"], rotation=30, ha="right")
    ax.set_title("Designed consequence")
    _style(ax)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "de_by_variant_class.pdf"), bbox_inches="tight")
    return fig, d


# ------------------------------------------------------------------- 14g

def source_14g() -> pd.DataFrame:
    """Guide-arms per arm, split by how many DE genes they called."""
    d = de_summary_stats()
    return d[["arm", "n_guide_arms", "active", "subthreshold", "none"]]


def _draw_14g(ax, letter=None):
    d = source_14g()
    x = np.arange(len(d))
    pct = d[["active", "subthreshold", "none"]].div(d["n_guide_arms"], axis=0) * 100
    bottom = np.zeros(len(d))
    for col, colour, label in [
            ("active", BLUE, "Active"),
            ("subthreshold", GREY_MID, f"1\u2013{ACTIVE_AT - 1} DE genes"),
            ("none", GREY_PALE, "No DE genes")]:
        ax.bar(x, pct[col], 0.65, bottom=bottom, color=colour, label=label)
        bottom += pct[col].values
    for i, row in d.iterrows():
        ax.text(i, 103, f"{row['active']:.0f}/{row['n_guide_arms']:.0f}",
                ha="center", va="bottom", fontsize=9, color=BLUE,
                fontweight="bold")
    ax.set_ylabel("% of guide-arms")
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylim(0, 112)          # room for the active/total labels
    _arm_labels(ax)
    ax.legend(frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(-0.02, 0.92))
    _title(ax, letter, f"DE calls ({d['active'].sum():.0f}/"
                       f"{d['n_guide_arms'].sum():.0f} active)")
    _style(ax)


def draw_14g(outdir=None, **kw):
    """How many guide-arms called differential expression, per arm."""
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    _draw_14g(ax)
    fig.tight_layout()
    return fig, source_14g()


# ------------------------------------------------------------------- 14h

def source_14h() -> pd.DataFrame:
    """Cells and DE genes for each of the 288 targeting guide-arms."""
    t = targeting()
    return t[["arm", "display_label", "num_cells", "num_de_genes", "is_active"]]


def rho_cells_vs_de() -> pd.Series:
    """Spearman rho of cells against DE genes, per arm -- 14h's title.

    Recomputed rather than read, and it lands on the delivered table's own
    values: +0.064, -0.025, -0.046, +0.051.
    """
    t = targeting()
    return pd.Series(
        {a: g["num_cells"].corr(g["num_de_genes"], method="spearman")
         for a, g in t.groupby("arm")}).reindex(ARMS)


def _draw_14h(ax, letter=None):
    d = source_14h()
    for arm in ARMS:
        g = d[d["arm"] == arm]
        ax.scatter(g["num_cells"], g["num_de_genes"], s=14,
                   color=ARM_COLOUR[arm], label=arm, linewidths=0, alpha=0.85)
    ax.axhline(ACTIVE_AT, color=RED, linestyle="--", linewidth=1.2)
    _depth_axis(ax, "x")
    ax.set_xlabel("Cells per guide")
    ax.set_ylabel("DE genes")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    # The published title quotes `rho_cells_vs_de` from the delivered table,
    # which is already rounded to three places. Formatting the full-precision
    # fit instead turns A549/CBE's -0.0247 into "-0.02" where the panel says
    # "-0.03". `rho_cells_vs_de()` recomputes it as a check; this prints it.
    rho = ", ".join(f"{v:+.2f}" for v in de_summary_stats()["rho_cells_vs_de"])
    _title(ax, letter, f"Not power-limited ($\\rho$ {rho})")
    _style(ax)


def draw_14h(outdir=None, **kw):
    """DE genes against sequencing depth: the calls are not depth-driven."""
    fig, ax = plt.subplots(figsize=(3.8, 3.2))
    _draw_14h(ax)
    fig.tight_layout()
    return fig, source_14h()


# ------------------------------------------------------------------- 14i

def source_14i() -> pd.DataFrame:
    """DE genes in each cell line for the guides tried in both.

    Paired on the protospacer, so a guide meets itself across the two lines:
    144 pairs, which is every targeting guide-arm counted once.
    """
    t = targeting()
    w = t.pivot_table(index=["guide", "editor", "display_label"],
                      columns="cell_line", values="num_de_genes").dropna()
    w = w.reset_index()
    w["active_either"] = w[["A549", "T47D"]].max(axis=1) >= ACTIVE_AT
    return w


def rho_between_lines() -> float:
    """Spearman rho between the two cell lines -- 14i's title. 0.264."""
    w = source_14i()
    return w["A549"].corr(w["T47D"], method="spearman")


def _draw_14i(ax, letter=None):
    w = source_14i()
    lim = max(w["A549"].max(), w["T47D"].max()) * 1.08
    ax.plot([0, lim], [0, lim], color="0.75", linewidth=1, zorder=1)
    rest = w[~w["active_either"]]
    hit = w[w["active_either"]]
    ax.scatter(rest["A549"], rest["T47D"], s=18, color="0.8", linewidths=0,
               zorder=2)
    ax.scatter(hit["A549"], hit["T47D"], s=28, color=BLUE, linewidths=0,
               zorder=3)
    for _, r in hit.iterrows():
        ax.annotate(r["display_label"], (r["A549"], r["T47D"]),
                    textcoords="offset points", xytext=(6, 2), fontsize=8)
    ax.set_xlabel("DE genes, A549")
    ax.set_ylabel("DE genes, T47D")
    _title(ax, letter,
           f"Same guide, both lines ($\\rho$ = {rho_between_lines():.2f})")
    _style(ax)


def draw_14i(outdir=None, **kw):
    """The same guide in both cell lines, DE genes against DE genes."""
    fig, ax = plt.subplots(figsize=(3.8, 3.2))
    _draw_14i(ax)
    fig.tight_layout()
    return fig, source_14i()


# ------------------------------------------------------------------- 14j

def source_14j() -> pd.DataFrame:
    """The eight active guide-arms, with cells and dominant edit rate."""
    return de_hits()


def _draw_14j(ax, letter=None):
    d = source_14j()
    y = np.arange(len(d))[::-1]          # longest bar at the top
    ax.barh(y, d["num_de_genes"], 0.62,
            color=[ARM_COLOUR[a] for a in d["arm"]])
    ax.axvline(ACTIVE_AT, color=RED, linestyle="--", linewidth=1.2)
    for yy, (_, r) in zip(y, d.iterrows()):
        ax.text(r["num_de_genes"] + 1, yy,
                f"{r['n_cells']:,.0f} cells · "
                f"{r['dominant_edit_pct']:.0f}% edit",
                va="center", fontsize=8, color=ANNOT)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['display_label']} ({r['arm']})"
                        for _, r in d.iterrows()], fontsize=8)
    ax.set_xlabel("DE genes")
    ax.set_xlim(0, d["num_de_genes"].max() * 1.45)
    _title(ax, letter, "Active guide-arms")
    _style(ax)


def draw_14j(outdir=None, **kw):
    """Every guide-arm that cleared the activity threshold."""
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    _draw_14j(ax)
    fig.tight_layout()
    return fig, source_14j()


# ------------------------------------------------------------------- 14k

def source_14k() -> dict:
    """The clustermap's inputs: its annotation tracks and its null bands.

    Not redrawn -- see this module's docstring. What is here is what the
    published panel annotates with, so the figure can be read against it.
    """
    return {"clustermap_null": io.read(FIGS / "clustermaps" /
                                       "clustermap_null.csv"),
            "guide_annotation": guide_annotation()}


# ------------------------------------------------------------ the composites

def draw_experiment_summary(outdir=None, **kw):
    """14a-14d on one row, as `experiment_summary.pdf` has them."""
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.4))
    for ax, (letter, fn) in zip(axes, [("a", _draw_14a), ("b", _draw_14b),
                                       ("c", _draw_14c), ("d", _draw_14d)]):
        fn(ax, letter)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "experiment_summary.pdf"), bbox_inches="tight")
    return fig


def draw_de_summary(outdir=None, **kw):
    """14g-14j on one row, as `de_summary.pdf` has them."""
    fig, axes = plt.subplots(1, 4, figsize=(17, 3.4))
    for ax, (letter, fn) in zip(axes, [("a", _draw_14g), ("b", _draw_14h),
                                       ("c", _draw_14i), ("d", _draw_14j)]):
        fn(ax, letter)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "de_summary.pdf"), bbox_inches="tight")
    return fig


#: the four composites that can be redrawn whole, for comparison with the
#: delivered PDFs. `pooled_corr` is absent for the reason in the docstring.
COMPOSITE_FIGURES = {"experiment_summary": draw_experiment_summary,
                     "cells_per_guide": draw_14e,
                     "de_by_variant_class": draw_14f,
                     "de_summary": draw_de_summary}

PANELS = {"S14a": source_14a, "S14b": source_14b, "S14c": source_14c,
          "S14d": source_14d, "S14e": source_14e, "S14f": source_14f,
          "S14g": source_14g, "S14h": source_14h, "S14i": source_14i,
          "S14j": source_14j, "S14k": source_14k}

#: 14k alone is not redrawn: the gene filter behind it is not in this repo
FIGURES = {"S14a": draw_14a, "S14b": draw_14b, "S14c": draw_14c,
           "S14d": draw_14d, "S14e": draw_14e, "S14f": draw_14f,
           "S14g": draw_14g, "S14h": draw_14h, "S14i": draw_14i,
           "S14j": draw_14j}
