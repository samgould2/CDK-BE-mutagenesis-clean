"""Cross-drug resistance heatmap and the hits-per-drug bar chart (figures 4a, 4b).

Ported from logo_plots_CDK2_4_6_canonical.ipynb
cells 29 and 33. The drawing logic is kept as it was -- its comments record
several decisions that are easy to undo by accident, in particular that the
row linkage is computed on the six drug columns *before* the T0 column is
appended, so baseline fitness never drives the drug clustering and a `_withT0`
figure keeps exactly the row order of its plain counterpart.

What changed: the function took its parameters from notebook globals, so they
are arguments here, and the file writing moved to the caller.
"""
from __future__ import annotations
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage

from .plotting import EDITOR_COLORS, GENE_COLORS_S3 as GENE_COLORS

#: the six CDK4/6- and CDK2-directed inhibitors, in the repo's display order
CONDITIONS = {"RIBO": "Ribociclib", "ABEMA": "Abemaciclib", "ATIRMO": "Atirmociclib",
              "PALBO": "Palbociclib", "TAGTO": "Tagtociclib", "INX": "INX-315"}
CLASS_CONDS = {"CDK4/6i": ["RIBO", "ABEMA", "ATIRMO", "PALBO"], "CDK2i": ["TAGTO", "INX"]}
GENES = ["CDK2", "CDK4", "CDK6"]

P_CUTOFF = 0.01          # empirical p; defines which conditions a guide is resistant in
T0_FDR_CUTOFF = 0.01     # on the T0 track an asterisk means FDR < this, either direction

#: every figure was written twice -- default text size and a '_large' variant
STYLES = {
    "": dict(ylabel=11, xtick=11, xlabel=12, logo=9, title=11, legend=10,
             panel_title=12, tick_len=3.5, tick_w=0.8),
    "_large": dict(ylabel=17, xtick=15, xlabel=18, logo=11, title=12, legend=14,
                   panel_title=16, tick_len=8.0, tick_w=1.8),
}


def cross_drug_heatmap(sub: pd.DataFrame, title: str, st: dict | None = None,
                       suffix: str = "", per_gene: bool = False, width: float = 9.5,
                       include_T0: bool = False):
    """Clustered LFC heatmap for a set of hits.

    Editor is carried by the colour of the row label (purple ABE, blue CBE); a
    gene colour block is drawn only when the figure mixes genes.
    """
    st = STYLES[suffix] if st is None else st
    # bbox_inches='tight' grows the saved figure to fit the title, so a long one-line
    # title silently un-narrows the narrow variants. Wrap to roughly the figure width
    # and size the top margin from the resulting line count.
    wrapped = "\n".join(textwrap.fill(line, max(40, int(width * 11.5)))
                        for line in title.split("\n"))
    title_h = (wrapped.count("\n") + 1) * st["title"] * 1.3 / 72        # inches

    m = sub[[f"LFC_{c}" for c in CONDITIONS]].copy()
    m.columns = list(CONDITIONS.values())
    # the gene is already in the title of a per-gene figure
    m.index = [r["Top_HGVSp"] if per_gene else f'{r["Gene"]} {r["Top_HGVSp"]}'
               for _, r in sub.iterrows()]
    sig = np.column_stack([(sub[f"p_{c}"] < P_CUTOFF) & (sub[f"LFC_{c}"] > 0)
                           for c in CONDITIONS])
    editors = list(sub["Editor"])

    # correlation distance groups variants by the SHAPE of their six-drug profile;
    # a row with zero variance has no defined correlation, so fall back to euclidean
    metric = "correlation" if (m.to_numpy().std(axis=1) > 0).all() else "euclidean"
    Z = linkage(m.to_numpy(), method="average", metric=metric)

    if include_T0:
        # appended AFTER the linkage: baseline fitness must not drive drug clustering
        m["T0 vs plasmid"] = sub["LFC_T0"].to_numpy()
        sig = np.column_stack([sig, (sub["FDR_T0"] < T0_FDR_CUTOFF).to_numpy()])

    v = np.nanmax(np.abs(m.to_numpy()))        # symmetric limits put LFC 0 at the midpoint

    kw = {}
    if not per_gene:
        kw["row_colors"] = pd.DataFrame(
            {"Gene": [GENE_COLORS[g] for g in sub["Gene"]]}, index=m.index)

    cg = sns.clustermap(m, row_linkage=Z, col_cluster=False, cmap="RdBu_r", center=0,
                        vmin=-v, vmax=v, figsize=(width, 0.20 * len(m) + 3.2),
                        dendrogram_ratio=(0.11, 0.015), colors_ratio=0.018,
                        cbar_pos=None, linewidths=0.3, linecolor="white", **kw)
    fig, ax = cg.figure, cg.ax_heatmap
    H = fig.get_figheight()
    # clustermap fills the figure, so title/legend/colorbar would land on the first
    # rows; shrink the grid and reflow. Offsets are inches over H so they do not
    # drift with the number of rows.
    cg.gs.update(top=1 - (title_h + 1.25) / H)
    for a in [ax, cg.ax_row_dendrogram] + ([cg.ax_row_colors] if not per_gene else []):
        a.set_position(a.get_subplotspec().get_position(fig))
    order = cg.dendrogram_row.reordered_ind

    for i, orig in enumerate(order):
        for j in range(m.shape[1]):
            if sig[orig, j]:
                ax.text(j + 0.5, i + 0.52, "*", ha="center", va="center",
                        fontsize=st["logo"], color="black")

    ax.axvline(len(CLASS_CONDS["CDK4/6i"]), color="black", linewidth=2)
    if include_T0:
        ax.axvline(len(CONDITIONS), color="black", linewidth=2)   # the T0 track is not a drug
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=st["xtick"])
    ax.tick_params(axis="y", labelsize=(5 if suffix == "" else 6.5) if len(m) > 55
                   else (7 if suffix == "" else 9))
    ax.set_ylabel("")
    for lbl, orig in zip(ax.get_yticklabels(), order):
        lbl.set_color(EDITOR_COLORS[editors[orig]])

    cax = fig.add_axes([0.33, 1 - (title_h + 0.85) / H, 0.26, 0.13 / H])
    cb = fig.colorbar(ax.collections[0], cax=cax, orientation="horizontal")
    cb.set_label("Median LFC (vs. DMSO)", fontsize=st["legend"])
    cb.ax.tick_params(labelsize=st["logo"])

    handles = [Patch(facecolor=c, edgecolor="none", label=e) for e, c in EDITOR_COLORS.items()]
    if not per_gene:
        handles = [Patch(facecolor=c, edgecolor="none", label=g)
                   for g, c in GENE_COLORS.items()] + handles
    fig.legend(handles=handles, ncol=5, frameon=False, fontsize=st["legend"],
               loc="upper center", bbox_to_anchor=(0.5, 1 - (title_h + 0.30) / H))
    fig.suptitle(wrapped, fontsize=st["title"], y=1 - 0.14 / H)
    return fig


def hits_per_drug(hit_counts: pd.DataFrame, n_hits: int, st: dict | None = None,
                  suffix: str = "", figsize=(9, 5)):
    """Resistance hits per drug by gene.

    The gene palette is low-contrast against white (yellow and orange sit near
    2.2:1), so every bar carries a black edge and a printed count -- identity is
    never colour-alone.
    """
    st = STYLES[suffix] if st is None else st
    order = [CONDITIONS[c] for c in CLASS_CONDS["CDK4/6i"] + CLASS_CONDS["CDK2i"]]
    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(data=hit_counts, x="Drug", y="n_hits", hue="Gene", order=order,
                hue_order=GENES, palette=GENE_COLORS, ax=ax, linewidth=1,
                edgecolor="black")
    for cont in ax.containers:
        ax.bar_label(cont, fmt="%d", fontsize=st["logo"], padding=2, color="black")
    ax.axvline(len(CLASS_CONDS["CDK4/6i"]) - 0.5, color="black", linestyle="dashed",
               linewidth=1)
    ax.set_xlabel("")
    ax.set_ylabel("Resistance hits (guides)", fontsize=st["xlabel"])
    ax.tick_params(axis="both", labelsize=st["xtick"], length=st["tick_len"],
                   width=st["tick_w"])
    ax.set_ylim(0, hit_counts["n_hits"].max() * 1.18)
    ax.spines[["right", "top"]].set_visible(False)
    ax.legend(title="Gene", fontsize=st["legend"], title_fontsize=st["legend"],
              frameon=False, loc="upper left")
    ax.set_title(
        f"Resistance hits per drug by gene | {n_hits} hits total, counted once per "
        f"drug they are resistant in\nresistant = empirical p < {P_CUTOFF} & LFC > 0  |  "
        f"CDK4/6 inhibitors left of the divider, CDK2 inhibitors right",
        fontsize=st["title"])
    return fig
