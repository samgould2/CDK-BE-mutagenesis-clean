"""LFC correlation clustermaps (supplementary figure 4).

Ported from `04_LFC_FDR_tables/LFC_correlation_clustermaps.ipynb`.

For each editor and subpool, every comparison the screen supports -- each drug
against DMSO, DMSO against T0, and T0 against Plasmid -- becomes a column of a
wide gRNA x comparison table, and the panel is the correlation matrix of those
columns, clustered on correlation distance.

Two choices from the original worth keeping in view:

* the table is reduced to complete cases (`dropna`), so every correlation in a
  panel is computed over exactly the same guides, and
* subpool 1 carries the SY-5609 retest as a second screen joined on gRNA_id, so
  its matrix has two `DMSO vs. T0` and two `T0 vs. Plasmid` columns -- the
  suffix is what tells them apart, and their correlation is a useful read on
  how reproducible the shared arms are.
"""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from . import io, paths

LFC_STAT = "LFC_median"
CORR_METHOD = "pearson"
GUIDE_CLASSES = ["targeting", "intron", "non-targeting control"]
EDITORS = ["ABE", "CBE"]

#: subpool -> {Screen_ID in the master table: label suffix for its shared arms}
SUBPOOLS = {
    "subpool1": {"Subpool1": "", "SY-5609": " (SY-5609 retest)"},
    "subpool2": {"Subpool2": ""},
    "subpool3": {"Subpool3": ""},
}


def master(editor: str) -> pd.DataFrame:
    return io.read(paths.LFC / f"{editor}_master_LFC_table.csv")


def pretty_condition(cond: str) -> str:
    """'KB-0742_1500nM' -> 'KB-0742 1,500 nM'; no concentration -> unchanged."""
    m = re.match(r"^(.*)_(\d+)nM$", cond)
    return f"{m.group(1)} {int(m.group(2)):,} nM" if m else cond


def screen_lfcs(df: pd.DataFrame, screen_id: str, suffix: str = "") -> pd.DataFrame:
    """Wide gRNA x comparison LFCs for one screen."""
    s = df[df["Screen_ID"] == screen_id]
    pre = f"{LFC_STAT}_"

    dmso = s[s["Base_Value"] == "DMSO"].set_index("gRNA_id")
    assert dmso.index.is_unique, f"{screen_id}: duplicate gRNA_id in the DMSO arm"
    drug_cols = [c for c in dmso.columns
                 if c.startswith(pre) and c.endswith("_DMSO")
                 and c not in (f"{pre}T0_DMSO", f"{pre}Plasmid_DMSO")
                 and dmso[c].notna().any()]
    drugs = dmso[drug_cols]
    drugs.columns = [pretty_condition(c[len(pre):-len("_DMSO")]) + " vs. DMSO"
                     for c in drug_cols]

    t0 = s[s["Base_Value"] == "T0"].set_index("gRNA_id")
    assert t0.index.is_unique, f"{screen_id}: duplicate gRNA_id in the T0 arm"
    dmso_t0 = t0[[f"{pre}DMSO_T0"]].rename(
        columns={f"{pre}DMSO_T0": f"DMSO vs. T0{suffix}"})

    pl = s[s["Base_Value"] == "Plasmid"].set_index("gRNA_id")
    assert pl.index.is_unique, f"{screen_id}: duplicate gRNA_id in the Plasmid arm"
    t0_pl = pl[[f"{pre}T0_Plasmid"]].rename(
        columns={f"{pre}T0_Plasmid": f"T0 vs. Plasmid{suffix}"})

    return pd.concat([drugs, dmso_t0, t0_pl], axis=1)


def subpool_lfcs(df: pd.DataFrame, screen_ids: dict) -> pd.DataFrame:
    wide = pd.concat([screen_lfcs(df, sid, suf) for sid, suf in screen_ids.items()],
                     axis=1, join="inner")
    classes = df.drop_duplicates("gRNA_id").set_index("gRNA_id")["classification"]
    wide = wide[classes.reindex(wide.index).isin(GUIDE_CLASSES)]
    # complete cases only, so every correlation in the panel uses one guide set
    return wide.dropna()


def correlation(editor: str, subpool: str):
    """-> (correlation matrix, the wide table it was computed from)."""
    wide = subpool_lfcs(master(editor), SUBPOOLS[subpool])
    return wide.corr(method=CORR_METHOD), wide


def clustermap(corr: pd.DataFrame, title: str):
    n = corr.shape[0]
    # correlation distance, average linkage; the same ordering on rows and columns
    Z = hierarchy.linkage(squareform(1 - corr.values, checks=False), method="average")
    size = 0.6 * n + 3
    g = sns.clustermap(
        corr, row_linkage=Z, col_linkage=Z,
        cmap="RdBu_r", vmin=-1, vmax=1, center=0,
        annot=True, fmt=".2f", annot_kws={"fontsize": 8 if n > 10 else 10},
        linewidths=1, linecolor="white", square=True,
        figsize=(size, size), dendrogram_ratio=0.12,
        cbar_kws={"label": f"{CORR_METHOD.capitalize()} r",
                  "ticks": [-1, -0.5, 0, 0.5, 1]},
        cbar_pos=(1.02, 0.35, 0.025, 0.3))
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(axis="both", length=0, labelsize=11)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.fig.suptitle(title, fontsize=14, y=1.02)
    return g


def panel(editor: str, subpool: str, outdir=None):
    """-> (figure, correlation matrix) for one editor x subpool clustermap."""
    corr, wide = correlation(editor, subpool)
    title = (f'{editor} {subpool.replace("subpool", "Subpool ")} '
             f'({LFC_STAT.replace("LFC_", "")} LFC, n = {len(wide):,} guides)')
    g = clustermap(corr, title)
    if outdir is not None:
        g.savefig(str(outdir / f"{editor}_{subpool}_LFC_correlation_clustermap.pdf"),
                  transparent=True, bbox_inches="tight")
    return g.fig, corr
