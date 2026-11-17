"""Supplementary figure 6e-6j: the CDK9 Boltz-2 panels.

From `manuscript_source_data/cdk9_boltz/cdk9_ifp.ipynb`. Each panel below names
the cell it came from and the file that cell writes, which is the file
`docs/source_data_mapping_supp.csv` lists for that panel:

    6e  cell 12   plots/model_confidence_vs_binding_probability.svg
    6f  cell 8    plots/cdk9_delta_binding_probabilities_l156f.svg
    6h  cell 23   plots/pca_WT_all_models.svg
    6i  cell 28   plots/PC1_features_tb_5.svg
    6j  cell 26   plots/PC2_feature_loadings_WT_all_models.png

6a, 6g and 6k-6m have no source file recorded, and 6c and 6d are western blots.
6b is still marked Missing.

The analysis itself is `cdkbe.boltz.ifp`; this module only draws. Two numbers
the README states independently of the figures -- 132 retained fingerprint
features and 58 retained WT ATP-pocket features -- come out of the port exactly,
which is the check that the fingerprint was carried over faithfully.

**6i and 6j are not the same kind of plot**, despite their names suggesting a
pair. 6i is the five most negative and five most positive PC1 loadings (cell 28,
`PC1_features_tb_5`); 6j is the ten largest-magnitude PC2 loadings (cell 26).
The notebook also has a PC1 version of 6j's style, cell 25, which the mapping
does not reference. Their colour rules differ too, and in opposite directions --
see `_loading_panel`.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..boltz import ifp

#: the notebook sets these at the top of every figure cell
RC = {"svg.fonttype": "none"}


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(1)
    ax.tick_params(axis="both", colors="black", labelsize=12)


def _palette():
    import seaborn as sns

    return dict(zip(ifp.DRUG_ORDER,
                    sns.color_palette("tab10", len(ifp.DRUG_ORDER))))


# ------------------------------------------------------------------- 6e

def source_6e() -> pd.DataFrame:
    """Per-seed mean iPTM and predicted binding probability, WT."""
    return ifp.seed_level("WT")


def draw_6e(outdir=None, **kw):
    """Model confidence against binding probability, WT, one point per seed."""
    s = source_6e()
    pal = _palette()
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.grid(True, linewidth=0.5)
        for drug in ifp.DRUG_ORDER:
            sub = s[s["drug"] == drug]
            c = pal[drug]
            ax.scatter(sub["iptm"], sub["binding_probability"], s=45, color=c,
                       alpha=0.45, edgecolor="none", zorder=2)
            mx, my = sub["iptm"].mean(), sub["binding_probability"].mean()
            ax.errorbar(mx, my, xerr=sub["iptm"].std(ddof=1),
                        yerr=sub["binding_probability"].std(ddof=1),
                        fmt="o", ms=11, color=c, ecolor=c, elinewidth=1.5,
                        capsize=4, zorder=3, markeredgecolor="white",
                        markeredgewidth=1.2)
            ax.annotate(ifp.DRUG_LABELS[drug], (mx, my), textcoords="offset points",
                        xytext=(10, 6), fontsize=11, color=c)
        ax.set_xlabel("iPTM (interface confidence)", fontsize=12)
        ax.set_ylabel("Binding probability (Boltz-2)", fontsize=12)
        ax.set_title("Model confidence vs binding probability", fontsize=14)
        ax.margins(0.1)
        for sp in ax.spines.values():
            sp.set_color("black")
            sp.set_linewidth(0.5)
        ax.tick_params(axis="both", colors="black", labelsize=11)
        fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "model_confidence_vs_binding_probability.svg"),
                    format="svg", bbox_inches="tight", transparent=True, dpi=600)
    return fig, s


# ------------------------------------------------------------------- 6f

def source_6f() -> pd.DataFrame:
    """L156F - WT binding probability, paired within a seed."""
    return ifp.delta_binding_probability("L156F")


def draw_6f(outdir=None, **kw):
    """Change in predicted binding probability, L156F, one panel per ligand."""
    d = source_6f()
    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, len(ifp.DRUG_ORDER), figsize=(12, 4), sharey=True)
        rng = np.random.default_rng(1)
        for ax, drug in zip(axes, ifp.DRUG_ORDER):
            ax.grid(False)
            vals = d.loc[d["drug"] == drug, "delta_binding_probability"].dropna().values
            ax.boxplot([vals], positions=[0], widths=0.4, showfliers=False,
                       medianprops={"color": "k", "linewidth": 1})
            ax.scatter(rng.normal(0, 0.025, len(vals)), vals, s=50, alpha=0.8,
                       color="red", zorder=3)
            ax.axhline(0, linestyle="--", linewidth=1.5, c="gray")
            ax.set_xticks([0])
            ax.set_xticklabels([ifp.DRUG_LABELS.get(drug, drug)], fontsize=12)
            ax.set_xlim(-0.5, 0.5)
            ax.set_ylim(-0.5, 0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        axes[0].set_ylabel(r"$\Delta$ Binding probability", fontsize=13)
        fig.suptitle("Change in predicted binding probability (L156F)", fontsize=16)
        fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "cdk9_delta_binding_probabilities_l156f.svg"),
                    format="svg", bbox_inches="tight", transparent=True, dpi=600)
    return fig, d


# ------------------------------------------------------------------- 6h

def source_6h() -> pd.DataFrame:
    """PC1/PC2 of the weighted WT ATP-pocket fingerprint PCA, all 36 poses."""
    X = ifp.pocket_matrix("WT")
    scores, var, _, _ = ifp.weighted_pca(X)
    out = pd.DataFrame(scores, columns=["PC1", "PC2"], index=X.index).reset_index()
    out.attrs["explained"] = tuple(var)
    return out


def draw_6h(outdir=None, **kw):
    """ATP-pocket fingerprints in PC space.

    The PCA is fitted on all 36 WT poses; only the model-0 pose of each seed is
    drawn. It is not refitted on the model-0 subset.
    """
    X = ifp.pocket_matrix("WT")
    scores, var, _, _ = ifp.weighted_pca(X)
    drug_all = X.index.get_level_values("drug").to_numpy()
    model_all = X.index.get_level_values("model").astype(int).to_numpy()
    pal = _palette()
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.grid(False)
        ax.axhline(0, color="0.7", lw=0.8, zorder=0, linestyle="--")
        ax.axvline(0, color="0.7", lw=0.8, zorder=0, linestyle="--")
        for dg in sorted(set(drug_all)):
            # The notebook scatters the non-model-0 poses at alpha=0. They are
            # invisible, but they still set the data limits, which is why the
            # published axes are wider than the drawn points need. Dropping
            # this as a no-op shrinks the axes and the legend lands on a point.
            other = (drug_all == dg) & (model_all != 0)
            ax.scatter(scores[other, 0], scores[other, 1], s=100, color=pal[dg],
                       alpha=0.0, edgecolor="none", zorder=2)
            m0 = (drug_all == dg) & (model_all == 0)
            ax.scatter(scores[m0, 0], scores[m0, 1], s=100, color=pal[dg],
                       alpha=0.9, edgecolor="black", linewidth=0.8,
                       label=ifp.DRUG_LABELS.get(dg, dg), zorder=3)
        ax.set_xlabel(f"PC1 ({var[0]:.0f}%)", fontsize=12)
        ax.set_ylabel(f"PC2 ({var[1]:.0f}%)", fontsize=12)
        ax.set_title("ATP-pocket fingerprints", fontsize=14)
        ax.legend(fontsize=11, frameon=False, loc="lower right")
        _style(ax)
        fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "pca_WT_all_models.svg"), format="svg",
                    bbox_inches="tight", transparent=True, dpi=600)
    return fig, source_6h()


# --------------------------------------------------------------- 6i and 6j

def _loading_panel(values: pd.Series, title: str, positive_colour: str):
    """One loading scatter, labels coloured to match their points.

    `positive_colour` differs between the two published panels -- cell 28 draws
    positive loadings blue and negative red, cell 26 the other way round -- so
    it is an argument rather than a constant.
    """
    negative_colour = "red" if positive_colour == "cornflowerblue" else "cornflowerblue"
    colours = [positive_colour if v > 0 else negative_colour for v in values.values]
    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.grid(False)
        ax.scatter(values.values, values.index, s=100, color=colours, alpha=0.85,
                   edgecolor="k")
        by_label = dict(zip(values.index, colours))
        for tick in ax.get_yticklabels():
            tick.set_color(by_label[tick.get_text()])
        ax.axvline(0, color="0.5", lw=0.8, linestyle="--")
        ax.set_xlabel(f"{values.name.split('_')[0]} loading", fontsize=15)
        ax.set_ylabel("")
        ax.set_title(title, fontsize=15)
        lim = max(0.5, np.ceil(np.abs(values.values).max() * 10) / 10)
        ax.set_xlim(-lim, lim)
        _style(ax)
        fig.tight_layout()
    return fig


def source_6i() -> pd.Series:
    """The five most negative and five most positive PC1 loadings."""
    pc1 = ifp.pc_loadings(0)
    return pd.concat([pc1.nsmallest(5), pc1.nlargest(5)]).sort_values()


def draw_6i(outdir=None, **kw):
    """ATP-pocket features at the two ends of PC1."""
    v = source_6i()
    fig = _loading_panel(v, "ATP-pocket IFP features along PC1", "cornflowerblue")
    if outdir is not None:
        fig.savefig(str(outdir / "PC1_features_tb_5.svg"), format="svg",
                    bbox_inches="tight", transparent=True, dpi=600)
    return fig, v.rename("PC1_loading").to_frame()


def source_6j() -> pd.Series:
    """The ten largest-magnitude PC2 loadings."""
    pc2 = ifp.pc_loadings(1)
    ordered = pc2.reindex(pc2.abs().sort_values(ascending=False).index)
    return ordered.head(10).sort_values()


def draw_6j(outdir=None, **kw):
    """The WT ATP-pocket features that dominate PC2."""
    v = source_6j()
    fig = _loading_panel(v, "WT ATP-pocket IFP features along PC2", "red")
    if outdir is not None:
        fig.savefig(str(outdir / "PC2_feature_loadings_WT_all_models.png"),
                    dpi=600, bbox_inches="tight")
    return fig, v.rename("PC2_loading").to_frame()


PANELS = {"S6e": source_6e, "S6f": source_6f, "S6h": source_6h,
          "S6i": source_6i, "S6j": source_6j}
FIGURES = {"S6e": draw_6e, "S6f": draw_6f, "S6h": draw_6h,
           "S6i": draw_6i, "S6j": draw_6j}
