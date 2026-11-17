"""Draw figure 1's panels.

1e — editing_analysis_v2.ipynb cell 23
1f — fig1_analysis.ipynb cell 13
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import counts, io, lollipop, paths
from ..plotting import CDK_COLORS
from . import library as lib

GENES = lib.GENES


# ------------------------------------------------------------------ 1e

def designed_coverage() -> pd.DataFrame:
    """Codon coverage of the designed library, before any editing is measured.

    Drawn as the dashed ceilings on 1e: the fraction of each protein's codons
    that a guide was even designed for, which is what the measured curves
    approach as the editing cutoff drops.
    """
    cdk_info = io.read(paths.REFERENCE / "library_generation" / "cdks_info_updated.csv")
    cdk_info = cdk_info[cdk_info["Gene"].isin(GENES)]
    rows = []
    for g in GENES:
        n_codons = len(cdk_info.loc[cdk_info["Gene"] == g, "Protein"].values[0])
        d = io.read(paths.REFERENCE / "library_generation" / "edit_libraries" / f"{g}_edits.csv")
        d = d.dropna(subset=["Codon"])
        rows.append(dict(
            Gene=g,
            combined=100 * d["Codon"].nunique() / n_codons,
            ABE=100 * d.loc[d["Editor"] == "ABE", "Codon"].nunique() / n_codons,
            CBE=100 * d.loc[d["Editor"] == "CBE", "Codon"].nunique() / n_codons))
    return pd.DataFrame(rows).set_index("Gene")


def figure_1e(coverage: pd.DataFrame | None = None, outdir=None, figsize=(12, 3.8)):
    """Codon coverage vs sensor-editing cutoff: CBE, ABE, and the two combined."""
    d = lib.figure_1e() if coverage is None else coverage
    designed = designed_coverage()
    fig, ax = plt.subplots(1, 3, figsize=figsize, sharey=True, sharex=True)

    for col, a, title, designed_col in (("CBE_codons", ax[0], "CBE", "CBE"),
                                        ("ABE_codons", ax[1], "ABE", "ABE"),
                                        ("total_codons", ax[2], "Combined ABE & CBE", "combined")):
        for g in GENES:
            s = d[d["Gene"] == g].sort_values("editing_cutoff")
            a.plot(s["editing_cutoff"], s[col], label=g, linewidth=2, alpha=0.8,
                   c=CDK_COLORS[g])
            # dashed: what the designed library could reach at best
            a.plot([0, 100], [designed.loc[g, designed_col]] * 2, linewidth=0.5,
                   linestyle="dashed", alpha=0.8, c=CDK_COLORS[g])
        a.set_xlabel("Sensor Editing Cutoff (%)", fontsize=15)
        a.set_title(title, fontsize=15)
        a.spines[["right", "top"]].set_visible(False)
        a.tick_params(axis="both", which="major", labelsize=15)

    ax[0].set_ylabel("Percentage of AAs Mutated", fontsize=15)
    ax[0].set_xlim(0, 101)
    ax[0].set_ylim(0)
    ax[0].set_xticks([0, 20, 40, 60, 80, 100])
    ax[2].legend(fontsize=10, ncol=2)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_1e_codon_coverage.pdf"), transparent=True)
    return fig, d


# ------------------------------------------------------------------ 1f

def protein_lengths() -> dict:
    """Canonical isoform length per CDK, from the per-residue disorder tables."""
    return {g: len(io.read(paths.REFERENCE / "disorder" / f"{g}_disorder_scores.csv"))
            for g in GENES}


def klifs_table() -> pd.DataFrame:
    return io.read(paths.REFERENCE / "KLIFS" / "KLIFS_compiled.csv")


def figure_1f(variants: pd.DataFrame | None = None, outdir=None,
              res_per_inch: float = 200, panel_h: float = 0.72, gap_h: float = 0.10,
              pad_l: float = 1.25, pad_r: float = 0.15, pad_t: float = 0.10,
              pad_b: float = 0.55, bar: int = 200, tick_w: float = 0.30,
              align_domains: bool = True):
    """Resistance variants on each CDK's own residue axis.

    Every panel is drawn in native residue coordinates at a shared
    inches-per-residue scale, so a shared axis is avoided (it would
    misrepresent the proteins' relative sizes).

    `align_domains` shifts each panel so the kinase domain starts at the same x
    across genes. **The published figure does not do this** — measured on
    `res_ABE_CBE_combined_length_scaled.pdf`, all nine panels share one left
    edge, while the notebook cell that names that file offsets them by up to
    260 pt. The figure on disk predates the alignment change. `True` follows
    the code as last saved; pass `False` to reproduce the published panel.
    """
    d = lib.figure_1f() if variants is None else variants
    klifs = klifs_table()
    length = protein_lengths()
    start = {g: lollipop.find_ranges(klifs[g].dropna())[0][0] for g in GENES}

    scale = 1.0 / res_per_inch
    max_n = max(start[g] for g in GENES)
    max_c = max(length[g] - start[g] for g in GENES)
    fig_w = pad_l + (max_n + max_c) * scale + pad_r
    fig_h = pad_t + len(GENES) * panel_h + (len(GENES) - 1) * gap_h + pad_b
    fig = plt.figure(figsize=(fig_w, fig_h))

    axes = {}
    for i, g in enumerate(GENES):
        L = length[g]
        left = pad_l + ((max_n - start[g]) * scale if align_domains else 0.0)
        bot = fig_h - pad_t - (i + 1) * panel_h - i * gap_h
        ax = fig.add_axes([left / fig_w, bot / fig_h, (L * scale) / fig_w, panel_h / fig_h],
                          sharey=axes[GENES[0]] if i else None)
        axes[g] = ax

        for lo, hi in lollipop.find_ranges(klifs[g].dropna()):
            ax.axvspan(lo, hi, alpha=0.2)

        for editor, colour in (("CBE", "tab:blue"), ("ABE", "tab:purple")):
            s = d[(d["gene"] == g) & (d["Editor"] == editor) & (d["Codon"] != "WT")]
            pos = np.array([int(c) for c in s["Codon"]])
            ax.vlines(pos, 0, s["max_LFC"], linewidth=0.4, color="black", alpha=0.2)
            ax.scatter(pos, s["max_LFC"], c=colour, s=20, edgecolors="black",
                       linewidths=0.5, clip_on=False)

        ax.axhline(0, linestyle="dashed", color="black", linewidth=0.5)
        ax.set_xlim(0, L)
        ax.set_xticks([])
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.tick_params(axis="y", left=True, labelleft=True, labelsize=7, length=2.5, pad=1.5)
        ax.text(1, 0, f"{L} aa", transform=ax.transAxes, ha="right", va="top",
                fontsize=6, color="0.45")
        fig.text((left - tick_w) / fig_w, (bot + panel_h / 2) / fig_h, g,
                 ha="right", va="center", fontsize=11)

    # a scale bar, since each panel now has its own x coordinates
    y0, x0 = (pad_b - 0.32) / fig_h, pad_l / fig_w
    fig.add_artist(plt.Line2D([x0, x0 + bar * scale / fig_w], [y0, y0], color="black", lw=1.2))
    fig.text(x0 + bar * scale / fig_w / 2, y0 - 0.008, f"{bar} aa",
             ha="center", va="top", fontsize=8)
    fig.text(0.008, (pad_b + (fig_h - pad_b - pad_t) / 2) / fig_h, "max LFC",
             rotation=90, va="center", ha="left", fontsize=10)

    if outdir is not None:
        fig.savefig(str(outdir / "Figure_1f_variants_by_residue.pdf"), transparent=True)
    return fig, d


FIGURES = {"1e": figure_1e, "1f": figure_1f}
