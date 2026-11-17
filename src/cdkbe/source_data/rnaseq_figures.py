"""PC3 RNA-seq panels 3f-3h, redrawn from the tables their own scripts saved.

The re-analysis scripts under
`data/processed/rnaseq/PC3_RNAseq/reanalysis/scripts`
each write the plotted frame next to the PDF, so these drawers read that frame
rather than re-deriving anything from the DESeq2 output:

    3f  03_figures.py          -> Fig2b_DEG_counts_with_baseline.csv
    3g  03_figures.py          -> Fig3_pct_response_per_replicate.csv
    3h  06_elongation_combined -> Fig10_replicate_decile_profiles.csv

Style constants are copied from `reanalysis/scripts/_style.py` so the panels
match the published ones; that file is not importable here because it sets a
non-interactive matplotlib backend at import and resolves its own paths.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba

from . import rnaseq
from .. import io, paths

RESULTS = (paths.PROCESSED / "rnaseq" / "PC3_RNAseq"
           / "reanalysis" / "results")

# --- from reanalysis/scripts/_style.py --------------------------------------
GENOS = ["NT", "I733V", "G822E", "G731K"]
VARS = GENOS[1:]
DRUGS = ["HQ461", "ZLY025"]
#: neutral ink for the NT reference + blue / gold / plum for the variants
GCOL = {"NT": "#3D3D3D", "I733V": "#1B6CA8", "G822E": "#E1A730", "G731K": "#8E4585"}
TLAB = {"Control": "DMSO", "HQ461": "HQ461", "ZLY025": "ZLY-025"}
INK, MUTED = "#0B0B0B", "#898781"
BAR_EDGE = dict(edgecolor="black", linewidth=1.0)

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def figure_3f(outdir=None):
    """DEG counts per variant and drug, with the DMSO-vs-NT baseline alongside.

    The y-axis is shared across all three: the baseline counts really are ~10x
    smaller than the drug effects, and that comparison is the point of the panel.
    """
    d = rnaseq.figure_3f()
    counts = {(r.drug, r.genotype): (int(r.up), int(r.down)) for r in d.itertuples()}

    def bars(ax, genos, drug, title):
        for i, g in enumerate(genos):
            up, dn = counts[(drug, g)]
            ax.bar(i, up, width=0.72, color=GCOL[g], **BAR_EDGE)
            ax.bar(i, -dn, width=0.72, color=to_rgba(GCOL[g], 0.45), **BAR_EDGE)
            ax.text(i, up + 90, f"{up:,}", ha="center", va="bottom", fontsize=5.5)
            ax.text(i, -dn - 90, f"{dn:,}", ha="center", va="top", fontsize=5.5)
        ax.axhline(0, color=INK, lw=0.6)
        ax.set_xticks(range(len(genos)))
        ax.set_xticklabels(genos, rotation=90)
        ax.set_xlim(-0.65, len(genos) - 0.35)
        ax.set_title(title)
        ax.set_ylim(-5300, 4400)
        ax.set_yticks([-4000, -2000, 0, 2000, 4000])
        ax.set_yticklabels(["4,000", "2,000", "0", "2,000", "4,000"])

    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, 3, figsize=(3.7, 2.2), sharey=True,
                                 gridspec_kw={"wspace": 0.14, "width_ratios": [4, 4, 3]})
        for ax, drug in zip(axes, DRUGS):
            bars(ax, GENOS, drug, f"{TLAB[drug]} vs DMSO")
        bars(axes[2], VARS, "Control_vs_NT", "DMSO, vs NT")
        axes[0].set_ylabel("DEGs (padj < 0.05, |log$_2$FC| > 1)\n← down    up →")
        for ax in axes[1:]:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", length=0)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_3f_DEG_counts.pdf"), transparent=True,
                    bbox_inches="tight")
    return fig, d


PADJ, LFC = 0.05, 1.0


def core_nt_signature():
    """Genes the NT drug response moves the same way under *both* drugs.

    03_figures.py line 34: `nt_sig` per drug, then intersect. The panel title
    carries these two counts, so they are recomputed here from the shipped
    contrast tables rather than hard-coded.
    """
    def sig(drug):
        r = io.read(RESULTS / "contrasts" / f"NT_{drug}_vs_NT_Control.csv",
                    index_col="gene_id")
        s = r[r.padj < PADJ]
        return (set(s[s.log2FC_shrunk > LFC].index),
                set(s[s.log2FC_shrunk < -LFC].index))

    (hu, hd), (zu, zd) = sig("HQ461"), sig("ZLY025")
    return sorted(hu & zu), sorted(hd & zd)


def figure_3g(outdir=None):
    """Fraction of the NT drug response each variant retains.

    Bars are the replicate mean, error bars SD (n = 3), points the replicates.
    The jitter uses the scripts' own seed so the points land where they do in
    the published panel.
    """
    pct = rnaseq.figure_3g()["drug_response_retained"]
    rng = np.random.default_rng(0)

    with plt.rc_context(RC):
        fig, ax = plt.subplots(figsize=(2.3, 2.2))
        for j, drug in enumerate(DRUGS):
            for i, g in enumerate(GENOS):
                x = j * 5 + i
                v = pct[(pct.genotype == g) & (pct.drug == drug)].pct_response
                ax.bar(x, v.mean(), width=0.78, color=GCOL[g], **BAR_EDGE)
                ax.errorbar(x, v.mean(), yerr=v.std(ddof=1), color=INK, lw=0.7,
                            capsize=1.5, capthick=0.7)
                ax.scatter(x + rng.uniform(-0.15, 0.15, len(v)), v, s=6, color="white",
                           edgecolor=INK, lw=0.4, zorder=3)
                ax.text(x, max(v.max(), v.mean()) + 4, f"{v.mean():.0f}",
                        ha="center", va="bottom", fontsize=5.5)
        ax.set_xticks([j * 5 + i for j in range(2) for i in range(4)])
        ax.set_xticklabels(GENOS * 2, rotation=90)
        for j, drug in enumerate(DRUGS):
            ax.text(j * 5 + 1.5, 1.02, TLAB[drug], transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=7)
        # 100% is the NT anchor, so ticks stop there; the limit runs a little
        # further only to fit NT's own replicate scatter rather than clip it.
        ax.set_ylim(min(-6, pct.pct_response.min() - 6), 112)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.axhline(0, color=INK, lw=0.6)
        # 100% is the NT anchor the percentages are expressed against
        ax.axhline(100, color=MUTED, lw=0.6, ls="--")
        ax.set_ylabel("Drug response retained\n(% of NT signature shift)")
        up, dn = core_nt_signature()
        ax.set_title(f"Core NT signature\n{len(up):,} up / {len(dn):,} down genes",
                     fontsize=6.5, pad=14)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_3g_response_retained.pdf"), transparent=True,
                    bbox_inches="tight")
    return fig, pct


def figure_3h(outdir=None):
    """Elongation profiles across gene-length deciles.

    Error bars are mean +/- SD across the three biological replicates, not the
    gene-to-gene IQR. The third panel keeps its own, much tighter y-scale: the
    no-drug variant-vs-NT differences are an order of magnitude smaller, and
    sharing the drug axis would flatten them to a line.
    """
    d = rnaseq.figure_3h()
    DRUG_YLIM, BASE_YLIM = (-0.45, 0.85), (-0.12, 0.16)

    def profile(ax, sub, g):
        sub = sub.sort_values("decile")
        ax.errorbar(sub.span_median_bp, sub.mean_log2ratio, yerr=sub.sd,
                    color=GCOL[g], lw=1.3, marker="o", ms=3.0, mec=INK, mew=0.5,
                    capsize=1.6, capthick=0.6, elinewidth=0.7, label=g,
                    zorder=3, clip_on=False)

    with plt.rc_context(RC):
        fig, axes = plt.subplots(1, 3, figsize=(6.7, 2.05),
                                 gridspec_kw={"width_ratios": [1, 1, 1], "wspace": 0.3})
        for c_i, drug in enumerate(DRUGS):
            ax = axes[c_i]
            for g in GENOS:
                sub = d[(d.group == "drug") & (d.contrast == f"{g}_{drug}")]
                if len(sub):
                    profile(ax, sub, g)
            ax.set_ylim(*DRUG_YLIM)
            # The original left the locator to matplotlib; pinning the ticks
            # keeps the panel looking like the published one across versions.
            ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
            ax.set_title(f"{TLAB[drug]} vs DMSO\n ", fontsize=6.5)
            if c_i == 0:
                ax.set_ylabel("median log$_2$ ratio per decile\n(mean ± SD, n = 3)")
            else:
                ax.tick_params(labelleft=False)

        ax = axes[2]
        for g in VARS:
            sub = d[(d.group == "baseline") & (d.contrast == f"{g}_Control")]
            if len(sub):
                profile(ax, sub, g)
        ax.set_ylim(*BASE_YLIM)
        ax.set_yticks([-0.10, -0.05, 0.00, 0.05, 0.10, 0.15])
        ax.set_title("No drug: variant vs NT\n(note expanded y-scale)", fontsize=6.5)
        ax.set_ylabel("median log$_2$ ratio per decile", fontsize=6.5)

        for ax in axes:
            ax.axhline(0, color=MUTED, lw=0.6, ls="--", zorder=1)
            ax.set_xscale("log")
            ax.set_xlabel("Gene length (bp)")
            ax.set_axisbelow(True)
        axes[0].legend(loc="upper right", handlelength=1.2, fontsize=5.5,
                       ncol=2, columnspacing=0.8)
        axes[2].legend(loc="lower left", handlelength=1.2, fontsize=5.5)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_3h_elongation_deciles.pdf"), transparent=True,
                    bbox_inches="tight")
    return fig, d


FIGURES = {"3f": figure_3f, "3g": figure_3g, "3h": figure_3h}
