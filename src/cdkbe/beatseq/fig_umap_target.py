"""KEY 20-21 — cell-level UMAPs coloured by the CDK each cell's gRNA targets,
before and after the cell cycle is regressed out.

Run `python -m analysis.umap_target` first; this module only draws.

WHAT THE TWO EMBEDDINGS ARE. Identical pipeline, identical 2,000 HVGs, identical
seed; they differ in exactly two regression terms (S_score, G2M_score). Both
already have a per-gem-group intercept and log-depth removed, so the raw panels
are not confounded by the batch structure either. See `analysis/umap_target.py`.

WHY THE COLOURS ARE A PROBLEM, AND WHAT IS DONE ABOUT IT. The parent manuscript's
nine per-CDK colours (`style.GENE_COLOR`) were chosen to label nine well-separated
structures in Figures 1-7, not nine interleaved point clouds. Three of them are
blues and two are greens. Run through a CVD/contrast validator as a nine-slot
categorical palette it fails four checks outright -- CDK6-orange against
CDK2-goldenrod is dE 10 even to normal colour vision (15 is the floor), and
CDK6-orange against CDK12-green collapses to dE 0.7 under protanopia. Nine-way
colour is therefore NOT a dependable identity channel in an overplotted scatter,
whatever palette is used.

The palette is kept anyway, because keying to Figures 1-7 matters more than
separating nine categories the reader should not be asked to separate by eye. The
fix is structural rather than chromatic:

    KEY 20  the nine-colour overlay, as asked for, with the cell-cycle and
            editor-arm panels beside it so the reader can see what the axes are
            actually made of
    KEY 21  small multiples -- one CDK per panel against the full cloud in grey.
            One series per panel, named in the panel title, so nine-way colour
            discrimination is never on the reader's critical path. This is the
            panel to read if you want to know where a given CDK's cells sit.

Cell-cycle phase gets a single-hue ordinal ramp rather than three categorical
hues, because G1 -> S -> G2M is ordered; the three steps pass the ordinal
lightness-monotonicity and light-end contrast checks. The editor arms keep the
manuscript's ABE-purple / CBE-blue, which fails the same validator on its own
(dE 1.7 under protanopia), so those panels carry the same-arm mixing RATIO in
writing: the claim there is "the two arms are not separated", and it is settled by
the number, not by telling the two hues apart.

THE DETACHED MINORITY CLOUD IN T47D is not a batch: 795 of 17,339 cells, split
583/212 between the ABE and CBE arms -- the same 71/29 proportion as the whole
line -- and almost exactly Leiden cluster 9 (660 cells) plus part of cluster 7
from the delivered pipeline. It carries every target in proportion to its share of
the population, NT included (12.3% there vs 10.5% elsewhere), and has LOWER
mitochondrial content (1.45% vs 2.83%), so it is a pre-existing cell state that
every gRNA samples, not a perturbation response.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from .. import io as _io  # noqa: E402
from . import paths, style  # noqa: E402

# Text tokens -- labels never wear a data colour.
INK, MUTED = "#202124", "#5F6368"
# Cells whose gRNA targets nothing, and the grey backdrop in the small multiples.
NT_GREY, BACKDROP = "#B9BDC1", "#E4E6E8"
# NT needs two colours, because it plays two roles. In the ten-way overlay it
# should recede (NT_GREY). When NT is the HIGHLIGHTED series in a small multiple it
# sits on the BACKDROP grey, where NT_GREY is nearly invisible -- 1.2:1 against it
# -- so the mark colour is a dark neutral instead: still obviously "not one of the
# nine CDK hues", but 4.4:1 on the backdrop and 5.4:1 on the surface.
NT_MARK = "#5F6368"
PHASE_COLOR = {"G1": "#6BAED6", "S": "#2171B5", "G2M": "#08306B"}
PHASE_ORDER = ["G1", "S", "G2M"]
TARGET_ORDER = style.CDK_ORDER + ["NT"]

EMBEDDINGS = [("raw", "cell cycle retained"), ("cc", "cell cycle regressed out")]
# What each column of KEY 20 is coloured by, and the statistic printed under it.
VIEWS = [
    ("phase", "cell-cycle phase", lambda m: f"phase = {m.phase_var_frac:.0%} of PC variance"),
    ("arm", "editor arm", lambda m: f"same-arm neighbours {m.arm_mixing_ratio:.2f}$\\times$ chance"),
    ("guide_target", "CDK targeted", lambda m: f"own-target neighbours {m.mixing_ratio:.2f}$\\times$ chance"),
]


def target_color(gene: str) -> str:
    """Colour for the ten-way overlay, where NT should recede."""
    return NT_GREY if gene == "NT" else style.GENE_COLOR[gene]


def mark_color(gene: str) -> str:
    """Colour for a small multiple, where the highlighted series must be visible."""
    return NT_MARK if gene == "NT" else style.GENE_COLOR[gene]


def arm_color(arm: str) -> str:
    return style.EDITOR_COLOR[arm.split("/")[1]]


def load(cell_line: str, with_gRNA_id: bool = False) -> pd.DataFrame:
    """Per-cell embeddings; `gRNA_id` back-filled from the singlet object on request.

    The delivered `umap_cells_*.csv.gz` carries `guide_identity` (allele|protospacer)
    but not `gRNA_id`. Back-filling it means opening the singlet `.h5ad`, which is one
    of the ~34 GB of objects deliberately not copied into this repository
    (`data/external/NOT_INCLUDED.csv`), so the back-fill is OFF by default and the
    caller has to ask. Nothing that draws a figure asks: every panel here keys on
    `guide_identity`, `guide_target`, `arm` or `phase`. The join is on the cell
    barcode index, so it cannot silently mismatch.
    """
    f = _io.resolve(paths.RESULTS / f"umap_cells_{cell_line}.csv.gz")
    d = pd.read_csv(f, index_col=0)
    if with_gRNA_id and "gRNA_id" not in d.columns:
        import anndata as ad
        obs = ad.read_h5ad(paths.SINGLETS[cell_line], backed="r").obs
        d["gRNA_id"] = obs["gRNA_id"].astype(str).reindex(d.index)
        assert d["gRNA_id"].notna().all(), "gRNA_id back-fill left gaps"
    return d


def natural_limits(d: pd.DataFrame, tag: str, pad: float = 0.04):
    """The embedding's own bounding box, padded.

    Used instead of a square window because a UMAP cloud is rarely square, and
    T47D's detached minority cloud stretches the x-range: forcing a square window
    then pays for that stretch in empty vertical space in every panel. `frame`
    keeps the aspect ratio honest by matching the axes box to the window, so
    distances are still undistorted.
    """
    x, y = d[f"umap1_{tag}"].to_numpy(), d[f"umap2_{tag}"].to_numpy()
    dx, dy = np.ptp(x), np.ptp(y)
    return ((x.min() - pad * dx, x.max() + pad * dx),
            (y.min() - pad * dy, y.max() + pad * dy))


def window_aspect(xlim, ylim) -> float:
    """Height/width of a data window — the shape every panel drawn on it takes."""
    return (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])


def plan_grid(fig_w: float, plot_w: float, aspect: float, nrow: int, *,
              head_in: float, foot_in: float, hgap_in: float):
    """Figure height and gridspec top/bottom/hspace for a grid of fixed-shape panels.

    Laying the grid out in inches and converting once is what stops a wide, short
    UMAP window from leaving an inch of dead space under every row: the row pitch
    is derived from the panel shape instead of being whatever is left over after
    the figure height was guessed.
    """
    box_h = plot_w * aspect
    fig_h = head_in + nrow * box_h + (nrow - 1) * hgap_in + foot_in
    return (fig_h, 1 - head_in / fig_h, foot_in / fig_h, hgap_in / box_h)


def frame(ax, xlim, ylim) -> None:
    """Strip the box and ticks, and match the axes box to the data window.

    A UMAP has no meaningful units to tick. Setting `box_aspect` from the window
    rather than calling `set_aspect("equal")` keeps x and y on the same scale
    without letting matplotlib resize the axes box to suit the data, which is
    what makes every panel in a grid come out the same physical size.
    """
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_box_aspect((ylim[1] - ylim[0]) / (xlim[1] - xlim[0]))


def axis_key(ax, frac: float = 0.16, lw: float = 0.9) -> None:
    """Small L-shaped orientation key, in place of a full axis frame."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    ox, oy = x0 + 0.015 * (x1 - x0), y0 + 0.015 * (y1 - y0)
    lx, ly = frac * (x1 - x0), frac * (y1 - y0)
    ax.plot([ox, ox, ox + lx], [oy + ly, oy, oy], color=MUTED, lw=lw,
            solid_capstyle="round", zorder=5)
    ax.text(ox + lx * 0.55, oy - 0.005 * (y1 - y0), "UMAP 1", ha="center", va="top",
            fontsize=6.0, color=MUTED)
    ax.text(ox - 0.012 * (x1 - x0), oy + ly * 0.55, "UMAP 2", ha="right", va="center",
            rotation=90, fontsize=6.0, color=MUTED)


def scatter(ax, d, tag, color_by, size, alpha, seed=0):
    """One rasterised scatter call, with point order randomised.

    Drawing category by category would put whichever category is drawn last on top
    of every other, which in an overplotted embedding is the difference between
    'CDK4 is everywhere' and 'CDK4 is nowhere'. Shuffling makes occlusion
    unbiased, so the visible colour mix at any spot is the real local mix.
    Rasterising the points keeps the PDF small while the text stays live vector.
    """
    rng = np.random.default_rng(seed)
    o = rng.permutation(len(d))
    cmap = {"guide_target": target_color, "phase": PHASE_COLOR.get,
            "arm": arm_color}[color_by]
    ax.scatter(d[f"umap1_{tag}"].to_numpy()[o], d[f"umap2_{tag}"].to_numpy()[o],
               s=size, c=[cmap(v) for v in d[color_by].to_numpy()[o]],
               alpha=alpha, linewidths=0, rasterized=True)


def note(ax, text: str) -> None:
    """Metric line just below the axes, right-aligned.

    Below rather than inside: these embeddings are wide and short, so there is no
    interior corner a caption can sit in without landing on cells.
    """
    ax.text(1.0, -0.045, text, transform=ax.transAxes, ha="right", va="top",
            fontsize=6.6, color=MUTED)


def legend_panel(ax, entries, title: str) -> None:
    """A legend alone on a blank axes, so it never displaces a plot."""
    ax.axis("off")
    h = [Line2D([], [], marker="o", ls="", markersize=5.4, markerfacecolor=c,
                markeredgecolor="none", label=l) for l, c in entries]
    lg = ax.legend(handles=h, loc="upper left", bbox_to_anchor=(-0.06, 1.0),
                   title=title, handletextpad=0.3, borderpad=0.0, labelspacing=0.5)
    lg.get_title().set_color(MUTED)
    lg.get_title().set_fontsize(mpl.rcParams["legend.fontsize"])
    for t in lg.get_texts():
        t.set_color(INK)


def detached_cloud(d: pd.DataFrame, tag: str):
    """Centroid and size of the minority cloud left of the embedding's midpoint.

    Found from the data rather than hard-coded, so the annotation still points at
    it if the embedding is recomputed and lands the cloud somewhere else.
    """
    x = d[f"umap1_{tag}"].to_numpy()
    y = d[f"umap2_{tag}"].to_numpy()
    s = x < (x.min() + x.max()) / 2
    if not s.any():
        return None
    return float(x[s].mean()), float(y[s].mean()), int(s.sum())


def point_spec(n: int):
    """Marker size / opacity for an overplotted cloud of `n` cells."""
    return (2.0, 0.40) if n > 12_000 else (3.2, 0.55)


# ---------------------------------------------------------------- KEY 20
def fig_overlay(data: dict, mix: pd.DataFrame) -> None:
    """Rows = cell line; left three columns = raw, right three = cell-cycle regressed.

    The two halves are the same cells through the same pipeline, so reading across
    a row shows what removing the cell cycle does to the embedding, and reading
    down a column shows whether the two lines agree.
    """
    style.apply(8.6)
    nv = len(VIEWS)
    nrow = len(paths.CELL_LINES)
    fig_w, left, right, leg = 15.8, 0.026, 0.995, 0.42
    plot_w = (right - left) * fig_w / (2 * nv + leg)
    aspect = max(window_aspect(*natural_limits(data[cl], t))
                 for cl in paths.CELL_LINES for t, _ in EMBEDDINGS)
    fig_h, gtop, gbot, hsp = plan_grid(fig_w, plot_w, aspect, nrow,
                                       head_in=1.55, foot_in=0.72, hgap_in=1.12)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(nrow, 2 * nv + 1,
                          width_ratios=[1] * (2 * nv) + [leg],
                          hspace=hsp, wspace=0.08,
                          left=left, right=right, top=gtop, bottom=gbot)
    letters = iter("abcdefghijkl")
    half_x: dict[int, tuple[float, float]] = {}

    for i, cl in enumerate(paths.CELL_LINES):
        d = data[cl]
        size, alpha = point_spec(len(d))
        for j, (tag, _) in enumerate(EMBEDDINGS):
            m = mix[(mix.cell_line == cl) & (mix.embedding == tag)].iloc[0]
            xlim, ylim = natural_limits(d, tag)
            for k, (cby, what, stat) in enumerate(VIEWS):
                ax = fig.add_subplot(gs[i, j * nv + k])
                scatter(ax, d, tag, cby, size, alpha)
                frame(ax, xlim, ylim)
                if j == 0 and k == 0:
                    axis_key(ax)
                style.panel(ax, next(letters), f"{cl} · {what}", pad=5.0, size=8.7)
                note(ax, stat(m))
                if i == 0:
                    pos = ax.get_position()
                    lo, hi = half_x.get(j, (pos.x0, pos.x1))
                    half_x[j] = (min(lo, pos.x0), max(hi, pos.x1))
                # The detached T47D cloud is the first thing a reader asks about;
                # answer it on the canvas rather than only in the legend.
                if cl == "T47D" and k == 0:
                    cloud = detached_cloud(d, tag)
                    if cloud:
                        cxx, cyy, cn = cloud
                        # A line below `note`, so the two captions never share a
                        # baseline; the arrow rises on the far left, clear of it.
                        ax.annotate(f"Leiden 9 cell state, {cn / len(d):.1%}",
                                    xy=(cxx, cyy), xytext=(0.0, -0.175),
                                    textcoords="axes fraction",
                                    fontsize=6.2, color=MUTED, ha="left", va="top",
                                    arrowprops=dict(arrowstyle="-", color=MUTED,
                                                    lw=0.6, shrinkA=2, shrinkB=3))

        # Legends live in their own column so they never displace a plot.
        legend_panel(fig.add_subplot(gs[i, 2 * nv]),
                     [(g, target_color(g)) for g in TARGET_ORDER] if i == 0 else
                     [(p_, PHASE_COLOR[p_]) for p_ in PHASE_ORDER]
                     + [(a, style.EDITOR_COLOR[a]) for a in paths.EDITORS],
                     "gRNA target" if i == 0 else "phase / editor arm")

    # Block headers spanning each half, drawn once rather than per panel, and a
    # hairline rule so the raw / regressed split is visible at a glance.
    for j, (_, tag_label) in enumerate(EMBEDDINGS):
        x0, x1 = half_x[j]
        fig.text(x0, 1 - 1.22 / fig_h, tag_label.upper(), ha="left", va="bottom",
                 fontsize=8.0, color=MUTED, weight="bold")
        fig.add_artist(plt.Line2D([x0, x1], [1 - 1.30 / fig_h] * 2, color="#D5D7DA",
                                  lw=0.8, transform=fig.transFigure))
    gap = (half_x[1][0] + half_x[0][1]) / 2
    fig.add_artist(plt.Line2D([gap, gap], [gbot * 0.55, 1 - 1.36 / fig_h],
                              color="#E8EAEC", lw=0.8, transform=fig.transFigure))

    n_txt = " · ".join(
        f"{cl} {int(mix[mix.cell_line == cl].n_cells.iloc[0]):,}"
        f"/{int(mix[mix.cell_line == cl].n_cells_total.iloc[0]):,} cells"
        for cl in paths.CELL_LINES)
    fig.text(left, 1 - 0.32 / fig_h,
             "Cells organise by cell-cycle phase, not by which CDK their gRNA targets",
             ha="left", va="top", fontsize=11.8, color=INK, weight="bold")
    fig.text(left, 1 - 0.62 / fig_h,
             "Both halves are the same cells, the same 2,000 HVGs and the same seed, "
             "with a per-gem-group intercept and log-depth regressed out. The right "
             "half adds S and G2M scores to that regression, so the cell cycle is the "
             "only difference between the halves.\n"
             "Phase falls from ~22% of leading-PC variance to ~1%, the editor arms "
             "were never separated to begin with (1.01–1.12$\\times$ chance), and "
             "target identity stays at chance throughout (1.02–1.03$\\times$). "
             + n_txt + ".",
             ha="left", va="top", fontsize=7.9, color=MUTED, linespacing=1.55)
    style.save(fig, "KEY20_umap_target")


# ---------------------------------------------------------------- KEY 21
def fig_facets(cell_line: str, d: pd.DataFrame, mix: pd.DataFrame) -> None:
    """Small multiples: one target per panel over the full cloud in grey."""
    style.apply(8.2)
    bg_s, fg_s = (1.4, 3.0) if len(d) > 12_000 else (2.0, 4.2)
    ncol, nrow = 5, 2
    fig_w, left, right = 12.0, 0.022, 0.995
    # Geometry in inches first: the panel shape follows the embedding's own aspect,
    # so the figure height has to follow the panels rather than the reverse.
    plot_w = (right - left) * fig_w / ncol
    box_h = max(plot_w * max(window_aspect(*natural_limits(d, t)) for t, _ in EMBEDDINGS),
                1.15)
    head_in, hgap_in, blockgap_in, foot_in = 1.52, 0.40, 0.80, 0.12
    block_in = nrow * box_h + (nrow - 1) * hgap_in
    fig_h = head_in + 2 * block_in + blockgap_in + foot_in
    blocks = [((fig_h - head_in - i * (block_in + blockgap_in)) / fig_h,
               (fig_h - head_in - i * (block_in + blockgap_in) - block_in) / fig_h)
              for i in range(2)]
    fig = plt.figure(figsize=(fig_w, fig_h))

    for b, (tag, tag_label) in enumerate(EMBEDDINGS):
        m = mix[(mix.cell_line == cell_line) & (mix.embedding == tag)].iloc[0]
        top, bottom = blocks[b]
        gs = fig.add_gridspec(nrow, ncol, top=top, bottom=bottom,
                              left=left, right=right, hspace=hgap_in / box_h,
                              wspace=0.07)
        xlim, ylim = natural_limits(d, tag)
        x, y = d[f"umap1_{tag}"].to_numpy(), d[f"umap2_{tag}"].to_numpy()
        for t, gene in enumerate(TARGET_ORDER):
            ax = fig.add_subplot(gs[t // ncol, t % ncol])
            s = (d.guide_target == gene).to_numpy()
            ax.scatter(x[~s], y[~s], s=bg_s, color=BACKDROP, linewidths=0,
                       rasterized=True, zorder=1)
            ax.scatter(x[s], y[s], s=fg_s, color=mark_color(gene), alpha=0.85,
                       linewidths=0, rasterized=True, zorder=2)
            frame(ax, xlim, ylim)
            # Identity above the axes: a colour dot beside ink-coloured text, so a
            # pale hue like CDK7's cyan never has to be legible as type.
            ax.plot([0.028], [1.055], marker="o", ms=4.4, ls="", clip_on=False,
                    color=mark_color(gene), transform=ax.transAxes)
            ax.text(0.075, 1.052, gene, transform=ax.transAxes, ha="left",
                    va="center", fontsize=8.8, color=INK, weight="bold")
            ax.text(0.99, 1.052, f"{int(s.sum()):,} cells", transform=ax.transAxes,
                    ha="right", va="center", fontsize=6.8, color=MUTED)
            if t == 0:
                axis_key(ax)

        fig.text(left, top + 0.30 / fig_h,
                 f"$\\bf{{{'ab'[b]}}}$   {cell_line}, {tag_label} — own-target "
                 f"neighbours {m.mixing_ratio:.2f}$\\times$ chance, same-guide "
                 f"{m.guide_mixing_ratio:.2f}$\\times$",
                 ha="left", va="bottom", fontsize=9.6, color=INK)

    n = int(mix[mix.cell_line == cell_line].n_cells.iloc[0])
    n_tot = int(mix[mix.cell_line == cell_line].n_cells_total.iloc[0])
    fig.text(left, 1 - 0.30 / fig_h,
             f"{cell_line}: every CDK's cells fill the whole embedding, with or "
             f"without the cell cycle",
             ha="left", va="top", fontsize=11.6, color=INK, weight="bold")
    fig.text(left, 1 - 0.60 / fig_h,
             "Each panel highlights one gRNA target over all other cells in grey. "
             "Colours are the parent manuscript's per-CDK colours, but identity here "
             f"comes from the panel, not the hue. {n:,} of {n_tot:,} cells.",
             ha="left", va="top", fontsize=7.9, color=MUTED)
    style.save(fig, f"KEY21_umap_facets_{cell_line}")


def main() -> None:
    mix = pd.read_csv(paths.RESULTS / "umap_mixing.csv")
    data = {cl: load(cl) for cl in paths.CELL_LINES}
    fig_overlay(data, mix)
    for cl in paths.CELL_LINES:
        fig_facets(cl, data[cl], mix)


if __name__ == "__main__":
    main()
