"""Main-figure panels, each written as its OWN file for assembly in Illustrator.

Everything lands in `analysis/figures/panels/` as `PANEL_<letter><n>_<slug>.{pdf,svg,png}`.
Nothing is composited here: each panel is a self-contained figure sized to be placed and
resized. PDF and SVG both carry live text (`pdf.fonttype = 42`, `svg.fonttype = "none"`),
point clouds are rasterised at 400 dpi and nothing else is, so a panel opens in Illustrator
as editable type over a single embedded raster.

    a  UMAPs of A549 and T47D, coloured by CDK targeted and by cell-cycle phase
    b  per-gRNA cell-cycle arrest in T47D — the nonsense and missense alleles that move G1
    c  editing evidence from the target transcript, and its relation to the arrest
    d  transcriptional severity per CDK (KEY 15c-d, redrawn one panel per cell line)
    e  severity against the parent screen's maximum drug log2FC

Every number is read from `analysis/results/`; nothing is hardcoded. MAIN_PANELS.md gives
the panel-by-panel statistics and the caveats each panel has to be shown with.

Run:
    python -m analysis.mainpanel_stats      # writes the tables this module needs
    python -m analysis.fig_panels           # all panels
    python -m analysis.fig_panels b c       # only the b and c blocks
"""

from __future__ import annotations

import contextlib
import re
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from .. import io as _io  # noqa: E402
from . import paths, style  # noqa: E402
from .fig_umap_target import (INK, MUTED, PHASE_COLOR, PHASE_ORDER,  # noqa: E402
                                      axis_key, frame, load, natural_limits,
                                      point_spec, scatter, target_color,
                                      window_aspect)

class _Results:
    """`R / "x.csv"` finds whichever spelling is on disk -- `x.csv` or `x.csv.gz`.

    The delivered tables were copied into this repository with the larger ones
    gzipped, so a literal `"x.csv"` here would miss `x.csv.gz`. Resolution happens
    at the one place the path is built rather than at the ~30 read sites. A name
    that resolves to nothing comes back unchanged, so the `(R / f"{n}.csv").exists()`
    guards below still read False instead of raising.
    """

    def __truediv__(self, rel):
        q = paths.RESULTS / rel
        try:
            return _io.resolve(q)
        except FileNotFoundError:
            return q


R = _Results()
OUT = paths.FIGURES / "panels"
GC = style.GENE_COLOR

# --- compact mode -----------------------------------------------------------------
# `python -m analysis.fig_panels clean` writes a second set to `panels_clean/` with the
# methods footnotes and the descriptive panel titles stripped and the canvas scaled down.
# Type stays at its original point size, so panels placed at 100% in Illustrator keep a
# consistent type size while the plot area shrinks. FIG_SCALE defaults to 1.0 because
# type is in points and does not scale with the canvas, so shrinking it crowds tick labels;
# most of the compaction comes from dropping prose, not from scaling. Pass `clean=0.85` to
# scale anyway. Everything that carries a NUMBER --
# p-values, effect sizes, counts, axis labels, legends, threshold labels -- is kept; only
# prose is removed. The full-caption set in `panels/` is the one to read the methods from.
CLEAN = False
FIG_SCALE = 1.0


def _sz(w: float, h: float) -> tuple:
    return (w * FIG_SCALE, h * FIG_SCALE) if CLEAN else (w, h)

# Variant class is the one categorical the b and c blocks turn on. Editor arm is carried
# by marker shape there, so these hues cannot collide with ABE-purple / CBE-blue.
VCLASS_COLOR = {
    "nonsense": "tab:red",
    "nonsense_competent": "tab:red",
    "nonsense_escaping": "#E8A33D",
    "start_loss": "tab:orange",
    "missense": "#35618F",
    "non-targeting": "#9AA0A6",
}
VCLASS_LABEL = {
    "nonsense": "nonsense (PTC)",
    "nonsense_competent": "nonsense, NMD-competent",
    "nonsense_escaping": "nonsense, NMD-escaping",
    "start_loss": "start loss",
    "missense": "missense",
    "non-targeting": "non-targeting",
}
ARM_MARKER = {"ABE": "o", "CBE": "s"}

# Cell-cycle phase: three DISTINCT hues rather than the single-hue blue ramp that
# `fig_umap_target.PHASE_COLOR` uses (KEY 20/24 keep that one). Phase is ordinal, so the
# assignment is by OKLCH lightness — G1 0.671 > S 0.575 > G2M 0.433 — and the cycle order
# still reads light-to-dark while the hues stay separable. Validated as a CATEGORICAL
# palette under all-pairs, which is the test a scatter needs (any two marks can be
# neighbours), with `dataviz/scripts/validate_palette.js`:
#   lightness band PASS · chroma floor PASS · contrast vs surface PASS (all >= 3:1)
#   CVD separation  worst all-pairs dE 13.0 deuteranopia / 17.4 tritanopia  (target >= 8)
#   normal vision   worst all-pairs dE 16.3                                (hard gate 15)
# It was the only trio of the 56 from the documented palette that cleared the contrast
# check on all three slots AND kept CVD separation out of the 6-8 floor band.
PHASE_COLOR_CB = {"G1": "#eb6834", "S": "#2a78d6", "G2M": "#4a3aa7"}


_LETTER = re.compile(r"^(\$\\bf\{[^}]*\}\$)")


def _strip_titles(fig) -> None:
    """Keep the bold panel letter, drop the sentence after it.

    Titles are set two ways in this module -- `style.panel` (an axes left-title) and a
    `fig.text` for panels whose sub-axes carry their own titles. Both are handled here
    rather than at the ~20 call sites, so the two output sets cannot drift apart.
    """
    for ax in fig.get_axes():
        t = getattr(ax, "_left_title", None)
        if t is not None and t.get_text():
            m = _LETTER.match(t.get_text())
            t.set_text(m.group(1) if m else "")
    # A figure-level title sits above the whole canvas; stripping it to a bare letter
    # would leave the gap it had reserved for a sentence. Move the letter onto the first
    # axes instead, where it sits flush, and drop the figure text.
    axs = fig.get_axes()
    for t in list(fig.texts):
        m = _LETTER.match(t.get_text())
        if not m:
            continue
        t.set_text("")
        if axs:
            a = axs[0]
            cur = getattr(a, "_left_title", None)
            if cur is None or not cur.get_text():
                a.set_title(m.group(1), loc="left", pad=4,
                            fontsize=t.get_fontsize())


#: set by `capture()`; while it is a dict, `save_panel` hands the figure over
#: instead of writing and closing it.
_CAPTURE: dict | None = None


@contextlib.contextmanager
def capture(clean: bool = True):
    """Collect the figures a panel block draws, instead of writing them to disk.

        with fig_panels.capture() as panels:
            fig_panels.panel_d()
        fig = panels["PANEL_d12w_severity_by_CDK_stacked"]

    The point is that a caller that wants ONE panel gets the figure the real
    pipeline drew, not a reimplementation of it: nothing here is duplicated, so
    the captured figure and the file in `figures/panels_clean/` cannot drift.
    `clean` selects the compact variant (titles stripped to the bare letter),
    which is the set the manuscript was assembled from.

    Blocks draw several panels each, so capturing one panel still draws its
    whole block. That costs a few seconds and buys exact fidelity.
    """
    global _CAPTURE, CLEAN
    prev_cap, prev_clean = _CAPTURE, CLEAN
    _CAPTURE, CLEAN = {}, clean
    try:
        yield _CAPTURE
    finally:
        _CAPTURE, CLEAN = prev_cap, prev_clean


def save_panel(fig, stem: str, dpi: int = 400) -> None:
    """PDF + SVG (live text, for Illustrator) + PNG (for looking at)."""
    if CLEAN:
        _strip_titles(fig)
    if _CAPTURE is not None:
        # handed to the caller, still open: closing it here would leave them a
        # figure matplotlib has already torn the canvas off.
        _CAPTURE[stem] = fig
        return
    out = OUT.parent / "panels_clean" if CLEAN else OUT
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=dpi, transparent=True,
                    bbox_inches="tight", pad_inches=0.01 if CLEAN else 0.02)
    plt.close(fig)
    print(f"  wrote {out.name}/{stem}.pdf / .svg / .png")


# ------------------------------------------------------------------ label placement
def _offset_ring(radii=(6, 11, 17, 24, 32, 42, 54), n_ang=16):
    """Candidate label offsets, in points, ordered near-to-far and top-first.

    A fixed hand-written list runs out in a dense cloud and the fallback then drops a
    label on top of another one. Sweeping rings guarantees a free slot exists unless the
    axes really are full, and ordering by radius keeps labels close to their point.
    """
    out = []
    for r in radii:
        ang = np.linspace(0, 2 * np.pi, n_ang, endpoint=False) + np.pi / 2
        pts = [(float(r * np.cos(a)), float(r * np.sin(a))) for a in ang]
        out += sorted(pts, key=lambda t: -t[1])          # prefer above the point
    return out


_OFFSETS = _offset_ring()


def place_labels(ax, xs, ys, texts, colors, fontsize=5.4, obstacles=None,
                 leader_colour="#B8BCBF"):
    """Annotate points, trying offsets until the label misses every other label/point.

    Matplotlib has no collision handling, and these panels label 8-15 points inside a
    dense cloud, so hand-tuned offsets stop being maintainable the moment a number
    changes. `obstacles` are display-space boxes (built from the plotted points) that a
    label must also avoid; a leader line is drawn whenever the label ends up far enough
    from its point that the association would otherwise be ambiguous.
    """
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    taken = list(obstacles or [])
    axbb = ax.get_window_extent(rend)
    for x, y, t, c in zip(xs, ys, texts, colors):
        chosen = None
        for dx, dy in _OFFSETS:
            ha = "left" if dx > 0 else ("right" if dx < 0 else "center")
            va = "bottom" if dy > 0 else ("top" if dy < 0 else "center")
            tx = ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points",
                             fontsize=fontsize, color=c, ha=ha, va=va, zorder=6)
            bb = tx.get_window_extent(rend).expanded(1.06, 1.22)
            if not any(bb.overlaps(o) for o in taken) and axbb.contains(
                    bb.x0, bb.y0) and axbb.contains(bb.x1, bb.y1):
                chosen = (tx, bb, dx, dy)
                break
            tx.remove()
        if chosen is None:
            # every candidate collided; place it anyway and say so, rather than
            # silently drawing a label on top of another one
            print(f"    note: no free slot for label {t!r}; placed at the default offset")
            tx = ax.annotate(t, (x, y), xytext=(6, 5), textcoords="offset points",
                             fontsize=fontsize, color=c, ha="left", va="bottom", zorder=6)
            chosen = (tx, tx.get_window_extent(rend), 6, 5)
        tx, bb, dx, dy = chosen
        if abs(dx) + abs(dy) > 14:
            # a distant label needs a leader; arrowprops is set at construction, so the
            # measured placement is rebuilt rather than mutated
            ha, va = tx.get_ha(), tx.get_va()
            tx.remove()
            tx = ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points",
                             fontsize=fontsize, color=c, ha=ha, va=va, zorder=6,
                             arrowprops=dict(arrowstyle="-", lw=0.35,
                                             color=leader_colour, shrinkA=0.5,
                                             shrinkB=1.5))
        taken.append(bb)


def point_boxes(ax, xs, ys, half=2.5):
    """Display-space boxes around plotted points, for `place_labels` to avoid."""
    from matplotlib.transforms import Bbox
    pts = ax.transData.transform(np.column_stack([np.asarray(xs, float),
                                                  np.asarray(ys, float)]))
    return [Bbox.from_extents(px - half, py - half, px + half, py + half)
            for px, py in pts]


def vlegend(ax, classes, loc="lower right", title=None, fontsize=5.6, **kw):
    h = [Line2D([], [], marker="o", ls="", ms=4.4, mfc=VCLASS_COLOR[c], mec="none",
                label=VCLASS_LABEL[c]) for c in classes]
    lg = ax.legend(handles=h, loc=loc, fontsize=fontsize, handletextpad=0.3,
                   borderpad=0.25, labelspacing=0.4, title=title, **kw)
    if title:
        lg.get_title().set_fontsize(fontsize)
        lg.get_title().set_color(MUTED)
    return lg


def footnote(ax, text: str, fontsize=5.6, y=-0.30, width: int | None = None):
    """A caption block under the axes, where it cannot land on data.

    `width` re-wraps to a fixed character count. Without it a long unbroken line widens
    the saved bounding box and the panel no longer matches its siblings' width.
    """
    if CLEAN:
        return
    if width:
        text = textwrap.fill(" ".join(text.split()), width=width)
    ax.text(0.0, y, text, transform=ax.transAxes, ha="left", va="top",
            fontsize=fontsize, color=MUTED)


STAR_KEY = "*** p < 0.001, ** p < 0.01, * p < 0.05, ns otherwise"


def stars(p: float) -> str:
    """Significance code. Printed WITH the p-value, never as a substitute for it."""
    if not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def bh_threshold(p: np.ndarray, q: np.ndarray, alpha: float = 0.05):
    ok = np.isfinite(p) & np.isfinite(q) & (q < alpha)
    return float(np.max(p[ok])) if ok.any() else None


def nmd_class_map() -> dict:
    """allele short name -> 'competent' / 'escaping' from the 50-nt rule."""
    n = pd.read_csv(R / "nmd_classification.csv")
    return dict(zip(n.id.str.split("|").str[0], n.nmd_class))


def refine_class(short: str, vclass: str, cls: dict) -> str:
    """Split `nonsense` into its NMD classes wherever the 50-nt rule has been applied."""
    if vclass == "nonsense" and short in cls:
        return f"nonsense_{cls[short]}"
    return vclass


# =========================================================================== a — UMAPs
def panel_a() -> None:
    """a1-a4: one UMAP per (cell line x colouring), cell-cycle-retained embedding.

    The `raw` embedding is used because these panels are a pair: the phase panel is only
    interpretable where phase structure is present, and colouring the SAME embedding by
    target is what makes "organised by phase, not by target" a within-panel comparison
    rather than a claim about two different pictures. `python -m analysis.fig_umap_target`
    draws the cell-cycle-regressed twin.
    """
    mix = pd.read_csv(R / "umap_mixing.csv")
    tag = "raw"
    for cl in paths.CELL_LINES:
        d = load(cl)
        m = mix[(mix.cell_line == cl) & (mix.embedding == tag)].iloc[0]
        xlim, ylim = natural_limits(d, tag)
        size, alpha = point_spec(len(d))
        for key, letter, note_txt in [
                ("guide_target", "a1" if cl == "A549" else "a2",
                 f"own-target neighbours {m.mixing_ratio:.2f}$\\times$ chance"),
                ("phase", "a3" if cl == "A549" else "a4",
                 f"phase = {m.phase_var_frac:.0%} of leading-PC variance"),
        ]:
            style.apply(8.0)
            w = 3.35
            fig = plt.figure(figsize=_sz(w, w * window_aspect(xlim, ylim) + 0.95))
            ax = fig.add_axes([0.01, 0.16, 0.98, 0.72])
            if key == "phase":
                # own scatter, so the phase panels use PHASE_COLOR_CB rather than the
                # single-hue ramp in fig_umap_target; alpha is raised because an
                # alpha-blended mark renders lighter than the hex the validator checked
                rng = np.random.default_rng(0)
                o = rng.permutation(len(d))
                ax.scatter(d[f"umap1_{tag}"].to_numpy()[o],
                           d[f"umap2_{tag}"].to_numpy()[o], s=size,
                           c=[PHASE_COLOR_CB[v] for v in d["phase"].to_numpy()[o]],
                           alpha=0.62, linewidths=0, rasterized=True)
            else:
                scatter(ax, d, tag, key, size, alpha)
            frame(ax, xlim, ylim)
            axis_key(ax)
            ttl = "CDK targeted" if key == "guide_target" else "cell-cycle phase"
            style.panel(ax, letter, f"{cl} — {ttl}", pad=4, size=8.4)
            ax.text(1.0, -0.03, f"{len(d):,} cells · {note_txt}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6.0, color=MUTED)
            if key == "guide_target":
                h = [Line2D([], [], marker="o", ls="", ms=4.2, mfc=target_color(g),
                            mec="none", label=g) for g in style.CDK_ORDER + ["NT"]]
                lg = fig.legend(handles=h, loc="lower center", ncol=5, frameon=False,
                                fontsize=6.0, handletextpad=0.25, columnspacing=0.9,
                                bbox_to_anchor=(0.5, -0.01))
            else:
                # direct labels as well as the legend, so identity is never colour-alone
                import matplotlib.patheffects as pe
                for ph in PHASE_ORDER:
                    m = d["phase"] == ph
                    if not m.any():
                        continue
                    ax.text(d.loc[m, f"umap1_{tag}"].median(),
                            d.loc[m, f"umap2_{tag}"].median(), ph, ha="center",
                            va="center", fontsize=7.2, weight="bold",
                            color=PHASE_COLOR_CB[ph], zorder=6,
                            path_effects=[pe.withStroke(linewidth=2.2,
                                                        foreground="white")])
                h = [Line2D([], [], marker="o", ls="", ms=4.6, mfc=PHASE_COLOR_CB[p],
                            mec="none", label=p) for p in PHASE_ORDER]
                lg = fig.legend(handles=h, loc="lower center", ncol=3, frameon=False,
                                fontsize=6.4, handletextpad=0.25, columnspacing=1.2,
                                bbox_to_anchor=(0.5, 0.0))
            for t in lg.get_texts():
                t.set_color(INK)
            slug = "target" if key == "guide_target" else "phase"
            save_panel(fig, f"PANEL_{letter}_umap_{slug}_{cl}")


# ================================================================ b — cell cycle, T47D
def panel_b() -> None:
    f = R / "g1_by_guide.csv"
    if not f.exists():
        raise SystemExit("missing g1_by_guide.csv — run `python -m analysis.mainpanel_stats`")
    g1 = pd.read_csv(f)
    cls = nmd_class_map()
    g1["vclass"] = [refine_class(s, v, cls) for s, v in zip(g1["short"], g1.variant_class)]
    t = g1[g1.cell_line == "T47D"].copy()
    tg, nt = t[~t.is_control], t[t.is_control]

    # ------------------------------------------------- b1: every T47D gRNA, effect vs p
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(4.1, 3.2))
    ax.axvspan(nt.delta_pp.min(), nt.delta_pp.max(), color=style.BAND, alpha=0.55, lw=0,
               zorder=0)
    ax.axvline(0, color=INK, lw=0.5, zorder=1)
    ax.axhline(-np.log10(0.05), color=MUTED, ls=":", lw=0.8, zorder=1)
    thr = bh_threshold(tg.p.to_numpy(), tg.q.to_numpy())
    if thr:
        ax.axhline(-np.log10(thr), color=style.HIT, ls="--", lw=0.8, zorder=1)
    for ed, sub in tg.groupby("editor"):
        for vc, s in sub.groupby("vclass"):
            ax.scatter(s.delta_pp, -np.log10(s.p), s=13, marker=ARM_MARKER[ed],
                       facecolor=VCLASS_COLOR.get(vc, MUTED), lw=0.3, edgecolor="white",
                       alpha=0.92, zorder=3)
    ax.scatter(nt.delta_pp, -np.log10(nt.p), s=16, marker="o", facecolor="none",
               edgecolor="#6E7377", lw=0.7, zorder=4)
    ax.set_ylim(-0.18, -np.log10(tg.p.min()) * 1.12)
    ax.set_xlim(tg.delta_pp.min() - 3, tg.delta_pp.max() + 6)
    ax.text(ax.get_xlim()[0] + 0.4, -np.log10(0.05), "raw p = 0.05", fontsize=5.4,
            color=MUTED, va="bottom", ha="left")
    if thr:
        ax.text(ax.get_xlim()[0] + 0.4, -np.log10(thr) - 0.06, "BH q = 0.05",
                fontsize=5.4, color=style.HIT, va="top", ha="left")
    hs = [Line2D([], [], marker=ARM_MARKER[e], ls="", ms=4.0, mfc=MUTED, mec="none",
                 label=e) for e in ["ABE", "CBE"]]
    hs += [Line2D([], [], marker="o", ls="", ms=4.4, mfc=VCLASS_COLOR[c], mec="none",
                  label=VCLASS_LABEL[c])
           for c in ["nonsense_competent", "nonsense_escaping", "start_loss", "missense"]]
    hs += [Line2D([], [], marker="o", ls="", ms=4.4, mfc="none", mec="#6E7377",
                  markeredgewidth=0.7, label="non-targeting")]
    ax.legend(handles=hs, loc="lower right", fontsize=5.2, handletextpad=0.3,
              borderpad=0.25, labelspacing=0.32)
    ax.set_xlabel("$\\Delta$G$_1$ vs held-out controls (percentage points)")
    ax.set_ylabel("$-$log$_{10}$ p   (Fisher exact, two-sided)")
    style.despine(ax)
    style.panel(ax, "b1", "Cell-cycle arrest is a CDK4 allele series in T47D", pad=4,
                size=8.2)
    lab = tg[tg.p < 0.05]
    place_labels(ax, lab.delta_pp, -np.log10(lab.p),
                 [s.replace("_", " ") for s in lab["short"]],
                 [VCLASS_COLOR.get(v, INK) for v in lab.vclass], fontsize=5.3,
                 obstacles=point_boxes(ax, t.delta_pp, -np.log10(t.p)))
    footnote(ax, f"T47D, {len(tg)} targeting gRNAs with $\\geq$20 cells, both arms. "
                 f"Shading = span of the {len(nt)} NT gRNA$\\times$arm units "
                 f"({int(nt.n_cells.min())}–{int(nt.n_cells.max())} cells each);\n"
                 f"none of them reaches raw p < 0.05. BH over the "
                 f"{len(tg)} targeting gRNAs of this cell line.", y=-0.215)
    save_panel(fig, "PANEL_b1_cellcycle_perguide_T47D")

    # ------------------------------------- b2: the CDK4 allele series, effect + 95% CI
    arm = "T47D/CBE"
    s = t[(t.arm == arm) & (t.Symbol == "CDK4")].sort_values("delta_pp")
    ntb = t[(t.arm == arm) & t.is_control]
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.6, 2.7))
    ax.axvspan(ntb.delta_pp.min(), ntb.delta_pp.max(), color=style.BAND, alpha=0.6, lw=0,
               zorder=0)
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    y = np.arange(len(s))
    for yi, (_, r) in zip(y, s.iterrows()):
        c = VCLASS_COLOR.get(r.vclass, INK)
        ax.plot([r.ci_lo_pp, r.ci_hi_pp], [yi, yi], color=c, lw=1.1, zorder=2)
        ax.scatter([r.delta_pp], [yi], s=26, color=c, zorder=3, lw=0.4, edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels([r["short"].split("_", 1)[1].replace("_", " ")
                        for _, r in s.iterrows()], fontsize=6.6)
    for tick, (_, r) in zip(ax.get_yticklabels(), s.iterrows()):
        tick.set_color(VCLASS_COLOR.get(r.vclass, INK))
    ax.set_ylim(-0.7, len(s) - 0.3)
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(y)
    ax2.set_yticklabels([f"{int(r.n_cells)} " + ("††" if r.q < 0.05 else
                                                 ("†" if r.p < 0.05 else ""))
                         for _, r in s.iterrows()], fontsize=5.8, color=MUTED)
    ax2.tick_params(length=0)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("$\\Delta$G$_1$ (pp), 95% Newcombe CI")
    vlegend(ax, ["nonsense_competent", "nonsense_escaping", "missense"], loc="lower right",
            fontsize=5.4)
    style.despine(ax)
    style.panel(ax, "b2", "Nonsense and missense alleles both arrest", pad=4, size=8.2)
    footnote(ax, f"{arm}, all eight CDK4 gRNAs. Right column = cells; "
                 "†† BH q < 0.05 over the 144 targeting gRNAs of the line, † raw p < 0.05 "
                 "(an asterisk in a name is a stop codon).\n"
                 "Shading = span of the nine NT gRNAs in this arm. D196N is missense and "
                 "its own transcript is intact (panel c3).", y=-0.26)
    save_panel(fig, "PANEL_b2_cellcycle_CDK4_forest_T47D")

    # --------------------------------------------- b3: phase composition, same gRNAs
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.6, 2.6))
    rows = [(f"NT pooled (n={int(ntb.n_cells.sum())})",
             float((ntb.g1 * ntb.n_cells).sum() / ntb.n_cells.sum()),
             float((ntb.s_frac * ntb.n_cells).sum() / ntb.n_cells.sum()),
             float((ntb.g2m_frac * ntb.n_cells).sum() / ntb.n_cells.sum()),
             "non-targeting")]
    for _, r in s.iloc[::-1].iterrows():
        rows.append((r["short"].split("_", 1)[1].replace("_", " "), r.g1, r.s_frac,
                     r.g2m_frac, r.vclass))
    yy = np.arange(len(rows))[::-1]
    left = np.zeros(len(rows))
    for k, ph in enumerate(PHASE_ORDER):
        vals = np.array([x[1 + k] for x in rows]) * 100
        ax.barh(yy, vals, left=left, height=0.72, color=PHASE_COLOR_CB[ph],
                edgecolor="white", lw=0.4, label=ph)
        left += vals
    ax.set_yticks(yy)
    ax.set_yticklabels([x[0] for x in rows], fontsize=6.2)
    for tick, x in zip(ax.get_yticklabels(), rows):
        tick.set_color(VCLASS_COLOR.get(x[4], INK))
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of cells")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=3, fontsize=6.0,
              handlelength=1.0, handletextpad=0.35, columnspacing=1.0)
    style.despine(ax)
    style.panel(ax, "b3", "Phase composition, T47D/CBE CDK4", pad=20, size=8.2)
    footnote(ax, "Phase calls from the delivered per-cell cell-cycle table. NT row pools "
                 "the nine non-targeting gRNAs of the arm.", y=-0.24)
    save_panel(fig, "PANEL_b3_phase_composition_T47D_CBE")


# ======================================================== c — editing / NMD evidence
def panel_c() -> None:
    tt = pd.read_csv(R / "target_transcript_all_guides.csv")
    tt["trunc"] = tt.variant_class.isin(["nonsense_competent", "start_loss"])
    meta = pd.read_csv(R / "editing_meta.csv")
    con = pd.read_csv(R / "editing_class_contrast.csv")
    order = ["missense", "start_loss", "nonsense_escaping", "nonsense_competent"]

    # --------------------------- c1: target-transcript loss by variant class, per line
    style.apply(8.0)
    fig, axes = plt.subplots(1, 2, figsize=_sz(4.7, 2.9), sharey=True)
    fig.subplots_adjust(top=0.78, wspace=0.10)
    rng = np.random.default_rng(0)
    for ax, cl in zip(axes, paths.CELL_LINES):
        d = tt[tt.cell_line == cl]
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        for i, vc in enumerate(order):
            v = d[d.variant_class == vc].z.dropna()
            if not len(v):
                continue
            ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=11,
                       color=VCLASS_COLOR[vc], lw=0.25, edgecolor="white", alpha=0.85,
                       zorder=3)
            ax.plot([i - 0.33, i + 0.33], [v.median()] * 2, color=INK, lw=1.4, zorder=4)
        c = con[(con.cell_line == cl) & con.is_headline].iloc[0]
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(["missense", "start\nloss", "PTC\nescaping", "PTC\ncompetent"],
                           fontsize=5.9)
        ax.set_title(f"{cl}\ntruncating < missense, p = {c.p:.3f}", fontsize=6.6, pad=3)
        ax.set_xlim(-0.6, len(order) - 0.4)
        style.despine(ax)
    axes[0].set_ylabel("target-transcript z\n(Poisson GLM vs controls)", fontsize=7.0)
    # the two sub-axes carry their own two-line titles, so the panel letter goes above
    # the figure rather than on an axes, where it would collide with them
    fig.text(0.0, 1.045, "$\\bf{c1}$   Truncating alleles lose their own transcript in "
                         "both lines", fontsize=8.2, ha="left", va="bottom")
    footnote(axes[0], "One-sided Mann–Whitney U on the per-gRNA z, truncating "
                      "(NMD-competent PTC + start loss, k = 26) vs missense (k = 116); "
                      "bar = median.\nThe two cell lines are tested separately: a gRNA "
                      "measured in both is one allele seen twice, not two observations.",
             y=-0.26)
    save_panel(fig, "PANEL_c1_editing_by_variant_class")

    # ------------- c1s: both cell lines on one narrow axes, with significance brackets
    # Classes are re-ordered so the two TRUNCATING classes are adjacent — the contrast
    # being tested is truncating vs missense, and a bracket cannot span a gap honestly.
    order_s = ["missense", "nonsense_escaping", "start_loss", "nonsense_competent"]
    tick_s = ["missense", "PTC\nescaping", "start\nloss", "PTC\ncompetent"]
    dodge = {"A549": -0.19, "T47D": 0.19}
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.3, 3.0))
    rng = np.random.default_rng(0)
    ax.axhline(0, color=INK, lw=0.5, zorder=1)
    for cl in paths.CELL_LINES:
        dd = tt[tt.cell_line == cl]
        for i, vc in enumerate(order_s):
            v = dd[dd.variant_class == vc].z.dropna()
            if not len(v):
                continue
            x0 = i + dodge[cl]
            ax.scatter(x0 + rng.uniform(-0.12, 0.12, len(v)), v, s=9,
                       marker=style.LINE_MARKER[cl], color=VCLASS_COLOR[vc], lw=0.2,
                       edgecolor="white", alpha=0.85, zorder=3)
            ax.plot([x0 - 0.15, x0 + 0.15], [v.median()] * 2, color=INK, lw=1.3, zorder=4)
    lo, hi = tt.z.min(), tt.z.max()
    # the two truncating classes, marked as one group under the axis
    ax.plot([2 - 0.34, 3 + 0.34], [lo - 0.45] * 2, color=INK, lw=0.9, clip_on=False,
            zorder=5)
    ax.text(2.5, lo - 0.62, "truncating", ha="center", va="top", fontsize=5.6, color=INK)
    # one significance bracket per cell line: missense vs the truncating pair
    y = hi + 0.30
    for cl in paths.CELL_LINES:
        c = con[(con.cell_line == cl) & con.is_headline].iloc[0]
        # The comparison is missense vs BOTH truncating classes. A bracket with one right
        # leg lands in the gap between them and reads as if only the nearer class were
        # tested, so the right side is a foot with a leg on each class.
        xa = 0 + dodge[cl]
        x1, x2 = 2 + dodge[cl], 3 + dodge[cl]
        xc = (x1 + x2) / 2
        top, mid = y + 0.30, y + 0.14
        ax.plot([xa, xa, xc], [y, top, top], color=MUTED, lw=0.7, clip_on=False, zorder=5)
        ax.plot([xc, xc], [top, mid], color=MUTED, lw=0.7, clip_on=False, zorder=5)
        ax.plot([x1, x1, x2, x2], [y, mid, mid, y], color=MUTED, lw=0.7, clip_on=False,
                zorder=5)
        ax.text((xa + xc) / 2, top + 0.08,
                f"{cl}  {stars(c.p)}  p = {c.p:.3f}", ha="center", va="bottom",
                fontsize=5.6, color=INK, clip_on=False)
        y += 1.55
    ax.set_xticks(range(len(order_s)))
    ax.set_xticklabels(tick_s, fontsize=5.9)
    ax.set_xlim(-0.6, len(order_s) - 0.4)
    ax.set_ylim(lo - 1.25, y + 0.25)
    ax.set_ylabel("target-transcript z  (Poisson GLM vs controls)", fontsize=7.0)
    hs = [Line2D([], [], marker=style.LINE_MARKER[c], ls="", ms=4.0, mfc=MUTED,
                 mec="none", label=c) for c in paths.CELL_LINES]
    ax.legend(handles=hs, loc="lower left", fontsize=5.6, handletextpad=0.3,
              borderpad=0.25, labelspacing=0.35)
    style.despine(ax)
    style.panel(ax, "c1s", "", pad=4, size=8.2)
    footnote(ax, f"Both cell lines on one axes; classes re-ordered so the two truncating "
                 f"classes are adjacent. The bracket forks so that it lands on BOTH "
                 f"truncating classes — the test pools them, it is not a comparison "
                 f"against either one alone. One-sided Mann–Whitney U on the per-gRNA z, "
                 f"truncating (NMD-competent PTC + start loss, k = 26) vs missense "
                 f"(k = 116); bar = median. The two lines are tested SEPARATELY and never "
                 f"pooled — a gRNA measured in both is one allele seen twice. "
                 f"Stars: {STAR_KEY}. Both contrasts land at a single star, so the p-value "
                 f"is printed beside it rather than left to the code.", y=-0.30, width=104)
    save_panel(fig, "PANEL_c1s_editing_by_variant_class_combined")

    # ------------------------------------------------ c2: the 50-nt rule, per allele
    nc = pd.read_csv(R / "nmd_readout_classified.csv")
    nc = nc[nc.on_target].copy()
    piv = nc.pivot_table(index=["gid", "nmd_class"], columns="cell_line",
                         values="pct_reduction").reset_index()
    piv = piv.sort_values(["nmd_class", "T47D"], ascending=[False, True])
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.6, 3.0))
    y = np.arange(len(piv))[::-1]
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    cols = [VCLASS_COLOR["nonsense_competent"] if c == "competent"
            else VCLASS_COLOR["nonsense_escaping"] for c in piv.nmd_class]
    for yi, (_, r) in zip(y, piv.iterrows()):
        ax.plot([r.A549, r.T47D], [yi, yi], color="#C8CCCF", lw=0.7, zorder=2)
    for mk, cl in [("o", "A549"), ("s", "T47D")]:
        ax.scatter(piv[cl], y, s=20, marker=mk, facecolor=cols, lw=0.4,
                   edgecolor="white", zorder=3, alpha=0.95)
    ax.set_yticks(y)
    ax.set_yticklabels([g.replace("_", " ") for g in piv.gid], fontsize=5.8)
    for tick, c in zip(ax.get_yticklabels(), cols):
        tick.set_color(c)
    ax.set_ylim(-0.8, len(piv) - 0.2)
    ax.set_xlabel("% reduction in the allele's own transcript")
    hs = [Line2D([], [], marker=m, ls="", ms=4.2, mfc=MUTED, mec="none", label=c)
          for m, c in [("o", "A549"), ("s", "T47D")]]
    hs += [Line2D([], [], marker="o", ls="", ms=4.4,
                  mfc=VCLASS_COLOR["nonsense_competent"], mec="none",
                  label="NMD-competent"),
           Line2D([], [], marker="o", ls="", ms=4.4,
                  mfc=VCLASS_COLOR["nonsense_escaping"], mec="none",
                  label="NMD-escaping")]
    ax.legend(handles=hs, loc="upper right", fontsize=5.4, handletextpad=0.3,
              borderpad=0.25, labelspacing=0.35)
    style.despine(ax)
    style.panel(ax, "c2", "Premature stops, by 50-nt-rule NMD class", pad=4, size=8.2)
    comp = nc[nc.nmd_class == "competent"].pct_reduction
    esc = nc[nc.nmd_class == "escaping"].pct_reduction
    footnote(ax, f"NMD class from the 50-nt rule on the MANE transcript "
                 f"(`nmd_rule.py`): a PTC >50 nt upstream of the last exon–exon junction "
                 f"is competent.\nMedian % reduction: competent {comp.median():+.1f}% "
                 f"(k = {len(comp)} guide$\\times$line), escaping {esc.median():+.1f}% "
                 f"(k = {len(esc)}). k = 4 is a qualitative internal control, not a test.",
             y=-0.17)
    save_panel(fig, "PANEL_c2_nmd_50nt_rule")

    # ------------------------- c3: transcript loss against the arrest, CDK4 in T47D/CBE
    e = pd.read_csv(R / "editing_target_transcript_phaseadj.csv")
    e = e.merge(tt[["arm", "guide", "variant_class"]], on=["arm", "guide"], how="left")
    arm = "T47D/CBE"
    d = e[e.arm == arm]
    cdk4 = d[d.target_gene == "CDK4"]
    other = d[d.target_gene != "CDK4"]
    rho = stats.spearmanr(cdk4.lfc_adj, cdk4.g1_shift_pp)
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.6, 2.95))
    ax.axhline(0, color=MUTED, lw=0.5, ls=":")
    ax.axvline(0, color=MUTED, lw=0.5, ls=":")
    ax.scatter(other.lfc_adj, other.g1_shift_pp, s=9, color="#D3D6D9", lw=0, zorder=2)
    for _, r in cdk4.iterrows():
        ax.scatter([r.lfc_adj], [r.g1_shift_pp], s=28,
                   color=VCLASS_COLOR.get(r.variant_class, INK), lw=0.4,
                   edgecolor="white", zorder=4)
    ax.set_xlabel("target-transcript log FC, phase-adjusted")
    ax.set_ylabel("$\\Delta$G$_1$ (percentage points)")
    ax.set_ylim(d.g1_shift_pp.min() - 3, d.g1_shift_pp.max() + 6)
    hs = [Line2D([], [], marker="o", ls="", ms=4.4, mfc=VCLASS_COLOR[c], mec="none",
                 label=VCLASS_LABEL[c])
          for c in ["nonsense_competent", "nonsense_escaping", "missense"]]
    hs += [Line2D([], [], marker="o", ls="", ms=4.0, mfc="#D3D6D9", mec="none",
                  label="other CDKs, same arm")]
    ax.legend(handles=hs, loc="lower left", fontsize=5.3, handletextpad=0.3,
              borderpad=0.25, labelspacing=0.35)
    style.despine(ax)
    style.panel(ax, "c3", "Arrest tracks transcript loss — except D196N", pad=4, size=8.2)
    place_labels(ax, cdk4.lfc_adj, cdk4.g1_shift_pp,
                 [g.split("_", 1)[1].replace("_", " ") for g in cdk4.guide],
                 [VCLASS_COLOR.get(v, INK) for v in cdk4.variant_class], fontsize=5.5,
                 obstacles=point_boxes(ax, d.lfc_adj, d.g1_shift_pp))
    footnote(ax, f"{arm}, the eight CDK4 gRNAs. Spearman $\\rho$ = {rho.statistic:+.2f}, "
                 f"p = {rho.pvalue:.3f} (n = {len(cdk4)}) — a description of eight points, "
                 f"not a test that survives\ncorrection. x is the Poisson-GLM coefficient "
                 f"with cell-cycle phase in the model, so it is transcript loss NOT "
                 f"attributable to the arrest itself.", y=-0.235)
    save_panel(fig, "PANEL_c3_editing_vs_cellcycle_CDK4_T47D")

    # --------------------------------------- c4: A549 — the same readout, forest by allele
    a = tt[(tt.cell_line == "A549") & tt.trunc].copy()
    a["key"] = a.guide.str.replace("_", " ") + "  (" + a.arm.str.split("/").str[1] + ")"
    a = a.sort_values("lfc")
    m = meta[(meta.cell_line == "A549") & (meta.group == "truncating")].iloc[0]
    ca = con[(con.cell_line == "A549") & con.is_headline].iloc[0]
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.7, 3.7))
    y = np.arange(len(a))[::-1]
    ax.axvspan(m.lo_re, m.hi_re, color=style.BAND, alpha=0.85, lw=0, zorder=0)
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    ax.axvline(m.mu_re, color=style.HIT, lw=0.9, ls="--", zorder=1)
    for yi, (_, r) in zip(y, a.iterrows()):
        c = VCLASS_COLOR.get(r.variant_class, INK)
        ax.plot([r.lfc - 1.96 * r.se, r.lfc + 1.96 * r.se], [yi, yi], color=c, lw=0.9,
                zorder=2)
        ax.scatter([r.lfc], [yi], s=17, color=c, zorder=3, lw=0.35, edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels(a.key, fontsize=5.4)
    for tick, vc in zip(ax.get_yticklabels(), a.variant_class):
        tick.set_color(VCLASS_COLOR.get(vc, INK))
    ax.set_ylim(-0.8, len(a) - 0.2)
    ax.set_xlabel("target-transcript log FC vs controls (95% CI)")
    vlegend(ax, ["nonsense_competent", "start_loss"], loc="upper right", fontsize=5.4)
    style.despine(ax)
    style.panel(ax, "c4", "The same readout in A549", pad=4, size=8.2)
    footnote(ax, f"A549, truncating alleles only, both arms (k = {int(m.k)}). Dashed line "
                 f"and band = random-effects pool (DerSimonian–Laird $\\tau^2$ with the\n"
                 f"Knapp–Hartung adjustment): {m.mu_re:+.3f} "
                 f"[{m.lo_re:+.3f}, {m.hi_re:+.3f}], p = {m.p_re:.2f}. The alleles are "
                 f"heterogeneous (I² = {m.I2:.0%}), so the fixed-effect value "
                 f"({m.mu_fe:+.3f},\np = {m.p_fe:.3f}) is an optimistic bound and is not "
                 f"the number to quote. The supported A549 statement is the "
                 f"class-level one in c1:\ntruncating alleles rank below missense, "
                 f"one-sided Mann–Whitney p = {ca.p:.3f}.", y=-0.13)
    save_panel(fig, "PANEL_c4_editing_A549_forest")


# ============================================================ d — severity per CDK
def _classify_severity() -> pd.DataFrame:
    """severity.csv with a `cls` column: nonsense split by the 50-nt NMD rule."""
    d = pd.read_csv(R / "severity.csv")
    d["cell_line"] = d.arm.str.split("/").str[0]
    d["short"] = d.guide_identity.str.split("|").str[0]
    cls = nmd_class_map()
    d["cls"] = [("non-targeting" if c else refine_class(s, v, cls))
                for s, v, c in zip(d["short"], d.variant_class, d.is_control)]
    return d


def panel_d() -> None:
    """d1-d2: severity per gRNA by CDK (KEY 15c-d). d3: severity by gRNA class."""
    sev = pd.read_csv(R / "severity.csv")
    sev["short"] = sev.guide_identity.str.split("|").str[0]
    sev["line"] = sev.arm.str.split("/").str[0]
    t = sev[~sev.is_control].copy()

    # (figsize, stem suffix, label font). The short variant trades height for width so the
    # label placer keeps somewhere to put the names it would otherwise have to drop.
    SHAPES = [((4.4, 2.9), "", 5.0), ((5.4, 2.0), "s", 4.5)]
    for k, ln in enumerate(paths.CELL_LINES):
      for shape, sfx, lab_fs in SHAPES:
          style.apply(8.0)
          fig, ax = plt.subplots(figsize=_sz(*shape))
          d = t[t.line == ln].copy()
          order = {g: i for i, g in enumerate(style.CDK_ORDER)}
          d["gi"] = d.Symbol.map(order)
          d = d.sort_values(["gi", "severity"]).reset_index(drop=True)
          rng = np.random.default_rng(0)
          xs = d.gi + rng.uniform(-0.30, 0.30, len(d))
          nt = sev[sev.is_control & (sev.line == ln)].severity
          ax.axhspan(nt.min(), nt.max(), color=style.BAND, zorder=0, lw=0)
          ax.axhline(0, color=INK, lw=0.5, zorder=1)
          ax.scatter(xs, d.severity, s=6 + 46 * (1 - d.q.clip(0, 1)),
                     c=[GC[g] for g in d.Symbol], lw=0.25, edgecolor="white", alpha=0.92,
                     zorder=3)
          # BH survivors get a dark ring rather than a text marker: "*" already means
          # "stop codon" in every variant name on this figure
          sig = d.q < 0.05
          ax.scatter(xs[sig.values], d.severity[sig], s=6 + 46 * (1 - d.q[sig].clip(0, 1)),
                     facecolor="none", edgecolor=INK, lw=0.7, zorder=4)
          ax.set_xticks(range(len(style.CDK_ORDER)))
          ax.set_xticklabels(style.CDK_ORDER, fontsize=6.2, rotation=35, ha="right")
          ax.set_xlim(-0.7, len(style.CDK_ORDER) - 0.3)
          ax.set_ylim(d.severity.min() - 0.5, d.severity.max() + 1.0)
          ax.set_ylabel("severity (SD of the size-matched null)", fontsize=7.2)
          hs = [Line2D([], [], marker="o", ls="", ms=np.sqrt(6 + 46 * (1 - q)) / 1.6,
                       color=style.CTRL, label=f"q = {q}") for q in [0.9, 0.5, 0.05]]
          ax.legend(handles=hs, fontsize=5.2, loc="upper left", handletextpad=0.3,
                    borderpad=0.25, labelspacing=0.5, title="point size", title_fontsize=5.2)
          style.despine(ax)
          style.panel(ax, f"d{k + 1}", f"Transcriptional severity per gRNA — {ln}", pad=4,
                      size=8.2)
          # label every gRNA that clears the NT band, not only the BH survivors: the band
          # is the calibration the reader is being asked to compare against, so the points
          # above it are the ones the panel is about. Ordered largest-first so that if the
          # placer runs out of room it is the smallest that loses its label.
          hit = d[d.severity > nt.max()].sort_values("severity", ascending=False)
          place_labels(ax, xs[hit.index], hit.severity,
                       [r["short"].replace("_", " ").replace(r.Symbol + " ", "")
                        for _, r in hit.iterrows()],
                       [GC.get(g, INK) for g in hit.Symbol], fontsize=lab_fs,
                       obstacles=point_boxes(ax, xs, d.severity, half=2.2))
          note = ("A549/ABE's high-severity gRNAs do not survive the variant-specific "
                  "decomposition (0/10 top-10 overlap, `guide_phenotypes.py`);\n"
                  "they are the shared targeting component, not allele biology."
                  if ln == "A549" else
                  "T47D/CBE is the one arm whose high-severity gRNAs are also "
                  "variant-specific (7/10 top-10 overlap, `guide_phenotypes.py`).")
          footnote(ax, f"{ln}, both editor arms, one point per gRNA. Severity = cross-fitted "
                       f"split-half $\\hat{{D}}$ standardised against a null of 9,000 "
                       f"pseudo-guides per arm (nine cell-count bands × 1,000 draws) built "
                       f"from control cells matched to each guide's size and gem-group "
                       f"composition; p is the one-sided empirical tail of that null and q is "
                       f"Benjamini–Hochberg within arm over all gRNAs including the nine "
                       f"NT. Shading = span of the nine NT gRNAs. Labelled points are the "
                       f"arm's {len(hit)} gRNAs above the NT band; ringed points "
                       f"are the {int((d.q < 0.05).sum())} at BH q < 0.05 "
                       f"— leads, not validated hits: per-variant "
                       f"effects are at or near chance in every arm (best arm T47D/CBE, "
                       f"exact p = 0.041, not significant over the four arms tested; "
                       f"STATS_AUDIT.md §1). Do not quote a count of q < 0.05 gRNAs. "
                       f"{note}",
                   y=-0.265, width=128)
          save_panel(fig, f"PANEL_d{k + 1}{sfx}_severity_by_CDK_{ln}")

    # ------------- d12: both cell lines stacked on one shared x axis -----------------
    # One figure, two rows, same CDK order and one x axis. The y scales are NOT shared:
    # severity is standardised within arm, so the two lines are on the same units but
    # their ranges differ (A549 reaches 10.0, T47D 8.7) and forcing a common y would
    # waste half of the T47D row.
    # The wide variant gives the label placer horizontal room instead of forcing it to
    # stack names vertically, so the labels can also come back up to a readable size.
    D12_SHAPES = [((5.2, 4.4), "", 4.5), ((7.6, 4.2), "w", 5.1)]
    for shape, sfx, lab_fs in D12_SHAPES:
        style.apply(8.0)
        fig, axes = plt.subplots(2, 1, figsize=_sz(*shape), sharex=True)
        fig.subplots_adjust(hspace=0.13)
        for ax, ln in zip(axes, paths.CELL_LINES):
            d = t[t.line == ln].copy()
            order = {g: i for i, g in enumerate(style.CDK_ORDER)}
            d["gi"] = d.Symbol.map(order)
            d = d.sort_values(["gi", "severity"]).reset_index(drop=True)
            rng = np.random.default_rng(0)
            xs = d.gi + rng.uniform(-0.30, 0.30, len(d))
            nt = sev[sev.is_control & (sev.line == ln)].severity
            ax.axhspan(nt.min(), nt.max(), color=style.BAND, zorder=0, lw=0)
            ax.axhline(0, color=INK, lw=0.5, zorder=1)
            ax.scatter(xs, d.severity, s=6 + 46 * (1 - d.q.clip(0, 1)),
                       c=[GC[g] for g in d.Symbol], lw=0.25, edgecolor="white", alpha=0.92,
                       zorder=3)
            sig = d.q < 0.05
            ax.scatter(xs[sig.values], d.severity[sig], s=6 + 46 * (1 - d.q[sig].clip(0, 1)),
                       facecolor="none", edgecolor=INK, lw=0.7, zorder=4)
            ax.set_xlim(-0.7, len(style.CDK_ORDER) - 0.3)
            ax.set_ylim(d.severity.min() - 0.5, d.severity.max() + 1.0)
            ax.set_ylabel(f"{ln}\nseverity", fontsize=7.0)
            style.despine(ax)
            hit = d[d.severity > nt.max()].sort_values("severity", ascending=False)
            place_labels(ax, xs[hit.index], hit.severity,
                         [r["short"].replace("_", " ").replace(r.Symbol + " ", "")
                          for _, r in hit.iterrows()],
                         [GC.get(g, INK) for g in hit.Symbol], fontsize=lab_fs,
                         obstacles=point_boxes(ax, xs, d.severity, half=2.2))
        axes[1].set_xticks(range(len(style.CDK_ORDER)))
        axes[1].set_xticklabels(style.CDK_ORDER, fontsize=6.2, rotation=35, ha="right")
        hs = [Line2D([], [], marker="o", ls="", ms=np.sqrt(6 + 46 * (1 - q)) / 1.6,
                     color=style.CTRL, label=f"q = {q}") for q in [0.9, 0.5, 0.05]]
        axes[0].legend(handles=hs, fontsize=5.2, loc="upper left", handletextpad=0.3,
                       borderpad=0.25, labelspacing=0.5, title="point size",
                       title_fontsize=5.2)
        style.panel(axes[0], "d12", "Transcriptional severity per gRNA", pad=4, size=8.2)
        footnote(axes[1],
                 "Both cell lines, both editor arms, one point per gRNA, on a shared CDK axis. "
                 "Severity = cross-fitted split-half D-hat standardised against a null of 9,000 "
                 "pseudo-guides per arm (nine cell-count bands x 1,000) matched to each guide's "
                 "size and gem-group composition; q is Benjamini-Hochberg within arm. Shading = "
                 "span of the nine NT gRNAs of that line; ringed points are at BH q < 0.05; "
                 "labelled points are those above the NT band. Leads, not validated hits: "
                 "per-variant effects are at or near chance in every arm (best arm T47D/CBE, "
                 "exact p = 0.041, not significant over the four arms tested; STATS_AUDIT.md "
                 "S1). Do not quote a count of q < 0.05 gRNAs. The same gRNA appears in both "
                 "rows -- the two cell lines are never pooled.", y=-0.40, width=128)
        save_panel(fig, f"PANEL_d12{sfx}_severity_by_CDK_stacked")

    # ------------------------------------------------- d3: severity by gRNA class
    for f in ["severity_by_class", "severity_class_contrast"]:
        if not (R / f"{f}.csv").exists():
            raise SystemExit(f"missing {f}.csv — run `python -m analysis.mainpanel_stats`")
    bc = pd.read_csv(R / "severity_by_class.csv")
    con = pd.read_csv(R / "severity_class_contrast.csv")
    cls_col = _classify_severity()
    order = ["non-targeting", "missense", "nonsense_escaping", "nonsense_competent",
             "start_loss"]
    tick = ["non-\ntargeting", "missense", "PTC\nescaping", "PTC\ncompetent",
            "start\nloss"]

    style.apply(8.0)
    fig, axes = plt.subplots(1, 2, figsize=_sz(5.2, 3.0), sharey=True)
    fig.subplots_adjust(wspace=0.08, top=0.78)
    rng = np.random.default_rng(0)
    for ax, ln in zip(axes, paths.CELL_LINES):
        d = cls_col[cls_col.cell_line == ln]
        nt = d[d.is_control].severity
        ax.axhspan(nt.min(), nt.max(), color=style.BAND, lw=0, zorder=0)
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        for i, c in enumerate(order):
            v = d[d.cls == c].severity
            if v.empty:
                continue
            ax.scatter(i + rng.uniform(-0.24, 0.24, len(v)), v, s=11,
                       color=VCLASS_COLOR.get(c, MUTED), lw=0.25, edgecolor="white",
                       alpha=0.85, zorder=3)
            ax.plot([i - 0.34, i + 0.34], [v.median()] * 2, color=INK, lw=1.5, zorder=4)
        r = con[(con.embedding == "raw") & (con.level == "cell line")
                & (con.group == ln)].iloc[0]
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(tick, fontsize=5.8)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_title(f"{ln}\ntruncating > missense  p = {r.p_one_sided:.3f}   "
                     f"($r_{{rb}}$ = {r.rank_biserial:+.2f})\n"
                     f"Kruskal–Wallis over classes  p = {r.kruskal_p:.2f}",
                     fontsize=6.2, pad=3)
        style.despine(ax)
    axes[0].set_ylabel("severity (SD of the size-matched null)", fontsize=7.0)
    fig.text(0.0, 1.055, "$\\bf{d3}$   Severity by gRNA class — a consistent direction, "
                         "no single significant test", fontsize=8.2, ha="left", va="bottom")
    per_arm = "; ".join(
        f"{r.group} p = {r.p_one_sided:.3f}"
        for _, r in con[(con.embedding == "raw") & (con.level == "arm")].iterrows())
    cc = con[(con.embedding == "cell-cycle-regressed") & (con.level == "cell line")]
    cc_txt = "; ".join(f"{r.group} p = {r.p_one_sided:.3f} "
                       f"($r_{{rb}}$ = {r.rank_biserial:+.2f})" for _, r in cc.iterrows())
    footnote(axes[0],
             f"One point per gRNA; bar = median; shading = span of the 18 NT gRNA × arm "
             f"units of the line. Nonsense is split by the 50-nt rule, because an "
             f"NMD-escaping premature stop is a known hypomorph and lumping it with the "
             f"competent stops dilutes the class the contrast is about. Within a cell line "
             f"each gRNA appears exactly once (144 rows, 144 distinct guide identities), so "
             f"pooling the two editor arms does not pseudo-replicate; the two CELL LINES are "
             f"never pooled, because there a gRNA is one allele measured twice. "
             f"Truncating = NMD-competent PTC + start loss (n = 26) vs missense (n = 116), "
             f"one-sided Mann–Whitney — one-sided because the direction is pre-specified. "
             f"Per arm: {per_arm}. On the cell-cycle-regressed severity: {cc_txt}. "
             f"The rank-biserial effect size is POSITIVE in all 12 line/arm × embedding "
             f"combinations (+0.10 to +0.28) and 2 of the 12 reach p < 0.05, with the "
             f"omnibus null every time — read this as a consistent direction across weak "
             f"tests, not as a detection.", y=-0.33, width=112)
    save_panel(fig, "PANEL_d3_severity_by_class")


# ============================================== e — severity vs the resistance axis
def panel_e() -> None:
    f = R / "severity_vs_resistance.csv"
    if not f.exists():
        raise SystemExit("missing severity_vs_resistance.csv — run "
                         "`python -m analysis.mainpanel_stats`")
    sr = pd.read_csv(f).dropna(subset=["max_drug_lfc", "severity"])

    # (figsize, suffix, label font, legend font). The narrow variant drops the legend
    # inside the axes to a smaller size and shrinks the allele labels, because the
    # placer has much less horizontal room to work with.
    E_SHAPES = [((4.0, 3.0), "", 5.1, 5.2), ((2.6, 3.2), "n", 4.4, 4.6)]
    for k, ln in enumerate(paths.CELL_LINES):
      for shape, sfx, lab_fs, leg_fs in E_SHAPES:
          d = sr[sr.cell_line == ln].copy()
          rho = stats.spearmanr(d.max_drug_lfc, d.severity)
          pear = stats.pearsonr(d.max_drug_lfc, d.severity)
          top = d.nlargest(10, "max_drug_lfc")
          rest = d.drop(top.index)
          mw = stats.mannwhitneyu(top.severity, rest.severity)

          style.apply(8.0)
          fig, ax = plt.subplots(figsize=_sz(*shape))
          ax.axhspan(d.nt_min.iloc[0], d.nt_max.iloc[0], color=style.BAND, lw=0, zorder=0)
          ax.axhline(0, color=INK, lw=0.5, zorder=1)
          ax.axvline(0, color=MUTED, lw=0.5, ls=":", zorder=1)
          for ed, sub in rest.groupby("editor"):
              ax.scatter(sub.max_drug_lfc, sub.severity, s=14, marker=ARM_MARKER[ed],
                         facecolor="#C3C7CA", lw=0.25, edgecolor="white", zorder=2)
          for ed, sub in top.groupby("editor"):
              ax.scatter(sub.max_drug_lfc, sub.severity, s=26, marker=ARM_MARKER[ed],
                         facecolor=style.HIT, lw=0.4, edgecolor="white", zorder=4)
          # the axis labels are the widest artists in the narrow variant, so they set
          # its bounding box; shorten them there rather than letting them define the width
          ax.set_xlabel("max drug log$_2$FC" if sfx
                        else "maximum drug log$_2$FC in the parent screen",
                        fontsize=7.2 if sfx else None)
          ax.set_ylabel("severity" if sfx else "severity (SD of the size-matched null)",
                        fontsize=7.2)
          ax.set_ylim(d.severity.min() - 0.6, d.severity.max() + 1.4)
          ax.set_xlim(d.max_drug_lfc.min() - 0.4, d.max_drug_lfc.max() + 1.3)
          hs = [Line2D([], [], marker=ARM_MARKER[e], ls="", ms=4.0, mfc=MUTED, mec="none",
                       label=e) for e in ["ABE", "CBE"]]
          hs += [Line2D([], [], marker="o", ls="", ms=4.4, mfc=style.HIT, mec="none",
                        label="10 strongest resistance alleles"),
                 Line2D([], [], marker="o", ls="", ms=4.0, mfc="#C3C7CA", mec="none",
                        label="all other gRNAs"),
                 Line2D([], [], marker="s", ls="", ms=5.0, mfc=style.BAND, mec="none",
                        label="span of the 9 NT gRNAs")]
          ax.legend(handles=hs, loc="upper right", fontsize=leg_fs,
                    handletextpad=0.3, borderpad=0.25, labelspacing=0.35)
          style.despine(ax)
          style.panel(ax, f"e{k + 1}", f"The strongest resistance alleles are "
                                       f"transcriptionally silent — {ln}", pad=4, size=8.2)
          place_labels(ax, top.max_drug_lfc, top.severity,
                       [s.replace("_", " ") for s in top["short"]],
                       [style.HIT] * len(top), fontsize=lab_fs,
                       obstacles=point_boxes(ax, d.max_drug_lfc, d.severity,
                                             half=3.0 if not sfx else 2.2))
          footnote(ax, f"{ln}, {len(d)} gRNAs. x = maximum log$_2$FC across every "
                       f"drug-vs-DMSO column of the parent screen — the most generous "
                       f"resistance coordinate\navailable, chosen so the negative result is "
                       f"argued against the strongest possible x. Spearman "
                       f"$\\rho$ = {rho.statistic:+.3f}, p = {rho.pvalue:.2f}; "
                       f"Pearson r = {pear.statistic:+.3f}, p = {pear.pvalue:.2f}.\n"
                       f"The 10 strongest resistance alleles: "
                       f"{int(top.within_nt_band.sum())}/10 inside the NT band, "
                       f"{int((top.q < 0.05).sum())}/10 at q < 0.05; vs all others "
                       f"Mann–Whitney p = {mw.pvalue:.2f}.", y=-0.235)
          save_panel(fig, f"PANEL_e{k + 1}{sfx}_severity_vs_drug_lfc_{ln}")


# ================================ f — is the NMD readout limited by transcript detection?
def panel_f() -> None:
    """f1-f4: what the target-transcript readout can and cannot see, and why.

    The short answer is two-part and the panels are built to keep the parts apart:
    per-ALLELE significance is detection-limited (f1-f3), but the CLASS-level absence of
    NMD outside CDK4 is not — it is bounded at a few percent (f4).
    """
    need = ["nmd_detection", "nmd_power", "nmd_bound", "nmd_thinning"]
    for n in need:
        if not (R / f"{n}.csv").exists():
            raise SystemExit(f"missing {n}.csv — run `python -m analysis.nmd_power`")
    det = pd.read_csv(R / "nmd_detection.csv")
    pw = pd.read_csv(R / "nmd_power.csv")
    bd = pd.read_csv(R / "nmd_bound.csv")
    th = pd.read_csv(R / "nmd_thinning.csv")
    tr = pw[pw.is_truncating].copy()
    tr["exposure"] = tr.n_cells * tr.lam

    # -------------------------------------------------- f1: the detection landscape
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.7, 2.8))
    y = {g: i for i, g in enumerate(style.CDK_ORDER[::-1])}
    for arm, sub in det.groupby("arm"):
        ed = arm.split("/")[1]
        for _, r in sub.iterrows():
            ax.scatter(r.mean_umi, y[r.gene] + (0.16 if r.cell_line == "T47D" else -0.16),
                       s=22, marker=ARM_MARKER[ed], facecolor=GC[r.gene],
                       edgecolor="white", lw=0.35, zorder=3,
                       alpha=0.95 if r.cell_line == "T47D" else 0.55)
    ax.axvline(1.0, color=MUTED, ls=":", lw=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_yticks(list(y.values()))
    ax.set_yticklabels(list(y.keys()), fontsize=6.6)
    for tick, g in zip(ax.get_yticklabels(), y.keys()):
        tick.set_color(GC[g])
    ax.set_xlabel("target transcript in control cells (UMI per cell)")
    hs = [Line2D([], [], marker=ARM_MARKER[e], ls="", ms=4.2, mfc=MUTED, mec="none",
                 label=e) for e in ["ABE", "CBE"]]
    hs += [Line2D([], [], marker="o", ls="", ms=4.4, mfc="#4D4D4D", mec="none",
                  label="T47D (solid)"),
           Line2D([], [], marker="o", ls="", ms=4.4, mfc="#4D4D4D", mec="none",
                  alpha=0.5, label="A549 (pale)")]
    ax.legend(handles=hs, loc="lower right", fontsize=5.3, handletextpad=0.3,
              borderpad=0.25, labelspacing=0.35)
    style.despine(ax)
    style.panel(ax, "f1", "CDK4 in T47D is the only well-detected target", pad=4, size=8.2)
    top = det.nlargest(1, "mean_umi").iloc[0]
    footnote(ax, f"Mean UMI per cell in each arm's clean control cells, 3' 10x. Range "
                 f"{det.mean_umi.min():.2f}–{det.mean_umi.max():.2f} "
                 f"({det.mean_umi.max() / det.mean_umi.min():.0f}-fold); "
                 f"{int((det.mean_umi < 1).sum())} of {len(det)} arm × gene combinations "
                 f"sit below 1 UMI/cell (dotted line). The maximum is {top.gene} in "
                 f"{top.arm} at {top.mean_umi:.1f}.", y=-0.24, width=104)
    save_panel(fig, "PANEL_f1_target_transcript_detection")

    # ------------------------- f2: minimum detectable effect against the observed effect
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.7, 3.0))
    lim = (0, max(tr.mde_pct.max(), tr.pct_reduction.max()) + 8)
    ax.fill_between(lim, lim, lim[1], color=style.BAND, alpha=0.5, lw=0, zorder=0)
    ax.plot(lim, lim, color=INK, lw=0.8, zorder=1)
    ax.axhline(0, color=MUTED, lw=0.5, ls=":", zorder=1)
    for ln, sub in tr.groupby("cell_line"):
        ax.scatter(sub.mde_pct, sub.pct_reduction, s=22,
                   marker=style.LINE_MARKER[ln],
                   facecolor=[GC[g] for g in sub.target_gene], lw=0.35,
                   edgecolor="white", alpha=0.92, zorder=3)
    ax.set_xlim(*lim)
    ax.set_ylim(min(tr.pct_reduction.min() - 6, -30), lim[1])
    ax.set_xlabel("minimum detectable reduction (%), 80% power")
    ax.set_ylabel("observed reduction (%)")
    res = tr[tr.pct_reduction > tr.mde_pct]
    hs = [Line2D([], [], marker=style.LINE_MARKER[l], ls="", ms=4.2, mfc=MUTED,
                 mec="none", label=l) for l in paths.CELL_LINES]
    hs += [Line2D([], [], marker="s", ls="", ms=5.0, mfc=style.BAND, mec="none",
                  label="resolvable")]
    ax.legend(handles=hs, loc="lower right", fontsize=5.3, handletextpad=0.3,
              borderpad=0.25, labelspacing=0.35)
    style.despine(ax)
    style.panel(ax, "f2", "Most truncating alleles fall below their own detection floor",
                pad=4, size=8.2)
    place_labels(ax, res.mde_pct, res.pct_reduction,
                 [f"{g.replace('_', ' ')} ({c})"
                  for g, c in zip(res.guide, res.cell_line)],
                 [GC.get(g, INK) for g in res.target_gene], fontsize=5.2,
                 obstacles=point_boxes(ax, tr.mde_pct, tr.pct_reduction, half=3.0))
    footnote(ax, f"One point per truncating allele × arm (k = {len(tr)}). x is the "
                 f"reduction detectable at 80% power and two-sided α = 0.05 given that "
                 f"allele's cells and its target's abundance, from "
                 f"Var(β̂) ≈ 1/(n_g·λ) + 1/(n_c·λ). Points in the shaded region exceed "
                 f"their own floor: {len(res)} of {len(tr)}, "
                 f"{int((res.target_gene == 'CDK4').sum())} of them CDK4. This panel is "
                 f"about which alleles can be CALLED individually — it is not evidence "
                 f"that the others have a large effect; see f4.", y=-0.235, width=104)
    save_panel(fig, "PANEL_f2_detection_floor")

    # ------------------------------------------------- f3: the thinning demonstration
    keep = ["CDK4_W198*", "CDK4_Q222*", "CDK4_R62*", "CDK4_Q270*"]
    t = th[th.guide.isin(keep) & ~th.separated]
    g = (t.groupby(["guide", "thin"])
           .agg(lam=("lam", "mean"), lfc=("lfc", "mean"), lfc_sd=("lfc", "std"),
                z=("z", "mean"), z_sd=("z", "std")).reset_index())
    other_lam = det[det.gene != "CDK4"].mean_umi
    style.apply(8.0)
    fig, axes = plt.subplots(2, 1, figsize=_sz(3.7, 4.0), sharex=True)
    fig.subplots_adjust(hspace=0.16)
    for ax, col, sd, lab in [
            (axes[0], "lfc", "lfc_sd", "effect estimate\n(target-transcript log FC)"),
            (axes[1], "z", "z_sd", "z  (effect / SE)")]:
        ax.axvspan(other_lam.min(), other_lam.max(), color=style.BAND, alpha=0.55, lw=0,
                   zorder=0)
        ax.axhline(0, color=MUTED, lw=0.5, ls=":", zorder=1)
        for gd, sub in g.groupby("guide"):
            sub = sub.sort_values("lam")
            c = (VCLASS_COLOR["nonsense_escaping"] if gd.endswith("Q270*")
                 else VCLASS_COLOR["nonsense_competent"])
            ax.plot(sub.lam, sub[col], color=c, lw=1.0, marker="o", ms=3.0, zorder=3,
                    alpha=0.95)
            ax.fill_between(sub.lam, sub[col] - sub[sd].fillna(0),
                            sub[col] + sub[sd].fillna(0), color=c, alpha=0.16, lw=0,
                            zorder=2)
            ax.annotate(gd.split("_", 1)[1], (sub.lam.iloc[-1], sub[col].iloc[-1]),
                        xytext=(4, 0), textcoords="offset points", fontsize=5.4,
                        color=c, va="center")
        ax.set_xscale("log")
        ax.set_ylabel(lab, fontsize=6.8)
        ax.margins(x=0.20)
        style.despine(ax)
    axes[1].axhspan(-1.96, 1.96, color="#EDEEEF", zorder=0, lw=0)
    axes[1].set_xlabel("target transcript after thinning (UMI per cell)")
    style.panel(axes[0], "f3", "Thinning the transcript destroys the signal, not the effect",
                pad=4, size=8.2)
    footnote(axes[1], f"T47D/CBE CDK4, the one well-powered series, with ONLY its target "
                      f"transcript's counts binomially thinned — same cells, same library "
                      f"sizes, same biology, {int(th.rep.max()) + 1} replicates per level. "
                      f"Thinning scales λ and leaves the true rate ratio unchanged, and "
                      f"the top panel confirms the estimate is flat. The bottom panel "
                      f"collapses into the grey |z| < 1.96 band. Shading = the range of "
                      f"the other eight CDK transcripts "
                      f"({other_lam.min():.2f}–{other_lam.max():.2f} UMI/cell): at their "
                      f"abundance this series would not be detected either. Q270* is the "
                      f"NMD-escaping allele and is null at every level.",
             y=-0.40, width=104)
    save_panel(fig, "PANEL_f3_thinning")

    # ---------------------------------------------- f4: what can actually be excluded
    want = ["all truncating", "CDK4 only", "excluding CDK4", "exposure high",
            "exposure mid", "exposure low"]
    b = bd[(bd.estimator == "robust HC0") & bd.group.isin(want)].copy()
    b["ord"] = b.group.map({g: i for i, g in enumerate(want)})
    b = b.sort_values("ord", ascending=False)
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.8, 2.5))
    y = np.arange(len(b))
    for v, lab, c in [(50, "50% (a typical NMD effect)", style.HIT),
                      (25, "25%", "#E8A33D")]:
        ax.axvline(v, color=c, ls="--", lw=0.9, zorder=1)
        ax.text(v, len(b) - 0.35, f" {lab}", fontsize=5.3, color=c, va="top", ha="left")
    ax.axvline(0, color=INK, lw=0.6, zorder=1)
    for yi, (_, r) in zip(y, b.iterrows()):
        c = style.HIT if r.group == "CDK4 only" else "#35618F"
        ax.plot([r.pct_lo, r.pct_hi], [yi, yi], color=c, lw=1.2, zorder=2)
        ax.scatter([r.pct_reduction], [yi], s=26, color=c, lw=0.4, edgecolor="white",
                   zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.group}  (k={int(r.k)})" for _, r in b.iterrows()],
                       fontsize=6.2)
    ax.set_ylim(-0.7, len(b) - 0.1)
    ax.set_xlim(-12, 58)
    ax.set_xlabel("pooled reduction in the allele's own transcript (%), 95% CI")
    style.despine(ax)
    style.panel(ax, "f4", "The absence of NMD is a bound, not a blind spot", pad=4,
                size=8.2)
    lowest = b[b.group == "exposure low"].iloc[0]
    footnote(ax, f"Inverse-variance pool of the per-allele log fold changes, "
                 f"heteroskedasticity-robust (HC0) SEs. Exposure = n_cells × λ; strata are "
                 f"tertiles. Even the LOWEST-exposure tertile pools to "
                 f"{lowest.pct_reduction:+.1f}% [{lowest.pct_lo:+.1f}, "
                 f"{lowest.pct_hi:+.1f}], excluding both a 25% and a 50% reduction. "
                 f"Identical conclusions under Poisson and quasi-Poisson SEs "
                 f"(`nmd_bound.csv`); overdispersion is mild (median quasi-Poisson scale "
                 f"1.08). So the readout is not simply blind outside CDK4 — the average "
                 f"effect there really is small.", y=-0.30, width=104)
    save_panel(fig, "PANEL_f4_nmd_bound")


# ======================== g — pre-registered gene panels: what should each CDK's loss do?
PANEL_ROW_ORDER = ["e2f_dream", "g2m_checkpoint", "short_half_life", "myc_targets",
                   "ddr_long", "dna_repair", "ifn_stat1", "sterol_srebp",
                   "translation_biosynth"]


def panel_g() -> None:
    """g1-g3: the panel x target specificity matrix, its calibration, and the one result.

    Run `python -m analysis.targeted_panels` first.
    """
    for n in ["targeted_panel_tests", "targeted_panel_composite", "targeted_panel_perm",
              "targeted_panel_definitions"]:
        if not (R / f"{n}.csv").exists():
            raise SystemExit(f"missing {n}.csv — run `python -m analysis.targeted_panels`")
    M = pd.read_csv(R / "targeted_panel_tests.csv")
    CP = pd.read_csv(R / "targeted_panel_composite.csv")
    PM = pd.read_csv(R / "targeted_panel_perm.csv")
    DEF = pd.read_csv(R / "targeted_panel_definitions.csv").set_index("panel")
    S = pd.read_csv(R / "targeted_panel_scores.csv.gz")

    rows = [k for k in PANEL_ROW_ORDER if k in set(M.panel)]
    lab = {k: f"{DEF.loc[k, 'label']}  ({int(DEF.loc[k, 'n_genes'])})" for k in rows}

    # ------------------------------------------------------------- g1: the matrix per arm
    for arm in paths.ARMS:
        style.apply(8.0)
        fig, axes = plt.subplots(1, 2, figsize=_sz(6.6, 2.9), sharey=True)
        fig.subplots_adjust(wspace=0.06, top=0.80)
        for ax, emb, ttl in zip(axes, ["raw", "cc"],
                                ["cell cycle retained", "cell cycle regressed out"]):
            d = M[(M.arm == arm) & (M.embedding == emb) & (M.unit == "LoF gRNAs")]
            Z = np.full((len(rows), len(style.CDK_ORDER)), np.nan)
            for i, pk in enumerate(rows):
                for j, g in enumerate(style.CDK_ORDER):
                    r = d[(d.panel == pk) & (d.target == g)]
                    if len(r):
                        r = r.iloc[0]
                        sgn = np.sign(r.median_z - r.median_z_rest) or 1.0
                        Z[i, j] = sgn * -np.log10(max(r.p, 1e-6))
            im = ax.imshow(Z, cmap=style.DIVERGING, vmin=-4, vmax=4, aspect="auto")
            for i, pk in enumerate(rows):
                for j, g in enumerate(style.CDK_ORDER):
                    r = d[(d.panel == pk) & (d.target == g)]
                    if not len(r):
                        ax.plot([j - 0.22, j + 0.22], [i, i], color="#B0B0B0",
                                lw=0.8, zorder=5)
                        continue
                    r = r.iloc[0]
                    if r.is_prediction:
                        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                                   edgecolor=INK, lw=1.1, zorder=4))
                    if r.q < 0.05:
                        # a drawn marker, not a glyph: Helvetica has no BLACK STAR and
                        # the PDF would carry a missing-glyph box into Illustrator
                        ax.scatter([j], [i], s=20, marker="o", facecolor="white",
                                   edgecolor=INK, lw=1.0, zorder=5)
            ax.set_xticks(range(len(style.CDK_ORDER)))
            ax.set_xticklabels(style.CDK_ORDER, fontsize=5.8, rotation=45, ha="right")
            for tick, g in zip(ax.get_xticklabels(), style.CDK_ORDER):
                tick.set_color(GC[g])
            ax.set_title(ttl, fontsize=6.6, pad=3)
            ax.tick_params(length=0)
            for sp in ax.spines.values():
                sp.set_visible(False)
        axes[0].set_yticks(range(len(rows)))
        axes[0].set_yticklabels([lab[k] for k in rows], fontsize=5.8)
        for tick, k in zip(axes[0].get_yticklabels(), rows):
            if DEF.loc[k, "source"].startswith("post-hoc"):
                tick.set_color(MUTED)
        cb = fig.colorbar(im, ax=axes, fraction=0.022, pad=0.015)
        cb.set_label("signed $-$log$_{10}$ p   (blue = panel down)", fontsize=5.8)
        cb.ax.tick_params(labelsize=5.4)
        cb.outline.set_visible(False)
        fig.text(0.0, 1.02, f"$\\bf{{g1}}$   Pre-registered panel $\\times$ target "
                            f"matrix — {arm}", fontsize=8.2, ha="left", va="bottom")
        pm = PM[(PM.arm == arm) & (PM.unit == "LoF gRNAs")]
        txt = "; ".join(f"{r.embedding} {int(r.observed_cells_p05)} observed vs "
                        f"{r.null_mean:.1f} expected (permutation p = {r.p_perm:.2f})"
                        for _, r in pm.iterrows())
        footnote(axes[0], f"Mann–Whitney U, two-sided, of {arm}'s loss-of-function gRNAs "
                          f"(NMD-competent PTC + start loss) for each target against every "
                          f"other targeting gRNA in the arm, on each panel's NT-standardised "
                          f"score. A black box marks the pre-registered prediction for that "
                          f"CDK; a white ring marks BH q < 0.05 over the 81 cells of that "
                          f"matrix; a short dash marks a CDK with no LoF gRNA in this arm. The bottom row is a post-hoc panel "
                          f"(grey label) carried as a comparator, not a prediction. "
                          f"Matrix-wide label-permutation calibration, cells at p < 0.05: "
                          f"{txt}.", y=-0.42, width=150)
        save_panel(fig, f"PANEL_g1_panel_matrix_{arm.replace('/', '_')}")

    # ------------------------------------------- g2: the composite pre-registered test
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(3.8, 2.6))
    c = CP.sort_values(["embedding", "arm"], ascending=[True, False]).reset_index(drop=True)
    y = np.arange(len(c))
    ax.axvline(0.5, color=INK, lw=0.7, zorder=2)
    for yi, (_, r) in zip(y, c.iterrows()):
        ax.plot([r.null_mean - 1.645 * r.null_sd, r.null_mean + 1.645 * r.null_sd],
                [yi, yi], color="#C8CCCF", lw=5.0, solid_capstyle="butt", zorder=1)
        col = style.HIT if r.p_perm_one_sided < 0.05 else "#35618F"
        ax.scatter([r.T_obs], [yi], s=30, color=col, lw=0.4, edgecolor="white", zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.arm}  ·  {r.embedding}" for _, r in c.iterrows()], fontsize=6.2)
    ax.set_ylim(-0.7, len(c) - 0.3)
    ax.set_xlabel("T — mean normalised rank of the LoF gRNAs,\noriented so higher = as "
                  "predicted", fontsize=6.8)
    for yi, (_, r) in zip(y, c.iterrows()):
        # outside the axes, so a p-value can never land on its own point
        ax.annotate(f"p = {r.p_perm_one_sided:.3f}", (1.0, yi),
                    xycoords=("axes fraction", "data"), xytext=(4, 0),
                    textcoords="offset points", fontsize=5.4, ha="left", va="center",
                    color=style.HIT if r.p_perm_one_sided < 0.05 else MUTED)
    style.despine(ax)
    style.panel(ax, "g2", "All predictions as one test per arm", pad=4, size=8.2)
    footnote(ax, f"Every pre-registered (target, panel) pair with a specified direction "
                 f"({int(c.n_cells.iloc[0])} pairs) combined into one statistic: the mean "
                 f"normalised rank of each target's LoF gRNAs on d × z, where d is the "
                 f"predicted direction. 0.5 is chance. Grey bar = central 90% of "
                 f"{int(c.n_perm.iloc[0])} gene-label permutations, which preserve the "
                 f"scores and the number of LoF gRNAs per gene and break only the "
                 f"correspondence being tested. One-sided, because every contributing "
                 f"direction was fixed in advance. Nothing survives correction for the "
                 f"{len(c)} arm × embedding tests shown.", y=-0.36, width=112)
    save_panel(fig, "PANEL_g2_prediction_composite")

    # -------------------------------------------- g3: the E2F/DREAM panel, all four arms
    e = S[(S.panel == "e2f_dream") & (S.embedding == "raw")].copy()
    style.apply(8.0)
    fig, axes = plt.subplots(2, 2, figsize=_sz(6.8, 4.6), sharey=True)
    fig.subplots_adjust(hspace=0.42, top=0.90)
    rng = np.random.default_rng(0)
    for ax, arm in zip(axes.ravel(), paths.ARMS):
        d = e[e.arm == arm]
        tgt, nt = d[~d.is_control], d[d.is_control]
        ax.axhspan(nt.z.min(), nt.z.max(), color=style.BAND, lw=0, zorder=0)
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        for j, g in enumerate(style.CDK_ORDER):
            sub = tgt[tgt.Symbol == g]
            if sub.empty:
                continue
            x = j + rng.uniform(-0.26, 0.26, len(sub))
            is_lof = sub.is_lof.to_numpy(bool)
            ax.scatter(x[~is_lof], sub.z[~is_lof], s=9, color="#C7CBCE", lw=0, zorder=3)
            ax.scatter(x[is_lof], sub.z[is_lof], s=26, color=GC[g], lw=0.4,
                       edgecolor="white", zorder=4)
        r = M[(M.arm == arm) & (M.embedding == "raw") & (M.unit == "LoF gRNAs")
              & (M.panel == "e2f_dream") & M.is_prediction].nsmallest(1, "p")
        note = (f"best pre-registered cell: {r.target.iloc[0]} p = {r.p.iloc[0]:.2g}"
                if len(r) else "")
        ax.set_xticks(range(len(style.CDK_ORDER)))
        if arm.startswith("T47D"):
            ax.set_xticklabels(style.CDK_ORDER, fontsize=5.6, rotation=45, ha="right")
            for tick, g in zip(ax.get_xticklabels(), style.CDK_ORDER):
                tick.set_color(GC[g])
        else:
            ax.set_xticklabels([])
        ax.set_title(f"{arm}   ·   {note}", fontsize=6.2, pad=3)
        ax.set_xlim(-0.7, len(style.CDK_ORDER) - 0.3)
        style.despine(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("E2F / DREAM panel z\n(vs the nine NT gRNAs)", fontsize=6.8)
    fig.text(0.0, 1.015, "$\\bf{g3}$   The one panel that carries a result: E2F / DREAM",
             fontsize=8.2, ha="left", va="bottom")
    footnote(axes[1, 0], "Mean control-normalised residual over the 30 E2F/DREAM genes per "
                         "gRNA, shared targeting axis removed leave-one-out, standardised "
                         "against the arm's nine NT gRNAs (shaded). Coloured points are the "
                         "pre-specified loss-of-function gRNAs (NMD-competent PTC + start "
                         "loss); grey are that CDK's other gRNAs. Predicted direction is "
                         "DOWN for CDK4, CDK6, CDK2 and CDK7. Cell-cycle-regressed twin in "
                         "`targeted_panel_scores.csv.gz`.", y=-0.42, width=126)
    save_panel(fig, "PANEL_g3_e2f_dream_by_target")


# ============================ h — the same figures drawn with energy distance instead
def panel_h() -> None:
    """h1-h5: what panels d1-d3 look like on the field-standard metric.

    Run `python -m analysis.edist_calibrate` first (which needs `benchmark.py`'s
    `edistance.csv` / `edistance_null.csv`).

    THE POINT OF THIS BLOCK is not that one metric wins. It is that energy distance and
    severity answer different questions on this dataset, and the difference is visible
    rather than argued: E-distance separates targeting from non-targeting almost perfectly
    and then has nothing left to rank variants with, while severity has dynamic range but
    a much weaker targeting-vs-control signal.
    """
    f = R / "edist_calibrated.csv"
    if not f.exists():
        raise SystemExit("missing edist_calibrated.csv — run "
                         "`python -m analysis.edist_calibrate`")
    E = pd.read_csv(f)
    con = pd.read_csv(R / "edist_class_contrast.csv")
    scon = pd.read_csv(R / "severity_class_contrast.csv")
    model = pd.read_csv(R / "edist_size_model.csv")
    null = pd.read_csv(R / "edistance_null.csv")

    # ------------------------- h1, h2: per-gRNA by CDK — the d1/d2 twin, same layout
    for k, ln in enumerate(paths.CELL_LINES):
        d = E[(E.cell_line == ln) & (~E.is_control)].reset_index(drop=True)
        nt = E[(E.cell_line == ln) & E.is_control].z_edist
        style.apply(8.0)
        fig, ax = plt.subplots(figsize=_sz(4.4, 2.9))
        order = {g: i for i, g in enumerate(style.CDK_ORDER)}
        d["gi"] = d.Symbol.map(order)
        d = d.sort_values(["gi", "z_edist"]).reset_index(drop=True)
        rng = np.random.default_rng(0)
        xs = d.gi + rng.uniform(-0.30, 0.30, len(d))
        ax.axhspan(nt.min(), nt.max(), color=style.BAND, zorder=0, lw=0)
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        ax.scatter(xs, d.z_edist, s=16, c=[GC[g] for g in d.Symbol], lw=0.25,
                   edgecolor="white", alpha=0.92, zorder=3)
        ax.set_xticks(range(len(style.CDK_ORDER)))
        ax.set_xticklabels(style.CDK_ORDER, fontsize=6.2, rotation=35, ha="right")
        ax.set_xlim(-0.7, len(style.CDK_ORDER) - 0.3)
        ax.set_ylabel("energy distance, size-calibrated\n(SD of the matched null)",
                      fontsize=7.0)
        above = int((d.z_edist > nt.max()).sum())
        sev_above = int((d.severity > E[(E.cell_line == ln)
                                        & E.is_control].severity.max()).sum())
        ax.text(0.015, 0.975, f"{above}/{len(d)} targeting gRNAs above the NT band\n"
                              f"(severity, same gRNAs: {sev_above}/{len(d)})",
                transform=ax.transAxes, va="top", fontsize=6.0, color=style.HIT)
        style.despine(ax)
        style.panel(ax, f"h{k + 1}", f"Energy distance per gRNA — {ln}", pad=4, size=8.2)
        footnote(ax, f"The panel d{k + 1} layout, same gRNAs, different metric. Unbiased "
                     f"energy distance E = 2E‖x−y‖ − E‖x−x′‖ − E‖y−y′‖ between the gRNA's "
                     f"cells and clean control cells in 30 PCs (PCA fitted on control "
                     f"cells only), both sets capped at 400 cells "
                     f"(`benchmark.py`, the scPerturb statistic). Standardised here "
                     f"against its own size-matched control-only null, because that "
                     f"null's SD falls {model[model.arm.str.startswith(ln)].sd.max() / model[model.arm.str.startswith(ln)].sd.min():.1f}-fold "
                     f"across the cell counts in this line while its mean stays at ~0 "
                     f"(`edist_calibrate.py`). Shading = span of the nine NT gRNAs per arm. "
                     f"No labels, because essentially every point clears the band — that "
                     f"is the result.", y=-0.265, width=128)
        save_panel(fig, f"PANEL_h{k + 1}_edistance_by_CDK_{ln}")

    # ------------------------------------------- h3: by gRNA class — the d3 twin
    order = ["non-targeting", "missense", "nonsense_escaping", "nonsense_competent",
             "start_loss"]
    tick = ["non-\ntargeting", "missense", "PTC\nescaping", "PTC\ncompetent",
            "start\nloss"]
    style.apply(8.0)
    fig, axes = plt.subplots(1, 2, figsize=_sz(5.2, 3.0), sharey=True)
    fig.subplots_adjust(wspace=0.08, top=0.78)
    rng = np.random.default_rng(0)
    for ax, ln in zip(axes, paths.CELL_LINES):
        d = E[E.cell_line == ln]
        nt = d[d.is_control].z_edist
        ax.axhspan(nt.min(), nt.max(), color=style.BAND, lw=0, zorder=0)
        ax.axhline(0, color=INK, lw=0.5, zorder=1)
        for i, c in enumerate(order):
            v = d[d.cls == c].z_edist
            if v.empty:
                continue
            ax.scatter(i + rng.uniform(-0.24, 0.24, len(v)), v, s=11,
                       color=VCLASS_COLOR.get(c, MUTED), lw=0.25, edgecolor="white",
                       alpha=0.85, zorder=3)
            ax.plot([i - 0.34, i + 0.34], [v.median()] * 2, color=INK, lw=1.5, zorder=4)
        r = con[(con.metric == "z_edist") & (con.cell_line == ln)].iloc[0]
        sr = scon[(scon.embedding == "raw") & (scon.level == "cell line")
                  & (scon.group == ln)].iloc[0]
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(tick, fontsize=5.8)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_title(f"{ln}\ntruncating > missense  p = {r.p_one_sided:.2f}   "
                     f"($r_{{rb}}$ = {r.rank_biserial:+.2f})\n"
                     f"severity gave  p = {sr.p_one_sided:.3f}   "
                     f"($r_{{rb}}$ = {sr.rank_biserial:+.2f})", fontsize=6.0, pad=3)
        style.despine(ax)
    axes[0].set_ylabel("energy distance, size-calibrated\n(SD of the matched null)",
                       fontsize=7.0)
    fig.text(0.0, 1.055, "$\\bf{h3}$   The panel d3 contrast on energy distance — it does "
                         "not reproduce", fontsize=8.2, ha="left", va="bottom")
    footnote(axes[0],
             f"Same gRNAs, same classes, same one-sided Mann–Whitney as d3. Energy distance "
             f"separates targeting from non-targeting almost completely — the gap between "
             f"the NT band and everything else is the largest effect on this figure — but "
             f"it does NOT order the loss-of-function classes above missense. In A549 the "
             f"contrast is reversed (r_rb = "
             f"{con[(con.metric == 'z_edist') & (con.cell_line == 'A549')].rank_biserial.iloc[0]:+.2f}) "
             f"and in T47D it is flat "
             f"({con[(con.metric == 'z_edist') & (con.cell_line == 'T47D')].rank_biserial.iloc[0]:+.3f}), "
             f"where severity gave +0.14 and +0.22. The same reversal is present on the RAW "
             f"E-distance (`edist_class_contrast.csv`), so it is not an artefact of the "
             f"calibration. Read this as: the truncating > missense ordering is "
             f"metric-dependent and should not be stated as a property of the data without "
             f"naming the statistic.", y=-0.33, width=112)
    save_panel(fig, "PANEL_h3_edistance_by_class")

    # ------------------------------------------------ h4: the two metrics against each other
    style.apply(8.0)
    fig, axes = plt.subplots(2, 2, figsize=_sz(5.4, 4.6))
    fig.subplots_adjust(hspace=0.46, wspace=0.30, top=0.90)
    for ax, arm in zip(axes.ravel(), paths.ARMS):
        d = E[(E.arm == arm) & (~E.is_control)]
        c = E[(E.arm == arm) & E.is_control]
        ax.scatter(d.severity, d.z_edist, s=12, c=[GC[g] for g in d.Symbol], lw=0.25,
                   edgecolor="white", alpha=0.9, zorder=3)
        ax.scatter(c.severity, c.z_edist, s=16, facecolor="none", edgecolor="#6E7377",
                   lw=0.7, zorder=4)
        r = stats.spearmanr(d.severity, d.z_edist)
        ax.set_title(f"{arm}   ·   Spearman $\\rho$ = {r.statistic:+.2f} "
                     f"(p = {r.pvalue:.1g})", fontsize=6.2, pad=3)
        ax.axhline(0, color=MUTED, lw=0.5, ls=":")
        ax.axvline(0, color=MUTED, lw=0.5, ls=":")
        ax.tick_params(labelsize=5.8)
        style.despine(ax)
    for ax in axes[1, :]:
        ax.set_xlabel("severity", fontsize=6.8)
    for ax in axes[:, 0]:
        ax.set_ylabel("energy distance (calibrated)", fontsize=6.8)
    fig.text(0.0, 1.015, "$\\bf{h4}$   The two metrics, gRNA by gRNA", fontsize=8.2,
             ha="left", va="bottom")
    footnote(axes[1, 0], "Open circles are the nine NT gRNAs. Both metrics come from the "
                         "SAME cells, so they share sampling noise and this is agreement, "
                         "not independent validation. The agreement is moderate in T47D and "
                         "near zero in A549, which is itself informative: the two statistics "
                         "are not measuring the same thing, and neither is a drop-in "
                         "substitute for the other.", y=-0.42, width=118)
    save_panel(fig, "PANEL_h4_severity_vs_edistance")

    # --------------------------------------- h5: why the published-style call is miscalibrated
    style.apply(8.0)
    fig, axes = plt.subplots(1, 4, figsize=_sz(7.4, 2.3), sharey=False)
    fig.subplots_adjust(wspace=0.34, top=0.80)
    raw = pd.read_csv(R / "edistance.csv").rename(columns={"guide": "guide_identity"})
    for ax, arm in zip(axes, paths.ARMS):
        d = raw[raw.arm == arm]
        m = model[model.arm == arm].sort_values("n")
        nl = null[null.arm == arm]
        ax.scatter(nl.n, nl.edist, s=3, color="#DDE0E2", lw=0, rasterized=True, zorder=1)
        t, c = d[~d.is_control], d[d.is_control]
        ax.scatter(np.clip(t.n_cells, 20, 400), t.edist, s=9, color="#8E9499", lw=0,
                   zorder=3)
        ax.scatter(np.clip(c.n_cells, 20, 400), c.edist, s=18, facecolor="none",
                   edgecolor=INK, lw=0.7, zorder=5)
        ax.axhline(d.null_p95.iloc[0], color=style.HIT, ls="--", lw=0.9, zorder=4)
        ax.plot(m.n, m["mean"] + 1.645 * m["sd"], color="#35618F", lw=1.1, zorder=4)
        # linear, not log: the cell-count range inside one arm is narrow, and a log
        # axis here puts unlabelled minor ticks on top of each other
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
        ax.set_title(arm, fontsize=6.4, pad=3)
        ax.set_xlabel("cells (capped at 400)", fontsize=6.2)
        ax.tick_params(labelsize=5.6)
        style.despine(ax)
    axes[0].set_ylabel("energy distance", fontsize=6.8)
    hs = [Line2D([], [], color=style.HIT, ls="--", lw=1.0,
                 label="pooled 95th percentile (as published)"),
          Line2D([], [], color="#35618F", lw=1.2, label="size-matched 95th percentile"),
          Line2D([], [], marker="o", ls="", ms=3.6, mfc="#DDE0E2", mec="none",
                 label="control-only null draws"),
          Line2D([], [], marker="o", ls="", ms=3.6, mfc="#8E9499", mec="none",
                 label="targeting gRNAs"),
          Line2D([], [], marker="o", ls="", ms=4.0, mfc="none", mec=INK,
                 markeredgewidth=0.7, label="the nine NT gRNAs")]
    fig.legend(handles=hs, loc="lower center", ncol=5, frameon=False, fontsize=5.4,
               handletextpad=0.3, columnspacing=1.1, bbox_to_anchor=(0.5, -0.10))
    fig.text(0.0, 1.02, "$\\bf{h5}$   Why the published-style E-distance call is "
                        "miscalibrated by cell count", fontsize=8.2, ha="left", va="bottom")
    footnote(axes[0], "`benchmark.py` calls a gRNA significant if its E-distance exceeds a "
                      "single pooled 95th percentile (red dashed). The null's MEAN is ~0 at "
                      "every size — the unbiased form of E does its job — but its SPREAD "
                      "falls 1.5–4.6× across the cell counts in an arm, so the correct "
                      "threshold is the blue curve. Recalibrating removes the 2 of 36 NT "
                      "gRNAs that the pooled threshold called significant and changes no "
                      "targeting call at all: 278 of 288 targeting gRNAs clear either "
                      "threshold. The calibration was worth fixing; it is not what limits "
                      "the metric here.", y=-0.62, width=150)
    save_panel(fig, "PANEL_h5_edistance_calibration")

    # ------------- h6: the control that decides whether any of h1-h3 means anything
    f6 = R / "edist_pca_control.csv"
    if not f6.exists():
        print("  skipping h6 — run `python -m analysis.edist_decompose pca`")
        return
    C = pd.read_csv(f6)
    Cd = pd.read_csv(R / "edist_pca_control_dispersion.csv")
    fits = ["PCA on ALL controls (as published)", "PCA on held-out half (symmetric)"]
    style.apply(8.0)
    fig, axes = plt.subplots(1, 2, figsize=_sz(5.6, 3.0), sharey=True)
    fig.subplots_adjust(wspace=0.08, top=0.80)
    rng = np.random.default_rng(0)
    for ax, fit in zip(axes, fits):
        d = C[C.fit == fit]
        ax.axhline(0, color=INK, lw=0.6, zorder=1)
        for i, (lab, sub) in enumerate([("targeting", d[~d.is_control]),
                                        ("non-targeting", d[d.is_control])]):
            col = "#8E9499" if i == 0 else VCLASS_COLOR["non-targeting"]
            col = style.HIT if i == 1 else "#4E5A63"
            ax.scatter(i + rng.uniform(-0.22, 0.22, len(sub)), sub.edist, s=10,
                       color=col, lw=0.2, edgecolor="white", alpha=0.8, zorder=3)
            ax.plot([i - 0.32, i + 0.32], [sub.edist.median()] * 2, color=INK, lw=1.6,
                    zorder=4)
            ax.annotate(f"{int(sub.sig.sum())}/{len(sub)}\nover null",
                        (i, ax.get_ylim()[1]), xytext=(0, -2),
                        textcoords="offset points", ha="center", va="top", fontsize=5.8,
                        color=col)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["targeting\ngRNAs", "the nine\nNT gRNAs"], fontsize=6.2)
        ax.set_xlim(-0.6, 1.6)
        disp = Cd[Cd.fit == fit]
        ax.set_title(f"{fit}\nwithin-set spread vs matched controls:  targeting "
                     f"{disp[~disp.is_control].ratio_vs_matched.median():.2f}   ·   NT "
                     f"{disp[disp.is_control].ratio_vs_matched.median():.2f}",
                     fontsize=5.8, pad=3)
        style.despine(ax)
    axes[0].set_ylabel("energy distance vs control cells", fontsize=7.0)
    fig.text(0.0, 1.055, "$\\bf{h6}$   The E-distance excess follows the PCA fit, not the "
                         "gRNA", fontsize=8.2, ha="left", va="bottom")
    footnote(axes[0],
             "All four arms pooled; each point is one gRNA, scored against a size-matched "
             "DISJOINT control-only null. LEFT is the published computation: the PCA is "
             "fitted on every clean-control cell, so the NT gRNAs are inside the fit and the "
             "targeting gRNAs are not. RIGHT fits the PCA on half the control cells and "
             "scores against the other half, so NEITHER side is in the fit. The excess "
             "swaps sides. A variance-maximising projection fitted on one set makes any "
             "set held out of it look tighter, and E-distance reads that through its "
             "−E‖x−x′‖ − E‖y−y′‖ terms. Under the symmetric fit the targeting gRNAs' median "
             "E-distance is ~0 and 8–19% clear the null, against 5% expected — which is "
             "the same weak signal severity reports, not a large missed effect. "
             "`analysis/edist_decompose.py pca`.", y=-0.33, width=118)
    save_panel(fig, "PANEL_h6_edistance_pca_control")


# ================================================== i — library composition donuts
LIB_CATEGORIES = ["non-targeting", "nonsense / start loss", "resistance hit",
                  "fitness hit"]
# Outer-ring hues are keyed to the variant-class palette used in blocks b-d, so the two
# rings are told apart by what they encode rather than by chroma alone.
CAT_COLOR = {
    "non-targeting": VCLASS_COLOR["non-targeting"],
    "nonsense / start loss": VCLASS_COLOR["nonsense"],
    "resistance hit": VCLASS_COLOR["missense"],
    "fitness hit": "#C08A2E",
}


def panel_i() -> None:
    """i1, i2: nested donut of the designed library, one panel per editor.

    Outer ring = the four library categories. Inner ring = the CDK each gRNA targets,
    nested inside its category so the two rings line up.
    """
    f = R / "library_composition.csv"
    if not f.exists():
        raise SystemExit("missing library_composition.csv — run "
                         "`python -m analysis.mainpanel_stats`")
    L = pd.read_csv(f)

    for k, ed in enumerate(["ABE", "CBE"]):
        d = L[L.editor == ed]
        outer_v, outer_c, inner_v, inner_c, inner_lab = [], [], [], [], []
        for cat in LIB_CATEGORIES:
            sub = d[d.category == cat]
            if sub.empty:
                continue
            outer_v.append(int(sub.n_guides.sum()))
            outer_c.append(CAT_COLOR[cat])
            genes = [g for g in style.CDK_ORDER if g in set(sub.Symbol)]
            genes += [g for g in sub.Symbol if g not in genes]      # NT falls out here
            for g in genes:
                n = int(sub[sub.Symbol == g].n_guides.sum())
                inner_v.append(n)
                inner_c.append("#C7CBCE" if g == "NT" else GC[g])
                inner_lab.append((g, n))
        total = sum(outer_v)

        style.apply(8.0)
        fig, ax = plt.subplots(figsize=_sz(3.5, 3.5))
        w = 0.26
        ax.pie(outer_v, radius=1.0, colors=outer_c, startangle=90, counterclock=False,
               wedgeprops=dict(width=w, edgecolor="white", linewidth=0.9))
        wedges, _ = ax.pie(inner_v, radius=1.0 - w - 0.02, colors=inner_c, startangle=90,
                           counterclock=False,
                           wedgeprops=dict(width=w, edgecolor="white", linewidth=0.7))
        ax.set(aspect="equal")
        ax.text(0, 0.06, ed, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(0, -0.10, f"{total} gRNAs", ha="center", va="center", fontsize=7.0,
                color=MUTED)

        # outer category labels, outside the ring, on the wedge bisector
        ang = 90.0
        for v, cat in zip(outer_v, [c for c in LIB_CATEGORIES
                                    if not d[d.category == c].empty]):
            span = 360.0 * v / total
            mid = np.deg2rad(ang - span / 2)
            x, y = np.cos(mid), np.sin(mid)
            ax.annotate(f"{cat}\n{v}  ({v / total:.0%})",
                        (1.02 * x, 1.02 * y), xytext=(1.30 * x, 1.30 * y),
                        ha="left" if x >= 0 else "right", va="center", fontsize=6.0,
                        color=CAT_COLOR[cat],
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="#C3C7CA"))
            ang -= span

        # inner gene labels, only where the wedge is wide enough to hold one
        for wd, (g, n) in zip(wedges, inner_lab):
            span = abs(wd.theta2 - wd.theta1)
            # 16 degrees is about the narrowest wedge that can hold two lines of type
            # without the text spilling over its own boundary; the rest are in the legend
            if span < 16:
                continue
            mid = np.deg2rad((wd.theta1 + wd.theta2) / 2)
            r = (1.0 - w - 0.02) - w / 2
            ax.text(r * np.cos(mid), r * np.sin(mid), f"{g}\n{n}", ha="center",
                    va="center", fontsize=5.0, color="white", linespacing=0.95)
        style.panel(ax, f"i{k + 1}", f"{ed} library composition", pad=8, size=8.4)
        hs = [Line2D([], [], marker="o", ls="", ms=4.2, mfc=GC[g], mec="none", label=g)
              for g in style.CDK_ORDER]
        hs += [Line2D([], [], marker="o", ls="", ms=4.2, mfc="#C7CBCE", mec="none",
                      label="NT")]
        ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=5,
                  fontsize=5.6, handletextpad=0.25, columnspacing=0.8, frameon=False,
                  title="inner ring — CDK targeted", title_fontsize=5.6)
        footnote(ax, f"Designed guide library ({total} gRNAs: "
                     f"{total - int(d[d.category == 'non-targeting'].n_guides.sum())} "
                     f"targeting + "
                     f"{int(d[d.category == 'non-targeting'].n_guides.sum())} "
                     f"non-targeting), not recovered cells. The same ten NT gRNAs are in "
                     f"both libraries. Categories are assigned by PRECEDENCE, not as a "
                     f"partition of the source annotation: molecular consequence first "
                     f"(a * in the mutation is a premature stop, M1x a start loss), then "
                     f"the parent screen's `Type` for the remaining missense gRNAs. "
                     f"{int(L.n_also_screen_hit.sum())} truncating gRNAs across both "
                     f"libraries were also listed as a resistance or fitness hit and are "
                     f"counted once, under their consequence. Counts in "
                     f"`results/library_composition.csv`.", y=-0.20, width=96)
        save_panel(fig, f"PANEL_i{k + 1}_library_donut_{ed}")

# ============================ j — what severity is, and how it differs from E-distance
def panel_j() -> None:
    """j1: a schematic of the two perturbation-magnitude statistics used in this set.

    The point clouds are SYNTHETIC — this is a schematic of the estimators, not a plot of
    the data. Every NUMBER quoted on it is read from `results/`, so the annotations cannot
    drift away from the figures they describe.
    """
    G = int(np.load(R / "gene_universe.npy", allow_pickle=True).size)
    sev = pd.read_csv(R / "severity.csv")
    sev["line"] = sev.arm.str.split("/").str[0]
    nmin, nmax = int(sev[~sev.is_control].n_cells.min()), int(sev[~sev.is_control].n_cells.max())
    null = pd.read_csv(R / "pseudo_guide_null.csv")
    n_draw = int(null.groupby("arm").size().iloc[0])
    n_band = int(null.groupby("arm").n_cells.nunique().iloc[0])
    above = {}
    for ln, g in sev.groupby("line"):
        nt, tg = g[g.is_control].severity, g[~g.is_control]
        above[ln] = (int((tg.severity > nt.max()).sum()), len(tg))
    C = pd.read_csv(R / "edist_pca_control.csv")
    D = pd.read_csv(R / "edist_pca_control_dispersion.csv")
    ed = {}
    for fit, g in C.groupby("fit"):
        tg = g[~g.is_control]
        d = D[(D.fit == fit) & (~D.is_control)]
        key = "pub" if "ALL" in fit else "sym"
        ed[key] = (int(tg.sig.sum()), len(tg), float(d.ratio_vs_matched.median()))

    SEV_C, EDI_C = "#2a78d6", "#eb6834"
    CTRL_C = "#C3C7CA"
    rng = np.random.default_rng(3)

    style.apply(8.0)
    fig = plt.figure(figsize=_sz(7.6, 7.0))

    def cloud(ax, n, cx, cy, sx, sy, c, s=5, alpha=0.75, marker="o"):
        x = rng.normal(cx, sx, n)
        y = rng.normal(cy, sy, n)
        ax.scatter(x, y, s=s, color=c, lw=0, alpha=alpha, marker=marker, zorder=3)
        return x, y

    def blank(ax):
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)

    # ---------------------------------------------------------- row 1: construction
    axA = fig.add_axes([0.035, 0.620, 0.42, 0.275])
    axB = fig.add_axes([0.545, 0.620, 0.42, 0.275])
    for ax, ttl, col in [(axA, "severity — cross-fitted magnitude of the MEAN shift", SEV_C),
                         (axB, "energy distance — a full DISTRIBUTIONAL distance", EDI_C)]:
        blank(ax)
        ax.set_xlim(0, 10); ax.set_ylim(0, 6)
        ax.text(0, 6.35, ttl, fontsize=7.4, weight="bold", color=col, va="bottom")

    # -- severity: one guide's cells, split in half
    axA.text(1.55, 5.5, "one gRNA's cells", fontsize=6.2, ha="center", color=INK)
    xa, ya = cloud(axA, 45, 1.2, 3.4, 0.42, 0.55, SEV_C, s=6)
    xb, yb = cloud(axA, 45, 1.9, 3.4, 0.42, 0.55, "#8FBCE8", s=6)
    axA.plot([1.55, 1.55], [2.0, 4.8], color=INK, ls="--", lw=0.7, zorder=4)
    axA.text(1.05, 1.75, "half A", fontsize=5.8, ha="center", color=SEV_C)
    axA.text(2.15, 1.75, "half B", fontsize=5.8, ha="center", color="#5A93CC")
    axA.annotate("", xy=(3.35, 3.4), xytext=(2.6, 3.4),
                 arrowprops=dict(arrowstyle="-|>", lw=0.9, color=MUTED))
    # two pseudobulk vectors as gene strips
    for i, (yy, c, lab) in enumerate([(4.1, SEV_C, "$m_A$"), (2.6, "#8FBCE8", "$m_B$")]):
        v = rng.normal(0, 1, 34)
        axA.imshow(v[None, :], extent=(3.6, 7.4, yy - 0.34, yy + 0.34), aspect="auto",
                   cmap="seismic", vmin=-2.6, vmax=2.6, zorder=3)
        axA.add_patch(plt.Rectangle((3.6, yy - 0.34), 3.8, 0.68, fill=False,
                                    edgecolor=INK, lw=0.6, zorder=4))
        axA.text(3.45, yy, lab, fontsize=7.0, ha="right", va="center", color=c)
    axA.text(5.5, 4.72, f"control-normalised pseudobulk mean\nover $G$ = {G:,} genes",
             fontsize=5.9, ha="center", color=MUTED)
    axA.annotate("", xy=(8.15, 3.4), xytext=(7.6, 3.4),
                 arrowprops=dict(arrowstyle="-|>", lw=0.9, color=MUTED))
    axA.text(8.3, 3.4, r"$\hat{D}=\dfrac{\langle m_A,\,m_B\rangle}{G}$", fontsize=8.2,
             va="center", color=INK)
    axA.text(3.6, 1.68, "the two halves are disjoint, so their sampling noise is\n"
                        "independent and the cross term has mean zero:", fontsize=5.9,
             color=MUTED)
    axA.text(3.6, 0.95, r"$\mathbb{E}\,\langle m_A,m_B\rangle=\|\mu\|^{2}$"
                        "     (no floor)", fontsize=7.0, color=SEV_C)
    axA.text(3.6, 0.25, r"$\mathbb{E}\,\|\hat m\|^{2}=\|\mu\|^{2}"
                        r"+\mathrm{tr}(\Sigma)/n$" "   (has one)", fontsize=7.0,
             color=style.HIT)

    # -- energy distance: two clouds and the three distance terms
    gx, gy = cloud(axB, 40, 3.1, 3.5, 0.62, 0.72, EDI_C, s=6)
    cx, cy = cloud(axB, 40, 5.6, 3.5, 0.62, 0.72, CTRL_C, s=6)
    axB.text(3.1, 5.15, "gRNA cells $X$", fontsize=6.2, ha="center", color=EDI_C)
    axB.text(5.6, 5.15, "control cells $Y$", fontsize=6.2, ha="center", color="#8A9095")
    for (x0, y0, x1, y1, c, ls) in [(3.3, 3.8, 5.4, 3.2, INK, "-"),
                                    (2.7, 3.1, 3.4, 3.9, EDI_C, "-"),
                                    (5.3, 3.9, 5.9, 3.1, "#8A9095", "-")]:
        axB.plot([x0, x1], [y0, y1], color=c, lw=0.8, ls=ls, zorder=5)
    axB.text(4.35, 2.55, r"$\|x-y\|$", fontsize=6.4, ha="center", color=INK)
    axB.text(2.45, 4.3, r"$\|x-x'\|$", fontsize=6.4, ha="center", color=EDI_C)
    axB.text(6.25, 4.3, r"$\|y-y'\|$", fontsize=6.4, ha="center", color="#8A9095")
    axB.text(0.1, 1.55, r"$\mathcal{E}=2\,\overline{\|x-y\|}"
                        r"-\overline{\|x-x'\|}-\overline{\|y-y'\|}$", fontsize=8.2,
             color=INK)
    axB.text(0.1, 0.72, "zero if and only if the two DISTRIBUTIONS are identical —\n"
                        "so it answers a strictly broader question than $\\hat{D}$ does.",
             fontsize=5.9, color=MUTED)

    # ---------------------------------------------- row 2: what each one responds to
    labs = ["mean shifts,\nspread equal", "spread differs,\nmeans equal"]
    verdict = {("sev", 0): ("detected", SEV_C), ("sev", 1): ("invisible", style.HIT),
               ("edi", 0): ("detected", EDI_C), ("edi", 1): ("detected", EDI_C)}
    for j, (tag, x0, col) in enumerate([("sev", 0.035, SEV_C), ("edi", 0.545, EDI_C)]):
        for i in range(2):
            ax = fig.add_axes([x0 + i * 0.215, 0.415, 0.175, 0.145])
            blank(ax)
            ax.set_xlim(0, 10); ax.set_ylim(0, 10)
            if i == 0:
                cloud(ax, 34, 3.4, 5, 0.75, 1.5, CTRL_C, s=5)
                cloud(ax, 34, 6.6, 5, 0.75, 1.5, col, s=5)
            else:
                cloud(ax, 34, 5, 5, 0.75, 2.5, CTRL_C, s=5)
                cloud(ax, 34, 5, 5, 0.42, 1.0, col, s=5)
            v, vc = verdict[(tag, i)]
            ax.set_title(labs[i], fontsize=5.8, color=MUTED, pad=2)
            ax.text(5, -1.4, v, fontsize=6.6, weight="bold", ha="center", color=vc)
    fig.text(0.035, 0.592, "what the statistic responds to", fontsize=7.0,
             weight="bold", color=INK)

    # --------------------------------------------------------- row 3: the consequences
    axT = fig.add_axes([0.035, 0.030, 0.93, 0.285])
    blank(axT)
    axT.set_xlim(0, 1); axT.set_ylim(0, 1)
    rows = [
        ("estimand",
         r"$\|\mu\|^{2}$ — squared magnitude of the mean shift",
         "any difference in distribution: location, spread or shape"),
        ("null it is scored against",
         f"{n_draw:,} pseudo-gRNAs per arm\n({n_band} cell-count bands "
         f"$\\times$ {n_draw // n_band:,}), matched to each\ngRNA's $n$ "
         f"({nmin}–{nmax}) and gem groups",
         "control-only pseudo-gRNAs, matched $n$;\nthe published call used ONE pooled "
         "threshold"),
        ("standardised form",
         r"severity $=\dfrac{\hat{D}-\mu_{\mathrm{null}}(n)}"
         r"{\sigma_{\mathrm{null}}(n)}$",
         r"$z=\dfrac{\mathcal{E}-\mu_{\mathrm{null}}(n)}"
         r"{\sigma_{\mathrm{null}}(n)}$   (`edist_calibrate.py`)"),
        ("on this dataset",
         f"above the NT band: {above['A549'][0]}/{above['A549'][1]} A549, "
         f"{above['T47D'][0]}/{above['T47D'][1]} T47D",
         f"over its null: {ed['pub'][0]}/{ed['pub'][1]} with the PCA fitted on controls,\n"
         f"{ed['sym'][0]}/{ed['sym'][1]} with a symmetric fit"),
        ("the trap",
         "none known; the cross-fit removes the $1/n$ floor\nthat "
         r"$\|\hat m\|^{2}$ would carry",
         f"reads spread, so a projection fitted on $Y$ makes anything held\nout of it look "
         f"tighter ({ed['pub'][2]:.2f}$\\times$ vs {ed['sym'][2]:.2f}$\\times$ symmetric) — "
         f"panel h6"),
    ]
    ys = [0.86, 0.66, 0.44, 0.24, 0.05]
    axT.text(0.005, 1.03, "side by side", fontsize=7.0, weight="bold", color=INK)
    axT.plot([0, 1], [0.99, 0.99], color=INK, lw=0.7)
    axT.text(0.34, 0.965, "severity", fontsize=7.2, weight="bold", color=SEV_C,
             ha="center", va="top")
    axT.text(0.78, 0.965, "energy distance", fontsize=7.2, weight="bold", color=EDI_C,
             ha="center", va="top")
    for (lab, a, b), y in zip(rows, ys):
        axT.text(0.0, y, lab, fontsize=6.0, color=MUTED, va="top")
        axT.text(0.34, y, a, fontsize=6.0, color=INK, va="top", ha="center")
        axT.text(0.78, y, b, fontsize=6.0, color=INK, va="top", ha="center")
        axT.plot([0, 1], [y - 0.045, y - 0.045], color="#E4E6E8", lw=0.5)

    fig.text(0.035, 0.955, "$\\bf{j1}$   What severity measures, and how energy distance "
                           "differs", fontsize=8.6, ha="left", va="bottom")
    footnote(axT, "Point clouds are SCHEMATIC; every number quoted is read from "
                  "`results/`. Both statistics are computed from the same cells and both "
                  "are direction-free magnitudes — neither says which genes moved or which "
                  "way, which is why the directed E2F/DREAM panel (block g) finds in 30 "
                  "genes what severity misses in 9,099. Severity is the narrower "
                  "question, and on this dataset that is the point: energy distance's "
                  "extra sensitivity to spread is what the control-fitted projection "
                  "exploited.", y=-0.10, width=128)
    save_panel(fig, "PANEL_j1_severity_vs_edistance_schematic")

# ================= k — was it spread, or was it means? The decomposition, per gRNA
def panel_k() -> None:
    """k1, k2: which term of the energy distance made it positive.

    E = 2*mean||x-y|| - mean||x-x'|| - mean||y-y'||, so E can rise either because the two
    clouds move apart (a mean shift, which shows up in mean||x-y|| and in the centroid
    distance) or because the gRNA's own cloud is tighter than the control cloud (a spread
    difference, which shows up in mean||x-x'||). These panels separate the two, per gRNA,
    with the nine NT gRNAs of each arm given identical treatment.
    """
    f = R / "edist_decomposition.csv"
    if not f.exists():
        raise SystemExit("missing edist_decomposition.csv — run "
                         "`python -m analysis.edist_decompose`")
    D = pd.read_csv(f)
    D["disp_ratio"] = D.d_xx_heldout / D.d_yy_heldout
    # centroid distance relative to that arm's NT median: 1.0 = no further apart than a
    # non-targeting gRNA already is, which is the only reference that controls for n
    nt_cent = (D[D.is_control].groupby("arm").centroid_heldout.median())
    D["cent_ratio"] = D.centroid_heldout / D.arm.map(nt_cent)

    TGT, NTC = "#eb6834", "#5F6368"
    rng = np.random.default_rng(0)

    # ------------------------------------------ k1: spread and centroid, side by side
    style.apply(8.0)
    fig, axes = plt.subplots(1, 2, figsize=_sz(5.6, 2.9))
    fig.subplots_adjust(wspace=0.30, top=0.80)
    specs = [(axes[0], "disp_ratio",
              "within-gRNA spread\n$\\overline{\\|x-x'\\|}$ / $\\overline{\\|y-y'\\|}$",
              "1.0 = same spread as the control pool"),
             (axes[1], "cent_ratio",
              "centroid distance from control\nrelative to that arm's NT median",
              "1.0 = no further out than an NT gRNA")]
    for ax, col, ylab, ref in specs:
        ax.axhline(1.0, color=INK, lw=0.8, ls="--", zorder=1)
        for i, arm in enumerate(paths.ARMS):
            g = D[D.arm == arm]
            for j, (lab, sub, c) in enumerate([("targeting", g[~g.is_control], TGT),
                                               ("NT", g[g.is_control], NTC)]):
                x = i + (-0.17 if j == 0 else 0.17)
                v = sub[col]
                ax.scatter(x + rng.uniform(-0.09, 0.09, len(v)), v, s=7,
                           facecolor=c if j == 0 else "none",
                           edgecolor="white" if j == 0 else c,
                           lw=0.2 if j == 0 else 0.6, alpha=0.85, zorder=3)
                ax.plot([x - 0.14, x + 0.14], [v.median()] * 2, color=INK, lw=1.4,
                        zorder=4)
        ax.set_xticks(range(len(paths.ARMS)))
        ax.set_xticklabels([a.replace("/", "\n") for a in paths.ARMS], fontsize=5.8)
        ax.set_xlim(-0.55, len(paths.ARMS) - 0.45)
        ax.set_ylabel(ylab, fontsize=6.8)
        ax.text(0.5, 1.015, ref, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=5.8, color=MUTED)
        style.despine(ax)
    hs = [Line2D([], [], marker="o", ls="", ms=4.0, mfc=TGT, mec="none",
                 label="targeting gRNAs"),
          Line2D([], [], marker="o", ls="", ms=4.0, mfc="none", mec=NTC,
                 markeredgewidth=0.7, label="the nine NT gRNAs")]
    axes[0].legend(handles=hs, loc="lower left", fontsize=5.6, handletextpad=0.3,
                   borderpad=0.25, labelspacing=0.35)
    fig.text(0.0, 1.03, "$\\bf{k1}$   The energy-distance excess is a spread difference, "
                        "not a mean shift", fontsize=8.2, ha="left", va="bottom")
    dmed = D[~D.is_control].groupby("arm").disp_ratio.median()
    cmed = D[~D.is_control].groupby("arm").cent_ratio.median()
    footnote(axes[0],
             f"Every gRNA scored against a control pool with its OWN cells removed. LEFT: "
             f"targeting gRNAs' cells are "
             f"{100 * (1 - dmed.max()):.0f}–{100 * (1 - dmed.min()):.0f}% TIGHTER than the "
             f"control pool (median {dmed.min():.2f}–{dmed.max():.2f}) while the nine NT "
             f"gRNAs sit on 1.00 — and E reads spread through its "
             f"$-\\overline{{\\|x-x'\\|}}$ term. RIGHT: targeting centroids are "
             f"{cmed.min():.2f}–{cmed.max():.2f}$\\times$ the NT centroid distance, i.e. "
             f"CLOSER to the control centroid than a non-targeting gRNA is, so there is no "
             f"mean shift to find. NT gRNAs match the pool's spread because they ARE "
             f"control cells and sit inside the PCA fit; targeting gRNAs do not — which is "
             f"what panel h6 tests and what makes the whole effect vanish under a "
             f"symmetric fit.", y=-0.34, width=118)
    save_panel(fig, "PANEL_k1_edistance_spread_not_means")

    # ---------------- k2: the two components of E, with NT beside targeting in both
    # This is the panel that settles "spread, not means": the LOCATION component is the
    # same for targeting and non-targeting gRNAs -- it is sampling noise in the centroid,
    # and the NT gRNAs actually carry slightly more of it. Only the SHAPE component
    # differs, and for NT gRNAs it is negative.
    style.apply(8.0)
    fig, axes = plt.subplots(1, 2, figsize=_sz(5.4, 2.8), sharex=True)
    fig.subplots_adjust(wspace=0.28, top=0.80)
    for ax, col, ttl in [(axes[0], "e_location_heldout", "location component\n"
                          "(centroids apart)"),
                         (axes[1], "e_shape_heldout", "shape component\n"
                          "(spread / shape differs)")]:
        ax.axhline(0, color=INK, lw=0.6, zorder=1)
        for i, arm in enumerate(paths.ARMS):
            g = D[D.arm == arm]
            for j, (sub, c) in enumerate([(g[~g.is_control], TGT),
                                          (g[g.is_control], NTC)]):
                x = i + (-0.18 if j == 0 else 0.18)
                ax.bar(x, sub[col].median(), 0.32,
                       color=c if j == 0 else "none", edgecolor=c,
                       lw=0.8, zorder=3,
                       label=("targeting" if j == 0 else "the nine NT gRNAs")
                       if i == 0 else None)
        ax.set_title(ttl, fontsize=6.6, pad=3)
        ax.set_xticks(range(len(paths.ARMS)))
        ax.set_xticklabels([a.replace("/", "\n") for a in paths.ARMS], fontsize=5.8)
        ax.set_xlim(-0.55, len(paths.ARMS) - 0.45)
        ax.set_ylabel("median energy distance", fontsize=6.8)
        style.despine(ax)
    axes[0].legend(fontsize=5.6, loc="upper left", handlelength=1.0, handletextpad=0.4,
                   borderpad=0.25, labelspacing=0.35)
    fig.text(0.0, 1.03, "$\\bf{k2}$   Only the shape component differs — the location "
                        "component is matched", fontsize=8.2, ha="left", va="bottom")
    L = D.groupby(["arm", "is_control"])[["e_location_heldout", "e_shape_heldout"]].median()
    lo_t = [L.loc[(a, False), "e_location_heldout"] for a in paths.ARMS]
    lo_n = [L.loc[(a, True), "e_location_heldout"] for a in paths.ARMS]
    sh_t = [L.loc[(a, False), "e_shape_heldout"] for a in paths.ARMS]
    sh_n = [L.loc[(a, True), "e_shape_heldout"] for a in paths.ARMS]
    footnote(axes[0],
             f"E is recomputed after translating the gRNA's cloud onto the control "
             f"centroid; what survives is SHAPE and the remainder is LOCATION. "
             f"Location is {min(lo_t):+.2f} to {max(lo_t):+.2f} for targeting gRNAs and "
             f"{min(lo_n):+.2f} to {max(lo_n):+.2f} for the NT gRNAs — the same, and if "
             f"anything larger for NT, which is what a centroid offset driven purely by "
             f"sampling noise looks like. Shape is {min(sh_t):+.2f} to {max(sh_t):+.2f} "
             f"for targeting and {min(sh_n):+.2f} to {max(sh_n):+.2f} for NT: the sign "
             f"flips. So the entire targeting-vs-control energy distance lives in the "
             f"spread term, and panel h6 shows that term follows the PCA fit rather than "
             f"the gRNA.", y=-0.34, width=118)
    save_panel(fig, "PANEL_k2_edistance_decomposition")

# =============== l — centroid distance as a simpler stand-in for severity: does it work?
def panel_l() -> None:
    """l1, l2: centroid distance against severity, and what each statistic owes to n.

    Both estimate the magnitude of the same mean shift, so they should agree -- and they
    do. The question these panels answer is what is lost by taking the simpler one.
    """
    f = R / "edist_decomposition.csv"
    if not f.exists():
        raise SystemExit("missing edist_decomposition.csv — run "
                         "`python -m analysis.edist_decompose`")
    sev = pd.read_csv(R / "severity.csv")
    dec = pd.read_csv(R / "edist_decomposition.csv")
    m = sev[["arm", "guide_identity", "n_cells", "D", "naive", "severity",
             "is_control", "Symbol"]].merge(
        dec[["arm", "guide_identity", "centroid_heldout"]],
        on=["arm", "guide_identity"], how="inner")

    # ------------------------------------------------- l1: the two against each other
    style.apply(8.0)
    fig, axes = plt.subplots(2, 2, figsize=_sz(5.4, 4.8))
    fig.subplots_adjust(hspace=0.46, wspace=0.30, top=0.89)
    for ax, arm in zip(axes.ravel(), paths.ARMS):
        g = m[m.arm == arm]
        t, c = g[~g.is_control], g[g.is_control]
        ax.scatter(t.severity, t.centroid_heldout, s=12,
                   c=[GC[s] for s in t.Symbol], lw=0.25, edgecolor="white", alpha=0.9,
                   zorder=3)
        ax.scatter(c.severity, c.centroid_heldout, s=16, facecolor="none",
                   edgecolor="#5F6368", lw=0.7, zorder=4)
        r = stats.spearmanr(t.severity, t.centroid_heldout)
        ax.set_title(f"{arm}   ·   Spearman $\\rho$ = {r.statistic:+.2f} "
                     f"(p = {r.pvalue:.1g})", fontsize=6.2, pad=3)
        ax.tick_params(labelsize=5.8)
        style.despine(ax)
    for ax in axes[1, :]:
        ax.set_xlabel("severity", fontsize=6.8)
    for ax in axes[:, 0]:
        ax.set_ylabel("centroid distance from control\n(30 PCs)", fontsize=6.8)
    hs = [Line2D([], [], marker="o", ls="", ms=4.0, mfc=MUTED, mec="none",
                 label="targeting gRNAs (coloured by CDK)"),
          Line2D([], [], marker="o", ls="", ms=4.2, mfc="none", mec="#5F6368",
                 markeredgewidth=0.7, label="the nine NT gRNAs")]
    fig.legend(handles=hs, loc="upper right", ncol=2, frameon=False, fontsize=5.8,
               handletextpad=0.3, columnspacing=1.0, bbox_to_anchor=(0.99, 1.005))
    fig.text(0.0, 1.015, "$\\bf{l1}$   Centroid distance tracks severity — they estimate "
                         "the same mean shift", fontsize=8.2, ha="left", va="bottom")
    rr = {a: stats.spearmanr(m[(m.arm == a) & (~m.is_control)].severity,
                             m[(m.arm == a) & (~m.is_control)].centroid_heldout).statistic
          for a in paths.ARMS}
    footnote(axes[1, 0],
             f"Both statistics are magnitudes of the same quantity — the shift in the mean "
             f"transcriptome — so the agreement is expected and is not independent "
             f"validation. It is much tighter in T47D "
             f"($\\rho$ = {rr['T47D/ABE']:+.2f}, {rr['T47D/CBE']:+.2f}) than in A549 "
             f"($\\rho$ = {rr['A549/ABE']:+.2f}, {rr['A549/CBE']:+.2f}), where ~500 cells "
             f"per gRNA puts both closer to their noise floors relative to the signal. "
             f"Note the NT gRNAs sit INSIDE the targeting cloud on both axes: neither "
             f"statistic separates targeting from control, which is the same bounded "
             f"negative the rest of this set reaches.", y=-0.40, width=118)
    save_panel(fig, "PANEL_l1_centroid_vs_severity")

    # ------------------------------- l2: what each candidate statistic owes to cell count
    cands = [("naive $\\|\\hat m\\|^2$\n(9,099 genes)", "naive"),
             ("centroid distance\n(30 PCs)", "centroid_heldout"),
             ("cross-fitted $\\hat D$\n(9,099 genes)", "D"),
             ("severity\n(standardised $\\hat D$)", "severity")]
    style.apply(8.0)
    fig, ax = plt.subplots(figsize=_sz(4.2, 2.7))
    w = 0.2
    for j, arm in enumerate(paths.ARMS):
        t = m[(m.arm == arm) & (~m.is_control)]
        inv = 1 / t.n_cells
        vals = [abs(stats.spearmanr(t[c], inv).statistic) for _, c in cands]
        ax.bar(np.arange(len(cands)) + (j - 1.5) * w, vals, w,
               color=plt.get_cmap("Blues")(0.35 + 0.18 * j), edgecolor="white", lw=0.4,
               label=arm, zorder=3)
    ax.axhline(0.3, color=style.HIT, ls="--", lw=0.8, zorder=2)
    ax.text(len(cands) - 0.45, 0.32, "a statistic above this line is\nranking gRNAs partly "
                                     "by their cell count", fontsize=5.4, color=style.HIT,
            ha="right", va="bottom")
    ax.set_xticks(range(len(cands)))
    ax.set_xticklabels([n for n, _ in cands], fontsize=5.8)
    ax.set_ylabel("|Spearman $\\rho$| with 1 / $n_{\\mathrm{cells}}$", fontsize=6.8)
    ax.set_ylim(0, 1.08)
    ax.legend(fontsize=5.6, loc="upper center", ncol=2, handlelength=1.0,
              handletextpad=0.4, borderpad=0.25, labelspacing=0.3, columnspacing=0.8)
    style.despine(ax)
    style.panel(ax, "l2", "What each statistic owes to cell count", pad=4, size=8.2)
    footnote(ax, "Cells per gRNA spans 22–836, so any statistic that carries a sampling "
                 "floor will rank gRNAs by n rather than by biology. The naive squared "
                 "norm is almost PURE cell count ($\\rho$ = 0.79–0.98) — that is what the "
                 "cross-fit exists to remove. Centroid distance escapes most of it "
                 "because 30 PCs is a far smaller space than 9,099 genes, but not all of "
                 "it: 0.24–0.28 in the two T47D arms, where n is smallest. So centroid "
                 "distance IS a usable simpler statistic, but it still has to be "
                 "standardised against a size-matched null before gRNAs are ranked on it, "
                 "and it can only see a shift that lies inside the leading 30 PCs.",
             y=-0.32, width=112)
    save_panel(fig, "PANEL_l2_statistic_vs_cellcount")


BLOCKS = {"a": panel_a, "b": panel_b, "c": panel_c, "d": panel_d, "e": panel_e,
          "f": panel_f, "g": panel_g, "h": panel_h, "i": panel_i, "j": panel_j, "k": panel_k, "l": panel_l}

if __name__ == "__main__":
    args = sys.argv[1:]
    clean = [a for a in args if a == "clean" or a.startswith("clean=")]
    if clean:
        CLEAN = True
        if "=" in clean[0]:
            FIG_SCALE = float(clean[0].split("=")[1])
        args = [a for a in args if a not in clean]
        print(f"compact mode: no footnotes, letters only, canvas x{FIG_SCALE}")
    want = [a for a in args if a in BLOCKS] or list(BLOCKS)
    for b in want:
        print(f"[{b}]")
        BLOCKS[b]()
    print(f"\npanels in {paths.rel(OUT)}")
