"""KLIFS pocket resistance alignment (supplementary 15c).

From `09_analysis_and_viz/CDK_resistance_patterns.ipynb`
cells 11-14. The 85 numbered KLIFS pocket positions as a sequence alignment,
one row per CDK, each cell shaded by how many drugs gave resistance at that
residue.

Everything derives from two shipped inputs -- `resistance_calls_canonical.csv`
for the analysis set and `KLIFS_compiled.csv` for the pocket numbering -- so
nothing re-runs the screen. The three shading states are distinct and worth
keeping straight:

* coloured  -- resistance here, in that many drugs
* `#f2f2f2` -- **assayed and negative**: a guide editing this residue passed
  filters in this CDK and gave no resistance
* white     -- not assayed in this CDK at all

Without the middle state an empty cell would be ambiguous, and the recurrence
markers below the alignment (hit in >= 2 or >= 3 CDKs) would be unreadable.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from .. import io, paths
from .structures import Renders


def _render_dir(name: str) -> Path:
    """Where a ChimeraX render subdirectory lives.

    `figures/chimerax/` is where `scripts/chimerax/render.py` writes and is not
    committed -- a full render is 119 MB of supersampled PNG. The copies that
    ship are under `data/reference/structure_renders/`, so a clone without
    ChimeraX can still assemble these figures. A fresh render wins.
    """
    for root in (paths.FIGURES / "chimerax", paths.REFERENCE / "structure_renders"):
        if (root / name).is_dir():
            return root / name
    return paths.FIGURES / "chimerax" / name

KLIFS = paths.REFERENCE / "KLIFS"
FDR_CUT = 0.01
SET = "core7"

GENES_ALL = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
#: CDK8/CDK19 are dropped from core7: their own inhibitors give almost no signal
EXCLUDE = {"core7": {"CDK8", "CDK19"}, "all9": set()}

#: subdomain blocks, matching the annotation used across the project
KLIFS_ELEMENTS = [
    ((1, 3), "I", "goldenrod"), ((4, 9), "g.l", "green"), ((10, 13), "II", "goldenrod"),
    ((14, 19), "III", "goldenrod"), ((20, 30), "aC", "red"), ((31, 37), "b.l", "green"),
    ((38, 41), "IV", "goldenrod"), ((42, 44), "V", "goldenrod"), ((45, 45), "GK", "orange"),
    ((46, 48), "Hinge", "purple"), ((49, 52), "Linker", "blue"), ((53, 59), "aD", "red"),
    ((60, 64), "aE", "red"), ((65, 67), "VI", "goldenrod"), ((68, 75), "c.l", "orange"),
    ((76, 78), "VII", "goldenrod"), ((79, 79), "VIII", "goldenrod"),
    ((80, 83), "xDFG", "blue"), ((84, 85), "a.l", "blue"),
]


def genes(aset=SET):
    return [g for g in GENES_ALL if g not in EXCLUDE[aset]]


def klifs_table(aset=SET):
    """KLIFS_compiled with a 1-based pocket position and its subdomain."""
    kl = io.read(KLIFS / "KLIFS_compiled.csv")
    kl["pos"] = np.arange(1, len(kl) + 1)
    kl["element"] = kl.KLIFS_idx.str.split(":").str[0]
    return kl


def matrices(aset=SET):
    """-> (mat_drugs, assayed, n_cdk, KL): the gene x position state of the pocket."""
    gs = genes(aset)
    kl = klifs_table(aset)
    # the shipped table is already one row per *enriched* guide x residue x drug
    # test, and already carries the KLIFS position, so the notebook's `enriched`
    # filter and its merge onto kmap are both baked in
    calls = io.read(paths.LFC / "resistance_patterns" / aset / "resistance_calls_canonical.csv")
    ck = calls.rename(columns={"KLIFS_pos": "pos"})
    ck = ck[ck["pos"].notna()].astype({"pos": int})
    hits = (ck.groupby(["Gene", "pos"])
            .agg(n_drugs=("Condition", "nunique")).reset_index())

    mat_drugs = pd.DataFrame(0, index=gs, columns=kl.pos)
    for r in hits.itertuples():
        mat_drugs.loc[r.Gene, r.pos] = r.n_drugs
    assayed = pd.DataFrame(False, index=gs, columns=kl.pos)
    for (g, p), _ in ck.groupby(["Gene", "pos"]):
        assayed.loc[g, p] = True
    return mat_drugs, assayed, (mat_drugs > 0).sum(axis=0), kl


def source(aset=SET) -> pd.DataFrame:
    """Long gene x position table: drugs with resistance, and whether assayed."""
    mat_drugs, assayed, n_cdk, kl = matrices(aset)
    long = (mat_drugs.stack().rename("n_drugs").reset_index()
            .rename(columns={"level_0": "Gene", "level_1": "pos"}))
    long["assayed"] = [assayed.loc[r.Gene, r.pos] for r in long.itertuples()]
    long["n_cdks_at_position"] = long["pos"].map(n_cdk)
    meta = kl[["pos", "KLIFS_idx", "element"]]
    return long.merge(meta, on="pos")


#: the OncoKB comparison kinases, in the order the published panel lists them
ONCOKB_GENES = ["ABL1", "ALK", "BTK", "EGFR", "KIT", "MET",
                "NTRK1", "NTRK3", "PDGFRA", "ROS1"]

#: the published lower panel's scale runs 1..5 distinct drugs
ONCOKB_MAX_DRUGS = 5


def oncokb_matrix():
    """-> (n_drugs per kinase x position, the 85-position sequences).

    The sequences come from each kinase's own `KLIFS2UniProtIdx` map in
    `data/reference/KLIFS/oncokb_klifs/`, the same route
    `klifs_resistance.annotated` uses to place a mutation, so a residue and the
    mutation reported at it cannot disagree.
    """
    from .klifs_resistance import klifs_index

    seqs = klifs_index(ONCOKB_GENES)
    kl = klifs_table()
    idx2pos = dict(zip(kl.KLIFS_idx, kl.pos))
    mat = pd.DataFrame(0, index=ONCOKB_GENES, columns=kl.pos)
    o = oncokb_by_position()
    for r in o.itertuples():
        pos = idx2pos.get(r.KLIFS_Index)
        if pos is not None:
            mat.loc[r.Gene, pos] = r.n_drugs
    return mat, seqs


def _write_fa(path, names, seqs):
    with open(path, "w") as fh:
        for n, sq in zip(names, seqs):
            fh.write(f">{n}\n{sq}\n")
    return path


def supp_15c(outdir=None, aset=SET, two_panel=True, **kw):
    """The published panel: the CDK screen above, the clinical kinases below.

    -> (figure, the long gene x position table).

    **The code that drew the published two-panel version was never found** --
    `CDK_resistance_patterns.ipynb` cell 14, the only cell that writes
    `KLIFS_resistance_MSA.pdf`, draws the CDK block alone under a different
    title. The lower block is rebuilt here from the repository's own inputs:
    the ten kinases' pocket sequences from `oncokb_klifs/`, and the drug counts
    from `oncokb_by_position()`, which already matched the published footnote
    (30 of 104 mutations outside the 85 numbered positions) before anything was
    drawn. The layout follows the published figure -- row order, the two colour
    scales, the recurrence markers, the diamonds on positions that appear in
    both blocks -- so it is a reconstruction of a lost figure, not a recovery of
    lost code. `two_panel=False` gives cell 14's single block verbatim.

    Both blocks are drawn by pymsaviz's own renderer onto axes of one figure,
    rather than by a second alignment painter written here, so the two cannot
    drift apart. That uses `MsaViz._plot_msa`, which is private: if a future
    pymsaviz drops it, this raises instead of quietly drawing one block.
    """
    from pymsaviz import MsaViz

    mat_drugs, assayed, n_cdk, kl = matrices(aset)
    gs = genes(aset)
    max_drugs = int(mat_drugs.values.max())
    # light-to-medium reds only, so the residue letters stay readable on the fill
    shade = {n: to_hex(plt.get_cmap("Reds")(0.18 + 0.20 * (n - 1)))
             for n in range(1, max_drugs + 1)}

    def colour(row_pos, col_pos, seq_char, msa):
        gene, pos = gs[row_pos], col_pos + 1       # pymsaviz is 0-indexed
        if seq_char == "-":
            return "#ffffff"
        n = mat_drugs.loc[gene, pos]
        if n > 0:
            return shade[n]
        return "#f2f2f2" if assayed.loc[gene, pos] else "#ffffff"

    # pymsaviz reads an alignment from disk, so the analysis-set rows are
    # written out first; it is derived, not an input
    fa_dir = paths.PROCESSED / "resistance_patterns"
    fa_dir.mkdir(parents=True, exist_ok=True)
    fa = _write_fa(fa_dir / f"KLIFS_analysis_{aset}.fa", gs,
                   ["".join(kl[f"{g}_seq"].fillna("-")) for g in gs])

    mv = MsaViz(str(fa), wrap_length=85, show_grid=True, show_consensus=True,
                show_count=False)
    mv.set_plot_params(ticks_interval=5, x_unit_size=0.165, y_unit_size=0.24)
    mv.set_custom_color_func(colour)
    for (lo, hi), label, col in KLIFS_ELEMENTS:
        mv.add_text_annotation((lo, hi), label, text_color=col, range_color=col)
    mv.add_markers(list(n_cdk.index[n_cdk >= 2]), marker="v", color="tab:grey", size=6)
    mv.add_markers(list(n_cdk.index[n_cdk >= 3]), marker="v", color="black", size=8)

    handles = [Patch(facecolor=shade[n], edgecolor="lightgrey",
                     label=f'{n} drug{"s" if n > 1 else ""}')
               for n in range(1, max_drugs + 1)]
    handles += [
        Patch(facecolor="#f2f2f2", edgecolor="lightgrey", label="assayed, no resistance"),
        Patch(facecolor="#ffffff", edgecolor="lightgrey", label="not assayed"),
        Line2D([], [], marker="v", color="none", markerfacecolor="tab:grey",
               markeredgecolor="none", markersize=7, label="hit in >=2 CDKs"),
        Line2D([], [], marker="v", color="none", markerfacecolor="black",
               markeredgecolor="none", markersize=8, label="hit in >=3 CDKs"),
    ]

    if not two_panel:
        fig = mv.plotfig()
        ax = fig.axes[0]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.5),
                  ncols=4, fontsize=9, frameon=False,
                  title="Resistance at residue", title_fontsize=10)
        fig.text(0.5, 1.13,
                 "KLIFS active-site resistance mutations across the CDKs  "
                 f"(FDR < {FDR_CUT}, LFC > 0)", ha="center", va="bottom", fontsize=13)
        if outdir is not None:
            fig.savefig(str(outdir / "Supp_15c_KLIFS_resistance_MSA.pdf"),
                        transparent=True, bbox_inches="tight")
        return fig, source(aset)

    # ---------------------------------------------------- the lower block
    omat, oseqs = oncokb_matrix()
    opurple = {n: to_hex(plt.get_cmap("Purples")(0.22 + 0.15 * (n - 1)))
               for n in range(1, ONCOKB_MAX_DRUGS + 1)}

    def ocolour(row_pos, col_pos, seq_char, msa):
        gene, pos = ONCOKB_GENES[row_pos], col_pos + 1
        n = int(omat.loc[gene, pos])
        return opurple[n] if n > 0 else "#f7f7f7"

    ofa = _write_fa(fa_dir / "KLIFS_oncokb_kinases.fa", ONCOKB_GENES,
                    ["".join(oseqs[f"{g}_seq"]) for g in ONCOKB_GENES])
    omv = MsaViz(str(ofa), wrap_length=85, show_grid=True, show_consensus=False,
                 show_count=False)
    omv.set_plot_params(ticks_interval=5, x_unit_size=0.165, y_unit_size=0.24)
    omv.set_custom_color_func(ocolour)
    # a position carrying resistance in both blocks -- the comparison the
    # published panel is making, so it is marked rather than left to the eye
    both = sorted(set(n_cdk.index[n_cdk > 0]) & set(omat.columns[(omat > 0).any()]))
    omv.add_markers(both, marker="D", color="black", size=5)

    for mv_obj in (mv, omv):
        if not hasattr(mv_obj, "_plot_msa"):            # pymsaviz changed
            raise RuntimeError(
                "this pymsaviz has no MsaViz._plot_msa; the two-panel layout "
                "composes both blocks onto one figure through it")

    u = 0.24                                            # y_unit_size, as set above
    heights = [len(gs) * u, u * 1.5, u * 2.0, u * 2.6, len(ONCOKB_GENES) * u]
    fig = plt.figure(figsize=(85 * 0.165, sum(heights) + 2.3))
    g = GridSpec(nrows=len(heights), ncols=1, height_ratios=heights)
    # the top block carries the subdomain annotation ABOVE its axes, so the
    # title needs headroom over that rather than over the alignment
    g.update(left=0, right=1, bottom=0.26, top=0.885, hspace=0)
    ax_cdk = fig.add_subplot(g[0])
    ax_con = fig.add_subplot(g[2])
    ax_onc = fig.add_subplot(g[4])
    mv._plot_msa(ax_cdk)
    mv._plot_consensus(ax_con)
    omv._plot_msa(ax_onc)

    ohandles = [Patch(facecolor=opurple[n], edgecolor="lightgrey",
                      label=f'{n} drug{"s" if n > 1 else ""}')
                for n in range(1, ONCOKB_MAX_DRUGS + 1)]
    ohandles += [
        Patch(facecolor="#f7f7f7", edgecolor="lightgrey",
              label="no reported resistance"),
        Line2D([], [], marker="D", color="none", markerfacecolor="black",
               markeredgecolor="none", markersize=6, label="position in both panels"),
    ]
    leg1 = fig.legend(handles=handles, loc="upper left",
                      bbox_to_anchor=(0.03, 0.175), ncols=7, fontsize=8.5,
                      frameon=False, title_fontsize=9,
                      title="CDK base-editing screen (top): distinct drugs with "
                            "resistance at the residue")
    leg1._legend_box.align = "left"
    leg2 = fig.legend(handles=ohandles, loc="upper left",
                      bbox_to_anchor=(0.03, 0.095), ncols=7, fontsize=8.5,
                      frameon=False, title_fontsize=9,
                      title="Other kinases (bottom), OncoKB R1/R2: distinct drugs "
                            "with reported clinical resistance")
    leg2._legend_box.align = "left"
    # the footnote's two numbers are counted here, not copied from the published
    # panel: the denominator is every OncoKB R1/R2 mutation in these kinases and
    # the numerator is those whose residue is not one of the 85 pocket positions
    from .klifs_resistance import annotated as _annotated
    _a = _annotated()
    n_tot, n_out = len(_a), int((~_a["KLIFS Catalytic Site"]).sum())
    fig.text(0.03, 0.015,
             f"{n_out} of {n_tot} OncoKB mutations fall outside the "
             "85 numbered pocket positions and have no column here.",
             fontsize=7.5, color="dimgrey", ha="left", va="bottom")
    fig.text(0.5, 0.955,
             "KLIFS active-site resistance: CDK base-editing screen "
             f"(FDR < {FDR_CUT}, LFC > 0) vs clinical resistance in other kinases",
             ha="center", va="bottom", fontsize=12)
    if outdir is not None:
        fig.savefig(str(outdir / "Supp_15c_KLIFS_resistance_MSA.pdf"),
                    transparent=True, bbox_inches="tight")
    return fig, source(aset)


def oncokb_by_position() -> pd.DataFrame:
    """The published figure's lower panel: OncoKB kinases on the same 85 positions.

    One row per (kinase, KLIFS position) with the number of distinct drugs that
    have a reported R1/R2 resistance mutation there. Verified against the
    published footnote: 30 of the 104 OncoKB mutations fall outside the 85
    numbered positions and have no column.
    """
    from .klifs_resistance import annotated

    a = annotated()
    inside = a[a["KLIFS Catalytic Site"]].copy()
    inside["KLIFS_Index"] = inside["KLIFS_Index"].astype(str)
    g = (inside.groupby(["Gene", "KLIFS_Index"])
         .agg(n_drugs=("Drug", "nunique"),
              mutations=("Mutation", lambda s: ",".join(sorted(set(s)))),
              drugs=("Drug", lambda s: ",".join(sorted(set(s)))))
         .reset_index())
    g.attrs["n_outside"] = int((~a["KLIFS Catalytic Site"]).sum())
    return g


def source_both(aset=SET) -> dict:
    """Both halves of the published panel, as two sheets."""
    return {"cdk_screen": source(aset), "oncokb_kinases": oncokb_by_position()}


def supp_15d(outdir=None, aset=SET, **kw):
    """Resistance enrichment along the CDK alignment, with peak detection.

    This is `canonical_resistance_figures.run()`, which draws the published
    panel and writes it to `figures/resistance_patterns/<set>/`. The regenerated
    PDF is pixel-identical to the published one.
    """
    from ..essentiality import canonical_resistance_figures as crf

    fig = crf.run(aset)
    return fig, peak_sites(aset)


def peak_sites(aset=SET) -> pd.DataFrame:
    """The detected peak columns and the variants under them -- 15d's data."""
    return io.read(paths.LFC / "resistance_patterns" / aset / "MSA_peak_sites.csv")


def peak_labels(aset=SET) -> pd.DataFrame:
    """One row per detected peak, with the variants it collects -- 15e's data.

    This is the table the leader-line blocks are written from: `P1  g.l:5
    CDK2-12  22 variants / 4 CDKs` and the gene lines under it all come from
    one row of this file.
    """
    return io.read(paths.LFC / "resistance_patterns" / aset
                   / "MSA_peak_sites_labels.csv")


def supp_15e(outdir=None, aset=SET, **kw):
    """-> ([image paths], the peak table). 15e is a ChimeraX render.

    Like the main-figure structure panels this is not drawn with matplotlib:
    the label positions are projected from the live ChimeraX camera, so the
    panel is produced by running the application. Two steps make it, and
    `scripts/chimerax/render.py` runs both:

        generate_peaks_cxc.py --set core7      # the spheres
        label_consensus.py core7 full peaks ...  # the leader-line blocks

    The regenerated SVG matches the published one label for label, with every
    coordinate inside half a pixel of a 1400-wide frame.
    """
    import shutil

    outputs = _render_dir(f"cxc_{aset}_peaks")
    found = Renders(outputs / n for n in ("peaks_leader_variants.svg",
                                          "peaks_leader_variants.png")
                    if (outputs / n).exists())
    if not found:
        raise FileNotFoundError(
            f"figure S15e needs peaks_leader_variants.* under {outputs}, which "
            "has not been rendered yet. Run:  python3 scripts/chimerax/render.py")
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for f in found:
            shutil.copy2(f, outdir / f.name)
    return found, peak_labels(aset)


# ------------------------------------------------------------- figure 6 ---

def consensus_sites(aset=SET) -> pd.DataFrame:
    """The alignment columns figure 6 marks, with their tier and variants."""
    return io.read(paths.LFC / "resistance_patterns" / aset / "consensus_sites.csv")


def homologous_resistance(aset=SET) -> pd.DataFrame:
    """Each consensus site's aligned position in kinases outside the CDKs.

    This is the table the mapping row for figure 6 asks after: the figure adds
    text naming the matching aligned positions in other kinases, and
    `homologous_clinical` is where that text comes from -- the clinical
    variants reported at the same KLIFS position in ABL1, ALK, BTK, EGFR, KIT,
    NTRK and ROS1, with their source (OncoKB or CIViC).
    `matching_sites_cdk_vs_clinical.csv` beside it is the same join the other
    way round, one row per KLIFS position rather than per site.
    """
    return io.read(paths.LFC / "resistance_patterns" / aset
                   / "consensus_sites_homologous_resistance.csv")


def matching_clinical_sites(aset=SET) -> pd.DataFrame:
    """The 11 KLIFS positions where a CDK hit meets a clinical variant."""
    return io.read(paths.LFC / "resistance_patterns" / aset
                   / "matching_sites_cdk_vs_clinical.csv")


def figure_6(outdir=None, aset="core7_consensus_strong", **kw):
    """-> ([image paths], the sites with their homologous clinical variants).

    Figure 6 is a whole figure rather than a lettered panel, and like
    supplementary 15e it is a ChimeraX render with leader-line labels rather
    than a matplotlib plot. `scripts/chimerax/render.py` makes it:

        generate_consensus_cxc.py --set core7      # the spheres
        label_consensus.py core7 variants ...      # the leader-line labels

    The background raster is copied beside the SVG, because the SVG is labels
    and leader lines over it and renders as labels floating on nothing without.
    """
    import shutil

    outputs = _render_dir(f"cxc_{aset}")
    names = ("consensus_leader_variants.svg", "consensus_leader_variants.png",
             "consensus_variants_plain.png")
    found = Renders(outputs / n for n in names if (outputs / n).exists())
    if not found:
        raise FileNotFoundError(
            f"figure 6 needs consensus_leader_variants.* under {outputs}, "
            "which has not been rendered yet. Run:  "
            "python3 scripts/chimerax/render.py")
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for f in found:
            shutil.copy2(f, outdir / f.name)
    return found, homologous_resistance()


FIGURES = {"S15c": supp_15c, "S15d": supp_15d, "S15e": supp_15e}
PANELS = {"S15c": lambda: source_both(),
          "S15d": lambda: peak_sites(),
          "S15e": lambda: peak_labels()}

#: figure 6 is a main-figure entry, so it is kept out of the supplementary
#: registries above and merged by `source_data.__init__` instead
MAIN_FIGURES = {"6": figure_6}
MAIN_PANELS = {"6": lambda: homologous_resistance()}

#: `scripts/render_figures.py` reads this to name the files a panel produces
PARTS = {"S15e": ["peaks_leader_variants_svg", "peaks_leader_variants"],
         "6": ["consensus_leader_variants_svg", "consensus_leader_variants",
               "consensus_variants_plain"]}
