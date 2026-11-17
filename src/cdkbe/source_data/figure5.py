"""Figure 5 -- which source panel each published panel is, and how to redraw it.

Figure 5 was assembled in Illustrator from panels drawn by two separate bodies of
code, neither of which letters its panels the way the manuscript does:

    * `cdkbe.beatseq.fig_panels`, whose blocks `panel_a` .. `panel_l` each write
      several files named `PANEL_<letter><n>_<slug>`. Its lettering is its own --
      its `panel_a` is the UMAPs, its `panel_i` is the library donut -- so
      `panel_f` is NOT figure 5f.
    * `cdkbe.source_data.rnaseq_t47d_figures`, which draws 5f-5i as one 2x2.

`SPEC` below is the mapping between the two, one entry per published panel. Each
part names the block that draws it, the stem it is written under, and the
evidence that it is the panel that was published. Nothing here is inferred from a
panel's name: the evidence is either an identity the file system can check or a
string in `figures/Main Figures/figure5_v2.ai` that only one candidate carries.

    from cdkbe.source_data import figure5
    fig = figure5.draw("5j")                 # one figure
    figs = figure5.draw("5b")                # {"target": fig, "phase": fig}

HOW THE MAPPING WAS ESTABLISHED. `tools/panel_source_finder.py` matched the four
UMAP rasters at figure level and ranked the rest by text similarity, which left
five panels wrong or unresolved (it put 5e on the severity panel and 5j on a
single-cell-line variant of it). The mapping below was settled instead by
rendering the published `.ai` and testing each candidate against something only
the right one satisfies -- an axis-label string, a tick-label set, or a plot
aspect ratio. `tests/test_figure5.py` re-runs every one of those checks, so a
wrong stem fails the suite rather than quietly drawing the wrong panel.

WHAT IS STILL NOT IN THE REPOSITORY. Figure 5a's experimental schematic -- the
cells, the arrows, the "+7 days" -- was drawn by hand in Illustrator. Only the
two donuts beside it come from code. `draw("5a")` returns the donuts.
"""
from __future__ import annotations

from . import rnaseq_t47d_figures as _rna

#: Published panel -> the parts that make it up, in the order they appear in the
#: figure (top to bottom, then left to right).
#:
#: Each part is (slug, source, evidence):
#:   slug      what the part is, used for the per-part file name
#:   source    ("beatseq", block letter, panel stem) or ("rnaseq", panel key)
#:   evidence  why that stem is this panel and not a near-twin of it
SPEC = {
    "5a": [
        ("donut_CBE", ("beatseq", "i", "PANEL_i2_library_donut_CBE"),
         "the published donut reads 'CBE / 72 gRNAs' with nonsense 17 (24%), "
         "fitness hit 16 (22%), resistance hit 29 (40%), non-targeting 10 (14%); "
         "every one of those numbers is in figure5_v2.ai and in this panel"),
        ("donut_ABE", ("beatseq", "i", "PANEL_i1_library_donut_ABE"),
         "likewise 'ABE / 92 gRNAs' with start loss 11 (12%), fitness hit 10 (11%), "
         "resistance hit 61 (66%), non-targeting 10 (11%)"),
    ],
    "5b": [
        ("umap_target", ("beatseq", "a", "PANEL_a1_umap_target_A549"),
         "pixel-identical raster embedded in figure5_v2.ai; the published header "
         "reads 'A549 (n=18,000 cells)', which is this panel's cell count"),
        ("umap_phase", ("beatseq", "a", "PANEL_a3_umap_phase_A549"),
         "pixel-identical raster; the lower half of the published 5b, carrying "
         "the G1/S/G2M direct labels"),
    ],
    "5c": [
        ("umap_target", ("beatseq", "a", "PANEL_a2_umap_target_T47D"),
         "pixel-identical raster; published header 'T47D (n=17,339 cells)'"),
        ("umap_phase", ("beatseq", "a", "PANEL_a4_umap_phase_T47D"),
         "pixel-identical raster, including the detached minority cloud"),
    ],
    "5d": [
        ("nmd_by_class", ("beatseq", "c", "PANEL_c1s_editing_by_variant_class_combined"),
         "the published panel is ONE axes carrying both cell lines and both "
         "bracket annotations (T47D p = 0.018, A549 p = 0.035). PANEL_c1 draws "
         "the same numbers as two side-by-side axes, one per cell line; the "
         "combined variant is the single-axes one"),
    ],
    "5e": [
        ("editing_vs_cellcycle", ("beatseq", "c", "PANEL_c3_editing_vs_cellcycle_CDK4_T47D"),
         "'percentage points' and 'phase-adjusted' both appear in figure5_v2.ai "
         "and only in this panel; it is the only candidate whose axes are dG1 "
         "against target-transcript log FC"),
    ],
    # 5f-5i are the four axes of one 2x2 figure, which exists in three layouts
    # that differ only in their tick locators. The published one is identified by
    # the tick labels it leaves in figure5_v2.ai -- {80, 160, 240} on the volcano
    # axes and 1.0 at the top of the E2F axis. The square layout would put
    # 60/120/180 there and the wide one 100/200/300; neither set appears.
    "5f": [("cdk4_tpm", ("rnaseq", "5f"),
            "panel a of pub_fig_square_4panel_mid; CDK4 TPM axis topping out at 400")],
    "5g": [("e2f_score", ("rnaseq", "5g"),
            "panel b; 'Hallmark E2F targets (mean z-score, n = 198)' with the "
            "mid layout's 1.0 top tick")],
    "5h": [("volcano_vs_NT", ("rnaseq", "5h"),
            "panel c; '343 up' / '477 down' and the mid layout's 80/160/240 ticks")],
    "5i": [("volcano_vs_D196N", ("rnaseq", "5i"),
            "panel d; '1 up' / '9 down', same gene universe (11,489)")],
    "5j": [
        ("severity_by_CDK", ("beatseq", "d", "PANEL_d12w_severity_by_CDK_stacked"),
         "the published panel is both cell lines stacked, which rules out the "
         "single-line d1/d1s/d2/d2s variants; and 'severity (SD of the "
         "size-matched null)', the y-label every non-stacked variant carries, "
         "does not appear anywhere in figure5_v2.ai. Between the two stacked "
         "variants it is the wide one: the published plot area is 1.73:1 and "
         "d12w's MediaBox is 464.3 x 267.7 pt = 1.734, against d12's 1.164"),
    ],
    "5k": [
        ("severity_vs_drug_T47D", ("beatseq", "e", "PANEL_e2n_severity_vs_drug_lfc_T47D"),
         "the published x-label is 'max drug log2FC'. The e1/e2 variants label it "
         "'maximum drug log2FC in the parent screen' and the y-axis 'severity (SD "
         "of the size-matched null)'; neither string is in figure5_v2.ai, and "
         "'max drug log' is"),
        ("severity_vs_drug_A549", ("beatseq", "e", "PANEL_e1n_severity_vs_drug_lfc_A549"),
         "same axis-label test; T47D is the left panel and A549 the right"),
    ],
}

PANELS = sorted(SPEC)

#: published panel -> the `fig_panels` block that has to be run to draw it
BLOCK = {k: {src[1] for _, src, _ in parts if src[0] == "beatseq"}
         for k, parts in SPEC.items()}


def parts(panel: str) -> list:
    """The ordered part slugs that make up a published panel."""
    return [slug for slug, _, _ in SPEC[panel]]


def draw(panel: str, outdir=None, clean: bool = True):
    """Redraw a published panel. -> a figure, or {slug: figure} if it has parts.

    With `outdir`, each part is also written there as `main_fig<panel>[_<slug>].pdf`
    -- the same names `scripts/render_figures.py` uses, so a notebook that saves
    its panels and a full render cannot disagree about what a file is called.

    Figures come from the code that drew them for the manuscript: the BEAT-seq
    parts are captured out of `fig_panels`, so they are the same objects that
    `figures/panels_clean/` is written from rather than a reimplementation.
    `clean=True` selects the compact variant the manuscript was assembled from
    (descriptive titles stripped to the bare panel letter).

    A block draws several panels and only some are wanted here, so the unwanted
    ones are closed rather than leaked.
    """
    import matplotlib.pyplot as plt

    spec = SPEC[panel]
    out = {}
    beatseq = [(slug, src[2]) for slug, src, _ in spec if src[0] == "beatseq"]
    if beatseq:
        from ..beatseq import fig_panels as _fp
        wanted = {stem for _, stem in beatseq}
        with _fp.capture(clean=clean) as got:
            for block in sorted(BLOCK[panel]):
                getattr(_fp, f"panel_{block}")()
        missing = wanted - set(got)
        if missing:
            raise RuntimeError(
                f"{panel}: block(s) {sorted(BLOCK[panel])} drew {sorted(got)}, "
                f"which does not include {sorted(missing)}")
        for stem, fig in got.items():
            if stem not in wanted:
                plt.close(fig)
        for slug, stem in beatseq:
            out[slug] = got[stem]
    for slug, src, _ in spec:
        if src[0] == "rnaseq":
            out[slug] = _rna.panel(src[1])
    ordered = {slug: out[slug] for slug in parts(panel)}
    if outdir is not None:
        from pathlib import Path
        d = Path(outdir)
        d.mkdir(parents=True, exist_ok=True)
        multi = len(ordered) > 1
        for slug, fig in ordered.items():
            name = f"main_fig{panel}" + (f"_{slug}" if multi else "")
            save(fig, d / f"{name}.pdf", panel, slug)
    return next(iter(ordered.values())) if len(ordered) == 1 else ordered


#: how a figure from each source is written, matching what that source's own
#: code does: the BEAT-seq panels are transparent with a hairline pad, the
#: RNA-seq ones opaque.
SAVE_KW = {
    "beatseq": dict(transparent=True, bbox_inches="tight", pad_inches=0.01),
    "rnaseq": dict(bbox_inches="tight"),
}

#: `pdf.fonttype` is read by the backend when the file is WRITTEN, not when the
#: figure is drawn, so a figure drawn under the right rcParams and saved under
#: the defaults comes out with Type 3 text. Every save here forces it.
_LIVE_TEXT = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}


def source_of(panel: str, slug: str) -> str:
    """"beatseq" or "rnaseq" -- which body of code drew this part."""
    for s, src, _ in SPEC[panel]:
        if s == slug:
            return src[0]
    raise KeyError(f"{panel} has no part {slug!r}; it has {parts(panel)}")


def save(fig, path, panel: str, slug: str | None = None, dpi: int = 400, **kw):
    """Write one part of a published panel, the way its own source writes it."""
    import matplotlib as mpl

    slug = slug or parts(panel)[0]
    opts = dict(SAVE_KW[source_of(panel, slug)], dpi=dpi, **kw)
    with mpl.rc_context(_LIVE_TEXT):
        fig.savefig(path, **opts)
    return path


#: one-line description of each published panel, for the generated figure list
WHAT = {
    "5a": "CROP-seq library composition: one donut per editor, gRNAs by class",
    "5b": "A549 UMAP (18,000 cells), coloured by CDK targeted and by cell-cycle phase",
    "5c": "T47D UMAP (17,339 cells), coloured by CDK targeted and by cell-cycle phase",
    "5d": "Target-transcript level by variant class, both cell lines, testing NMD",
    "5e": "CDK4 in T47D: change in G1 fraction against target-transcript knockdown",
    "5f": "CDK4 transcript level (TPM) across the CBE NT / E56K / D196N / W198* series",
    "5g": "Hallmark E2F target score across the same series",
    "5h": "Differential expression, CDK4 W198* against CBE NT",
    "5i": "Differential expression, CDK4 W198* against D196N (matched background)",
    "5j": "Transcriptional severity of every targeting gRNA, per CDK, both cell lines",
    "5k": "Severity against the parent screen's maximum drug log2FC, per cell line",
}


def _wrap(panel: str):
    f = lambda outdir=None, **kw: draw(panel, outdir=outdir, **kw)    # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    parts_txt = ", ".join(parts(panel))
    f.__doc__ = f"{WHAT[panel]}. Parts: {parts_txt}."
    return f


#: the registry `cdkbe.source_data` merges. One entry per published panel.
FIGURES = {k: _wrap(k) for k in SPEC}
