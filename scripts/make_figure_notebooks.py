#!/usr/bin/env python3
"""Generate notebooks/08_figures/figure_0N.ipynb."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _lib.nbwrite import md, code, write, HEADER

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
_fignum = None

ROOT = Path(__file__).resolve().parents[1]
#: Panels are cited against the notebook they were first computed in. Those
#: notebooks are not shipped -- see docs/REPRODUCE.md -- so the citation is a
#: name and a cell number, not a path into this checkout.
LEG = ""

# panel -> (description, original source, status)
#   "figure"  source data regenerates AND the panel is drawn here
#   "data"    source data regenerates; plotting not ported yet
#   "unported"  neither yet (kept as a status name only)
FIGS = {
 1: ("Library design and screen overview", {
   "e": ("Codon and variant coverage vs sensor-editing cutoff", "09_analysis_and_viz/editing_analysis_v2.ipynb cell 16", "figure"),
   "f": ("Resistance variants by residue, per CDK", "09_analysis_and_viz/fig1_analysis.ipynb cells 4-13", "figure"),
 }),
 2: ("CDK9 and CDK7 resistance", {
   "a": ("Lollipop, KI-CDK9d-32 1,000 nM, subpool 1", "09_analysis_and_viz/lollipop_plots_canonical.ipynb cell 15", "figure"),
   "b": ("CDK9 top variants, per-replicate LFC", "09_analysis_and_viz/logo_plots.ipynb cells 60-62", "figure"),
   "d": ("MTT dose-response, CDK9 L156F", "source_data/04_MTT_data/FACS_dose_response_2024_8_22.ipynb cells 126-127", "unported"),
   "h": ("Lollipop, SY-5609 10 nM, subpool 1", "09_analysis_and_viz/lollipop_plots_canonical.ipynb cell 15", "figure"),
   "i": ("CDK7 SY-5609 resistance hits", "09_analysis_and_viz/CDK7_SY5609_hits.ipynb cells 8/10/14", "figure"),
   "k": ("MTT dose-response, CDK7 L18F", "source_data/04_MTT_data/FACS_dose_response_2024_8_22.ipynb cell 130", "unported"),
   "m": ("Xenograft tumour volume fold change", "11_in_vivo_analysis/invivo_xenograft_CDK9_CDK7_combined.ipynb cell 5", "figure"),
   "n": ("Xenograft Ki67", "11_in_vivo_analysis/invivo_xenograft_CDK9_CDK7_combined.ipynb cell 7", "figure"),
 }),
 3: ("CDK12/13 and the PC3 transcriptional response", {
   "a": ("Lollipop, HQ461, subpool 2", "09_analysis_and_viz/lollipop_plots_canonical.ipynb cell 15", "figure"),
   "b": ("CDK12/13 top hits plus HQ461-specific", "09_analysis_and_viz/logo_plots_CDK12_13_canonical.ipynb cells 19-20", "figure"),
   "d": ("CDK12/13 validation MTT", "source_data/04_MTT_data/FACS_dose_response_2024_8_22.ipynb cells 136, 142-143", "figure"),
   "f": ("PC3 DEG counts", "reanalysis/scripts/03_figures.py", "data"),
   "g": ("PC3 drug response retained", "reanalysis/scripts/03_figures.py", "data"),
   "h": ("PC3 elongation deciles", "reanalysis/scripts/06_elongation_combined.py", "data"),
 }),
 4: ("CDK2/4/6 and clinical variants", {
   "a": ("Cross-drug resistance heatmap", "09_analysis_and_viz/logo_plots_CDK2_4_6_canonical.ipynb cell 29", "figure"),
   "b": ("Resistance hits per drug by gene", "09_analysis_and_viz/logo_plots_CDK2_4_6_canonical.ipynb cell 33", "figure"),
   "f": ("CDK4 validation allele frequency", "09_analysis_and_viz/analysis.ipynb cell 137", "unported"),
   # Included as the image that was published; this repository does not
   # regenerate these. The files are under data/raw/collaborator_panels/.
   "g": ("T47D pre- and post-selection plot",
         "data/raw/collaborator_panels/T47D_MCF7/cdk_variant_T47D_plot_t0_tf.svg", "artwork"),
   "i": ("T47D competition fold-change heatmap",
         "data/raw/collaborator_panels/T47D_MCF7/C010_T47D_comp_foldchange_heatmap_horizontal.svg", "artwork"),
   "j": ("C014",
         "data/raw/collaborator_panels/T47D_MCF7/20261011_C014.png", "artwork"),
   "k": ("T47D %S phase, ABE-derived lines",
         "data/raw/collaborator_panels/T47D_cell_cycle", "figure"),
   "l": ("T47D %S phase, CBE-derived lines",
         "data/raw/collaborator_panels/T47D_cell_cycle", "figure"),
   "m": ("RNA-seq marker-gene panel, DMSO and palbociclib",
         "data/raw/collaborator_panels/T47D_MCF7/marker_gene_panel_DMSO_palbo_stacked.svg", "artwork"),
   "n": ("Xenograft tumour volume, CDK6", "11_in_vivo_analysis/invivo_xenograft_CDK6_palbociclib.ipynb cell 20", "figure"),
   "o": ("Xenograft Ki67, CDK6", "11_in_vivo_analysis/invivo_xenograft_CDK6_palbociclib.ipynb cell 22", "figure"),
   "p": ("CDK6 clinical variants", "12_clinical_analysis/clinical.ipynb cells 19-23", "figure"),
 }),
 5: ("CROP-seq single-cell readout", {k: (v, "cdkbe.beatseq.fig_panels", "data") for k, v in {
   "a": "Library composition", "b": "UMAP, A549", "c": "UMAP, T47D",
   "d": "Editing by variant class", "e": "Editing vs cell cycle",
   "f": "CDK4 transcript level", "g": "Hallmark E2F score",
   "h": "Volcano, W198* vs NT", "i": "Volcano, W198* vs D196N",
   "j": "Severity by CDK", "k": "Severity vs drug LFC"}.items()}),
 6: ("KLIFS consensus resistance sites", {
   "(whole)": ("Consensus sites and homologous clinical resistance", "resistance-patterns pipeline (core7)", "unported"),
 }),
}

MARK = {"figure": "**figure**", "data": "source data", "unported": "not ported",
        "blocked": "blocked", "artwork": "published image"}


def reconcile(panels):
    """Replace the declared status with what the code actually provides.

    The tags below are easy to forget when a drawer lands, so they are only a
    floor: "blocked" stays whatever it was declared (the reason is in
    docs/DISCREPANCIES.md), and everything else is re-derived from the
    registries, so the notebook can never claim a panel it cannot draw --
    or omit one it can.
    """
    from cdkbe.source_data import DRAWN, PANELS

    # A panel can be registered and still be missing from FIGS above, which is
    # a hand-written list and falls behind. The notebook then silently omits a
    # panel it could draw, so anything the registries have for this figure and
    # FIGS does not mention is added here, described by its own drawer.
    n = str(_fignum)
    panels = dict(panels)
    for key in sorted(set(DRAWN) | set(PANELS)):
        letter = key[len(n):]
        if not key.startswith(n) or not (letter == "" or letter.isalpha()):
            continue
        if letter in panels or (letter == "" and any(
                k.startswith("(") for k in panels)):
            continue
        drawer = DRAWN.get(key)
        doc = (drawer.__doc__ or "").strip().splitlines()
        where = getattr(drawer, "__module__", "") if drawer else ""
        panels[letter or "(whole figure)"] = (
            doc[0].rstrip(".") if doc else f"panel {key}",
            where or "cdkbe.source_data", "figure" if key in DRAWN else "data")

    out = {}
    for letter, (desc, origin, status) in panels.items():
        if status == "blocked":
            out[letter] = (desc, origin, status)
            continue
        # figure 6 is a whole figure; the registries key it on the number
        key = panel_key(_fignum, letter)
        if key in DRAWN:
            status = "figure"
        elif key in PANELS:
            status = "data"
        out[letter] = (desc, origin, status)
    return out



def panel_key(n, letter: str) -> str:
    """The registry key for a panel.

    A whole-figure entry -- figure 6 -- is keyed on the number alone, and its
    letter is a placeholder like `(whole)`, so it must not be concatenated.
    """
    return f"{n}" + ("" if letter.startswith("(") else letter)


def build(n, title, panels):
    global _fignum
    _fignum = n
    panels = reconcile(panels)
    rows = "\n".join(
        f"| {p} | {d} | `{o}` | {MARK[s]} |" for p, (d, o, s) in sorted(panels.items()))
    drawn = [p for p, (_, _, s) in sorted(panels.items()) if s == "figure"]
    data = [p for p, (_, _, s) in sorted(panels.items()) if s in ("figure", "data")]

    cells = [md(f"""# Figure {n} — {title}

Reads only `data/processed/` and `data/reference/`. Nothing here recomputes a
screen statistic; if a number looks wrong the fix belongs upstream.

| panel | content | from | status |
|---|---|---|---|
{rows}

**figure** — this notebook regenerates both the source data and the panel.
**source data** — the table regenerates and is checked by `tests/`; the
plotting is still in the original notebook named above.
**not ported** — the plotting is not in this repository yet.
See `docs/COVERAGE.md` for the status of every panel."""),
             code(HEADER + """import matplotlib.pyplot as plt

from cdkbe import paths, io, editing, lfc, plotting
from cdkbe.source_data import PANELS

plotting.use_style()
paths.ensure(paths.FIG_MAIN)""")]

    if data:
        keys = ", ".join(f'"{panel_key(n, p)}"' for p in data)
        cells += [md("## Source data\n\nThe frames these panels plot."),
                  code(f"""for key in [{keys}]:
    result = PANELS[key]()
    frames = result if isinstance(result, dict) else {{"": result}}
    for name, df in frames.items():
        label = f"{{key}}/{{name}}" if name else key
        print(f"{{label:<28}} {{len(df):>8,}} rows x {{len(df.columns)}} cols")""")]

    if drawn:
        cells.append(md(f"""## Panels

Each panel below is drawn from `data/processed/` and rendered inline, and also
written to `figures/main/`. The numbers printed on the in vivo panels are
checked against the manuscript in `tests/test_invivo_figures.py`."""))
        cells.append(code("""from matplotlib.figure import Figure

from cdkbe.source_data import DRAWN
from cdkbe.source_data import invivo_figures as ivf

# Drawing functions predate this notebook and return four different shapes: a
# bare figure, (figure, source rows), a dict of either keyed by gene or by
# sub-panel, or -- for the in vivo cohorts -- a list of (figure, stats) pairs.
# `show` takes them all rather than each cell knowing which it is getting.
DRAW = {k: ("invivo" if k in ivf.FIGURES else "screen", v)
        for k, v in DRAWN.items()}


def show(fig, note=None):
    if note:
        print(note)
    display(fig)
    plt.close(fig)


def describe(obj):
    if isinstance(obj, dict):
        return ", ".join(f"{e}: {len(g):,} guides" for e, g in obj.items())
    return f"{len(obj):,} rows"


def draw(key):
    \"\"\"Draw one panel, show it, and report what it is made of.\"\"\"
    kind, fn = DRAW[key]
    result = fn(outdir=paths.FIG_MAIN)
    if kind == "invivo":
        for fig, stats in result:
            cols = [c for c in ["Cohort", "comparison", "group",
                                "n_WT", "n_mutant", "p_shown"] if c in stats.columns]
            display(stats[cols])
            show(fig)
        return
    if isinstance(result, Figure):
        show(result)
        return
    if isinstance(result, dict):              # one figure per gene or sub-panel
        for name, v in result.items():
            if isinstance(v, Figure):
                show(v, f"{key} {name}")
            else:
                fig, sub = v
                show(fig, f"{key} {name}: {describe(sub)}")
        return
    fig, obj = result
    if not isinstance(fig, Figure):           # a ChimeraX panel: rendered images
        print(describe(obj))
        display(fig)                          # tiled preview; see structures.Renders
        return
    show(fig, describe(obj))"""))
        for p in drawn:
            key = panel_key(n, p)
            desc = panels[p][0]
            cells.append(md(f"### {key} — {desc}"))
            cells.append(code(f'draw("{key}")'))

    artwork = [p for p, (_, _, s) in sorted(panels.items()) if s == "artwork"]
    if artwork:
        cells.append(md("## Included as published images\n\n" + "\n".join(
            f"- **{n}{p}** — {panels[p][0]}" for p in artwork) + """

These panels are included as the images that were published, under
`data/raw/collaborator_panels/`, and are shown below. This repository does not
regenerate them."""))
        cells.append(code("from IPython.display import Image, SVG"))
        for pl in artwork:
            rel = panels[pl][1]
            cells.append(code(
                f'print("{n}{pl}:", "{rel}")\n'
                f'_p = paths.ROOT / "{rel}"\n'
                f'display(SVG(filename=str(_p)) if _p.suffix == ".svg" '
                f'else Image(filename=str(_p)))'))

    remaining = [p for p, (_, _, s) in sorted(panels.items())
                 if s not in ("figure", "artwork")]
    if remaining:
        cells.append(md("## Still to port\n\n" + "\n".join(
            f"- **{p}** — {panels[p][0]} — `{panels[p][1]}`" for p in remaining)))
    return cells


for n, (title, panels) in FIGS.items():
    p = write(ROOT / "notebooks" / "08_figures" / f"figure_{n:02d}.ipynb", build(n, title, panels))
    panels = reconcile(panels)   # report what the notebook actually got
    drawn = sum(1 for _, _, s in panels.values() if s == "figure")
    dat = sum(1 for _, _, s in panels.values() if s in ("figure", "data"))
    print(f"  figure_{n:02d}.ipynb  {len(panels)} panels: {drawn} drawn, {dat} with source data")
