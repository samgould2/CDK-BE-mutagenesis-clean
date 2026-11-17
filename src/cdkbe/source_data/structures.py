"""The ChimeraX structure panels: figures 1f, 2c, 2j, 3c, 4c, 4d, 4e.

These are renders, not plots, so they are the one group of panels this package
does not draw with matplotlib. `scripts/chimerax/render.py` runs the generated
`.cxc` command files through ChimeraX and writes the images; the functions here
name the result for its place in the manuscript and carry the residue table the
render marks up.

    python3 scripts/chimerax/render.py        # produces the images
    python scripts/render_figures.py 2j       # copies them to main_fig2j.png

**Rendering needs the ChimeraX application**, because saving an image needs an
OpenGL context. It is reproducible, but by running ChimeraX rather than by
running the test suite, so `FIGURES` here returns image paths rather than
figures and a panel whose image has not been rendered yet says so instead of
drawing something misleading.

The mapping from panel to image is `docs/source_data_mapping_main.csv`, the same
file `render.py` reads to decide what to render, so the two cannot disagree
about which image belongs to which panel.
"""
from __future__ import annotations

import csv
import re
from io import BytesIO
from pathlib import Path

import pandas as pd

from .. import io, paths

CHIMERAX = paths.ROOT / "scripts" / "chimerax"
#: where `scripts/chimerax/render.py` writes. Not committed: it is build
#: output, and a full render is 119 MB of supersampled PNG.
RENDERS = paths.FIGURES / "chimerax"
#: the copies that ship, so a clone can assemble the figures without having
#: ChimeraX installed. `scripts/sync_structure_renders.py` keeps them equal to
#: what `RENDERS` holds; a fresh render takes precedence over them.
SHIPPED = paths.REFERENCE / "structure_renders"
RENDER_ROOTS = (RENDERS, SHIPPED)

#: panel -> the residue table the render marks up. Figure 1f's per-gene panels
#: are the library's own coverage, which `library_figures.figure_1f` already
#: provides, so only the hit-marked panels carry one here.
HIT_TABLES = {
    "2c": CHIMERAX / "cdk9_fig2b_hit_residues.csv",
    "2j": CHIMERAX / "cdk7_sy5609_hit_residues.csv",
    "3c": CHIMERAX / "cdk12_13_ddb1_hit_residues.csv",
    "4c": CHIMERAX / "cdk2_4_6_hit_residues.csv",
    "4d": CHIMERAX / "cdk2_4_6_hit_residues.csv",
    "4e": CHIMERAX / "cdk2_4_6_hit_residues.csv",
}

#: the gene each active-site panel is about, for slicing the shared hit table
PANEL_GENE = {"4c": "CDK4", "4d": "CDK6", "4e": "CDK2"}


def _mapping() -> dict:
    """{panel: [image stem, ...]} for every panel that uses a render."""
    out = {}
    with open(paths.ROOT / "docs" / "source_data_mapping_main.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            fig, panel = str(row["figure"]), str(row["panel"]).strip()
            m = re.search(r"\d+", fig)
            if not m:
                continue
            key = m.group() + ("" if panel.startswith("(") else panel)
            stems = [Path(t).stem for t in str(row.get("source_file", "")).split()
                     if "/chimerax/" in t]
            if stems:
                out.setdefault(key, []).extend(stems)
    return out


MAPPING = _mapping()
#: the panels this module is responsible for. Figure 6 is the whole figure
#: rather than a panel, and figure 1f's render is one part of a panel whose
#: other part `library_figures` draws, so neither is registered here.
PANELS_WITH_RENDERS = [k for k in MAPPING if k in HIT_TABLES]


def find_render(stem: str) -> list:
    """Every image file for one render stem, from the first root that has it.

    A freshly rendered file wins over the shipped copy, so re-rendering is
    always what you get; with no ChimeraX on the machine, the shipped copy is
    what you get. The two are the same panel either way -- ChimeraX is not
    deterministic to the byte, so they are rarely the same bytes.
    """
    for root in RENDER_ROOTS:
        if not root.is_dir():
            continue
        hit = sorted((p for p in root.rglob(f"{stem}.*")
                      if p.suffix.lower() in {".png", ".svg", ".pdf"}),
                     key=lambda p: (len(p.relative_to(root).parts), str(p)))
        # `structure_renders/` holds a second set of per-gene renders under
        # `essentiality/`, drawn with a different colouring, so a bare stem
        # like `CDK2` matches twice. The shallowest match is the panel's.
        best, seen = [], set()
        for p in hit:
            if p.suffix.lower() in seen:
                continue
            seen.add(p.suffix.lower())
            best.append(p)
        if best:
            return best
    return []


class Renders(list):
    """The rendered images for a structure panel: a list of paths, with two
    additions.

    It prints repository-relative paths, because an absolute one would be baked
    into the output of any notebook that displayed it. And in a notebook it
    shows the images themselves, tiled and scaled down to a preview -- the
    files on disk stay the full-resolution originals, and embedding those in a
    notebook would cost tens of megabytes.
    """

    #: preview geometry: inches across, and the dpi it is saved at
    WIDTH_IN = 7.5
    DPI = 110

    @staticmethod
    def _rel(p: Path) -> str:
        try:
            return str(Path(p).relative_to(paths.ROOT))
        except ValueError:
            return str(p)

    def __repr__(self) -> str:
        return "[" + ", ".join(self._rel(p) for p in self) + "]"

    def _repr_png_(self):
        """A tiled preview of the rasters, for `display()` in a notebook."""
        raster = [p for p in self if p.suffix.lower() in {".png", ".jpg", ".jpeg"}]
        if not raster:
            return None
        import matplotlib.pyplot as plt
        from matplotlib.image import imread

        imgs = [imread(p) for p in raster]
        n = len(imgs)
        cols = 1 if n == 1 else 2 if n <= 4 else 3
        rows = -(-n // cols)
        aspect = max(i.shape[0] / i.shape[1] for i in imgs)
        fig, axes = plt.subplots(
            rows, cols, squeeze=False,
            figsize=(self.WIDTH_IN, self.WIDTH_IN / cols * aspect * rows))
        for ax in axes.ravel():
            ax.axis("off")
        for ax, img, path in zip(axes.ravel(), imgs, raster):
            ax.imshow(img)
            ax.set_title(path.name, fontsize=6)
        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=self.DPI, bbox_inches="tight")
        plt.close(fig)
        return buf.getvalue()


def images(panel: str) -> Renders:
    """The rendered files for a panel, in the order the mapping lists them."""
    out = Renders()
    for stem in MAPPING.get(panel, []):
        out += find_render(stem)
    return out


def source(panel: str) -> pd.DataFrame:
    """The residues the render marks, for the panel's source-data sheet."""
    table = HIT_TABLES[panel]
    d = io.read(table)
    gene = PANEL_GENE.get(panel)
    if gene is not None and "gene" in d.columns:
        d = d[d["gene"] == gene].reset_index(drop=True)
    return d


def render(panel: str, outdir=None, **kw):
    """-> ([image paths], the residue table).

    Does not run ChimeraX itself: a figure render takes minutes and needs the
    application, so it stays an explicit step. If the images are not there yet
    this says how to make them rather than returning an empty list.
    """
    found = images(panel)
    if not found:
        stems = ", ".join(MAPPING.get(panel, [])) or "?"
        raise FileNotFoundError(
            f"figure {panel} needs {stems} under {RENDERS} or {SHIPPED}, and "
            "neither has it. Run:  python3 scripts/chimerax/render.py")
    if outdir is not None:
        import shutil

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for p in found:
            shutil.copy2(p, outdir / p.name)
    return found, source(panel)


def _wrap(panel: str):
    def drawer(outdir=None, **kw):
        return render(panel, outdir=outdir, **kw)
    drawer.__name__ = drawer.__qualname__ = f"render_{panel}"
    n = len(MAPPING.get(panel, []))
    drawer.__doc__ = (f"ChimeraX render{'s' if n > 1 else ''} for figure {panel}"
                      f" ({n} image{'s' if n > 1 else ''}), with the residues marked.")
    return drawer


PANELS = {k: (lambda _p=k: source(_p)) for k in PANELS_WITH_RENDERS}
FIGURES = {k: _wrap(k) for k in PANELS_WITH_RENDERS}

#: `scripts/render_figures.py` reads this to name the files a panel produces
PARTS = {k: [s for s in MAPPING[k]] for k in PANELS_WITH_RENDERS}
