#!/usr/bin/env python3
"""Generate notebooks/09_supplementary/supp_figure_NN.ipynb from the mapping.

Reads `docs/source_data_mapping_supp.xlsx` — the panel-to-source-file list —
and writes one notebook per supplementary figure, each carrying its panel
table, which sources resolve, and what is still outstanding.

Rerun after updating the mapping:

    python3 scripts/make_supplementary_notebooks.py
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd                                        # noqa: E402
from _lib.nbwrite import md, code, write, HEADER           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MAPPING = ROOT / "docs" / "source_data_mapping_supp.csv"

sys.path.insert(0, str(ROOT / "src"))
from cdkbe.source_data import SUPP_DRAWN, SUPP_PANELS, SUPP_PANELS      # noqa: E402

#: values in the source column that are statements, not paths
NO_DATA = re.compile(r"^\s*none\b", re.I)
TO_COME = re.compile(r"^\s*missing\??\s*$", re.I)


def classify(cell) -> tuple[str, list]:
    """-> (status, [paths]). status is one of: sources, no_data, awaited, blank."""
    s = "" if cell is None else str(cell).strip()
    if not s or s.lower() == "nan":
        return "blank", []
    if TO_COME.match(s):
        return "awaited", []
    if NO_DATA.match(s):
        return "no_data", [s]
    return "sources", [p.strip().rstrip("'\"") for p in re.split(r"\s+(?=/)", s) if p.strip()]


#: a mapping note that points at a panel included as its published image
VENDORED = re.compile(r"vendored to (\S[^;]*?)\s*(?:;|$)")


def vendored_artwork(note: str):
    """-> the repo path of a panel included as its published image, or None.

    These are the panels this repository does not regenerate. The notebook
    shows the published image rather than redrawing it.
    """
    m = VENDORED.search(note or "")
    if not m:
        return None
    path = ROOT / m.group(1).strip()
    if path.suffix.lower() in {".svg", ".png"} and path.is_file():
        return path
    return None


def main():
    d = pd.read_csv(MAPPING).dropna(subset=["figure"])
    d["panel"] = d["panel"].astype(str).str.strip()
    out_dir = ROOT / "notebooks" / "09_supplementary"

    for fig, grp in d.groupby("figure", sort=True):
        num = int(re.search(r"(\d+)", fig).group(1))
        rows, artwork = [], []
        counts = dict(sources=0, no_data=0, awaited=0, blank=0)
        for _, r in grp.iterrows():
            status, paths = classify(r["source_file"])
            counts[status] += 1
            if status == "sources":
                shown = "<br>".join(f"`{Path(p).name}`" for p in paths)
                where = "<br>".join(
                    ("in repo" if Path(p).exists() else "**not found**") for p in paths)
            elif status == "no_data":
                shown, where = paths[0], "no underlying data"
            elif status == "awaited":
                shown, where = "—", "**to be supplied**"
            else:
                shown, where = "—", "—"
            note = "" if pd.isna(r.get("notes")) else str(r["notes"])
            key = f"S{num}{str(r['panel']).strip()}"
            vendored = vendored_artwork(note)
            if key in SUPP_DRAWN:
                where = "**drawn here**"
            elif key in SUPP_PANELS:
                where = "source data here"
            elif vendored is not None:
                where = "**published image** -- included, not regenerated"
                shown = f"`{vendored.name}`"
                artwork.append((str(r["panel"]), vendored))
            elif status == "sources" and paths and all(
                    p.endswith((".png", ".svg")) for p in paths):
                where = "render (ships in data/reference)"
            rows.append(f"| {r['panel']} | {shown} | {where} | {note} |")

        total = len(grp)
        drawn = [str(r["panel"]) for _, r in grp.iterrows()
                 if f"S{num}{str(r['panel']).strip()}" in SUPP_DRAWN]
        status = (f"{len(drawn)} of them drawn here" if drawn
                  else "None are ported yet")
        cells = [md(f"""# {fig}

{total} panel(s): {counts['sources']} with a source file,
{counts['awaited']} awaited, {counts['no_data']} with no underlying data
(blots, Sanger traces), {counts['blank']} unassigned. {status}.

| panel | source | status | notes |
|---|---|---|---|
{chr(10).join(rows)}

Regenerate with `python3 scripts/make_supplementary_notebooks.py` after editing
`docs/source_data_mapping_supp.xlsx`."""),
                 code(HEADER + """import matplotlib.pyplot as plt

from pathlib import Path

from IPython.display import Image, SVG

from cdkbe import paths, io, editing, lfc, plotting
from cdkbe.source_data import SUPP_DRAWN, SUPP_PANELS

plotting.use_style()
paths.ensure(paths.FIG_SUPP)


def files_of(result):
    \"\"\"Paths, if this panel is a ChimeraX render rather than a plot.\"\"\"
    first = result[0] if isinstance(result, tuple) and result else result
    if isinstance(first, (list, tuple)) and first and all(
            isinstance(x, Path) for x in first):
        return list(first)
    return None


def draw(key):
    \"\"\"Draw one supplementary panel, show it, and say what it is made of.\"\"\"
    result = SUPP_DRAWN[key](outdir=paths.FIG_SUPP)
    rendered = files_of(result)
    if rendered is not None:                  # a ChimeraX panel: show the image
        print(f"{key}: {len(result[1]):,} rows; "
              + ", ".join(f.name for f in rendered))
        display(result[0])                # tiled preview; see structures.Renders
        return
    if isinstance(result, dict):              # one figure per editor
        for name, (fig, obj) in result.items():
            print(f"{key} {name}: {len(obj):,} rows")
            display(fig)
            plt.close(fig)
        return
    fig, obj = result
    if isinstance(obj, dict):
        print(", ".join(f"{e}: {len(g):,} guides" for e, g in obj.items()))
    else:
        print(f"{len(obj):,} rows")
    display(fig)
    plt.close(fig)""")]
        if drawn:
            cells.append(md("## Panels"))
            for letter in drawn:
                cells.append(md(f"### {num}{letter}"))
                cells.append(code(f'draw("S{num}{letter}")'))
        sourced = [str(r["panel"]) for _, r in grp.iterrows()
                   if f"S{num}{str(r['panel']).strip()}" in SUPP_PANELS]
        if sourced and len(sourced) > len(drawn):
            cells.append(md(
                "## Source data without a drawing\n\n"
                + ", ".join(f"**{num}{p}**" for p in sourced if p not in drawn)
                + "\n\nThese reproduce their numbers but are not redrawn; "
                  "`docs/DISCREPANCIES.md` says why for each."))
            for letter in sourced:
                if letter in drawn:
                    continue
                cells.append(code(
                    f'for _sheet, _df in '
                    f'(lambda r: r if isinstance(r, dict) else {{"": r}})'
                    f'(SUPP_PANELS["S{num}{letter}"]()).items():\n'
                    f'    print(f"{num}{letter}" + (f" [{{_sheet}}]" if _sheet else ""),\n'
                    f'          f"-> {{len(_df):,}} rows x {{len(_df.columns)}} columns")\n'
                    f'    display(_df.head())'))
        if artwork:
            cells.append(md(
                "## Included as published images\n\n" + ", ".join(
                    f"**{num}{letter}**" for letter, _ in artwork)
                + """

These panels are included as the images that were published, under
`data/raw/collaborator_panels/`, and are shown below. This repository does not
regenerate them."""))
            for letter, path in artwork:
                rel = path.relative_to(ROOT)
                cells.append(code(
                    f'print("{num}{letter}:", "{rel}")\n'
                    f'_p = paths.ROOT / "{rel}"\n'
                    f'display(SVG(filename=str(_p)) if _p.suffix == ".svg" '
                    f'else Image(filename=str(_p)))'))

        remaining = [str(r["panel"]) for _, r in grp.iterrows()
                     if str(r["panel"]) not in drawn and str(r["panel"]) not in sourced
                     and str(r["panel"]) not in {a for a, _ in artwork}]
        if remaining:
            cells.append(md("## Still to port\n\n" + ", ".join(
                f"**{num}{p}**" for p in remaining) + """

To port one: find the notebook or script that wrote its source file, move the
reduction into `cdkbe` and the drawing into `cdkbe.source_data`, then add it to
that module's `FIGURES`. `cdkbe.lollipop`, `cdkbe.heatmap`, `cdkbe.logo_panel`
and `cdkbe.dose_response` already cover most of the layouts."""))
        p = write(out_dir / f"supp_figure_{num:02d}.ipynb", cells)
        print(f"  {p.name}  {total:>2} panels  "
              f"({len(drawn)} drawn, {counts['sources']} sourced, "
              f"{counts['awaited']} awaited)")


if __name__ == "__main__":
    main()
