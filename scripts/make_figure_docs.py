#!/usr/bin/env python3
"""Generate `docs/FIGURES.md` -- what every panel is, and what draws it.

    python scripts/make_figure_docs.py

Nothing here is written by hand. Each row is assembled from the registries in
`cdkbe.source_data`: the description is the drawing function's docstring (or the
source-data function's, if the panel is not drawn here), and the module, function
and output filename are read off the objects themselves. A panel renamed, moved
or newly drawn therefore changes this file on the next run instead of leaving a
stale sentence behind, which is the failure mode a hand-written figure list has.

Panels with source data but no plotting code are listed in their own section
with the reason from `docs/COVERAGE.md` where there is one.
"""
from __future__ import annotations

import csv
import inspect
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cdkbe import source_data as S  # noqa: E402

sys.path.insert(0, str(ROOT / "scripts"))
from render_figures import stem  # noqa: E402

DOC = ROOT / "docs" / "FIGURES.md"
MODULES = ROOT / "docs" / "MODULES.md"


#: docstring first lines sometimes end with the function's return shape
#: ("-> (figure, stats)"), which describes the API, not the panel. Dropped here
#: rather than edited out of the docstrings, where it is useful.
_API_NOTE = re.compile(r"\s*(?:--?>|Returns)\s.*$")


def first_line(fn) -> str:
    d = (inspect.getdoc(fn) or "").strip()
    if not d:
        return ""
    line = _API_NOTE.sub("", d.splitlines()[0].strip()).strip()
    if line.endswith("."):
        line = line[:-1]
    return line.replace("|", "\\|")      # a bare pipe would split the table cell


def where(fn) -> str:
    mod = (fn.__module__ or "").replace("cdkbe.", "")
    name = getattr(fn, "__qualname__", getattr(fn, "__name__", "?"))
    return f"`{mod}.{name}`"


def describe(key, drawn, panels) -> str:
    for fn in (drawn.get(key), panels.get(key)):
        if fn is not None:
            d = first_line(fn)
            if d:
                return d
    return "--"


def written() -> dict:
    """{panel: [file, ...]} from the last render, if there was one.

    A panel drawn once per gene or cohort only reveals how many files it makes
    by being run, so the authoritative list is the index the renderer wrote.
    Without it this falls back to the single predictable name.
    """
    idx = ROOT / "figures" / "manuscript" / "INDEX.csv"
    if not idx.exists():
        return {}
    out = {}
    with open(idx, newline="") as fh:
        for row in csv.DictReader(fh):
            out.setdefault(row["panel"], []).append(row["file"])
    return out


def files_for(key, drawn, index) -> str:
    """The file names `render_figures.py` writes for this panel."""
    from cdkbe.source_data import figure5
    base = stem(key)
    if key in index:
        return ", ".join(f"`{f}`" for f in sorted(index[key]))
    if key in figure5.SPEC and len(figure5.SPEC[key]) > 1:
        return ", ".join(f"`{base}_{s}.pdf`" for s in figure5.parts(key))
    return f"`{base}.pdf`"


def table(keys, drawn, panels, index) -> list:
    rows = ["| panel | what it shows | drawn by | file |",
            "|---|---|---|---|"]
    for k in keys:
        fn = drawn.get(k)
        rows.append(f"| **{k}** | {describe(k, drawn, panels)} | "
                    f"{where(fn) if fn else '--'} | "
                    f"{files_for(k, drawn, index) if fn else '--'} |")
    return rows


def sort_key(k: str):
    body = k[1:] if k.startswith("S") else k
    num = body.rstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    return (k.startswith("S"), int(num or 0), body[len(num):].lower())


def summary(mod) -> str:
    """A module's first sentence, not its first physical line.

    Docstring opening sentences routinely wrap, so taking `splitlines()[0]`
    truncates them mid-clause. The opening paragraph is rejoined and cut at the
    first sentence end instead.
    """
    d = (inspect.getdoc(mod) or "").strip()
    if not d:
        return "--"
    para = []
    for line in d.splitlines():
        if not line.strip():
            break
        para.append(line.strip())
    text = " ".join(para)
    m = re.search(r"(?<=[a-z0-9\)])\.(?:\s|$)", text)
    if m:
        text = text[:m.start()]
    return text.rstrip(".")


def write_modules() -> int:
    """`docs/MODULES.md`: every module in `cdkbe`, with its own first line.

    An inventory drifts the moment a module is added, so it is read off the
    package rather than kept by hand. The prose about what the methods DO lives
    in `docs/METHODS.md`, which is written by hand and does not change per module.
    """
    import importlib
    import pkgutil

    import cdkbe

    rows = []
    for m in sorted(x.name for x in pkgutil.walk_packages(cdkbe.__path__, "cdkbe.")):
        short = m[len("cdkbe."):]
        if short.startswith("_"):
            continue
        try:
            mod = importlib.import_module(m)
        except Exception as exc:                      # a missing optional dep
            rows.append((short, f"*(not importable: {type(exc).__name__})*"))
            continue
        rows.append((short, summary(mod).replace("|", "\\|")))

    body = [
        "# Modules",
        "",
        "Generated by `scripts/make_figure_docs.py` -- do not edit by hand.",
        "The first line of each module's own docstring. `docs/METHODS.md` is the",
        "prose account of what these do and why; this is the index.",
        "",
        "| module | |",
        "|---|---|",
        *[f"| `cdkbe.{n}` | {d} |" for n, d in rows],
        "",
    ]
    MODULES.write_text("\n".join(body))
    return len(rows)


def main() -> int:
    drawn = {**S.DRAWN, **S.SUPP_DRAWN}
    panels = {**S.PANELS, **S.SUPP_PANELS}

    index = written()
    main_k = sorted((k for k in S.PANELS), key=sort_key)
    supp_k = sorted((k for k in S.SUPP_PANELS), key=sort_key)
    undrawn = sorted((k for k in panels if k not in drawn), key=sort_key)

    out = [
        "# Figures",
        "",
        "Generated by `scripts/make_figure_docs.py` -- do not edit by hand.",
        "",
        f"{len(S.PANELS)} main-figure panels and {len(S.SUPP_PANELS)} supplementary "
        f"panels have source data in this repository; {len(drawn)} of them are also "
        "redrawn here.",
        "",
        "Write every drawable panel, named for its place in the manuscript:",
        "",
        "```bash",
        "python scripts/render_figures.py",
        "```",
        "",
        "Each panel's source-data table is available without drawing anything:",
        "",
        "```python",
        "from cdkbe.source_data import PANELS, SUPP_PANELS",
        'df = PANELS["3a"]()',
        "```",
        "",
        "`figures/manuscript/INDEX.csv` lists every file written and the function",
        "that drew it. A panel assembled in Illustrator from more than one drawn",
        "piece gets one file per piece.",
        "",
        "## Main figures",
        "",
        *table(main_k, drawn, panels, index),
        "",
        "## Supplementary figures",
        "",
        "Keys are the mapping's own spelling, `S` + figure + panel.",
        "",
        *table(supp_k, drawn, panels, index),
        "",
    ]
    if undrawn:
        out += [
            "## Source data only",
            "",
            "These reproduce as tables but have no plotting code in this "
            "repository. `docs/DISCREPANCIES.md` item 13 says why for each.",
            "",
            "| panel | what it shows |",
            "|---|---|",
            *[f"| **{k}** | {describe(k, drawn, panels)} |" for k in undrawn],
            "",
        ]
    DOC.write_text("\n".join(out))
    print(f"wrote {DOC.relative_to(ROOT)} "
          f"({len(main_k)} main, {len(supp_k)} supplementary, {len(undrawn)} undrawn)")
    n = write_modules()
    print(f"wrote {MODULES.relative_to(ROOT)} ({n} modules)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
