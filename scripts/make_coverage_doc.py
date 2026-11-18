#!/usr/bin/env python3
"""Generate `docs/COVERAGE.md`: every published panel, and what it is waiting on.

    python scripts/make_coverage_doc.py

The two mapping files kept beside the manuscript list every panel of every
published figure:

    docs/source_data_mapping_main.csv    57 main-figure rows
    docs/source_data_mapping_supp.csv   103 supplementary rows

Both are exports of the workbooks kept beside the manuscript, with the
leading path stripped -- the home directory, and the name of the tree the file
sat in. The mapping is about which file made which panel, not where that file
sat; `data/MANIFEST.csv` carries the full original path of everything that was
vendored.

This walks both against the registries in `cdkbe.source_data` and classifies
each panel:

    drawn        source data reproduces AND the panel is redrawn here
    source only  the table reproduces; no plotting code (see DISCREPANCIES 13)
    no data      the mapping says the panel has no underlying data -- a western
                 blot, a Sanger trace, a schematic, a structure render
    published image
                 the panel is included as the image that was published; this
                 repository does not regenerate it
    blocked      the mapping's source_file column is empty or says "Missing",
                 so the file that drew it has not been located

The point of generating it is that "blocked" shrinks on its own as panels are
added, instead of a hand-written list going stale and overstating what is left.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd  # noqa: E402

from cdkbe import source_data as S  # noqa: E402

DOC = ROOT / "docs" / "COVERAGE.md"
MAIN = ROOT / "docs" / "source_data_mapping_main.csv"
SUPP = ROOT / "docs" / "source_data_mapping_supp.csv"

#: how the mapping's source_file column says "there is nothing to port"
NO_FILE = {"nan", ""}
#: the mapping's own word for a panel that has no underlying data file
EXPLICIT_NONE = {"none"}
#: the mapping's own word for a file it could not find
NOT_FOUND = {"missing", "missing?"}

#: phrases in the mapping that say what kind of panel it is
NO_DATA = ("western", "sanger")
STRUCTURE = ("chimerax/", "pymol", "_ribbon", "_cartoon")


def classify(key: str, source_file: str, drawn: set, panels: set,
             notes: str = "") -> tuple:
    """-> (status, note).

    The mapping distinguishes three kinds of empty cell and they mean different
    things: "None" is the author saying the panel has no data file (a schematic,
    a micrograph), "Missing" is the author saying they could not find the file,
    and a blank cell is a row nobody filled in. Collapsing them would put a
    hand-drawn schematic and a lost analysis script in the same bucket.
    """
    s = str(source_file).strip()
    low = s.lower()
    if key in drawn:
        return "drawn", ""
    if key in panels:
        return "source only", "no plotting code"
    # "None (Sanger trace)" is the author saying there is no data; a file whose
    # NAME happens to contain "sanger" is a figure drawn FROM Sanger reads, and
    # conflating the two lost supplementary 13b for a while.
    if low.startswith("none") and any(w in low for w in NO_DATA):
        return "no data", s
    if "data/raw/collaborator_panels" in str(notes) and key not in panels:
        return "published image", ("included as the published image, under "
                                   "data/raw/collaborator_panels")
    if any(w in low for w in STRUCTURE):
        return "structure render", "scripts/chimerax generates these"
    if low in EXPLICIT_NONE or low.startswith("none ") or low.startswith("none("):
        return "no data", "mapping says None -- schematic or illustration"
    if low in NOT_FOUND:
        return "file not found", "the mapping says Missing"
    if low in NO_FILE:
        return "unfilled", "the mapping row is blank"
    return "not ported", s


def _mapping(path) -> pd.DataFrame:
    """The mapping as it is written, with no value read as a missing value.

    `pd.read_csv` treats the four-letter string `None` as NA by default, and
    the mapping uses that word deliberately -- it is the author saying a panel
    has no underlying data file. Coerced to NaN it arrived at `classify` as the
    string `"nan"`, which is what `NO_FILE` means, so nine main-figure panels
    whose row says `None` were counted as rows nobody had filled in. Panel 2o
    was the only one spelled `none`, which pandas leaves alone, and it was the
    only one of the ten classified correctly.
    """
    x = pd.read_csv(path, keep_default_na=False)
    return x[x["figure"].astype(str).str.strip() != ""]


def rows_main() -> list:
    x = _mapping(MAIN)
    x["n"] = x["figure"].astype(str).str.extract(r"(\d+)")
    x = x.dropna(subset=["n"])
    out = []
    for _, r in x.iterrows():
        panel = str(r["panel"]).strip()
        # figure 6 is a whole figure, and its "panel" reads "(whole figure)";
        # the registries key it on the number alone
        key = f"{int(r['n'])}" + ("" if panel.startswith("(") else panel)
        out.append((int(r["n"]), panel, key, r.get("source_file", ""),
                    r.get("notes", "")))
    return out


def rows_supp() -> list:
    x = _mapping(SUPP)
    x["n"] = x["figure"].astype(str).str.extract(r"(\d+)")
    x = x.dropna(subset=["n"])
    out = []
    for _, r in x.iterrows():
        panel = str(r["panel"]).strip()
        key = f"S{int(r['n'])}{panel}"
        out.append((int(r["n"]), panel, key, r.get("source_file", ""),
                    r.get("notes", "")))
    return out


def resolve(key: str, panels: set, drawn: set) -> str:
    """The registry key for a published panel, allowing for letter case."""
    for cand in (key, key[:-1] + key[-1].upper(), key[:-1] + key[-1].lower()):
        if cand in panels or cand in drawn:
            return cand
    return key


def section(title, rows, drawn, panels, prefix) -> list:
    body = [f"## {title}", "", "| panel | status | notes |", "|---|---|---|"]
    counts = {}
    for n, panel, key, src, notes in rows:
        k = resolve(key, panels, drawn)
        status, note = classify(k, src, drawn, panels, notes)
        counts[status] = counts.get(status, 0) + 1
        label = f"{prefix}{n}{panel}"
        note = re.sub(r"\s+", " ", str(note))[:110].replace("|", r"\|")
        body.append(f"| {label} | {status} | {note} |")
    head = ", ".join(f"{v} {k}" for k, v in sorted(counts.items()))
    body.insert(1, f"\n{sum(counts.values())} published panels: {head}.\n")
    return body, counts


def main() -> int:
    drawn = set(S.DRAWN) | set(S.SUPP_DRAWN)
    panels = set(S.PANELS) | set(S.SUPP_PANELS)

    main_body, mc = section("Main figures", rows_main(), drawn, panels, "")
    supp_body, sc = section("Supplementary figures", rows_supp(), drawn, panels, "S")

    total = sum(mc.values()) + sum(sc.values())
    out = [
        "# Coverage",
        "",
        "Generated by `scripts/make_coverage_doc.py` -- do not edit by hand.",
        "",
        f"{total} published panels in all: {sum(mc.values())} main, "
        f"{sum(sc.values())} supplementary.",
        "",
        "| status | main | supplementary | meaning |",
        "|---|---:|---:|---|",
        f"| drawn | {mc.get('drawn', 0)} | {sc.get('drawn', 0)} | "
        "source data reproduces and the panel is redrawn here |",
        f"| source only | {mc.get('source only', 0)} | {sc.get('source only', 0)} | "
        "the table reproduces; no plotting code in this repository |",
        f"| published image | {mc.get('published image', 0)} | "
        f"{sc.get('published image', 0)} | included as the image that was "
        "published; not regenerated here |",
        f"| no data | {mc.get('no data', 0)} | {sc.get('no data', 0)} | "
        "a western blot, a Sanger trace, a schematic or a structure render |",
        f"| not ported | {mc.get('not ported', 0)} | {sc.get('not ported', 0)} | "
        "a source file is named but nothing here reads it yet |",
        f"| structure render | {mc.get('structure render', 0)} | "
        f"{sc.get('structure render', 0)} | a ChimeraX/PyMOL image; "
        "`scripts/chimerax` generates these |",
        f"| file not found | {mc.get('file not found', 0)} | "
        f"{sc.get('file not found', 0)} | the mapping itself says Missing |",
        f"| unfilled | {mc.get('unfilled', 0)} | {sc.get('unfilled', 0)} | "
        "the mapping row is blank |",
        "",
        "`docs/FIGURES.md` describes the panels that draw. "
        "`docs/DISCREPANCIES.md` has the open questions behind every "
        "**source only**, **published image**, **not ported** and "
        "**file not found** row.",
        "",
        "A `source_file` in the two mapping CSVs is the path the panel's "
        "source had inside the tree it came from, with that tree's own name "
        "removed; it is not a path in this repository. Where that file ships "
        "here, `data/MANIFEST.csv` gives both its original path and where it "
        "now sits.",
        "",
        *main_body,
        "",
        *supp_body,
        "",
    ]
    DOC.write_text("\n".join(out))
    both = {k: mc.get(k, 0) + sc.get(k, 0) for k in set(mc) | set(sc)}
    print(f"wrote {DOC.relative_to(ROOT)}: {total} panels -- "
          + ", ".join(f"{v} {k}" for k, v in sorted(both.items(), key=lambda kv: -kv[1])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
