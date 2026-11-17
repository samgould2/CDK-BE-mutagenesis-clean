#!/usr/bin/env python3
"""Run the generated `.cxc` files through ChimeraX and save the images.

    python3 scripts/chimerax/render.py                 # every panel a figure needs
    python3 scripts/chimerax/render.py --all           # every .cxc under cxc/
    python3 scripts/chimerax/render.py CDK7_SY5609_labeled_zoom
    python3 scripts/chimerax/render.py --list          # what would run, and why

**ChimeraX cannot do this with `--nogui`.** It reports

    LimitationError: Unable to save images because OpenGL rendering is not available

because saving an image needs a GL context and `--nogui` has none. The app's own
binary, `ChimeraX.app/Contents/MacOS/ChimeraX`, does have one and renders fine
headless-ish: it opens a window, runs the commands, and exits. That is the
binary this script drives. `--nogui` is still right for anything that does not
save an image, such as `validate_structures.py`.

Each `.cxc` already carries an absolute `cd` into this directory and absolute
`save` paths -- ChimeraX silently aborts a render on a relative path -- so the
images land in `figures/chimerax/` wherever this is run from.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CXC = HERE / "cxc"

#: the app bundles that have a GL context, newest first
CANDIDATES = [
    "/Applications/ChimeraX-1.10.app/Contents/MacOS/ChimeraX",
    "/Applications/ChimeraX.app/Contents/MacOS/ChimeraX",
]

SAVE = re.compile(r"^save\s+(\S+)", re.M)

#: The leader-line panels are not command files. `label_consensus.py` projects
#: every label from the live ChimeraX camera, writes the 2dlabels commands, and
#: saves the panel and its SVG itself, so it is run with `runscript` and keyed
#: here by the image it produces. The arguments are the ones that reproduce the
#: published panels; `zoom` pulls the structure back to clear the gutters, and
#: without it the wide blocks of the peak panel overlap the spheres.
LEADER = {
    "consensus_leader_variants":
        "core7 variants consensus_strong recurrent,strong_singleton",
    "peaks_leader_variants":
        "core7 full peaks recurrent panel=peaks.cxc "
        "sites=MSA_peak_sites_labels.csv stem=peaks_leader zoom=0.66",
}


def chimerax_binary() -> str:
    """The GUI binary. `CHIMERAX` overrides it; `bin/ChimeraX` will not do."""
    env = os.environ.get("CHIMERAX")
    if env:
        return env
    for c in CANDIDATES:
        if Path(c).exists():
            return c
    for app in sorted(Path("/Applications").glob("ChimeraX*.app"), reverse=True):
        exe = app / "Contents" / "MacOS" / "ChimeraX"
        if exe.exists():
            return str(exe)
    raise SystemExit(
        "no ChimeraX found. Install it, or set CHIMERAX to the binary inside "
        "the app bundle -- ChimeraX.app/Contents/MacOS/ChimeraX, not "
        "Contents/bin/ChimeraX, which cannot render.")


def outputs_of(cxc: Path) -> list:
    """The image paths a command file writes."""
    return [Path(m) for m in SAVE.findall(cxc.read_text())]


def run(cxc: Path, binary: str, timeout: int = 900) -> list:
    """Render one command file. -> the images it produced."""
    want = outputs_of(cxc)
    for p in want:
        p.parent.mkdir(parents=True, exist_ok=True)
    before = {p: (p.stat().st_mtime if p.exists() else 0) for p in want}
    t0 = time.time()
    proc = subprocess.run(
        [binary, "--silent", "--exit", "--cmd", f"open {cxc}"],
        capture_output=True, text=True, timeout=timeout)
    made = [p for p in want if p.exists() and p.stat().st_mtime > before[p]]
    if not made:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-4:]
        raise RuntimeError(f"{cxc.name} wrote nothing\n    "
                           + "\n    ".join(tail))
    print(f"  {cxc.name:<46} {len(made)} image(s)  {time.time() - t0:.0f}s")
    return made


def run_leader(stem: str, binary: str, timeout: int = 900) -> list:
    """Render one leader-line panel. -> the files it produced."""
    args = LEADER[stem]
    t0 = time.time()
    proc = subprocess.run(
        [binary, "--silent", "--exit", "--cmd",
         f"runscript {HERE / 'label_consensus.py'} {args}"],
        capture_output=True, text=True, timeout=timeout, cwd=HERE)
    made = sorted((ROOT / "figures" / "chimerax").rglob(f"{stem}.*"))
    fresh = [p for p in made if time.time() - p.stat().st_mtime < time.time() - t0 + 60]
    if not fresh:
        log = HERE / "label_consensus.log"
        tail = (log.read_text().strip().splitlines()[-4:] if log.exists()
                else (proc.stderr or proc.stdout or "").strip().splitlines()[-4:])
        raise RuntimeError(f"{stem} wrote nothing\n    " + "\n    ".join(tail))
    print(f"  {stem:<46} {len(fresh)} file(s)  {time.time() - t0:.0f}s")
    return fresh


def needed() -> dict:
    """{cxc stem: [panel, ...]} for every image a published figure uses.

    Read from `docs/source_data_mapping_main.csv` rather than listed here, so a
    panel that starts or stops using a render changes this on its own.
    """
    import csv as _csv

    want = {}
    for which in ("main", "supp"):
        path = ROOT / "docs" / f"source_data_mapping_{which}.csv"
        if not path.exists():
            continue
        for row in _csv.DictReader(open(path, newline="")):
            fig, panel = str(row["figure"]), str(row["panel"])
            for token in str(row.get("source_file", "")).split():
                if "/chimerax/" not in token:
                    continue
                stem = Path(token).stem
                m = re.search(r"\d+", fig)
                label = ("S" if which == "supp" else "") + (
                    m.group() if m else fig) + (
                    "" if panel.startswith("(") else panel)
                want.setdefault(stem, []).append(label)
    return want


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("only", nargs="*", help="cxc stems to render")
    ap.add_argument("--all", action="store_true", help="every .cxc under cxc/")
    ap.add_argument("--list", action="store_true", help="say what would run")
    ap.add_argument("--timeout", type=int, default=900)
    a = ap.parse_args()

    # A command file often saves several images, so the index is by IMAGE, not
    # by file stem: `CDK7_SY5609_labeled_zoom.png` is written by
    # `CDK7_SY5609_labeled.cxc`, and looking for a cxc of that name finds
    # nothing. `render_all.cxc` only chains the others, so it is skipped.
    # A command file often saves several images, and the per-gene panels save
    # none at all -- `cxc/CDK2.cxc` only sets the view, and `render_all.cxc`
    # opens each of them in turn and does the saving. So the index is by IMAGE,
    # and the files to run are chosen to cover the wanted images in as few
    # ChimeraX launches as possible.
    files = sorted(CXC.glob("*.cxc"))
    for d in sorted(HERE.glob("cxc_*")):
        files += sorted(d.glob("*.cxc"))

    saves = {p: {img.stem for img in outputs_of(p)} for p in files}
    by_stem = {p.stem: p for p in files}

    want = needed()
    if a.only:
        images = set(a.only)
    elif a.all:
        images = set().union(*saves.values()) if saves else set()
    else:
        images = set(want)

    leaders = sorted(images & set(LEADER))
    images = images - set(LEADER)
    missing = sorted(images - set().union(*saves.values())) if saves else sorted(images)
    todo = images - set(missing)

    resolved = []
    while todo:
        best = max(saves, key=lambda p: (len(saves[p] & todo), -len(saves[p])))
        covers = saves[best] & todo
        if not covers:
            break
        resolved.append((sorted(covers)[0], best))
        todo -= covers
    # a stem asked for by name that saves nothing still gets run, so that
    # `render.py CDK2` does something visible rather than silently passing
    for s in (a.only or []):
        if s in by_stem and by_stem[s] not in {p for _, p in resolved} \
                and not saves[by_stem[s]]:
            resolved.append((s, by_stem[s]))

    if a.list:
        print(f"{len(resolved) + len(leaders)} job(s) would render:\n")
        for stem in leaders:
            covered = sorted(want.get(stem, []))
            print(f"  {'label_consensus.py ' + stem:<42} -> figure "
                  f"{', '.join(covered) or '-'}")
        for s, p in resolved:
            covered = sorted({w for img in outputs_of(p)
                              for w in want.get(img.stem, [])})
            print(f"  {p.name:<42} -> figure {', '.join(covered) or '-'}")
        if missing:
            print(f"\nno .cxc for: {', '.join(missing)}"
                  "\n  run the generate_*.py scripts first")
        return 0

    if missing:
        print(f"no .cxc for: {', '.join(missing)}\n"
              "  run the generate_*.py scripts in this directory first",
              file=sys.stderr)

    binary = chimerax_binary()
    print(f"rendering {len(resolved) + len(leaders)} job(s) with {binary}\n")
    made, failed = [], []
    for s, p in resolved:
        try:
            made += run(p, binary, a.timeout)
        except Exception as exc:
            failed.append((s, exc))
            print(f"  {s:<46} FAILED: {exc}", file=sys.stderr)
    for stem in leaders:
        try:
            made += run_leader(stem, binary, a.timeout)
        except Exception as exc:
            failed.append((stem, exc))
            print(f"  {stem:<46} FAILED: {exc}", file=sys.stderr)

    print(f"\n{len(made)} image(s) written under {ROOT / 'figures' / 'chimerax'}")
    if failed:
        print(f"{len(failed)} failed: {', '.join(s for s, _ in failed)}",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
