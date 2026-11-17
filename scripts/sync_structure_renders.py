#!/usr/bin/env python3
"""Copy the ChimeraX renders the figures need into the tier that ships.

    python3 scripts/sync_structure_renders.py           # copy and update the manifest
    python3 scripts/sync_structure_renders.py --check   # report drift, change nothing

`scripts/chimerax/render.py` writes into `figures/chimerax/`, which is build
output and is not committed: a full render is 119 MB of supersampled PNG, most
of it alternates and working views that no panel uses. But the structure
panels -- 1f, 2c, 2j, 3c, 4c, 4d, 4e, figure 6 and supplementary 15e -- cannot
be redrawn from a table, and ChimeraX is a GUI application that not everyone
reading this will have. So the renders those panels actually use are copied
into `data/reference/structure_renders/` and committed, and the modules look
there when `figures/chimerax/` has nothing.

What to copy is read from the registries rather than listed here, so a panel
that starts using a different render changes this on its own:

* `structures.MAPPING` names the per-panel stems,
* figure 6 and supplementary 15e name their own files, since they are whole
  figures rather than lettered panels and are not in that mapping,
* the essentiality montages under `structure_renders/essentiality/` are
  already shipped and are left alone.

The 2c zoom variants are copied too while which framing the published panel
used is still open (`docs/DISCREPANCIES.md` item 30).

ChimeraX is not deterministic to the byte, so a re-render differs from the
shipped copy even when nothing about the panel changed. `--check` therefore
reports which files are missing or differ, and it is for you to decide whether
a difference matters; it is not a test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cdkbe import paths                                        # noqa: E402
from cdkbe.source_data import structures                       # noqa: E402

RENDERS = paths.FIGURES / "chimerax"
SHIPPED = paths.REFERENCE / "structure_renders"
MANIFEST = paths.DATA / "MANIFEST.csv"
NOTE = "ChimeraX structure renders placed in the figures"

#: whole figures, which `structures.MAPPING` does not cover
EXTRA = [
    "cxc_core7_consensus_strong/consensus_leader_variants.svg",
    "cxc_core7_consensus_strong/consensus_leader_variants.png",
    "cxc_core7_consensus_strong/consensus_variants_plain.png",
    "cxc_core7_peaks/peaks_leader_variants.svg",
    "cxc_core7_peaks/peaks_leader_variants.png",
    # the alternative framing for figure 2c, while that is undecided
    "CDK9_fig2b_zoom.png",
    "CDK9_fig2b_labeled_zoom.png",
    # what each panel's structure, chains and ligand are
    "structure_table.csv",
    "structure_table.md",
]


def wanted() -> list:
    """Paths under `figures/chimerax/`, relative, that have to ship."""
    out = []
    for stems in structures.MAPPING.values():
        for stem in stems:
            hits = sorted(p for p in RENDERS.rglob(f"{stem}.*")
                          if p.suffix.lower() in {".png", ".svg", ".pdf"})
            out += [p.relative_to(RENDERS) for p in hits]
    out += [Path(name) for name in EXTRA]
    seen, unique = set(), []
    for rel in out:
        if str(rel) not in seen:
            seen.add(str(rel))
            unique.append(rel)
    return sorted(unique, key=str)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report what would change and write nothing")
    args = ap.parse_args()

    copied, same, differs, absent = [], [], [], []
    for rel in wanted():
        src, dst = RENDERS / rel, SHIPPED / rel
        if not src.exists():
            absent.append(rel)
            continue
        if dst.exists() and sha256(src) == sha256(dst):
            same.append(rel)
            continue
        (differs if dst.exists() else copied).append(rel)
        if not args.check:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

    verb = "would copy" if args.check else "copied"
    print(f"{len(same)} already identical, {len(differs)} {verb} (re-rendered), "
          f"{len(copied)} {verb} (new)")
    for rel in differs + copied:
        print(f"  {rel}")
    if absent:
        print(f"\n{len(absent)} not rendered yet, so not shipped:")
        for rel in absent:
            print(f"  {rel}")
        print("  run:  python3 scripts/chimerax/render.py")
    if args.check or not (differs or copied):
        return 1 if absent else 0

    # keep data/MANIFEST.csv true: one row per shipped file, with its source
    rows = list(csv.DictReader(open(MANIFEST, newline="")))
    fields = list(rows[0].keys())
    fresh = {}
    for rel in sorted(set(differs) | set(copied) | set(same), key=str):
        dst = SHIPPED / rel
        fresh[str(dst.relative_to(paths.ROOT))] = dict(
            source=str((RENDERS / rel).relative_to(paths.ROOT)),
            dest=str(dst.relative_to(paths.ROOT)),
            bytes_in=dst.stat().st_size, bytes_out=dst.stat().st_size,
            gzipped=False, note=NOTE, sha256=sha256(dst))
    kept = [r for r in rows if r["dest"] not in fresh]
    out = sorted(kept + list(fresh.values()), key=lambda r: r["dest"])
    with open(MANIFEST, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    print(f"\nmanifest: {len(rows)} -> {len(out)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
