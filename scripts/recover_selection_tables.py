#!/usr/bin/env python3
"""Recover the allele-frequency tables behind figure 4g and supplementary 13a.

    python3 scripts/recover_selection_tables.py          # write the tables
    python3 scripts/recover_selection_tables.py --check  # compare, write nothing

This repository has the published figure for both panels but not the table
behind it.

The numbers are still in the figure. Both are matplotlib line plots against a
labelled linear y axis, so `cdkbe.svg_values` inverts the axis and reads each
line's two vertices straight back out. That is an exact inversion rather than
an estimate, and this script will not write anything until it has demonstrated
as much: it first runs the same inversion over the two T47D cell-cycle panels,
whose published SVGs *and* underlying FlowJo export are both in this
repository, and stops unless the recovered bar heights equal the measured gate
frequencies exactly.

What comes back is 15 whole-number percentages out of 16 -- the same form the
hand-quantified figure 4f table takes -- plus one value, MCF7 CDK2 S261G after
selection, that is genuinely fractional at 4.264871%. Every value is within
3e-7 of what is written here, which is far inside the resolution of the source
coordinates, so the tables round to six decimals.

**These are recovered values, not the primary quantification.** If the Sanger
traces behind them ever turn up, that quantification replaces this table; what
this gives the repository meanwhile is a redrawable panel and a number a reader
can check, instead of a picture.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import pandas as pd                                           # noqa: E402

from cdkbe import io, paths, svg_values as sv                 # noqa: E402
from cdkbe.source_data import t47d_cell_cycle as cc           # noqa: E402

ARTWORK = paths.RAW / "collaborator_panels" / "T47D_MCF7"
OUT = paths.PROCESSED / "validation"

#: panel -> (the vendored artwork, the table written for it)
PANELS = {
    "4g": ("cdk_variant_T47D_plot_t0_tf.svg", "allele_frequency_4g_T47D.csv"),
    "S13a": ("cdk_variant_MCF7_plot_t0_tf.svg", "allele_frequency_S13a_MCF7.csv"),
}

#: how far a recovered number may sit from what is written, in percent
TOLERANCE = 1e-6


def self_check() -> None:
    """Invert two panels whose answer this repository already knows.

    The T47D cell-cycle panels came from the same author in the same style,
    and their FlowJo export is vendored, so the bar heights have a right
    answer. If the inversion cannot reproduce it there is no reason to trust
    it on a panel where nothing can be checked.
    """
    means = cc.means().set_index(["CellLine", "Treatment"])
    for stem, lines in (("T47D_cell_cycle_stacked_Sphase_ABE.svg", cc.ABE_LINES),
                        ("T47D_cell_cycle_stacked_Sphase_CBE.svg", cc.CBE_LINES)):
        got = [v for _, v in sv.bars(sv.load(ARTWORK / stem))]
        want = [float(means.loc[(line, treatment), "S"])
                for treatment in cc.TREATMENT_ORDER for line in lines]
        if len(got) != len(want):
            raise SystemExit(f"{stem}: read {len(got)} bars, expected {len(want)}")
        worst = max(abs(a - b) for a, b in zip(got, want))
        if worst > TOLERANCE:
            raise SystemExit(f"{stem}: recovered bars differ from the FlowJo "
                             f"export by up to {worst:.3g}%")
        print(f"  self-check {stem:<44} {len(got)} bars exact")


def recover(panel: str) -> pd.DataFrame:
    """-> one row per variant: its allele frequency before and after selection."""
    artwork, _ = PANELS[panel]
    rows = []
    for label, points in sv.labelled_series(ARTWORK / artwork):
        if len(points) != 2:
            raise SystemExit(f"{artwork}: {label} has {len(points)} points, expected 2")
        (_, pre), (_, post) = points
        rows.append({"Allele": label,
                     "Pre-Selection": round(pre, 6),
                     "Post-Selection": round(post, 6)})
    if not rows:
        raise SystemExit(f"{artwork}: no labelled series found")
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="compare with the tables on disk and write nothing")
    args = ap.parse_args()

    self_check()
    OUT.mkdir(parents=True, exist_ok=True)
    failed = False
    for panel, (artwork, name) in PANELS.items():
        table = recover(panel)
        path = OUT / name
        if args.check:
            if not path.exists():
                print(f"  {panel:<5} {name} is not there yet")
                failed = True
                continue
            have = io.read(path)
            same = (list(have["Allele"]) == list(table["Allele"])
                    and (have[["Pre-Selection", "Post-Selection"]].to_numpy()
                         - table[["Pre-Selection", "Post-Selection"]].to_numpy()
                         ).__abs__().max() <= TOLERANCE)
            print(f"  {panel:<5} {name:<36} {'matches' if same else 'DIFFERS'}")
            failed = failed or not same
        else:
            io.write(table, path, index=False)
            print(f"  {panel:<5} {len(table)} variants -> {path.relative_to(ROOT)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
