"""
Residue table for figure 2c: the CDK9 variants figure 2b plots, and only those.

Figure 2c is the CDK9/cyclin T1 structure with 2b's variants marked. The two
panels therefore have to agree about which variants those are, so this script
does not keep its own list: it calls the same selection figure 2b draws from,
`cdkbe.source_data.screen_figures._cdk9_top_hits`, and splits each guide's
HGVSp call into one row per substituted residue. Change 2b's rule and 2c
follows on the next run.

That is a narrower set than `hit_residues.csv`, which is every CDK9 residue
passing FDR in any arm (16 positions) and is what the figure 1f panel marks.
2b keeps the ten strongest guides at the ranking dose, and those ten name five
residues.

Input : the subpool-1 barcode LFC tables and the canonical HGVSp tables,
        by way of figure 2b's own selection
Output: cdk9_fig2b_hit_residues.csv  (one row per guide / residue)

Usage:  python3 scripts/chimerax/build_cdk9_fig2b_hits.py
"""
from __future__ import annotations

import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

OUT = os.path.join(HERE, "cdk9_fig2b_hit_residues.csv")
SUBST = re.compile(r"^([A-Z])(\d+)([A-Z]|\*)$")
FIELDS = ["gene", "editor", "position", "wt", "mut", "variant", "max_LFC",
          "min_FDR", "gRNA_id", "HGVSp", "screen", "variants_at_position"]


def parse_hgvsp(hgvsp: str):
    """'D104N_C106Y' -> [('D', 104, 'N'), ('C', 106, 'Y')]"""
    out = []
    for part in str(hgvsp).split("_"):
        m = SUBST.match(part.strip())
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def build(canonical: bool = True):
    from cdkbe.source_data import screen_figures as sf

    hits = sf._cdk9_top_hits(canonical=canonical)
    #: the five arms figure 2b plots; the LFC and FDR recorded per guide are the
    #: strongest and the smallest across exactly those arms, so a guide is
    #: described by the panel it appears in rather than by the whole screen
    arms = list(sf.CDK9_CONDITIONS)
    lfc = [f"LFC_median_{a}" for a in arms]
    fdr = [f"FDR_{a}" for a in arms]

    rows = []
    for _, v in hits.iterrows():
        subs = parse_hgvsp(v["label"])
        if not subs:
            continue
        at_position = v["label"].replace("_", ",")
        for wt, pos, mut in subs:
            rows.append(dict(
                gene="CDK9", editor=v["Editor"], position=pos, wt=wt, mut=mut,
                variant=f"{wt}{pos}{mut}",
                max_LFC=float(v[lfc].max()), min_FDR=float(v[fdr].min()),
                gRNA_id=v["gRNA_id"], HGVSp=v["label"], screen="Subpool1",
                variants_at_position=at_position))

    rows.sort(key=lambda r: (-r["max_LFC"], r["position"]))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    return rows


if __name__ == "__main__":
    rows = build()
    positions = sorted({(r["position"], r["wt"]) for r in rows})
    print(f"wrote {OUT}")
    print(f"{len(rows)} rows over {len(positions)} positions: "
          + ", ".join(f"{wt}{p}" for p, wt in positions))
