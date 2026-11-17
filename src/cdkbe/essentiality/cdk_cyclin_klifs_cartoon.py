#!/usr/bin/env python3
"""Schematic CDK/cyclin cartoon annotated with the 85-residue KLIFS pocket regions.

Draws a generic bilobal CDK kinase domain with a cyclin docked against the aC
helix. Secondary-structure elements belonging to the KLIFS pocket definition are
coloured by functional family; everything else is neutral scaffold. Region names
and residue counts follow source_data/16_KLIFS/KLIFS_compiled.csv.

Output is a layered SVG with named groups (cyclin, kinase-lobes,
kinase-scaffold, klifs-subpockets, klifs-regions, ligand, annotations, legend,
variant-anchors, variant-markers) so the figure can be recoloured or stripped
down in Illustrator. Every KLIFS element is its own `region-<name>` group, and
`variant-anchors` holds an invisible point per region for snapping resistance
variants onto the cartoon later.

    python cdk_cyclin_klifs_cartoon.py [-o OUT.svg]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from xml.sax.saxutils import escape

from . import _paths
from .. import paths

# ---------------------------------------------------------------- palette ---
# Categorical slots validated all-pairs on a light surface with the dataviz
# validator: worst CVD dE 9.2, worst normal-vision dE 16.3.
FAMILY = {
    "roof":  ("#2a78d6", "Glycine-rich loop / pocket roof"),
    "back":  ("#1baf7a", "β3–αC & back pocket"),
    "front": ("#eb6834", "Gatekeeper, hinge & front pocket"),
    "cat":   ("#4a3aa7", "Catalytic loop & DFG / activation loop"),
}

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a887f"

SCAFFOLD_FILL = "#dedcd5"
SCAFFOLD_LINE = "#a9a69d"
LOBE_FILL = "#f3f2ee"
LOBE_LINE = "#d2d0c7"
CYCLIN_FILL = "#cfc7b6"
CYCLIN_LINE = "#8e8676"
CYCLIN_BODY = "#eae4d8"

_uid = [0]


def _nid(prefix: str) -> str:
    _uid[0] += 1
    return f"{prefix}{_uid[0]}"


def mix(hex_a: str, hex_b: str, t: float) -> str:
    a = [int(hex_a[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(hex_b[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b))


def tint(c: str, amount: float) -> str:
    """Positive -> toward white, negative -> toward black."""
    return mix(c, "#ffffff" if amount >= 0 else "#000000", abs(amount))


# ------------------------------------------------------------- primitives ---
def _unit(p0, p1):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    n = math.hypot(dx, dy) or 1.0
    return (dx / n, dy / n), n


def _catmull(points, closed=False):
    """Catmull-Rom through points -> cubic bezier path data."""
    p = list(points)
    pts = ([p[-1]] + p + [p[0], p[1]]) if closed else ([p[0]] + p + [p[-1]])
    d = f"M {pts[1][0]:.1f},{pts[1][1]:.1f}"
    for i in range(1, len(pts) - 2):
        p0, p1, p2, p3 = pts[i - 1], pts[i], pts[i + 1], pts[i + 2]
        c1 = (p1[0] + (p2[0] - p0[0]) / 6.0, p1[1] + (p2[1] - p0[1]) / 6.0)
        c2 = (p2[0] - (p3[0] - p1[0]) / 6.0, p2[1] - (p3[1] - p1[1]) / 6.0)
        d += (f" C {c1[0]:.1f},{c1[1]:.1f} {c2[0]:.1f},{c2[1]:.1f} "
              f"{p2[0]:.1f},{p2[1]:.1f}")
    return d + (" Z" if closed else "")


def helix(p0, p1, color, width=19.0, pitch=17.0, halo=False):
    """Schematic helix: a capsule with slanted coil ridges clipped inside it."""
    (ux, uy), length = _unit(p0, p1)
    vx, vy = -uy, ux
    r = width / 2.0
    axis = f"M {p0[0]:.1f},{p0[1]:.1f} L {p1[0]:.1f},{p1[1]:.1f}"
    cid = _nid("hx")
    out = ""
    if halo:
        out += (f'<path d="{axis}" fill="none" stroke="{SURFACE}" '
                f'stroke-width="{width + 10:.1f}" stroke-linecap="round"/>')
    # capsule body, dark rim drawn first
    out += (f'<path d="{axis}" fill="none" stroke="{tint(color, -0.38)}" '
            f'stroke-width="{width + 3.4:.1f}" stroke-linecap="round"/>'
            f'<path d="{axis}" fill="none" stroke="{color}" '
            f'stroke-width="{width:.1f}" stroke-linecap="round"/>')
    # coil ridges, clipped to the capsule
    out += (f'<clipPath id="{cid}"><path d="{axis}" fill="none" '
            f'stroke="#000" stroke-width="{width:.1f}" stroke-linecap="round"/>'
            f'</clipPath><g clip-path="url(#{cid})">')
    n = max(2, int(length / pitch))
    slant = width * 0.55
    for i in range(n + 1):
        t = (i + 0.5) / (n + 1)
        cx, cy = p0[0] + ux * length * t, p0[1] + uy * length * t
        ax, ay = cx + vx * r - ux * slant, cy + vy * r - uy * slant
        bx, by = cx - vx * r + ux * slant, cy - vy * r + uy * slant
        out += (f'<path d="M {ax:.1f},{ay:.1f} L {bx:.1f},{by:.1f}" '
                f'stroke="{tint(color, -0.34)}" stroke-width="2" '
                f'stroke-linecap="round" fill="none"/>')
    return out + "</g>"


def strand(p0, p1, color, width=12.0, head=18.0, halo=True):
    """Beta-strand drawn as a flat arrow pointing at p1."""
    (ux, uy), length = _unit(p0, p1)
    vx, vy = -uy, ux
    head = min(head, length * 0.45)
    bx, by = p1[0] - ux * head, p1[1] - uy * head
    hw, hh = width / 2.0, width * 1.05
    pts = [
        (p0[0] + vx * hw, p0[1] + vy * hw), (bx + vx * hw, by + vy * hw),
        (bx + vx * hh, by + vy * hh), p1,
        (bx - vx * hh, by - vy * hh), (bx - vx * hw, by - vy * hw),
        (p0[0] - vx * hw, p0[1] - vy * hw),
    ]
    d = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in pts) + " Z"
    out = ""
    if halo:
        out += (f'<path d="{d}" fill="{SURFACE}" stroke="{SURFACE}" '
                f'stroke-width="9" stroke-linejoin="round"/>')
    return out + (f'<path d="{d}" fill="{color}" stroke="{tint(color, -0.38)}" '
                  f'stroke-width="1.8" stroke-linejoin="round"/>')


def loop(points, color, width=11.0, halo=True):
    """Coil / loop region drawn as a smooth thick stroke."""
    d = _catmull(points)
    out = ""
    if halo:
        out += (f'<path d="{d}" fill="none" stroke="{SURFACE}" '
                f'stroke-width="{width + 10:.1f}" stroke-linecap="round" '
                f'stroke-linejoin="round"/>')
    return out + (
        f'<path d="{d}" fill="none" stroke="{tint(color, -0.38)}" '
        f'stroke-width="{width + 3.4:.1f}" stroke-linecap="round" '
        f'stroke-linejoin="round"/>'
        f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width:.1f}" '
        f'stroke-linecap="round" stroke-linejoin="round"/>')


def blob(points, fill, stroke, width=2.0, eid=None, opacity=1.0, dash=None):
    idattr = f' id="{eid}"' if eid else ""
    dashattr = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<path{idattr} d="{_catmull(points, closed=True)}" fill="{fill}" '
            f'fill-opacity="{opacity}" stroke="{stroke}" '
            f'stroke-width="{width}"{dashattr} stroke-linejoin="round"/>')


def connector(points, color=SCAFFOLD_LINE, width=5.0):
    return (f'<path d="{_catmull(points)}" fill="none" stroke="{color}" '
            f'stroke-width="{width}" stroke-linecap="round" '
            f'stroke-linejoin="round"/>')


def label(x, y, text, size=12.5, color=INK, anchor="middle", weight="500",
          style="normal"):
    return (f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" '
            f'text-anchor="{anchor}" font-weight="{weight}" '
            f'font-style="{style}">{escape(text)}</text>')


def leader(points, color=INK_3, width=1.1):
    return (f'<path d="{_catmull(points)}" fill="none" stroke="{color}" '
            f'stroke-width="{width}" stroke-linecap="round"/>')


# ------------------------------------------------------------ region table --
# name, family, KLIFS residue count, (label x, label y), optional leader points
REGIONS = [
    ("I",       "roof",  3, (680, 116),  None),
    ("II",      "roof",  4, (594, 118),  None),
    ("g.l",     "roof",  6, (754, 288),  [(740, 290), (726, 288)]),
    ("III",     "back",  6, (512, 116),  None),
    ("αC", "back", 11, (398, 180),  [(412, 184), (442, 186), (460, 188)]),
    ("b.l",     "back",  7, (646, 250),  [(646, 240), (644, 226)]),
    ("IV",      "back",  4, (782, 196),  None),
    ("V",       "front", 3, (850, 190),  None),
    ("GK",      "front", 1, (878, 258),  None),
    ("hinge",   "front", 3, (908, 316),  None),
    ("linker",  "front", 4, (920, 382),  None),
    ("αD", "front", 7, (906, 452),  [(894, 454), (880, 458)]),
    ("αE", "cat",   5, (890, 584),  [(878, 588), (850, 596)]),
    ("VI",      "cat",   3, (806, 672),  [(794, 670), (778, 666)]),
    ("c.l",     "cat",   8, (676, 506),  [(678, 516), (684, 526)]),
    ("VII",     "cat",   3, (650, 598),  None),
    ("VIII",    "cat",   1, (552, 684),  [(562, 676), (580, 666)]),
    ("xDFG",    "cat",   4, (524, 536),  [(540, 538), (562, 540)]),
    ("a.l",     "cat",   2, (638, 438),  [(634, 446), (628, 462)]),
]


def build() -> str:
    W, H = 1160, 800
    roof = FAMILY["roof"][0]
    back = FAMILY["back"][0]
    front = FAMILY["front"][0]
    cat = FAMILY["cat"][0]
    s: list[str] = []
    add = s.append

    add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="{W}" height="{H}" '
        f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">')
    add(f'<rect id="background" width="{W}" height="{H}" fill="{SURFACE}"/>')

    # ------------------------------------------------------------- cyclin --
    add('<g id="cyclin">')
    add(blob([(120, 330), (130, 232), (182, 168), (270, 142), (358, 166),
              (414, 240), (420, 344), (372, 436), (302, 484), (210, 478),
              (144, 420)], CYCLIN_BODY, CYCLIN_LINE, 2.2))
    for a, b in [((172, 236), (208, 434)), ((228, 196), (250, 420)),
                 ((286, 182), (298, 406)), ((342, 222), (346, 382)),
                 ((382, 282), (374, 364))]:
        add(helix(a, b, CYCLIN_FILL, width=20, pitch=14))
    add(label(262, 126, "Cyclin", size=19, color="#6f6758", weight="700"))
    add('</g>')

    # --------------------------------------------------------- lobe bodies --
    add('<g id="kinase-lobes">')
    add(blob([(424, 330), (420, 196), (466, 118), (558, 72), (706, 60),
              (830, 80), (896, 136), (906, 238), (884, 312), (800, 344),
              (672, 354), (530, 350), (452, 344)],
             LOBE_FILL, LOBE_LINE, 2.2, eid="n-lobe"))
    add(blob([(414, 408), (466, 372), (560, 392), (660, 420), (762, 424),
              (858, 392), (926, 426), (942, 546), (912, 660), (816, 732),
              (664, 754), (514, 728), (430, 640), (400, 506)],
             LOBE_FILL, LOBE_LINE, 2.2, eid="c-lobe"))
    add('</g>')

    # ----------------------------------------------------- non-KLIFS helices --
    add('<g id="kinase-scaffold">')
    for a, b in [((462, 520), (448, 608)), ((560, 700), (700, 716))]:
        add(helix(a, b, SCAFFOLD_FILL, width=19, pitch=16))
    add('</g>')

    # --------------------------------------------- ATP cleft / subpockets ---
    cleft = [(642, 348), (684, 322), (746, 310), (808, 306), (850, 328),
             (854, 372), (812, 404), (752, 418), (692, 418), (658, 390)]
    add('<g id="klifs-subpockets">')
    add(blob(cleft, "#ffffff", INK_3, 1.4, eid="atp-cleft", opacity=0.6,
             dash="5 4"))
    for x, y, t in [(676, 356, "BP"), (740, 352, "GA"), (800, 348, "AP"),
                    (806, 382, "FP")]:
        add(label(x, y, t, size=11, color=INK_2, weight="700", style="italic"))
    add('</g>')

    # ------------------------------------------------------- KLIFS regions --
    add('<g id="klifs-regions">')

    # b.l runs behind the N-lobe sheet -> drawn first, thin, no halo
    add('<g id="region-b.l">'
        + loop([(486, 134), (556, 190), (638, 222), (712, 234), (750, 240)],
               tint(back, 0.3), width=7, halo=False) + '</g>')

    # --- N-lobe -----------------------------------------------------------
    add(f'<g id="region-I">{strand((666, 136), (678, 254), roof)}</g>')
    add('<g id="region-g.l">'
        + loop([(678, 254), (720, 280), (718, 320), (674, 328), (630, 308),
                (610, 274)], tint(roof, 0.14)) + '</g>')
    add(f'<g id="region-II">{strand((610, 270), (598, 136), tint(roof, -0.14))}</g>')
    add(f'<g id="region-III">{strand((518, 134), (532, 274), back)}</g>')
    add('<g id="region-aC">'
        + helix((432, 296), (484, 132), tint(back, 0.1), width=21, pitch=15)
        + '</g>')
    add(f'<g id="region-IV">{strand((752, 240), (744, 148), tint(back, -0.24))}</g>')
    add(connector([(744, 148), (780, 134), (812, 152)]))
    add(f'<g id="region-V">{strand((814, 154), (824, 234), front)}</g>')

    # gatekeeper: a single residue, drawn as a marked side chain
    add('<g id="region-GK">'
        f'<path d="M 828,252 L 806,288" stroke="{tint(front, -0.38)}" '
        f'stroke-width="6.5" stroke-linecap="round"/>'
        f'<circle cx="828" cy="252" r="10.5" fill="{front}" '
        f'stroke="{tint(front, -0.42)}" stroke-width="2"/>'
        f'<circle cx="804" cy="291" r="7.5" fill="{tint(front, 0.28)}" '
        f'stroke="{tint(front, -0.42)}" stroke-width="2"/></g>')

    add('<g id="region-hinge">'
        + loop([(834, 268), (864, 298), (866, 336)], tint(front, -0.2), width=12)
        + '</g>')
    add('<g id="region-linker">'
        + loop([(866, 336), (876, 374), (860, 404)], tint(front, 0.16)) + '</g>')

    # --- C-lobe -----------------------------------------------------------
    add('<g id="region-aD">'
        + helix((856, 410), (870, 500), tint(front, -0.3), width=19, pitch=15)
        + '</g>')
    add(connector([(870, 500), (862, 528), (852, 542)]))
    add('<g id="region-aE">'
        + helix((848, 548), (790, 632), tint(cat, 0.34), width=19, pitch=15)
        + '</g>')
    add(connector([(790, 632), (780, 656), (768, 668)]))
    add(f'<g id="region-VI">{strand((766, 666), (746, 606), cat)}</g>')
    add('<g id="region-c.l">'
        + loop([(746, 604), (738, 556), (700, 528), (656, 528), (624, 556)],
               tint(cat, 0.18), width=12) + '</g>')
    add(f'<g id="region-VII">{strand((620, 560), (608, 620), tint(cat, -0.18))}</g>')
    add('<g id="region-VIII">'
        f'<circle cx="598" cy="656" r="9.5" fill="{tint(cat, 0.06)}" '
        f'stroke="{tint(cat, -0.42)}" stroke-width="2"/></g>')
    add('<g id="region-xDFG">'
        + loop([(592, 664), (568, 606), (566, 542), (588, 496), (624, 474)],
               tint(cat, -0.34), width=12) + '</g>')
    add('<g id="region-a.l">'
        + loop([(624, 474), (652, 462), (680, 458)], tint(cat, 0.5), width=12)
        + '</g>')
    add('</g>')

    # activation loop continuation, drawn in front of the C-lobe core
    cont = _catmull([(680, 458), (738, 506), (800, 566), (858, 624), (894, 668)])
    add('<g id="activation-loop-continuation">'
        f'<path d="{cont}" fill="none" stroke="{SURFACE}" stroke-width="17" '
        f'stroke-linecap="round"/>'
        f'<path d="{cont}" fill="none" stroke="{SCAFFOLD_LINE}" '
        f'stroke-width="8" stroke-linecap="round" stroke-dasharray="1 14"/>'
        + label(926, 700, "activation loop", size=11.5, color=INK_2,
                style="italic")
        + leader([(912, 692), (900, 678)]) + '</g>')

    # ---------------------------------------------------------- the ligand --
    add('<g id="ligand">')
    add(f'<path d="M 772,330 l 20,-11 l 20,11 l 0,23 l -20,11 l -20,-11 Z" '
        f'fill="#ffffff" stroke="{INK}" stroke-width="2.2" '
        f'stroke-linejoin="round"/>')
    add(f'<circle cx="758" cy="366" r="12" fill="#ffffff" stroke="{INK}" '
        f'stroke-width="2.2"/>')
    add(f'<path d="M 772,353 L 766,357" stroke="{INK}" stroke-width="2.2"/>')
    add(f'<path d="M 816,336 L 846,328 M 816,348 L 846,354" stroke="{INK_2}" '
        f'stroke-width="1.6" stroke-dasharray="3 3"/>')
    add(label(766, 452, "ATP-competitive inhibitor", size=11.5, color=INK_2,
              style="italic"))
    add(leader([(757, 442), (756, 402), (756, 382)]))
    add('</g>')

    # ------------------------------------------------------- K-E salt bridge --
    add('<g id="k-e-pair">'
        f'<path d="M 526,250 L 464,226" stroke="{INK_2}" stroke-width="1.8" '
        f'stroke-dasharray="4 4"/>'
        f'<circle cx="526" cy="250" r="4.5" fill="{INK}"/>'
        f'<circle cx="464" cy="226" r="4.5" fill="{INK}"/>'
        + label(498, 312, "K–E pair", size=11, color=INK_2, style="italic")
        + '</g>')

    # ------------------------------------------------------------- labels --
    add('<g id="annotations">')
    for name, _fam, n, (lx, ly), lead in REGIONS:
        if lead:
            add(leader(lead))
        add(label(lx, ly, name, size=13.5, color=INK, weight="700"))
        add(label(lx, ly + 13, str(n), size=9.5, color=INK_3, weight="400"))
    add(label(652, 112, "N", size=12, color=INK_3, weight="700"))
    add(label(858, 108, "N-lobe", size=15, color="#a7a49a", weight="700"))
    add(label(470, 452, "C-lobe", size=15, color="#a7a49a", weight="700"))
    add('</g>')

    # ------------------------------------------------------------- legend --
    add('<g id="legend">')
    lx, ly = 96, 520
    add(label(lx, ly, "KLIFS pocket regions", size=14.5, color=INK,
              anchor="start", weight="700"))
    add(label(lx, ly + 18, "coloured by functional family · small "
              "numeral = residues in region", size=10.5, color=INK_3,
              anchor="start", weight="400"))
    y = ly + 46
    rows = [(FAMILY[k][0], tint(FAMILY[k][0], -0.38), FAMILY[k][1])
            for k in ("roof", "back", "front", "cat")]
    rows += [(SCAFFOLD_FILL, SCAFFOLD_LINE, "Kinase domain outside KLIFS pocket"),
             (CYCLIN_FILL, CYCLIN_LINE, "Cyclin")]
    for fill, line, desc in rows:
        add(f'<rect x="{lx}" y="{y - 10}" width="17" height="13" rx="3" '
            f'fill="{fill}" stroke="{line}" stroke-width="1.5"/>')
        add(label(lx + 27, y, desc, size=11.5, color=INK_2, anchor="start",
                  weight="500"))
        y += 24
    add(label(lx, y + 16, "Subpockets  AP adenine · FP front · "
              "GA gate area · BP back", size=10.5, color=INK_3,
              anchor="start", weight="400"))
    add('</g>')

    # ------------------------------------------- placeholders for variants --
    add('<g id="variant-anchors" opacity="0">')
    for name, _fam, _n, (ax, ay), _lead in REGIONS:
        safe = name.replace(".", "-")
        add(f'<circle id="anchor-{safe}" cx="{ax}" cy="{ay}" r="3" fill="none"/>')
    add('</g>')
    add('<g id="variant-markers"><!-- resistance variants go here --></g>')

    add('</svg>')
    return "\n".join(s)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    default = Path(_paths.FIG_DIR) / "CDK_cyclin_KLIFS_cartoon.svg"
    ap.add_argument("-o", "--out", type=Path, default=default)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(build(), encoding="utf-8")
    print(f"wrote {paths.rel(args.out)}")


if __name__ == "__main__":
    main()
