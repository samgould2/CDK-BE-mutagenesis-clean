#!/usr/bin/env python3
"""Vector cartoon of a real CDK/cyclin structure coloured by KLIFS pocket region.

Builds a publication-style ribbon figure directly from backbone coordinates
rather than drawing a schematic by hand: ribbon frames come from the peptide
planes, secondary structure from the mmCIF `_struct_conf` / `_struct_sheet_range`
records, and the image is assembled as depth-sorted polygons (painter's
algorithm) with Lambert shading and depth fog. Output is layered SVG, so it stays
editable in Illustrator.

KLIFS pocket regions are taken from source_data/16_KLIFS/KLIFS_compiled.csv and
coloured by functional family; the rest of the kinase is neutral scaffold and the
cyclin is a muted second tone.

    python cdk_cyclin_klifs_ribbon.py                      # CDK9/cyclin T1, 8k5r
    python cdk_cyclin_klifs_ribbon.py --structure 7nxj.cif --cdk CDK13 \
        --kinase-chain A --cyclin-chain B
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from . import _paths
from .. import paths

# ---------------------------------------------------------------- palette ---
# Categorical slots validated all-pairs on a light surface (dataviz validator):
# worst CVD dE 9.2, worst normal-vision dE 16.3.
FAMILY = {
    "roof":  ("#2a78d6", "Glycine-rich loop / pocket roof"),
    "back":  ("#1baf7a", "β3–αC & back pocket"),
    "front": ("#eb6834", "Gatekeeper, hinge & front pocket"),
    "cat":   ("#4a3aa7", "Catalytic loop & DFG / activation loop"),
}

# KLIFS region -> (family, lightness step within the family)
REGION_STYLE = {
    "I":     ("roof", 0.00), "g.l": ("roof", 0.20), "II": ("roof", -0.16),
    "III":   ("back", 0.00), "αC": ("back", 0.20), "b.l": ("back", -0.10),
    "IV":    ("back", -0.26),
    "V":     ("front", 0.18), "GK": ("front", -0.30), "hinge": ("front", 0.00),
    "linker": ("front", 0.30), "αD": ("front", -0.16),
    "αE": ("cat", 0.38), "VI": ("cat", 0.00), "c.l": ("cat", 0.20),
    "VII":   ("cat", -0.16), "VIII": ("cat", 0.08), "xDFG": ("cat", -0.34),
    "a.l":   ("cat", 0.52),
}

# order the regions run along the chain, for the legend / label ring
REGION_ORDER = ["I", "g.l", "II", "III", "αC", "b.l", "IV", "V", "GK",
                "hinge", "linker", "αD", "αE", "VI", "c.l", "VII",
                "VIII", "xDFG", "a.l"]

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a887f"
SCAFFOLD = "#cdd1d2"
CYCLIN = "#c0b095"
LIGAND = "#2e2c28"


def _hex2rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], float)


def _rgb2hex(c):
    c = np.clip(c, 0, 255).round().astype(int)
    return "#%02x%02x%02x" % tuple(c)


def tint(h, amount):
    """Positive -> toward white, negative -> toward black."""
    base = _hex2rgb(h)
    target = np.array([255.0, 255, 255]) if amount >= 0 else np.zeros(3)
    return _rgb2hex(base + (target - base) * abs(amount))


# ------------------------------------------------------------- structure ---
def load_structure(path: Path, kinase_chain: str, cyclin_chain: str):
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("s", str(path))[0]
    doc = MMCIF2Dict(str(path))

    def ranges(prefix):
        asym = doc.get(f"{prefix}.beg_auth_asym_id")
        beg = doc.get(f"{prefix}.beg_auth_seq_id")
        end = doc.get(f"{prefix}.end_auth_seq_id")
        if not asym:
            return []
        if isinstance(asym, str):
            asym, beg, end = [asym], [beg], [end]
        return [(a, int(b), int(e)) for a, b, e in zip(asym, beg, end)]

    helices = ranges("_struct_conf")
    sheets = ranges("_struct_sheet_range")

    chains = {}
    for cid in (kinase_chain, cyclin_chain):
        res = [r for r in model[cid] if r.id[0] == " " and "CA" in r]
        entries = []
        for r in res:
            ca = r["CA"].get_coord().astype(float)
            o = r["O"].get_coord().astype(float) if "O" in r else None
            entries.append({"resi": r.id[1], "ca": ca, "o": o, "ss": "C"})
        for a, b, e in helices:
            if a != cid:
                continue
            for en in entries:
                if b <= en["resi"] <= e:
                    en["ss"] = "H"
        for a, b, e in sheets:
            if a != cid:
                continue
            for en in entries:
                if b <= en["resi"] <= e:
                    en["ss"] = "E"
        chains[cid] = entries

    hetatms = []
    for ch in model:
        for r in ch:
            if r.id[0].startswith("H_") and r.id[0] not in ("H_HOH",):
                if len(r) >= 12:  # skip ions, glycerol, buffer fragments
                    hetatms.append((r.get_resname(),
                                    np.array([a.get_coord() for a in r], float)))
    return chains, hetatms


def klifs_map(csv: Path, cdk: str):
    """resi -> KLIFS region name, for the given CDK column."""
    df = pd.read_csv(csv)
    df = df[["KLIFS_idx", cdk]].dropna()
    out = {}
    for idx, resi in zip(df.KLIFS_idx, df[cdk]):
        out[int(resi)] = idx.split(":")[0]
    return out


# ---------------------------------------------------------- ribbon frames ---
def ribbon_frames(entries):
    """Guide points + side vectors, one per residue, with twist removed."""
    ca = np.array([e["ca"] for e in entries])
    n = len(ca)
    side = np.zeros((n, 3))
    prev = None
    for i in range(n):
        a = ca[min(i + 1, n - 1)] - ca[max(i - 1, 0)]
        o = entries[i]["o"]
        b = (o - ca[i]) if o is not None else np.array([0.0, 0.0, 1.0])
        c = np.cross(a, b)
        d = np.cross(c, a)
        nrm = np.linalg.norm(d)
        d = d / nrm if nrm > 1e-6 else (prev if prev is not None
                                        else np.array([0.0, 0.0, 1.0]))
        if prev is not None and float(np.dot(d, prev)) < 0:
            d = -d          # remove the 180-degree flip between residues
        side[i] = d
        prev = d

    # smooth the beta pleat and gentle coil noise; keep the helical spiral
    guide = ca.copy()
    for i in range(1, n - 1):
        if entries[i]["ss"] == "E":
            guide[i] = 0.25 * ca[i - 1] + 0.5 * ca[i] + 0.25 * ca[i + 1]
        elif entries[i]["ss"] == "C":
            guide[i] = 0.12 * ca[i - 1] + 0.76 * ca[i] + 0.12 * ca[i + 1]
    return guide, side


def half_widths(entries):
    """Ribbon half-width per residue, with arrowheads on strand ends."""
    n = len(entries)
    w = np.zeros(n)
    base = {"H": 2.05, "E": 1.35, "C": 0.42}
    for i, e in enumerate(entries):
        w[i] = base[e["ss"]]
    for i, e in enumerate(entries):
        if e["ss"] != "E":
            continue
        last = (i == n - 1) or entries[i + 1]["ss"] != "E"
        if last:
            w[i] = 0.05            # tip
            if i - 1 >= 0 and entries[i - 1]["ss"] == "E":
                w[i - 1] = 2.45    # flare
    # ease the width across SS boundaries so edges do not jump
    sm = w.copy()
    for i in range(1, n - 1):
        if entries[i]["ss"] == "C" and entries[i - 1]["ss"] != "C":
            sm[i] = 0.55 * w[i] + 0.45 * w[i - 1]
        if entries[i]["ss"] == "C" and entries[i + 1]["ss"] != "C":
            sm[i] = 0.55 * w[i] + 0.45 * w[i + 1]
    return sm


def resample(guide, side, width, colors, ss, per_res=9):
    """Catmull-Rom subdivision of the whole strip, splitting at chain breaks."""
    segs, cur = [], []
    for i in range(len(guide) - 1):
        cur.append(i)
        if np.linalg.norm(guide[i + 1] - guide[i]) > 4.6:
            segs.append(cur + [i])
            cur = []
    cur.append(len(guide) - 1)
    segs.append(cur)

    out = []
    for seg in segs:
        idx = sorted(set(seg))
        if len(idx) < 2:
            continue
        g, s, w = guide[idx], side[idx], width[idx]
        cols = [colors[i] for i in idx]
        sss = [ss[i] for i in idx]
        pts, sds, wds, cls, sst = [], [], [], [], []
        m = len(idx)
        for j in range(m - 1):
            p0 = g[max(j - 1, 0)]
            p1, p2 = g[j], g[j + 1]
            p3 = g[min(j + 2, m - 1)]
            for t in np.linspace(0, 1, per_res, endpoint=False):
                t2, t3 = t * t, t * t * t
                pts.append(0.5 * ((2 * p1) + (-p0 + p2) * t +
                                  (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 +
                                  (-p0 + 3 * p1 - 3 * p2 + p3) * t3))
                v = s[j] * (1 - t) + s[j + 1] * t
                nv = np.linalg.norm(v)
                sds.append(v / nv if nv > 1e-6 else s[j])
                wds.append(w[j] * (1 - t) + w[j + 1] * t)
                cls.append(cols[j])
                sst.append(sss[j])
        pts.append(g[-1]); sds.append(s[-1]); wds.append(w[-1])
        cls.append(cols[-1]); sst.append(sss[-1])
        out.append((np.array(pts), np.array(sds), np.array(wds), cls, sst))
    return out


# ------------------------------------------------------------- rendering ---
def shade(hex_color, lum, fog):
    c = _hex2rgb(hex_color)
    if lum >= 0.74:
        c = c + (255 - c) * min(1.0, (lum - 0.74) * 1.5)
    else:
        c = c * (1.0 - min(0.55, (0.74 - lum) * 0.95))
    bg = _hex2rgb(SURFACE)
    c = c + (bg - c) * fog
    return _rgb2hex(c)


def build(args) -> str:
    cif = Path(args.structure)
    if not cif.is_absolute():
        cif = Path(_paths.STRUCTURES) / cif
    chains, hets = load_structure(cif, args.kinase_chain, args.cyclin_chain)
    kl = klifs_map(Path(_paths.KLIFS_DIR) / "KLIFS_compiled.csv", args.cdk)

    kin = chains[args.kinase_chain]
    cyc = chains[args.cyclin_chain]
    # Clip the construct back to the kinase domain: the disordered N-terminal
    # extension loops across the N-lobe and reads as a hole in the fold.
    _lo, _hi = min(kl), max(kl)
    kin = [e for e in kin if _lo - 4 <= e["resi"] <= _hi + 170]

    # ---- view frame: N-lobe up, ligand (hence the cleft) toward the viewer --
    def centroid(entries, names):
        pts = [e["ca"] for e in entries if kl.get(e["resi"]) in names]
        return np.mean(pts, axis=0)

    kin_ca = np.array([e["ca"] for e in kin])

    # True lobe axis: the hinge splits the kinase domain into N- and C-lobe,
    # so take the bulk centroid of each side rather than pocket-lining regions
    # (which only span the cleft and give a badly tilted "up").
    hinge_resi = [r for r, g in kl.items() if g == "hinge"]
    h_lo, h_hi = min(hinge_resi), max(hinge_resi)
    kl_lo, kl_hi = min(kl), max(kl)
    n_lobe = np.array([e["ca"] for e in kin
                       if kl_lo - 12 <= e["resi"] < h_lo])
    c_lobe = np.array([e["ca"] for e in kin
                       if h_hi < e["resi"] <= kl_hi + 150])
    up = n_lobe.mean(axis=0) - c_lobe.mean(axis=0)
    up /= np.linalg.norm(up)

    # the inhibitor is whichever het sits in the ATP cleft (nearest the hinge)
    hinge_mid = centroid(kin, {"hinge", "GK"})
    ligs = [c for _n, c in hets
            if np.linalg.norm(c.mean(axis=0) - hinge_mid) < 16.0]
    lig_mid = np.vstack(ligs).mean(axis=0) if ligs else hinge_mid

    # Canonical kinase view: lay the long axis of the ATP cleft (back pocket ->
    # hinge) along the screen horizontal. With the lobe axis vertical that
    # fully determines the frame, and looking along the remaining axis looks
    # into the cleft rather than edge-on to it.
    back_a = centroid(kin, {"αC", "III"})
    front_a = centroid(kin, {"hinge", "linker"})
    right = front_a - back_a
    right -= up * float(np.dot(right, up))
    right /= np.linalg.norm(right)
    z = np.cross(right, up)                  # (right, up, z) right-handed
    # negating `right` negates z too, so this stays a rotation (no mirroring)
    if float(np.dot(lig_mid - kin_ca.mean(axis=0), z)) < 0:
        right, z = -right, -z
    R = np.vstack([right, up, z])            # world -> view

    # ---- per-residue colours ---------------------------------------------
    def colour_of(resi):
        reg = kl.get(resi)
        if reg is None:
            return SCAFFOLD
        fam, step = REGION_STYLE[reg]
        return tint(FAMILY[fam][0], step)

    strips = []
    for entries, base, back in ((kin, None, False), (cyc, CYCLIN, True)):
        guide, side = ribbon_frames(entries)
        width = half_widths(entries)
        cols = [base or colour_of(e["resi"]) for e in entries]
        ss = [e["ss"] for e in entries]
        for st in resample(guide, side, width, cols, ss):
            strips.append(st + (back,))

    # ---- project ----------------------------------------------------------
    allv = np.vstack([R @ p for s in strips for p in s[0]])
    # Frame on the kinase, not the complex: the cyclin is nearly as large and
    # sits off-axis, so including it whole shrinks the pocket to nothing. It
    # bleeds off the corner instead, which still reads as "bound to cyclin".
    kinv = np.vstack([R @ p for s in strips if not s[5] for p in s[0]])
    klo, khi = kinv[:, :2].min(axis=0), kinv[:, :2].max(axis=0)
    marg = np.array([0.30, 0.09]) * (khi - klo)
    lo, hi = klo - marg * np.array([0.25, 1.0]), khi + marg
    W, GUT, RPAD, VPAD = 1740, 250, 250, 76    # GUT/RPAD = label gutters
    scale = (W - GUT - RPAD) / (hi[0] - lo[0])
    Hc = int((hi[1] - lo[1]) * scale + 2 * VPAD)
    zmin, zmax = allv[:, 2].min(), allv[:, 2].max()

    def proj(p):
        v = R @ p
        return (GUT + (v[0] - lo[0]) * scale,
                Hc - VPAD - (v[1] - lo[1]) * scale, v[2])

    light = np.array([-0.35, 0.45, 0.82])
    light /= np.linalg.norm(light)

    prims = []     # (depth, svg)
    for pts, sds, wds, cls, sst, back in strips:
        P = np.array([proj(p) for p in pts])
        for i in range(len(pts) - 1):
            col = cls[i]
            zmid = 0.5 * (P[i, 2] + P[i + 1, 2])
            fog = 0.30 * (1.0 - (zmid - zmin) / max(1e-6, zmax - zmin))
            if back:
                fog = min(0.62, fog + 0.30)
            if sst[i] == "C" and wds[i] < 0.8 and wds[i + 1] < 0.8:
                lum = 0.80
                c = shade(col, lum, fog)
                prims.append((zmid,
                    f'<path d="M {P[i,0]:.1f},{P[i,1]:.1f} L '
                    f'{P[i+1,0]:.1f},{P[i+1,1]:.1f}" stroke="{c}" '
                    f'stroke-width="{0.86 * scale:.2f}" stroke-linecap="round" '
                    f'fill="none"/>'))
                continue
            a0 = pts[i] + sds[i] * wds[i]
            b0 = pts[i] - sds[i] * wds[i]
            a1 = pts[i + 1] + sds[i + 1] * wds[i + 1]
            b1 = pts[i + 1] - sds[i + 1] * wds[i + 1]
            q = [proj(a0), proj(a1), proj(b1), proj(b0)]
            nrm = np.cross(pts[i + 1] - pts[i], sds[i])
            nn = np.linalg.norm(nrm)
            lum = 0.74
            if nn > 1e-6:
                nv = R @ (nrm / nn)
                lum = 0.40 + 0.60 * abs(float(np.dot(nv, light)))
            c = shade(col, lum, fog)
            edge = shade(col, lum * (0.72 if back else 0.42), fog)
            d = "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y, _ in q) + " Z"
            prims.append((zmid,
                f'<path d="{d}" fill="{c}" stroke="{c}" stroke-width="0.7"/>'
                f'<path d="M {q[0][0]:.1f},{q[0][1]:.1f} L {q[1][0]:.1f},{q[1][1]:.1f}" '
                f'stroke="{edge}" stroke-width="1.0" fill="none" stroke-linecap="round"/>'
                f'<path d="M {q[3][0]:.1f},{q[3][1]:.1f} L {q[2][0]:.1f},{q[2][1]:.1f}" '
                f'stroke="{edge}" stroke-width="1.0" fill="none" stroke-linecap="round"/>'))

    # ligand sticks
    for het in ligs:
        pj = [proj(a) for a in het]
        for i in range(len(het)):
            for j in range(i + 1, len(het)):
                if np.linalg.norm(het[i] - het[j]) < 1.95:
                    zmid = 0.5 * (pj[i][2] + pj[j][2])
                    prims.append((zmid + 0.01,
                        f'<path d="M {pj[i][0]:.1f},{pj[i][1]:.1f} L '
                        f'{pj[j][0]:.1f},{pj[j][1]:.1f}" stroke="{LIGAND}" '
                        f'stroke-width="{0.55 * scale:.2f}" '
                        f'stroke-linecap="round" fill="none"/>'))

    prims.sort(key=lambda t: t[0])

    # ---- region label anchors (projected centroid of each region) ---------
    anchors = {}
    for reg in REGION_ORDER:
        pts = [proj(e["ca"])[:2] for e in kin if kl.get(e["resi"]) == reg]
        if pts:
            anchors[reg] = np.array(pts, float).mean(axis=0)
    # Labels live in a left-hand gutter, ordered by the vertical position of
    # their anchor and spread to a minimum spacing - a radial ring reads badly
    # here because every KLIFS region clusters around one cleft.
    keys = [r for r in REGION_ORDER if r in anchors]
    split_x = float(np.median([a[0] for a in anchors.values()]))
    side = {r: ("L" if anchors[r][0] <= split_x else "R") for r in keys}
    GAPY, TOP = 34.0, 46.0
    BOT = {"L": Hc - 200.0, "R": Hc - 70.0}

    def spread(group, bot):
        """Order by anchor height, then separate to a minimum spacing."""
        group = sorted(group, key=lambda r: anchors[r][1])
        ys = [anchors[r][1] for r in group]
        for i in range(1, len(ys)):
            ys[i] = max(ys[i], ys[i - 1] + GAPY)
        over = ys[-1] - bot
        if over > 0:
            ys = [y - over for y in ys]
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - GAPY)
        return group, [max(TOP, y) for y in ys]

    labels = {}
    for sd in ("L", "R"):
        grp = [r for r in keys if side[r] == sd]
        if not grp:
            continue
        grp, ys = spread(grp, BOT[sd])
        x = GUT - 108.0 if sd == "L" else W - RPAD + 118.0
        for r, y in zip(grp, ys):
            labels[r] = np.array([x, y])

    # ---------------------------------------------------------------- SVG --
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {Hc}" '
         f'width="{W}" height="{Hc}" '
         f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">',
         f'<rect id="background" width="{W}" height="{Hc}" fill="{SURFACE}"/>',
         '<g id="structure">']
    s += [p for _, p in prims]
    s.append('</g>')

    s.append('<g id="annotations">')
    for reg in keys:
        lp, ap = labels[reg], anchors[reg]
        fam = FAMILY[REGION_STYLE[reg][0]][0]
        n_res = sum(1 for v in kl.values() if v == reg)
        txt = f"{reg}"
        pw = 34 + 8.0 * len(txt)
        px, py = lp[0] - pw / 2, lp[1] - 12
        if side[reg] == "L":
            stub, elbow = px + pw, GUT - 52   # leaders turn here, then run in
        else:
            stub, elbow = px, W - RPAD + 56
        # pill -> horizontal stub -> straight run to the anchor dot
        s.append(f'<path d="M {stub:.1f},{lp[1]:.1f} L {elbow:.1f},{lp[1]:.1f} '
                 f'L {ap[0]:.1f},{ap[1]:.1f}" stroke="{INK_3}" '
                 f'stroke-width="1.1" fill="none" stroke-linejoin="round"/>')
        s.append(f'<circle cx="{ap[0]:.1f}" cy="{ap[1]:.1f}" r="3.4" '
                 f'fill="{INK}" stroke="{SURFACE}" stroke-width="1.5"/>')
        s.append(f'<rect x="{px:.1f}" y="{py:.1f}" width="{pw:.1f}" height="24" '
                 f'rx="12" fill="{SURFACE}" stroke="{fam}" stroke-width="2.2"/>')
        s.append(f'<text x="{lp[0] - 5:.1f}" y="{lp[1] + 4.5:.1f}" '
                 f'font-size="12.5" fill="{INK}" text-anchor="middle" '
                 f'font-weight="700">{escape(txt)}</text>')
        s.append(f'<text x="{px + pw - 8:.1f}" y="{lp[1] + 4:.1f}" '
                 f'font-size="9.5" fill="{INK_3}" text-anchor="end">'
                 f'{n_res}</text>')
    s.append('</g>')

    # legend
    s.append('<g id="legend">')
    lx, ly = 46, Hc - 160
    s.append(f'<text x="{lx}" y="{ly}" font-size="15" fill="{INK}" '
             f'font-weight="700">KLIFS pocket regions</text>')
    y = ly + 26
    rows = [(FAMILY[k][0], FAMILY[k][1]) for k in ("roof", "back", "front", "cat")]
    rows += [(SCAFFOLD, f"{args.cdk} outside KLIFS pocket"),
             (CYCLIN, args.cyclin_name)]
    for col, desc in rows:
        s.append(f'<rect x="{lx}" y="{y - 10}" width="17" height="13" rx="3" '
                 f'fill="{col}" stroke="{tint(col, -0.38)}" stroke-width="1.4"/>')
        s.append(f'<text x="{lx + 27}" y="{y}" font-size="12" fill="{INK_2}">'
                 f'{escape(desc)}</text>')
        y += 22
    s.append('</g>')

    s.append(f'<g id="caption"><text x="{W - 46}" y="{Hc - 34}" font-size="12" '
             f'fill="{INK_3}" text-anchor="end">{escape(args.cdk)}/'
             f'{escape(args.cyclin_name)} · PDB {escape(cif.stem.upper())}'
             f'</text></g>')
    # Invisible per-region anchors, plus a per-residue screen-coordinate table
    # written alongside the SVG, so variants can be placed on this exact view
    # without re-deriving the projection.
    s.append('<g id="variant-anchors" opacity="0">')
    for reg, ap in anchors.items():
        s.append(f'<circle id="anchor-{escape(reg).replace(".", "-")}" '
                 f'cx="{ap[0]:.1f}" cy="{ap[1]:.1f}" r="3" fill="none"/>')
    s.append('</g>')
    s.append('<g id="variant-markers"><!-- resistance variants go here --></g>')
    s.append('</svg>')

    rows = []
    for e in kin:
        x, y, zz = proj(e["ca"])
        rows.append({"resi": e["resi"], "klifs_region": kl.get(e["resi"], ""),
                     "x": round(x, 1), "y": round(y, 1),
                     "depth": round(float(zz), 2), "ss": e["ss"]})
    build.resi_xy = pd.DataFrame(rows)
    return "\n".join(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure", default="8k5r.cif")
    ap.add_argument("--cdk", default="CDK9")
    ap.add_argument("--kinase-chain", default="A")
    ap.add_argument("--cyclin-chain", default="B")
    ap.add_argument("--cyclin-name", default="cyclin T1")
    ap.add_argument("-o", "--out", type=Path)
    a = ap.parse_args()
    out = a.out or Path(_paths.FIG_DIR) / f"{a.cdk}_cyclin_KLIFS_ribbon.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build(a), encoding="utf-8")
    side = out.with_name(out.stem + "_residue_xy.csv")
    build.resi_xy.to_csv(side, index=False)
    print(f"wrote {paths.rel(out)}")
    print(f"wrote {side}  ({len(build.resi_xy)} residues, depth: +ve = toward viewer)")


if __name__ == "__main__":
    main()
