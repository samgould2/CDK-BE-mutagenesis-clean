"""
Label the consensus panel with leader lines instead of labels pinned to the atoms.

    /Applications/ChimeraX-1.10.app/Contents/bin/ChimeraX --exit \
        --cmd "runscript label_consensus.py core7"

Residue-attached labels (the `label ... offset` approach used by
generate_cdk2_4_6_active_site.py + deconflict_labels.py) work when a handful of
residues are spread around the cleft.  The consensus panel stacks 25-56 positions
into one view, most of them inside the ATP cleft, so no offset keeps them legible
and the ones near the frame edge get clipped.

Here every label is parked in a left or right gutter outside the structure and
joined to its residue by a thin leader line, so nothing overlaps and nothing is
cut off.  Projection to normalized screen coordinates uses the same camera maths
as deconflict_labels.py: a point at depth d maps to cc[0]/(th*d), cc[1]/(tv*d),
with th = tan(fov/2) and tv = th * height/width.

Writes, next to the panel:
    consensus_leader.cxc   the 2dlabels commands, so the figure is reproducible
    and renders consensus_leader.png into figures/chimerax/<tag>/
"""
import collections, csv, math, os, re, sys
from xml.sax.saxutils import escape

import numpy as np
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
# The frame the published leader panels were composed in. Every label position
# is a fraction of this frame, so the shape matters as much as the size: the
# same run in a square window moves all of them and stretches the gaps by 1.76.
# ChimeraX honours a taller request on a taller display, which is how the
# published set came out 1400x796 while the main-figure renders are square --
# so this is set rather than left to the machine. `window=WxH` overrides it.
WINDOW = (1400, 796)

# gutter geometry, in fractions of the window
GUTTER_L, GUTTER_R = 0.035, 0.965
TOP, BOT = 0.965, 0.035
MIN_GAP = 0.034            # minimum vertical spacing between stacked labels
FONT = 13                  # upper bound; the actual size adapts to the longest label
CHAR_W_PER_PT = 0.00047    # width of one character per point of font size, window fractions
LINE_H_PER_PT = 0.0020     # height of one line of text, likewise per point
BLOCK_GAP_PER_PT = 0.0009  # blank space under a multi-line block, likewise per point
MAX_W = 0.30               # widest a label may get, as a fraction of the window
MAX_W_FULL = 0.26          # ... and in `full` mode, where the blocks are much wider
WRAP = 52                  # characters per line before a gene's variant list wraps

POCKET_RAMP = {1: "#fee5d9", 2: "#fcae91", 3: "#fb6a4a", 4: "#cb181d", 5: "#67000d"}
GAP_RAMP = {1: "#eff3ff", 2: "#c6dbef", 3: "#6baed6", 4: "#2171b5", 5: "#08306b"}


SAV_RE = re.compile(r"^([A-Z])(\d+)([A-Z*])$")

# the alignment's own gene order, so a variant block always reads in the same order
GENE_ORDER = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]


def group_variants(raw):
    """'CDK12-D817G,CDK12-D817N,CDK4-D97N' -> {'CDK12': ['D817G/N'], 'CDK4': ['D97N']}

    Several substitutions at one residue collapse onto a shared stem, which shortens
    the list without dropping any variant.  Genes come back in GENE_ORDER and each
    gene's residues in sequence order.
    """
    stems = collections.OrderedDict()           # (gene, wt, pos) -> [alt, ...]
    for v in raw.split(","):
        if not v:
            continue
        gene, _, sav = v.partition("-")
        m = SAV_RE.match(sav)
        key = (gene, m.group(1), int(m.group(2))) if m else (gene, sav, 0)
        stems.setdefault(key, []).append(m.group(3) if m else "")
    out = collections.OrderedDict()
    order = sorted({k[0] for k in stems},
                   key=lambda g: (GENE_ORDER.index(g) if g in GENE_ORDER else 99, g))
    for gene in order:
        out[gene] = [f'{wt}{pos}{"/".join(alts)}'
                     for (g, wt, pos), alts in sorted(stems.items(),
                                                      key=lambda kv: kv[0][2])
                     if g == gene]
    return out


def compress_variants(raw):
    """One-line form: 'CDK12 D817G/N, D819G, CDK4 D97N'."""
    return ", ".join(f"{g} {', '.join(svs)}"
                     for g, svs in group_variants(raw).items())


def variant_block(head, raw, width=WRAP):
    """Header, then every contributing variant, one gene per line (wrapped).

    Returns (indent_in_characters, text) pairs rather than pre-padded strings: the
    indent is applied as an x offset in both renderers, because neither ChimeraX
    2dlabels nor SVG <text> can be relied on to keep leading spaces.
    """
    lines = [(0, head)]
    for gene, svs in group_variants(raw).items():
        hang = len(gene) + 2            # continuation lines sit under the first variant
        indent, cur = 1, f"{gene} {svs[0]}"
        for sv in svs[1:]:
            if indent + len(cur) + 2 + len(sv) > width:
                lines.append((indent, cur + ","))
                indent, cur = hang, sv
            else:
                cur = f"{cur}, {sv}"
        lines.append((indent, cur))
    return lines


def load_sites(aset, reference="CDK2", min_cdks=2, sites_csv=None,
               tiers=("recurrent",)):
    src = sites_csv or os.path.join(ROOT, "data", "processed", "lfc_fdr_tables",
                                    "resistance_patterns", aset, "MSA_recurrent_columns.csv")
    rows = []
    for r in csv.DictReader(open(src)):
        tier = r.get("tier", "recurrent")
        if tier not in tiers:
            continue
        if tier == "recurrent" and int(r["n_cdks"]) < min_cdks:
            continue
        resi = None
        for bit in r["residues"].split(";"):
            g, _, aa = bit.partition(":")
            if g == reference:
                resi = int(aa[1:])
                break
        if resi is None:
            continue
        in_k = r["in_KLIFS"].lower() == "true"
        seg = (r.get("segment") or "").strip()
        reg = (r.get("region") or "").strip()
        label = (r["KLIFS_idx"] if in_k else
                 seg if seg else "C-lobe" if reg == "C-lobe post-pocket" else "gap")
        raw = r.get("variants") or ""
        rows.append(dict(order=len(rows) + 1,
                         resi=resi, n=int(r["n_cdks"]), in_klifs=in_k, label=label,
                         cdks=r["cdks"], variants=compress_variants(raw), raw=raw,
                         n_vars=len([v for v in raw.split(",") if v]),
                         rank=r.get("rank") or "", tier=tier))
    return rows


def spread(items, key_y, key_h=None):
    """Order by projected height, then push apart without leaving the frame.

    key_h gives each item's height in window fractions, so multi-line labels
    reserve the room they actually need.
    """
    key_h = key_h or (lambda i: MIN_GAP)
    items = sorted(items, key=key_y, reverse=True)
    ys = [min(TOP, max(BOT, key_y(i))) for i in items]
    for i in range(1, len(ys)):
        ys[i] = min(ys[i], ys[i - 1] - max(MIN_GAP, key_h(items[i - 1])))
    # if we ran off the bottom, push the whole stack back up (the last block's own
    # height counts, otherwise a tall block hangs out of the frame)
    if ys and ys[-1] - key_h(items[-1]) < BOT:
        shift = BOT - (ys[-1] - key_h(items[-1]))
        ys = [min(TOP, y + shift) for y in ys]
        for i in range(1, len(ys)):
            ys[i] = min(ys[i], ys[i - 1] - max(MIN_GAP, key_h(items[i - 1])))
    return list(zip(items, ys))


def main(aset="core7", reference="CDK2", min_cdks=2, mode="klifs",
         tag="consensus", sites_csv=None, tiers=("recurrent",),
         panel_name="consensus.cxc", out_stem="consensus_leader", zoom=1.0):
    panel = os.path.join(HERE, f"cxc_{aset}_{tag}", panel_name)
    if not os.path.exists(panel):
        _log(f"no panel at {panel}; run the generator for it first")
        return
    run(session, f"cd {HERE}")
    run(session, f"open {panel}")
    run(session, f"windowsize {WINDOW[0]} {WINDOW[1]}")
    # `full` blocks are far wider than a one-line label, so the structure is pulled
    # back first to free the gutters.  This has to happen before the projection below,
    # which reads the live camera.
    if zoom and abs(zoom - 1.0) > 1e-9:
        run(session, f"zoom {zoom}")

    sites = load_sites(aset, reference, min_cdks, sites_csv, tiers)
    model = [m for m in session.models if getattr(m, "residues", None) is not None][0]
    kc = "A"
    byres = {r.number: r for r in model.residues if r.chain_id == kc}

    cam = session.main_view.camera
    inv = cam.position.inverse()
    w, h = session.main_view.window_size or WINDOW
    th = math.tan(math.radians(cam.field_of_view) / 2)
    tv = th * (h / w)

    placed = []
    for s in sites:
        r = byres.get(s["resi"])
        if r is None:
            _log(f"  {reference} {s['resi']} not modelled - skipped")
            continue
        c = np.array([a.scene_coord for a in r.atoms]).mean(axis=0)
        cc = np.array(inv * tuple(c))
        d = max(-cc[2], 1e-6)
        x = (cc[0] / (th * d) + 1) / 2          # normalized half-frame -> 0..1
        y = (cc[1] / (tv * d) + 1) / 2
        placed.append(dict(**s, x=x, y=y))

    # ---- label text -------------------------------------------------------
    # Three shapes of label. `klifs` and `variants` are each one line; `full` is
    # a block -- a header, then every contributing variant, one gene per line,
    # wrapped at WRAP -- which is what the peak panel (supplementary 15e) uses
    # and what `variant_block`, MAX_W_FULL and BLOCK_GAP_PER_PT are here for.
    for p in placed:
        if mode == "klifs":
            p["lines"] = [(0, f'{p["label"]} {reference}-{p["resi"]} ({p["n"]})')]
        elif mode == "full":
            head = (f'P{p["order"]}  {p["label"]}  {reference}-{p["resi"]}  '
                    f'{p["n_vars"]} variants / {p["n"]} CDKs')
            p["lines"] = variant_block(head, p["raw"])
        else:
            p["lines"] = [(0, p["variants"]
                           or f'{p["label"]} {reference}-{p["resi"]}')]
        p["text"] = p["lines"][0][1]
    # the font is set from the longest line of any label: keep it inside the
    # mode's own width budget, and never go above FONT
    longest = max(ind + len(t) for p in placed for ind, t in p["lines"])
    budget = MAX_W_FULL if mode == "full" else MAX_W
    font = max(8, min(FONT, int(budget / (longest * CHAR_W_PER_PT))))
    char_w = CHAR_W_PER_PT * font
    line_h = LINE_H_PER_PT * font
    for p in placed:
        p["w"] = char_w * max(ind + len(t) for ind, t in p["lines"])
        # a label reserves one line of height per line of text, plus a blank
        # line's worth beneath it so two of them never touch
        p["h"] = len(p["lines"]) * line_h + BLOCK_GAP_PER_PT * font

    mid = float(np.median([p["x"] for p in placed]))
    left = [p for p in placed if p["x"] <= mid]
    right = [p for p in placed if p["x"] > mid]

    cmds = ["# leader-line labels for the consensus panel",
            "# generated by label_consensus.py -- edit that script, not this file",
            "2dlabels delete all"]
    placed_boxes = []        # (text, label_x, label_y, anchor_x, anchor_y, colour)
    for group, side in ((left, "L"), (right, "R")):
        for p, ly in spread(group, lambda i: i["y"], lambda i: i["h"]):
            shade = 1 if p["tier"] == "strong_singleton" else min(max(p["n"], 1), 5)
            col = (POCKET_RAMP if p["in_klifs"] else GAP_RAMP)[shade]
            lx = GUTTER_L if side == "L" else GUTTER_R - p["w"]
            lines = []
            for k, (ind, text) in enumerate(p["lines"]):
                tx, ty = lx + ind * char_w, ly - k * line_h
                cmds.append(f'2dlabels text "{text}" xpos {tx:.4f} ypos {ty:.4f} '
                            f'size {font} color black bold false')
                lines.append((text, tx, ty))
            # the arrow leaves the inner edge of the label, level with the middle
            # of the block -- for a one-line label that is the line itself
            ax = lx + p["w"] + 0.004 if side == "L" else lx - 0.004
            ay = ly + 0.008 - (len(p["lines"]) - 1) / 2 * line_h
            cmds.append(f'2dlabels arrow start {ax:.4f},{ay:.4f} '
                        f'end {p["x"]:.4f},{p["y"]:.4f} color {col} weight 0.06 '
                        f'headStyle blocky')
            placed_boxes.append((lines, ax, ay, p["x"], p["y"], col))
    suffix = "" if mode == "klifs" else "_variants"
    out_cxc = os.path.join(HERE, f"cxc_{aset}_{tag}", f"{out_stem}{suffix}.cxc")
    open(out_cxc, "w").write("\n".join(cmds) + "\n")
    for c in cmds:
        if not c.startswith("#"):
            run(session, c)

    fig = os.path.join(ROOT, "figures", "chimerax", f"cxc_{aset}_{tag}")
    os.makedirs(fig, exist_ok=True)
    dest = os.path.join(fig, f"{out_stem}{suffix}.png")
    run(session, f"save {dest} width 3000 supersample 3 transparentBackground true")

    # ---- layered vector version -------------------------------------------------
    # ChimeraX renders the scene with OpenGL and has no vector export, so the
    # structure stays a raster layer; the labels and leader lines are written as real
    # SVG <text> and <path>, which keeps them editable in Illustrator.
    run(session, "2dlabels delete all")
    # consensus_leader + "_variants" -> consensus_variants_plain.png. The suffix
    # has to be in the name: the two modes pull the camera back by different
    # amounts, so one plain render is not the other's background.
    plain_stem = out_stem.rsplit("_", 1)[0]       # consensus_leader -> consensus
    plain = os.path.join(fig, f"{plain_stem}{suffix}_plain.png")
    run(session, f"save {plain} width 3000 supersample 3 transparentBackground true")
    svg = write_svg(os.path.join(fig, f"{out_stem}{suffix}.svg"),
                    os.path.basename(plain), placed_boxes, font, (w, h))
    _log(f"wrote {svg}")
    _log(f"labelled {len(placed)} residues; wrote {out_cxc} and {dest}")
    print(f"labelled {len(placed)} residues")


SVG_SIDE = 1400           # layout space; the raster layer is scaled into it


def write_svg(path, image_name, boxes, font_pt, window=None):
    """Raster structure + vector labels and leader lines.

    `window` is the live ChimeraX window, which is where the viewBox gets its
    shape. It matters: label positions are fractions of the frame, so writing
    them into a square viewBox when the frame was 1400x796 stretches every
    vertical gap by 1.76 and leaves the structure image distorted underneath.
    """
    W = SVG_SIDE
    H = SVG_SIDE if not window else int(round(SVG_SIDE * window[1] / window[0]))
    # ChimeraX 2dlabels y runs from the bottom of the frame, SVG y from the top
    def X(v): return v * W
    def Y(v): return (1.0 - v) * H
    fs = 1.2 * font_pt      # point size -> px in the 1400-wide viewBox

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" '
           f'xmlns:xlink="http://www.w3.org/1999/xlink" viewBox="0 0 {W} {H}" '
           f'width="{W}" height="{H}" '
           f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">',
           '<g id="structure">',
           f'  <image xlink:href="{escape(image_name)}" x="0" y="0" '
           f'width="{W}" height="{H}"/>',
           '</g>',
           '<g id="leader-lines" fill="none" stroke-linecap="round">']
    for _lines, ax, ay, px, py, col in boxes:
        out.append(f'  <path d="M {X(ax):.1f},{Y(ay):.1f} L {X(px):.1f},{Y(py):.1f}" '
                   f'stroke="{col}" stroke-width="1.6"/>')
    out.append('</g>')
    out.append('<g id="labels" fill="#1a1a1a">')
    for lines, *_rest in boxes:
        for t, tx, ty in lines:
            out.append(f'  <text x="{X(tx):.1f}" y="{Y(ty):.1f}" '
                       f'font-size="{fs:.1f}">{escape(t)}</text>')
    out.append('</g>')
    out.append('</svg>')
    open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")
    return path


def _log(msg):
    # ChimeraX sends print() to its own log window, so progress and tracebacks are
    # also written here, where a terminal caller can actually read them
    with open(os.path.join(HERE, "label_consensus.log"), "a") as fh:
        fh.write(str(msg) + "\n")


# ChimeraX `runscript` executes this file without setting __name__ == "__main__",
# so the entry point runs at module level, as in deconflict_labels.py.
_args = [a for a in sys.argv[1:] if not a.endswith(".py")]
open(os.path.join(HERE, "label_consensus.log"), "w").close()
try:
    # args: <set> [mode] [tag] [tiers] [panel=<file.cxc>] [sites=<csv>]
    #              [stem=<output name>] [mincdks=<n>]
    #   mode  = klifs | variants
    #   tag   = output dir suffix, e.g. consensus, consensus_strong or peaks
    #   tiers = comma-separated, e.g. recurrent,strong_singleton
    # The keyword extras let the same leader-line machinery label the peak panel
    # (generate_peaks_cxc.py), whose files are named peaks.cxc rather than
    # consensus.cxc and whose site table has one row per detected peak.
    _kw = dict(a.split("=", 1) for a in _args if "=" in a)
    _pos = [a for a in _args if "=" not in a]
    _set = _pos[0] if _pos else "core7"
    mode = _pos[1] if len(_pos) > 1 else "klifs"
    tag = _pos[2] if len(_pos) > 2 else "consensus"
    tiers = tuple((_pos[3] if len(_pos) > 3 else "recurrent").split(","))
    sites_csv = _kw.get("sites") or (
        None if tiers == ("recurrent",) else
        os.path.join(ROOT, "data", "processed", "lfc_fdr_tables",
                     "resistance_patterns", _set, "consensus_sites.csv"))
    if sites_csv and not os.path.isabs(sites_csv):
        sites_csv = os.path.join(ROOT, "data", "processed", "lfc_fdr_tables",
                                 "resistance_patterns", _set, sites_csv)
    if "window" in _kw:
        WINDOW = tuple(int(v) for v in _kw["window"].lower().split("x"))
    main(_set, "CDK2", int(_kw.get("mincdks", 2)), mode, tag, sites_csv, tiers,
         _kw.get("panel", "consensus.cxc"), _kw.get("stem", "consensus_leader"),
         float(_kw.get("zoom", 1.0)))
except Exception:
    import traceback
    _log(traceback.format_exc())
    raise
