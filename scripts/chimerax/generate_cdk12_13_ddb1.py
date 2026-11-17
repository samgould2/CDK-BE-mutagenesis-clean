"""
Generate the CDK12 and CDK13 DDB1 / molecular-glue panels.

Two separate figures, both anchored on the same DDB1 frame:

    cxc/CDK12_DDB1.cxc   8BUG as deposited -- DDB1 + CDK12 + cyclin K + HQ461
    cxc/CDK13_DDB1.cxc   7NXJ (CDK13 + cyclin K) superposed onto 8BUG's CDK12,
                         so CDK13 is placed against the same DDB1 and glue

8BUG is "Structure of DDB1 bound to HQ461-engaged CDK12-cyclin K" (X-ray,
3.53 A), and its ligand RPW *is* HQ461 -- the glue this screen's HQ461 arm used.
(cdk_config.py's note for 8BUG calls the ligand CR8; that note is wrong, and
this script does not rely on it.)

No CDK12-DDB1-free equivalent exists for CDK13, so the CDK13 panel is a
superposition, not an observed complex:  CDK13 from 7NXJ is matched onto CDK12
in 8BUG, and 8BUG's DDB1 and HQ461 are kept in place.  Everything about where
CDK13 sits relative to DDB1 is therefore inferred from the CDK12 complex.
8BUG's own CDK12 and cyclin K are hidden in that panel so only one kinase and
one cyclin are drawn.

Hit coloring
------------
red   a resistance hit from CDK12_13_resistance_variants.csv
gold  an HQ461-specific hit (resistant in HQ461 and in neither other condition)

A residue is gold only when *every* variant seen at it comes from an
HQ461-specific guide.  CDK12 Y742 is the one position where that matters: it
carries Y742C (HQ461-specific, from the Y742C_K743E guide) and also Y742H
(BSJ-4-116-specific), so it is drawn red while its neighbour K743 is gold.

Usage:  python3 generate_cdk12_13_ddb1.py
        then, from the ChimeraX GUI:  open cxc/render_cdk12_13_ddb1.cxc
"""
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cdk_config import (STRUCTURES, CDK_COLORS, CYCLIN_COLOR, HIT_COLOR, HIT_STYLE,
                        NONCATALYTIC_TRANSPARENCY, RENDER_PREAMBLE)

HITS_CSV = os.path.normpath(os.path.join(
    HERE, "..", "..", "data", "processed", "lfc_fdr_tables", "CDK12_13_hits",
    "CDK12_13_resistance_variants.csv"))
RESIDUE_CSV = os.path.join(HERE, "cdk12_13_ddb1_hit_residues.csv")
OUT = os.path.join(HERE, "cxc")

# the glue scaffold: DDB1 and HQ461 both come from 8BUG chain A / B:1101
SCAFFOLD = dict(pdb="8BUG", file="8bug.cif", ddb1_chain="A", cdk_chain="B",
                cyclin_chain="C", glue="RPW", glue_chain="B")

# DDB1 is context, not the subject -- mostly transparent white so the kinase and
# the glue site read through it
# DDB1 is context, not the subject -- but a *white* surface on a white background
# is invisible whatever its alpha: it reads as a bare silhouette outline, and
# changing the transparency percentage does essentially nothing to it. Grey is
# what makes the subunit read as a body. 70% transparency keeps it well
# subordinate to the kinase -- a readable shape rather than a competing object --
# and lets the interface hits that sit inside its envelope (CDK13 E803,
# CDK12 K743/G822) show through.
DDB1_COLOR = "#b8b8b8"
DDB1_TRANSPARENCY = 70

GLUE_CARBON = "#909090"
GLUE_STYLE = "stick"
HQ_HIT_COLOR = "#ffd700"          # gold, for the HQ461-specific positions

# Orientation from compute_ddb1_camera.py: DDB1 left, kinase centre, cyclin K
# right, viewed side-on.  Rotation only -- the panels call `view` right after to
# set the distance, which is far more robust than computing one here.
CAMERA_ROTATION = ("0.2461,-0.4802,0.8419,0.0000,"
                   "0.2931,0.8648,0.4076,0.0000,"
                   "-0.9239,0.1465,0.3536,0.0000")

# DDB1 is byte-identical between the two panels, so framing on it and then
# pulling back gives both figures the same camera. Fitting "everything displayed"
# instead would let the two kinases' different sizes shift the framing.
# Frame on 8BUG's DDB1 *and* cyclin K (chains A and C).  Both are present in
# both panels -- the CDK13 panel opens 8BUG too and merely hides its kinase and
# cyclin -- so the framing is identical, and spanning the two ends of the
# assembly centres the composition instead of leaving DDB1 in the middle with
# empty space to its left.
FRAME_SPEC = "#1/A,C"
FRAME_ZOOM = 0.95

LABEL_HEIGHT = 1.6
LABEL_COLOR = "black"

PANELS = {
    "CDK12": dict(
        gene="CDK12", pdb="8BUG", file="8bug.cif", chain="B", cyclin="C",
        modelled=True,
        note="8BUG as deposited: DDB1 (A) + CDK12 (B) + cyclin K (C) + HQ461 (RPW B/1101). "
             "Copy 1 of 3; HQ461 is modelled in this copy."),
    "CDK13": dict(
        gene="CDK13", pdb="7NXJ", file="7nxj.cif", chain="A", cyclin="B",
        modelled=False,
        note="7NXJ (CDK13 + cyclin K) superposed onto 8BUG chain B. The DDB1 and HQ461 "
             "shown are 8BUG's; CDK13's placement against them is inferred, not observed."),
}


def load_hits():
    """{gene: {pos: {'wt','variants','gold'}}} from the resistance-variant table."""
    by_gene = {}
    for r in csv.DictReader(open(HITS_CSV)):
        hq = r["which_conditions"] == "HQ461"
        for var in r["Top_HGVSp"].split("_"):
            m = re.match(r"^([A-Z*])(\d+)([A-Z*])$", var)
            if m is None:
                continue
            pos = int(m.group(2))
            d = by_gene.setdefault(r["Gene"], {}).setdefault(
                pos, dict(wt=m.group(1), variants=set(), hq=[]))
            d["variants"].add(var)
            d["hq"].append(hq)

    out = {}
    for gene, positions in by_gene.items():
        out[gene] = {p: dict(wt=d["wt"],
                             variants=",".join(sorted(d["variants"])),
                             gold=all(d["hq"]))
                     for p, d in sorted(positions.items())}
    return out


def write_residue_table(hits):
    """Provenance table next to the scripts, same idea as cdk7_sy5609_hit_residues.csv."""
    with open(RESIDUE_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "position", "wt", "variants", "class"])
        for gene in ("CDK12", "CDK13"):
            for p, h in hits[gene].items():
                w.writerow([gene, p, h["wt"], h["variants"],
                            "HQ461_specific" if h["gold"] else "resistance"])
    print("wrote", RESIDUE_CSV)


def build(gene, hits, labeled=False, opaque=False):
    cfg = PANELS[gene]
    kd_lo, kd_hi = STRUCTURES[gene]["kd"]
    color = CDK_COLORS[gene]
    borrowed = not cfg["modelled"]

    positions = hits[gene]
    red = [p for p, h in positions.items() if not h["gold"]]
    gold = [p for p, h in positions.items() if h["gold"]]

    # the kinase model is #1 for CDK12 (8BUG itself) and #2 for CDK13 (7NXJ on top)
    km = "#2" if borrowed else "#1"
    kc = cfg["chain"]
    ks = f"{km}/{kc}"
    cs = f"{km}/{cfg['cyclin']}"

    name = gene + "_DDB1" + ("_opaque" if opaque else "") + ("_labeled" if labeled else "")

    L = []
    L.append("# " + "=" * 74)
    L.append(f"# {gene} resistance variants on DDB1 + HQ461"
             + (", whole kinase opaque" if opaque else "")
             + (", residues labeled" if labeled else ""))
    L.append("# Generated by generate_cdk12_13_ddb1.py -- edit that script, not this file.")
    L.append(f"# {cfg['note']}")
    L.append("# red = resistance hit; gold = HQ461-specific hit")
    L.append("# " + "=" * 74)
    L.append("close session")
    L.append("cd " + HERE)
    L.append("")

    # ---- glue scaffold: DDB1 + HQ461, always model #1
    L.append("# ---------------------------------------- DDB1 + HQ461 scaffold (8BUG)")
    L.append(f"open structures/{SCAFFOLD['file']} name DDB1_scaffold")
    L.append("hide #1 target acsb")
    L.append("hide #1 models")

    if borrowed:
        L.append("")
        L.append(f"# ------------------------------------------------- {gene} ({cfg['pdb']})")
        L.append(f"open structures/{cfg['file']} name {gene}")
        L.append("hide #2 target acsb")
        L.append("hide #2 models")
        L.append("")
        L.append("# place CDK13 where CDK12 sits in the DDB1 complex. matchmaker moves the")
        L.append("# first spec onto the second, so 8BUG (and therefore DDB1 and HQ461)")
        L.append("# stays put and only 7NXJ is transformed.")
        L.append(f"matchmaker {ks} to #1/{SCAFFOLD['cdk_chain']}")

    # ---- DDB1 surface, mostly transparent white
    ds = f"#1/{SCAFFOLD['ddb1_chain']}"
    L.append("")
    L.append("# DDB1 as a mostly transparent white surface: context for the glue interface")
    L.append(f"show {ds} atoms")
    L.append(f"surface {ds}")
    L.append(f"hide {ds} atoms")
    L.append(f"color {ds} {DDB1_COLOR} target s")
    L.append(f"transparency {ds} {DDB1_TRANSPARENCY} target s")
    L.append(f"name DDB1 {ds}")

    # ---- HQ461 sticks
    gs = f"#1/{SCAFFOLD['glue_chain']} & :{SCAFFOLD['glue']}"
    L.append("")
    L.append("# HQ461 (ligand RPW) -- the molecular glue, at the CDK12/DDB1 interface.")
    L.append("# Element-coloured sticks on grey carbons, as in the per-kinase panels.")
    L.append(f"show {gs} target ab")
    L.append(f"style {gs} {GLUE_STYLE}")
    L.append(f"color {gs} {GLUE_CARBON} target ab")
    L.append(f"color {gs} byhetero")
    L.append(f"name HQ461 {gs}")

    # ---- kinase cartoon
    L.append("")
    if opaque:
        L.append(f"# {gene} cartoon, fully opaque -- no catalytic-domain fade")
        L.append(f"show {ks} cartoon")
        L.append(f"color {ks} {color} target c")
    else:
        L.append(f"# {gene} cartoon, catalytic core opaque and the rest faded back")
        L.append(f"show {ks} cartoon")
        L.append(f"color {ks} {color} target c")
        L.append(f"transparency {ks} & ~:{kd_lo}-{kd_hi} "
                 f"{NONCATALYTIC_TRANSPARENCY} cartoons")
    L.append(f"name {gene}_kinase {ks}")

    # ---- cyclin K surface
    L.append("")
    L.append("# cyclin K as a wheat surface, as in every other panel")
    L.append(f"show {cs} atoms")
    L.append(f"surface {cs}")
    L.append(f"hide {cs} atoms")
    L.append(f"color {cs} {CYCLIN_COLOR} target s")
    L.append(f"name {gene}_cyclin {cs}")

    if borrowed:
        L.append("")
        L.append("# 7NXJ's own inhibitor (5I1, THZ531) is an ATP-site covalent ligand, not a")
        L.append("# glue; it is left hidden so HQ461 is the only ligand drawn.")

    # ---- hits
    L.append("")
    L.append(f"# {len(red)} resistance positions in red")
    if red:
        hs = f"{ks}:{','.join(str(p) for p in red)}"
        L.append("# " + ", ".join(f"{positions[p]['wt']}{p} {positions[p]['variants']}"
                                  for p in red))
        L.append(f"show {hs} atoms")
        L.append("# ChimeraX hides backbone atoms under a cartoon, so a glycine hit would")
        L.append("# otherwise show nothing but a recolored ribbon")
        L.append(f"cartoon {hs} suppressBackboneDisplay false")
        L.append(f"style {hs} {HIT_STYLE}")
        L.append(f"color {hs} {HIT_COLOR} target ac")
        L.append(f"name {gene}_hits {hs}")

    L.append("")
    L.append(f"# {len(gold)} HQ461-specific positions in gold")
    if gold:
        gh = f"{ks}:{','.join(str(p) for p in gold)}"
        L.append("# " + ", ".join(f"{positions[p]['wt']}{p} {positions[p]['variants']}"
                                  for p in gold))
        L.append(f"show {gh} atoms")
        L.append(f"cartoon {gh} suppressBackboneDisplay false")
        L.append(f"style {gh} {HIT_STYLE}")
        L.append(f"color {gh} {HQ_HIT_COLOR} target ac")
        L.append(f"name {gene}_HQ461_hits {gh}")

    # ---- shared look and camera
    L.append("")
    L.append("# ---------------------------------------------------- look and camera")
    L.extend(RENDER_PREAMBLE.strip().splitlines())
    L.append("# the models were opened hidden; make them visible again. The earlier")
    L.append("# `hide #1 target acsb` still applies, so for the CDK13 panel only DDB1 and")
    L.append("# HQ461 come back -- 8BUG's own CDK12 and cyclin K stay hidden.")
    L.append("show #1 models")
    if borrowed:
        L.append("show #2 models")
    L.append("# mmCIF chain breaks otherwise draw dashed pseudobonds labelled")
    L.append('# "N residues", which land on top of the kinase in the render')
    L.append("label delete")
    L.append("hide #* target p")
    L.append("hide @@element=H target ab")

    if labeled:
        L.append("")
        L.append("# --------------------------------------------- residue labels")
        L.append("# emitted after `label delete`, which would otherwise remove them")
        for p, h in positions.items():
            col = LABEL_COLOR
            L.append(f'label {ks}:{p} residues text "{h["wt"]}{p}" '
                     f'height {LABEL_HEIGHT} color {col} bgColor none')

    L.append("")
    L.append("# Orientation first, then let `view` set the distance: `view` moves the")
    L.append("# camera along its own axis and leaves the rotation alone.")
    L.append(f"view matrix camera {CAMERA_ROTATION}")
    L.append("# `view` measures displayed *atoms*, and DDB1's are hidden under its surface,")
    L.append("# so show them for the one command and hide them again straight after.")
    L.append(f"show {FRAME_SPEC} atoms")
    L.append(f"view {FRAME_SPEC}")
    L.append(f"hide {FRAME_SPEC} atoms")
    L.append("# pull back so cyclin K, off to the right of DDB1, stays in frame")
    L.append(f"zoom {FRAME_ZOOM}")

    L.append("")
    L.append("# `view <spec>` sets the near/far clipping planes to that spec's own bounds,")
    L.append("# which slices any geometry nearer the camera than the target -- CDK13 T715")
    L.append("# was being cut in half by the near plane in the HQ461 zoom. Reset the slab")
    L.append("# to the whole scene. This has to come AFTER every view/zoom command, so if")
    L.append("# you re-aim the camera later, re-issue `clip off` afterwards.")
    L.append("clip off")
    L.append("")
    L.append(f"# zoom on the glue:  view HQ461 ; zoom 0.30 ; clip off")
    L.append("")
    L.append(f"# save figures/chimerax/{name}.png width 3000 supersample 3 "
             "transparentBackground true")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{name}.cxc")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote", path)
    return name


def write_render_script(names):
    fig = os.path.normpath(os.path.join(HERE, "..", "..", "figures", "chimerax"))
    R = ["# " + "=" * 74,
         "# Render the CDK12 / CDK13 DDB1 + HQ461 panels at publication size",
         "# Generated by generate_cdk12_13_ddb1.py -- edit that script, not this file.",
         "# ChimeraX needs a real GPU context for `save`, so run this from the ChimeraX GUI",
         '#   (or `ChimeraX --exit --cmd "open cxc/render_cdk12_13_ddb1.cxc"`), not --nogui.',
         "# " + "=" * 74,
         "close session",
         "cd " + HERE,
         "windowsize 1400 1400",
         ""]
    for name in names:
        gene = name.split("_")[0]
        R += [f"open cxc/{name}.cxc",
              f"save {fig}/{name}.png width 3000 supersample 3 transparentBackground true",
              "",
              "# Zoom on the glue site. Framing the hits instead does not work: CDK12's",
              "# hits are a tight cluster but CDK13's span 711-1004, so `view <gene>_hits`",
              "# frames most of the kinase and ends up wider than the panel above.",
              "# HQ461 is the same object in both panels, so this crop is comparable.",
              "view HQ461",
              "zoom 0.30",
              "# re-aiming the camera re-tightens the clipping slab; reset it again",
              "clip off",
              f"save {fig}/{name}_glue_zoom.png width 3000 supersample 3 "
              "transparentBackground true",
              ""]
    R += ["close session", ""]

    rpath = os.path.join(OUT, "render_cdk12_13_ddb1.cxc")
    with open(rpath, "w") as fh:
        fh.write("\n".join(R))
    print("wrote", rpath)


if __name__ == "__main__":
    hits = load_hits()
    write_residue_table(hits)
    names = []
    for gene in ("CDK12", "CDK13"):
        for opq in (False, True):
            for lab in (False, True):
                names.append(build(gene, hits, labeled=lab, opaque=opq))
        n = hits[gene]
        g = [f"{h['wt']}{p}" for p, h in n.items() if h["gold"]]
        print(f"  {gene}: {len(n)} positions, {len(g)} HQ461-specific (gold): {', '.join(g)}")
    write_render_script(names)
