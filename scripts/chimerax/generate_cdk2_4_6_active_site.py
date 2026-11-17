"""
Generate the zoomed ATP-site panels for the CDK2/4/6 resistance variants.

    python3 build_cdk2_4_6_hits.py                       # hit table
    ChimeraX --nogui --exit --cmd "runscript $PWD/analyze_cdk2_4_6_active_site.py"
    ChimeraX --nogui --exit --cmd "runscript $PWD/compute_active_site_camera.py"
    python3 generate_cdk2_4_6_active_site.py             # -> cxc/*_active_site*.cxc

then, from the ChimeraX GUI (`save` needs a GPU context):

    open cxc/render_cdk2_4_6_active_site.cxc

Writes, per gene, `{GENE}_active_site.cxc` and `{GENE}_active_site_labeled.cxc`.

Design notes
------------
* **Framing is on the kinase chain** (FRAME_ON), which puts the whole kinase and its
  cyclin in frame with margin, seen straight-on down the cleft.  Framing on the ATP
  pocket instead gives a true close-up but crops the surrounding kinase away on all
  four sides, so the view stops being placeable; that option is one constant away.
  Framing on the *hits* is not usable either -- CDK2's top hits span 2.7-30 A from the
  ATP, so `view <hits>` takes in everything.  Same reasoning the CDK12/13 panels use
  for framing on the glue rather than the hits.
* **The camera is the occlusion-derived one** from compute_active_site_camera.py,
  not cdk_config.CANONICAL_VIEW -- the canonical camera frames the whole kinase
  and does not look into the cleft.  Each panel carries the canonical matrix on a
  commented line so it can be swapped back with one edit.
* **The cyclin is drawn as a semi-transparent wheat surface** (opaque in the
  whole-kinase panels).  At this crop an opaque surface covers the cleft, so the
  transparency keeps the interface visible without hiding the active site.
* **CDK4 (2W96) is apo** -- it has no ligand to mark its pocket.  The superposed
  1QMZ ATP is drawn there instead, at ATP_BORROWED_TRANSPARENCY, following the
  same convention `cdk_config` uses for CDK19's borrowed cyclin: an inferred
  component is shown semi-transparent so it reads as inferred.  CDK2's ATP and
  CDK6's aminopurvalanol are native and drawn opaque.
* **Every modelled hit is labeled.**  An earlier version labeled only the hits
  within 16 A of the ATP, which left hits just past that cutoff (CDK6 L33, 16.5 A)
  drawn in frame as unlabeled spheres.  A label on a hit outside the frame does not
  render, so labeling all of them cannot under-label and cannot over-draw.  The
  distance and crop/off-crop columns stay in cdk2_4_6_hit_distances.csv.
"""
import csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# Optional variant tag.  Set CDK246_TAG to run the whole chain a second time against a
# different hit set without overwriting the first: every file this script reads or writes
# picks up the suffix.  Used for the stringent-FDR variant (CDK246_TAG=_strictFDR).
TAG = os.environ.get("CDK246_TAG", "")

sys.path.insert(0, HERE)
from cdk_config import (STRUCTURES, ALIGN_REFERENCE, CDK_COLORS, CYCLIN_COLOR,
                        HIT_COLOR, HIT_STYLE, NONCATALYTIC_TRANSPARENCY,
                        RENDER_PREAMBLE, CANONICAL_VIEW, ZOOM_OUT)

GENES = ["CDK2", "CDK4", "CDK6"]
OUT = os.path.join(HERE, "cxc")

LABEL_HEIGHT = 3.2          # scene units (A), so labels scale with the structure.
# The zoom camera sits 74.6 A out against the wide panel's ~270 A, i.e. ~3.6x the
# magnification, and scene-unit labels scale with it -- 3.2 A renders enormous and
# clipped there.  3.2 / 3.6 ~= 0.9 keeps the labels the same apparent size.
LABEL_HEIGHT_ZOOM = 0.9
# The shared camera sits 248 A back, further than the per-gene framing it replaced, so
# scene-unit labels render smaller: 1.8 was unreadable, 2.6 was marginal.
LABEL_COLOR = "black"
# Zoom applied after framing the pocket.  A sweep (0.8 -> 2.0, checking the camera
# frustum against every labeled hit) found the tightest value that keeps all of them
# in frame is 1.6 for CDK2, 1.8 for CDK6 and 2.0 for CDK4; 1.5 is the common value
# with margin for the labels, which stick out past the atoms they sit on.
# The first version used 0.80 -- a pull-back, which made this barely a zoom at all.
# 1.5 proved tight enough to crop away the surrounding kinase, so 1.25 is the
# compromise: the cleft fills the frame but the lobes around it are still readable.
# Per-gene zoom, measured by fit_active_site_zoom.py so each panel's complex fills
# ~88% of the frame.  A shared value cannot do this: at identical camera and framing
# CDK2's complex spans 0.90 of the half-frame but CDK4's and CDK6's only 0.53-0.55,
# because 1QMZ's cyclin A is far larger in projection than cyclin D1 or V-cyclin.
# TRADE-OFF: the three panels are therefore NOT at a common magnification, unlike the
# whole-kinase panels in cxc/{GENE}.cxc which share one camera and pull-back on purpose.
# Sphere and label sizes are comparable within a panel, not between them.  Set
# ZOOM_SHARED instead if cross-panel scale matters more than filling the frame.
PER_GENE_ZOOM = {"CDK2": 1.22, "CDK4": 2.07, "CDK6": 2.01}
# Set to a number to give all three panels ONE camera distance.  This is what makes the
# panels genuinely equivalent: the rotation is already shared (every kinase is
# matchmaker-superposed onto 1QMZ chain A), but a per-gene zoom changes the camera
# DISTANCE, and with a 30 deg field of view that changes the perspective foreshortening
# and the scale too -- so the three stopped being comparable even though the orientation
# was identical.  With FRAME_ON="kd" the camera also centres on the same superposed
# object in all three, so the cameras are identical, not merely parallel.
ZOOM_SHARED = 0.85
# What the camera frames.  Three options, and the choice matters more than it looks:
#   "kd" (used)       the catalytic domain, cdk_config's kd range.  These ranges were
#                     derived from the superposition (map_kd_ranges.py), so they are
#                     the same structural object in all three genes -- framing on them
#                     puts the camera in an identical place for every panel.
#   "complex"         kinase + cyclin together.  Fills the frame nicely per panel, but
#                     each gene's complex has a different size and centroid, so every
#                     panel gets a different camera distance and centre.
#   "kinase"          frames the kinase alone.  The cyclin then spills off to one side,
#                     which both pushes the composition off-centre and leaves a wide
#                     white margin, because the camera is sized to the kinase while the
#                     scene is much wider than it.
#   "pocket"          the ATP-site residues: a true close-up, but the surrounding kinase
#                     runs off all four edges and the view stops being placeable.
FRAME_ON = "kd"
ATP_BORROWED_TRANSPARENCY = 55   # CDK4 only; flags the ATP as superposed, not observed
# The cyclin is wheat and opaque in the whole-kinase panels.  At this zoom an opaque
# surface covers the cleft, so it is drawn semi-transparent here: still reads as the
# cyclin interface, without hiding what the panel is about.
CYCLIN_TRANSPARENCY = 0
LIGAND_STICK_RADIUS = 0.4   # default sticks vanish at complex scale


def read_camera():
    """The repo's house camera, so these panels line up with every other figure.

    cdk_config.CANONICAL_VIEW comes from viz2_cdk2.cxs and is what cxc/{GENE}.cxc,
    CDK7_SY5609 and the rest already use, with the same ZOOM_OUT pull-back.  Because
    every structure is superposed onto 1QMZ chain A, one matrix orients all of them.

    compute_active_site_camera.py can derive a bespoke active-site camera instead
    (occlusion-optimised, or constrained to put the hinge up and the cyclin right).
    Those are kept for reference in active_site_camera.txt but are NOT used: a panel
    in its own private orientation cannot be compared with the rest of the figure set,
    which is the whole point of the canonical camera.
    """
    return CANONICAL_VIEW


def read_zoom_cameras():
    """per-gene ligand-centred zoom cameras from compute_active_site_camera.py"""
    txt = open(os.path.join(HERE, "active_site_camera.txt")).read()
    out = {}
    for g in GENES:
        m = re.search(rf'ZOOM_CAMERA_{g} = "([^"]+)"', txt)
        if m:
            out[g] = m.group(1)
    return out


def read_csv(name):
    return list(csv.DictReader(open(os.path.join(HERE, name))))


def build(gene, camera, hits, pocket, labeled, zoom_cam=None):
    """zoom_cam set -> the ligand-centred zoom panel, with its own label height and
    its own offsets (the zoom_* columns), and only the hits that are in frame there."""
    zoom = ZOOM_SHARED if ZOOM_SHARED else PER_GENE_ZOOM[gene]
    cfg = STRUCTURES[gene]
    ref = STRUCTURES[ALIGN_REFERENCE]
    kc, color = cfg["kinase_chain"], CDK_COLORS[gene]
    lo, hi = cfg["kd"]
    is_ref = (gene == ALIGN_REFERENCE)
    gid = 1 if is_ref else 2          # model number of the gene structure
    name = f"{gene}_active_site{TAG}" + ("_zoom_labeled" if zoom_cam else
                                        ("_labeled" if labeled else ""))
    lab_h = LABEL_HEIGHT_ZOOM if zoom_cam else LABEL_HEIGHT
    pre = "zoom_" if zoom_cam else "label_"

    gh = [h for h in hits if h["gene"] == gene and h["modelled"] == "True"]
    inside = [h for h in gh if h["in_crop"] == "True"]
    resis = ",".join(sorted({h["position"] for h in gh}, key=int))
    pock = [p for p in pocket if p["gene"] == gene][0]["residues"]

    L = ["# " + "=" * 74,
         f"# {gene} ATP-site zoom with the top resistance variants ({cfg['pdb']})"
         + (", residues labeled" if labeled else ""),
         "# Generated by generate_cdk2_4_6_active_site.py -- edit that script, not this file.",
         f"# {cfg['note']}",
         "# " + "=" * 74,
         "close session",
         "cd " + HERE,
         ""]

    if not is_ref:
        L += ["# reference loaded hidden: supplies the common frame, and its ATP marks the pocket",
              f"open structures/{ref['file']} name {ALIGN_REFERENCE}_ref",
              "hide #1 target acsb",
              "hide #1 models",
              ""]

    L += [f"# ---------------------------------------------------------------- {gene}",
          f"open structures/{cfg['file']} name {gene}",
          f"hide #{gid} target acsb",
          f"hide #{gid} models",
          f"show #{gid}/{kc} cartoon",
          f"color #{gid}/{kc} {color} target c",
          f"transparency #{gid}/{kc} & ~:{lo}-{hi} {NONCATALYTIC_TRANSPARENCY} cartoons"]

    if not is_ref:
        L += ["", f"matchmaker #{gid}/{kc} to #1/{ref['kinase_chain']}"]

    # ---- cyclin partner
    cs = f"#{gid}/{','.join(cfg['cyclin_chains'])}"
    L += ["",
          f"# {cfg['cyclin_name']} as a semi-transparent surface -- context for the",
          "# interface without covering the cleft at this zoom",
          f"show {cs} atoms", f"surface {cs}", f"hide {cs} atoms",
          f"color {cs} {CYCLIN_COLOR} target s",
          f"transparency {cs} {CYCLIN_TRANSPARENCY} target s",
          f"name {gene}_cyclin {cs}"]

    # ---- what marks the pocket.  Collected into `lig_cmds` because the framing step
    # below hides every atom of the framing spec, and for FRAME_ON="complex" that spec
    # is the whole chain -- which includes the ligand's HET residue.  These have to be
    # re-issued afterwards or the ligand silently disappears from the figure.
    lig_cmds = []
    L.append("")
    if cfg["ligand"]:
        ls = f"#{gid}/{kc} & :{cfg['ligand']}"
        lig_cmds += [f"show {ls} target ab", f"style {ls} stick",
                     f"size {ls} stickRadius {LIGAND_STICK_RADIUS}",
                     f"color {ls} #909090 target ab", f"color {ls} byhetero"]
        if gene == "CDK2":
            lig_cmds += [f"show #{gid}/{kc} & :MG atoms",
                         f"style #{gid}/{kc} & :MG sphere",
                         f"size #{gid}/{kc} & :MG atomRadius 0.8"]
        L += [f"# {cfg['ligand']} is native to {cfg['pdb']}"] + lig_cmds
    else:
        L += [f"# {cfg['pdb']} is apo -- no ligand marks its pocket, so the superposed",
              "# 1QMZ ATP is drawn instead, semi-transparent because it is inferred",
              "# rather than observed (same convention as CDK19's borrowed cyclin)",
              "show #1 models"]
        lig_cmds += ["show #1/A & :ATP target ab", "style #1/A & :ATP stick",
                     f"size #1/A & :ATP stickRadius {LIGAND_STICK_RADIUS}",
                     "color #1/A & :ATP #909090 target ab",
                     "color #1/A & :ATP byhetero",
                     f"transparency #1/A & :ATP {ATP_BORROWED_TRANSPARENCY} target ab"]
        L += lig_cmds

    # ---- hits
    hs = f"#{gid}/{kc}:{resis}"
    L += ["",
          f"# {len(gh)} modelled hit positions: "
          + ", ".join(f'{h["wt"]}{h["position"]}' for h in gh),
          f"show {hs} atoms",
          "# ChimeraX hides backbone atoms under a cartoon, so a glycine hit would",
          "# otherwise draw nothing but a recolored ribbon",
          f"cartoon {hs} suppressBackboneDisplay false",
          f"style {hs} {HIT_STYLE}",
          f"color {hs} {HIT_COLOR} target ac",
          f"name {gene}_hits {hs}",
          f"name {gene}_pocket #{gid}/{kc}:{pock}",
          f"name {gene}_kinase #{gid}/{kc}",
          f"name {gene}_complex #{gid}/{kc},{','.join(cfg['cyclin_chains'])}",
          f"name {gene}_kd #{gid}/{kc}:{lo}-{hi}"]

    # ---- shared look
    L += ["", "# ---------------------------------------------------- look and camera"]
    L += RENDER_PREAMBLE.strip().splitlines()
    L += [f"show #{gid} models",
          "# mmCIF chain breaks otherwise draw dashed pseudobonds labelled",
          '# "N residues", which land on top of the kinase in the render',
          "label delete",
          "hide #* target p",
          "hide @@element=H target ab"]

    # ---- labels, after `label delete` so they survive it
    if labeled:
        L += ["", "# --------------------------------------------- residue labels",
              f"# EVERY modelled hit is labeled ({len(gh)} of them), not just the",
              f"# {len(inside)} inside the nominal 16 A crop.  Labeling only the crop set left"
              if not zoom_cam else
              "# Only the hits in frame at the zoom camera are labeled; the rest are off-crop.",
              "# hits that sit just past the cutoff (CDK6 L33 at 16.5 A) visible in frame",
              "# as unlabeled spheres, which reads as a bug.  A label on a hit that falls",
              "# outside the frame simply does not render, so this cannot under-label.",
              "# The crop/off-crop split is kept in cdk2_4_6_hit_distances.csv for reference."]
        lab_rows = [h for h in gh
                    if not zoom_cam or h.get("zoom_in_frame") == "True"]
        for h in sorted(lab_rows, key=lambda h: float(h["min_dist_to_ATP"])):
            off = ""
            if h.get(pre + "dx"):
                off = f' offset {h[pre+"dx"]},{h[pre+"dy"]},{h[pre+"dz"]}'
            L.append(f'label #{gid}/{kc}:{h["position"]} residues '
                     f'text "{h["wt"]}{h["position"]}" height {lab_h} '
                     f'color {LABEL_COLOR} bgColor none{off}')

    # ---- camera: the repo's canonical view plus the standard pull-back, exactly as
    # cxc/{GENE}.cxc and the CDK7_SY5609 panels do, so these line up with the rest of
    # the figure set.  Nothing per-gene is applied -- every structure is superposed
    # onto 1QMZ chain A, so one matrix orients all three identically.
    L += ["", f"view matrix camera {zoom_cam if zoom_cam else camera}"]
    if not zoom_cam:
        L += ["# uniform pull-back so no panel clips at the frame edge", f"zoom {ZOOM_OUT}"]
    L += [
          "# `view matrix camera` leaves the clipping planes alone, but reset them anyway",
          "# so this panel does not inherit a slab from whatever ran before it.",
          "clip off",
          "",
          f"# save figures/chimerax/{name}.png width 3000 supersample 3 transparentBackground true"]

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"{name}.cxc")
    open(path, "w").write("\n".join(L) + "\n")
    return path, len(gh), len(inside)


def write_render_script(zoom_cams):
    fig = os.path.normpath(os.path.join(HERE, "..", "..", "figures", "chimerax"))
    R = ["# " + "=" * 74,
         "# Render the CDK2/4/6 ATP-site panels to figures/chimerax/",
         "# Generated by generate_cdk2_4_6_active_site.py -- edit that script, not this file.",
         "# ChimeraX needs a real GPU context for `save`, so run this from the ChimeraX GUI",
         '#   (or `ChimeraX --exit --cmd "open cxc/render_cdk2_4_6_active_site.cxc"`), not --nogui.',
         "# " + "=" * 74,
         "close session",
         "cd " + HERE,
         "# square frame, as requested; matches the 1400x1400 the other panels use",
         "windowsize 1400 1400",
         ""]
    for gene in GENES:
        for name in (f"{gene}_active_site{TAG}", f"{gene}_active_site{TAG}_labeled"):
            R += [f"open cxc/{name}.cxc",
                  f"save {fig}/{name}.png width 3000 height 3000 supersample 3 "
                  "transparentBackground true",
                  ""]
            # Zoom only from the UNLABELED panel.  Labels are sized in scene units, so at
            # the zoom camera's 74.6 A they render ~3.6x larger than in the wide panel --
            # enormous and clipped at the frame edges.  A labeled zoom would need its own
            # label height and its own deconfliction pass against the zoom camera.
            if gene in zoom_cams and not name.endswith("_labeled"):
                R += ["# ligand-centred zoom: same rotation, moved in front of the ligand.",
                      "# `clip off` MUST be re-issued after re-aiming the camera -- see the",
                      "# README: a view onto a small ligand leaves a thin slab that slices",
                      "# any residue nearer the camera than it.",
                      f"view matrix camera {zoom_cams[gene]}",
                      "clip off",
                      f"save {fig}/{name}_zoom.png width 3000 height 3000 supersample 3 "
                      "transparentBackground true",
                      ""]
    for gene in GENES:
        if gene in zoom_cams:
            zn = f"{gene}_active_site{TAG}_zoom_labeled"
            R += [f"open cxc/{zn}.cxc",
                  f"save {fig}/{zn}.png width 3000 height 3000 supersample 3 "
                  "transparentBackground true",
                  ""]
    R += ["close session", ""]
    path = os.path.join(OUT, f"render_cdk2_4_6_active_site{TAG}.cxc")
    open(path, "w").write("\n".join(R))
    return path


if __name__ == "__main__":
    camera = read_camera()
    hits = read_csv(f"cdk2_4_6_hit_distances{TAG}.csv")
    pocket = read_csv(f"cdk2_4_6_pocket_residues{TAG}.csv")
    zoom_cams = read_zoom_cameras()
    for gene in GENES:
        for labeled in (False, True):
            path, n_hits, n_lab = build(gene, camera, hits, pocket, labeled)
            print(f"wrote {os.path.relpath(path, HERE)}  ({n_hits} hits, {n_lab} labeled)")
        if gene in zoom_cams:
            path, n_hits, _ = build(gene, camera, hits, pocket, True, zoom_cam=zoom_cams[gene])
            n_in = sum(1 for h in hits if h["gene"] == gene and h.get("zoom_in_frame") == "True")
            print(f"wrote {os.path.relpath(path, HERE)}  ({n_in} of {n_hits} hits in the zoom frame)")
        off = [h for h in hits if h["gene"] == gene and h["modelled"] == "True"
               and h["in_crop"] != "True"]
        un = [h for h in hits if h["gene"] == gene and h["modelled"] != "True"]
        if off:
            print(f"   off-crop (drawn, outside the zoom): "
                  + ", ".join(f'{h["wt"]}{h["position"]} {float(h["min_dist_to_ATP"]):.0f}A' for h in off))
        if un:
            print(f"   not modelled: " + ", ".join(f'{h["wt"]}{h["position"]}' for h in un))
    print(f"zoom cameras: {', '.join(sorted(zoom_cams)) or 'none -- run compute_active_site_camera.py'}")
    print("wrote", os.path.relpath(write_render_script(zoom_cams), HERE))
