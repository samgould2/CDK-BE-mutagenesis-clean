"""
Spread the active-site residue labels so they stop overlapping.

    ChimeraX --nogui --exit --cmd "runscript deconflict_labels.py"
    python3 generate_cdk2_4_6_active_site.py     # picks up the new offsets

Rewrites the label_dx/dy/dz columns of cdk2_4_6_hit_distances.csv, which
generate_cdk2_4_6_active_site.py emits as each label's `offset`.  Run it after the
panels exist: it opens each one to read the camera that panel actually ends up with.
The camera does not depend on the offsets, so one extra generate pass settles it.

Two earlier attempts and why they failed
----------------------------------------
1. **Radial offsets in 3D.**  Push each label away from the pocket centroid.  Overlap
   happens in the *projection*, and a 3D radial direction can project to nearly
   nothing, so two residues stacked along the view axis were pushed apart in a
   direction the camera cannot see.
2. **Force-directed relaxation in the camera plane, in Angstroms.**  Better, but
   measuring label boxes in scene units assumes an orthographic camera.  ChimeraX
   renders with perspective (fov 30 deg here) over a complex ~80 A deep, so a label
   near the camera draws much larger than one at the back and collides with
   neighbours the model thought were clear.  The geometric check reported 1
   overlapping pair for CDK6 while the render plainly showed six.

What this does
--------------
Works in **normalized screen coordinates**, taken from the panel's real camera, so
perspective is included by construction:

  * a scene-unit size s at depth d spans s / (d * tan(fov/2)) of the half-frame,
  * each label gets a box in those units from LABEL_HEIGHT and its text length,
  * labels are placed greedily, strongest hit first, trying the zero offset and then
    rings of candidates, keeping the first that hits no placed label and covers no
    hit sphere,
  * the chosen screen displacement is converted back to a scene offset along the
    camera's x/y axes (labels are billboarded, so that moves them in the screen plane).
"""
import csv, os, math, re, sys
import numpy as np
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
# Optional variant tag.  Set CDK246_TAG to run the whole chain a second time against a
# different hit set without overwriting the first: every file this script reads or writes
# picks up the suffix.  Used for the stringent-FDR variant (CDK246_TAG=_strictFDR).
TAG = os.environ.get("CDK246_TAG", "")

sys.path.insert(0, HERE)
from cdk_config import STRUCTURES
from generate_cdk2_4_6_active_site import LABEL_HEIGHT, LABEL_HEIGHT_ZOOM

# CDK246_ZOOM=1 solves the labels against the ligand-centred zoom camera instead,
# at the zoom label height, and writes the zoom_* columns.
ZOOM = bool(os.environ.get('CDK246_ZOOM'))

GENES = ["CDK2", "CDK4", "CDK6"]
WINDOW = (1400, 1400)
WIDTH_PER_CHAR = 0.62     # label width in units of LABEL_HEIGHT, per character
PAD = 0.30                # extra clearance between boxes, in units of LABEL_HEIGHT
N_DIRS = 24
RADII = [0.0, 2.0, 3.5, 5.0, 6.5, 8.0, 9.5, 11.0, 13.0]   # A, converted to screen per label
# Bigger labels need more room to escape.  Past ~13 A a label is far enough from its
# sphere to be ambiguous, and ChimeraX residue labels have no leader lines, so that
# is the practical ceiling -- a few residual overlaps are preferred to detached labels.
SPHERE_R = 2.0            # A; approximate hit-sphere radius, so labels avoid the spheres

dist = list(csv.DictReader(open(os.path.join(HERE, f"cdk2_4_6_hit_distances{TAG}.csv"))))
LH = LABEL_HEIGHT_ZOOM if ZOOM else LABEL_HEIGHT
PREFIX = "zoom_" if ZOOM else "label_"
for _r in dist:
    for _c in (PREFIX + "dx", PREFIX + "dy", PREFIX + "dz"):
        _r.setdefault(_c, "")
    if ZOOM:
        _r.setdefault("zoom_in_frame", "")
zoom_cams = {}
if ZOOM:
    _t = open(os.path.join(HERE, "active_site_camera.txt")).read()
    for _g in GENES:
        _m = re.search(rf'ZOOM_CAMERA_{_g} = "([^"]+)"', _t)
        if _m:
            zoom_cams[_g] = _m.group(1)
imp = {(r["gene"], r["position"]): float(r["active_median_LFC"])
       for r in csv.DictReader(open(os.path.join(HERE, f"cdk2_4_6_hit_residues{TAG}.csv")))}

for gene in GENES:
    cfg = STRUCTURES[gene]
    kc = cfg["kinase_chain"]
    gid = 1 if gene == "CDK2" else 2
    run(session, "close session")
    run(session, f"windowsize {WINDOW[0]} {WINDOW[1]}")
    run(session, f"open {os.path.join(HERE, 'cxc', gene + '_active_site' + TAG + '.cxc')}")
    if ZOOM:
        if gene not in zoom_cams:
            continue
        run(session, f"view matrix camera {zoom_cams[gene]}")
        run(session, "clip off")

    m = [x for x in session.models if getattr(x, "id", None) == (gid,)][0]
    cam = session.main_view.camera
    inv = cam.position.inverse()
    cx = np.array(cam.position.axes()[0]); cx /= np.linalg.norm(cx)
    cy = np.array(cam.position.axes()[1]); cy /= np.linalg.norm(cy)
    w, h = session.main_view.window_size or WINDOW
    th = math.tan(math.radians(cam.field_of_view) / 2)
    tv = th * (h / w)

    byres = {}
    for r in m.residues:
        if r.chain_id == kc:
            byres.setdefault(r.number, r)

    rows = [d for d in dist if d["gene"] == gene and d["modelled"] == "True"
            and int(d["position"]) in byres]
    if not rows:
        continue

    anchors, boxes, depths = [], [], []
    for d in rows:
        c = np.array([a.scene_coord for a in byres[int(d["position"])].atoms]).mean(axis=0)
        cc = np.array(inv * tuple(c))
        depth = max(-cc[2], 1e-6)
        # normalized half-frame coords
        anchors.append([cc[0] / (th * depth), cc[1] / (tv * depth)])
        depths.append(depth)
        t = f'{d["wt"]}{d["position"]}'
        boxes.append([LH * WIDTH_PER_CHAR * len(t) / (th * depth),
                      LH / (tv * depth)])
    anchors = np.array(anchors); boxes = np.array(boxes); depths = np.array(depths)
    spheres = anchors.copy()
    pad = np.array([PAD * LH / (th * depths), PAD * LH / (tv * depths)]).T
    sph = np.array([SPHERE_R / (th * depths), SPHERE_R / (tv * depths)]).T

    def overlap(pi, bi, pj, bj, pd):
        need = (bi + bj) / 2 + pd
        ov = need - np.abs(pi - pj)
        return float(min(ov)) if (ov > 0).all() else 0.0

    order = sorted(range(len(rows)), key=lambda i: -imp[(gene, rows[i]["position"])])
    pos = anchors.copy()
    placed = []
    for i in order:
        best, best_cost = None, None
        for r_ in RADII:
            for k in range(N_DIRS if r_ > 0 else 1):
                a = 2 * math.pi * k / N_DIRS
                # candidate radius is in A at this label's depth -> screen units
                off = np.array([r_ * math.cos(a) / (th * depths[i]),
                                r_ * math.sin(a) / (tv * depths[i])])
                p = anchors[i] + off
                cost = sum(overlap(p, boxes[i], pos[j], boxes[j], (pad[i] + pad[j]) / 2) * 50
                           for j in placed)
                half = boxes[i] / 2 + sph
                dd = np.abs(spheres - p)
                cost += float(((dd[:, 0] < half[:, 0]) & (dd[:, 1] < half[:, 1])).sum())
                cost += r_ * 0.02
                if best_cost is None or cost < best_cost:
                    best, best_cost = p, cost
            if best_cost is not None and best_cost <= r_ * 0.02 + 1e-9:
                break
        pos[i] = best
        placed.append(i)

    def n_over(p):
        return sum(1 for i in range(len(p)) for j in range(i + 1, len(p))
                   if overlap(p[i], boxes[i], p[j], boxes[j], (pad[i] + pad[j]) / 2) > 0)

    before, after = n_over(anchors), n_over(pos)
    # screen displacement -> scene offset along the camera axes
    shifts = []
    if ZOOM:
        for i, d in enumerate(rows):
            d["zoom_in_frame"] = str(bool(abs(anchors[i][0]) <= 1 and abs(anchors[i][1]) <= 1))
    for i, d in enumerate(rows):
        dx = (pos[i][0] - anchors[i][0]) * th * depths[i]
        dy = (pos[i][1] - anchors[i][1]) * tv * depths[i]
        v = dx * cx + dy * cy
        shifts.append(math.hypot(dx, dy))
        d[PREFIX + "dx"], d[PREFIX + "dy"], d[PREFIX + "dz"] = (round(float(v[0]), 2),
                                                                round(float(v[1]), 2),
                                                                round(float(v[2]), 2))
    tag = " [zoom]" if ZOOM else ""
    n_if = sum(1 for d in rows if d.get("zoom_in_frame") == "True") if ZOOM else len(rows)
    print(f"{gene}{tag}: {len(rows)} labels ({n_if} in frame) | overlapping pairs "
          f"{before} -> {after} | shift mean {np.mean(shifts):.1f} A "
          f"max {np.max(shifts):.1f} A (depth {depths.min():.0f}-{depths.max():.0f} A)")

with open(os.path.join(HERE, f"cdk2_4_6_hit_distances{TAG}.csv"), "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(dist[0].keys()))
    wr.writeheader()
    wr.writerows(dist)
print("\nrewrote label offsets in cdk2_4_6_hit_distances.csv")
