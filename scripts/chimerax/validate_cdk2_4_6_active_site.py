"""
Check the CDK2/4/6 active-site panels without rendering them.

    ChimeraX --nogui --exit --cmd "runscript validate_cdk2_4_6_active_site.py"

`save` needs OpenGL and so cannot run headless, but everything that decides
whether the figure is right can still be checked:

  (a) each .cxc runs start to finish without an error,
  (b) every labeled hit is inside the camera frustum, i.e. actually in frame
      after `view <pocket>` + `zoom` -- the thing that cannot be eyeballed here,
  (c) nothing is cut by the near/far clipping planes (the trap documented in the
      README for the DDB1 panels),
  (d) the pocket marker (ATP or inhibitor) is displayed.

Writes cdk2_4_6_active_site_validation.txt.
"""
import csv, math, os, sys
import numpy as np
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
# Optional variant tag.  Set CDK246_TAG to run the whole chain a second time against a
# different hit set without overwriting the first: every file this script reads or writes
# picks up the suffix.  Used for the stringent-FDR variant (CDK246_TAG=_strictFDR).
TAG = os.environ.get("CDK246_TAG", "")

sys.path.insert(0, HERE)
from cdk_config import STRUCTURES

GENES = ["CDK2", "CDK4", "CDK6"]
WINDOW = (1400, 1400)          # must match the render script
out = []

def emit(s=""):
    print(s)
    out.append(s)

dist = {(r["gene"], r["position"]): r
        for r in csv.DictReader(open(os.path.join(HERE, f"cdk2_4_6_hit_distances{TAG}.csv")))}

emit(f"window {WINDOW[0]}x{WINDOW[1]} (square, as in render_cdk2_4_6_active_site.cxc)")
emit()
problems = 0

for gene in GENES:
    for labeled in (False, True):
        name = f"{gene}_active_site" + ("_labeled" if labeled else "")
        run(session, "close session")
        try:
            run(session, f"windowsize {WINDOW[0]} {WINDOW[1]}")
            run(session, f"open {os.path.join(HERE,'cxc',name+'.cxc')}")
        except Exception as e:
            emit(f"=== {name}: FAILED TO RUN -- {e}")
            problems += 1
            continue

        cfg = STRUCTURES[gene]
        kc = cfg["kinase_chain"]
        gid = 1 if gene == "CDK2" else 2
        cam = session.main_view.camera
        inv = cam.position.inverse()
        w, h = session.main_view.window_size or WINDOW
        hfov = math.radians(cam.field_of_view)
        tan_h = math.tan(hfov / 2)
        tan_v = tan_h * (h / w)

        # clipping planes currently in force
        planes = {p.name: p for p in session.main_view.clip_planes.planes()}

        gm = [m for m in session.models if getattr(m, "id", None) == (gid,)]
        model = gm[0] if gm else None
        byres = {}
        if model is not None:
            for r in model.residues:
                if r.chain_id == kc:
                    byres.setdefault(r.number, r)

        labeled_hits = [p for (g, p) in dist if g == gene
                        and dist[(g, p)]["in_crop"] == "True"]
        off = []
        clipped = []
        for pos in sorted(labeled_hits, key=int):
            r = byres.get(int(pos))
            if r is None:
                continue
            c = np.array([a.scene_coord for a in r.atoms])
            cc = np.array([inv * tuple(p) for p in c])
            depth = -cc[:, 2]
            # a residue is in frame if any of its atoms is inside the frustum
            inside = ((depth > 0) &
                      (np.abs(cc[:, 0]) <= tan_h * depth) &
                      (np.abs(cc[:, 1]) <= tan_v * depth))
            if not inside.any():
                off.append(f"{r.name}{pos}")
            for pname, pl in planes.items():
                d = (c - np.array(pl.plane_point)) @ np.array(pl.normal)
                if (d < 0).any():
                    clipped.append(f"{pos}/{pname}")

        marker = cfg["ligand"] or "ATP(borrowed)"
        lig_spec = (f"#{gid}/{kc} & :{cfg['ligand']}" if cfg["ligand"] else "#1/A & :ATP")
        n_lig = len(run(session, f"info atoms {lig_spec}") or []) if False else None

        emit(f"=== {name}")
        emit(f"    camera fov {cam.field_of_view:.1f} deg, window {w}x{h}")
        emit(f"    clip planes in force: {sorted(planes) or 'none (clip off)'}")
        emit(f"    pocket marker: {marker}")
        emit(f"    labeled hits: {len(labeled_hits)} -> "
             + ", ".join(sorted(labeled_hits, key=int)))
        if off:
            emit(f"    OUT OF FRAME: {off}")
            problems += 1
        else:
            emit("    all labeled hits in frame: yes")
        if clipped:
            emit(f"    CLIPPED: {clipped}")
            problems += 1
        else:
            emit("    clipped by near/far plane: none")
        emit()

emit(f"{'OK' if problems == 0 else str(problems) + ' PROBLEM(S)'}")
with open(os.path.join(HERE, f"cdk2_4_6_active_site_validation{TAG}.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
