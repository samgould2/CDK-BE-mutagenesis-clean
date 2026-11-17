"""
Check that no hit residue is cut by the near/far clipping planes.

    ChimeraX --nogui --exit --cmd "runscript validate_ddb1_clipping.py"

`view <spec>` sets the clipping slab to that spec's own bounds, so a zoom onto
the small HQ461 ligand will slice any residue nearer the camera than the ligand
(CDK13 T715 was being cut in half this way).  `clip off` resets the slab, but it
has to be issued *after* the last view/zoom.  This checks both panel cameras for
both genes and reports any hit atom outside the slab.
"""
import csv, os
import numpy as np
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
KINASE_MODEL = {"CDK12": ("1", "B"), "CDK13": ("2", "A")}

rows = list(csv.DictReader(open(os.path.join(HERE, "cdk12_13_ddb1_hit_residues.csv"))))
report = []
bad = 0
for gene in ("CDK12", "CDK13"):
    mid, chain = KINASE_MODEL[gene]
    positions = sorted({int(r["position"]) for r in rows if r["gene"] == gene})
    for view_name, cmds in (("wide panel", []),
                            ("glue zoom", ["view HQ461", "zoom 0.30", "clip off"])):
        run(session, "close session")
        run(session, f"open {HERE}/cxc/{gene}_DDB1.cxc")
        run(session, "windowsize 1200 1200")
        for c in cmds:
            run(session, c)
        v = session.main_view
        near, far = v.near_far_distances(v.camera, 0)
        P = np.array(v.camera.position.matrix)
        origin, zax = P[:, 3], P[:, 2]
        clipped = []
        for m in session.models:
            if getattr(m, "residues", None) is None or m.id_string != mid:
                continue
            for r in m.residues:
                if r.chain_id != chain or r.number not in positions:
                    continue
                d = [np.dot(np.array(a.scene_coord) - origin, -zax) for a in r.atoms]
                if min(d) < near or max(d) > far:
                    clipped.append(f"{r.name}{r.number} ({min(d):.1f}-{max(d):.1f})")
        bad += len(clipped)
        report.append(f"{gene:6s} {view_name:10s} slab {near:7.2f}-{far:7.2f}  "
                      f"{len(positions)} hit positions  clipped: {clipped if clipped else 'none'}")

text = "\n".join(report) + f"\n\nTOTAL CLIPPED HIT RESIDUES: {bad}\n"
print(text)
with open(os.path.join(HERE, "ddb1_clipping_report.txt"), "w") as fh:
    fh.write(text)
