"""
The catalytic-domain ranges for CDK7/8/9/12/13/19 come from align_catalytic_v2.py
(KLIFS catalytic site).  CDK2, CDK4 and CDK6 were not in that script, so derive
their equivalent ranges from the superposition: for each known range endpoint,
find the closest CA in the aligned CDK2 reference, and vice versa.
"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cdk_config import STRUCTURES

def cas(m, ch):
    a = m.atoms
    sel = [x for x in a if x.name == "CA" and x.residue.chain_id == ch]
    return sel, np.array([x.scene_coord for x in sel])

ref_m = [m for m in session.models if getattr(m, "name", "") == "CDK2"][0]
ref_atoms, ref_xyz = cas(ref_m, STRUCTURES["CDK2"]["kinase_chain"])

def nearest(xyz, atoms, point):
    d = np.linalg.norm(xyz - point, axis=1)
    i = int(np.argmin(d))
    return atoms[i].residue.number, float(d[i])

print("known range -> equivalent CDK2 residues")
starts, ends = [], []
for gene in ("CDK7", "CDK8", "CDK9", "CDK12", "CDK13", "CDK19"):
    cfg = STRUCTURES[gene]
    m = [x for x in session.models if getattr(x, "name", "") == gene][0]
    atoms, xyz = cas(m, cfg["kinase_chain"])
    bynum = {a.residue.number: a for a in atoms}
    for label, num, bucket in (("start", cfg["kd"][0], starts), ("end", cfg["kd"][1], ends)):
        a = bynum.get(num)
        if a is None:
            print(f"  {gene} {label} {num}: not modelled"); continue
        r, d = nearest(ref_xyz, ref_atoms, a.scene_coord)
        bucket.append(r)
        print(f"  {gene:6s} {label:5s} {num:5d} -> CDK2 {r:4d}  ({d:.1f} A)")

lo, hi = int(np.median(starts)), int(np.median(ends))
print(f"\nconsensus CDK2 catalytic range: {lo}-{hi}")

ref_bynum = {a.residue.number: a for a in ref_atoms}
for gene in ("CDK4", "CDK6"):
    cfg = STRUCTURES[gene]
    m = [x for x in session.models if getattr(x, "name", "") == gene][0]
    atoms, xyz = cas(m, cfg["kinase_chain"])
    out = []
    for num in (lo, hi):
        a = ref_bynum.get(num)
        r, d = nearest(xyz, atoms, a.scene_coord)
        out.append(r)
        print(f"  CDK2 {num} -> {gene} {r} ({d:.1f} A)")
    print(f"  ==> {gene} kd = {tuple(out)}")
