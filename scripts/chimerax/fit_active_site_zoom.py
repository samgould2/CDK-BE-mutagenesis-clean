"""
Measure how much of the frame each active-site panel's complex actually fills.

    ChimeraX --nogui --exit --cmd "runscript fit_active_site_zoom.py"

Run this after changing the camera, FRAME_ON or any zoom in
generate_cdk2_4_6_active_site.py, then update PER_GENE_ZOOM there from the
numbers it prints:

    zoom_new = zoom_current * TARGET / span

A single shared zoom does not work: with the same camera and framing, CDK2's
complex spans 0.90 of the half-frame while CDK4's and CDK6's span 0.53-0.55,
because 1QMZ's cyclin A is much larger in projection than cyclin D1 or V-cyclin.
One value therefore either overflows CDK2 or leaves CDK4/CDK6 swimming in white.

NOTE on the ligand-occlusion number below: select the ligand by CHAIN, not by
residue name alone.  1QMZ contains two copies of ATP, and averaging both puts the
test centroid between them, which reported 44 blocking atoms for a ligand that is
actually in the open (the panel draws only chain A's copy).
"""
import math, os, sys
import numpy as np
from chimerax.core.commands import run
# This script is run by ChimeraX (`runscript`), whose working directory is
# wherever the user launched it, so resolve everything beside the script itself.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cdk_config import STRUCTURES
W = 1400
TARGET = 0.88     # fraction of the half-frame the complex should span

for gene in ["CDK2", "CDK4", "CDK6"]:
    cfg = STRUCTURES[gene]; kc = cfg["kinase_chain"]
    gid = 1 if gene == "CDK2" else 2
    run(session, "close session")
    run(session, f"windowsize {W} {W}")
    run(session, f"open {HERE}/cxc/{gene}_active_site.cxc")
    m = [x for x in session.models if getattr(x,'id',None)==(gid,)][0]
    cam = session.main_view.camera
    inv = cam.position.inverse()
    th = math.tan(math.radians(cam.field_of_view)/2)
    w,h = session.main_view.window_size or (W,W)
    tv = th*(h/w)

    chains = [kc] + cfg["cyclin_chains"]
    pts = np.array([a.scene_coord for a in m.atoms if a.residue.chain_id in chains])
    cc = np.array([inv*tuple(p) for p in pts]); depth = -cc[:,2]
    fx = np.abs(cc[:,0])/(th*depth); fy = np.abs(cc[:,1])/(tv*depth)
    span = max(fx.max(), fy.max())
    print(f"{gene}: complex spans {span:.3f} of the half-frame -> "
          f"zoom x{span/TARGET:.2f} would make it {TARGET:.2f}")

    # is the ligand actually unoccluded from this camera?
    lig_spec = (lambda: [a for a in m.atoms if a.residue.name == cfg["ligand"]]) if cfg["ligand"] else None
    if lig_spec:
        lig = np.array([a.scene_coord for a in lig_spec()])
        prot = np.array([a.scene_coord for a in m.atoms
                         if a.residue.chain_id == kc and a.residue.polymer_type != 0
                         and a.element.name != "H"])
        zc = cam.position.z_axis()
        zc = np.array(zc)/np.linalg.norm(zc)
        c_l = lig.mean(axis=0)
        v = prot - c_l; t = v @ zc
        perp = np.linalg.norm(v - t[:,None]*zc[None,:], axis=1)
        blockers = int(((t > 1.0) & (t < 25.0) & (perp < 4.0)).sum())
        print(f"   ligand {cfg['ligand']}: {len(lig)} atoms, {blockers} protein atoms in front of it")
