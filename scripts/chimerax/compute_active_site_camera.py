"""
Derive the shared camera orientation for the CDK2/4/6 active-site panels.

    ChimeraX --nogui --exit --cmd "runscript compute_active_site_camera.py"

The camera is fixed by two layout constraints, which together leave no freedom:

    camera x (right)          kinase centroid -> cyclin centroid, so the cyclin
                              sits on the right
    camera y (up)             ATP centroid -> hinge centroid, so the hinge is at
                              the top of the pocket rather than directly behind
                              the ligand
    camera z (toward viewer)  x cross y -- forced, once x and y are chosen

The hinge is CDK2 80-86 (gatekeeper F80 through D86); 80, 81, 82, 83 and 86 are
all in the independently computed ATP pocket, so the range is the right one.

With both constraints applied the hinge lands +7.8 A above the ATP in screen y at
a depth offset of 0.0 -- beside the ligand, not behind it -- and the cyclin +34 A
to the right.  The opposite handedness is rejected: it would put the cyclin on the
left and is more occluded (90 atoms vs 50).

An earlier version instead picked the **least-occluded** axis by sampling 800
directions (4 blocking atoms, against 57 for the naive kinase-centroid -> ATP
axis).  That gives the most open view into the pocket but puts the hinge directly
behind the ligand.  The occlusion of whatever axis is in use is still reported
below, because it is worth knowing -- though note the count assumes space-filling
protein, and these panels draw the kinase as a cartoon, so a ribbon in the ray
blocks far less than the number suggests.
    camera x (right)          kinase centroid -> cyclin centroid, orthogonalized
                              against z, so the kinase sits left and the cyclin
                              right and the pair reads side-on ("straight-on")
    camera y (up)             z x x

An earlier version took y from the N-lobe -> C-lobe axis instead.  That keeps the
same open view of the cleft but rolls the scene so the cyclin lands on the left
and the kinase is seen obliquely; aligning x with the kinase -> cyclin axis is
what makes the pair read straight-on.

Computed once from 1QMZ (the repo's superposition reference, which never moves),
so all three panels share one orientation and stay directly comparable -- the
same property cdk_config.CANONICAL_VIEW gives the whole-kinase panels.

The FULL camera matrix is emitted -- rotation *and* origin -- so every panel gets a
byte-identical camera and the three are genuinely equivalent in orientation,
perspective and scale.

This replaces an earlier approach of emitting only the rotation and letting each
panel run `view <spec>` to frame itself.  That does not give equal cameras:
re-issued identically, with the same number of superposed catalytic-domain atoms
displayed, `view <gene>_kd` put the camera 136 A from the domain for CDK2 but
259 A for CDK4 and 256 A for CDK6.  `view` is not simply framing the spec's
bounds -- its baseline shifts once a second model (the hidden 1QMZ reference) is
in the session, and the CDK2 panel has only one model because it *is* the
reference.  Computing the distance here removes that dependency entirely.

The distance is set so a single sphere, centred on the superposed catalytic domain
and large enough to contain the kinase+cyclin of ALL THREE genes, fits the frustum.
CDK2's cyclin A is the largest, so it sets the radius and the other two panels
carry more empty space -- that is the unavoidable cost of one shared camera, and
the alternative (per-gene framing) is what made the panels non-comparable.

Writes active_site_camera.txt.
"""
import os
import sys
import numpy as np
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cdk_config import STRUCTURES, ALIGN_REFERENCE, CANONICAL_VIEW
CHAIN = "A"            # 1QMZ kinase chain
CYCLIN_CHAIN = "B"     # 1QMZ cyclin A chain
HINGE = (80, 86)       # CDK2 gatekeeper F80 .. D86
N_LOBE = (8, 80)       # within the CDK2 catalytic range 8-149; hinge sits ~80-83
C_LOBE = (84, 149)
RAY_RADIUS = 5.0       # A; cylinder around the view ray
RAY_REACH = 25.0       # A; how far out from the ATP to look for blockers
N_SAMPLE = 800
FOV = 30.0             # ChimeraX default vertical field of view, degrees
FIT_MARGIN = 1.06      # a little air around the fitted sphere
ZOOM_HALF_WIDTH = 20.0 # A; half the field of the ligand-centred zoom cameras.  At 20 A
                       # the frame spans 40 A -- the pocket plus its surroundings, and
                       # 12/20, 25/29 and 26/34 of the drawn hits for CDK2/4/6.

run(session, "close session")
m = run(session, f"open {HERE}/structures/1qmz.cif")[0]

sel = lambda f: np.array([a.scene_coord for a in m.atoms if f(a)])
atp = sel(lambda a: a.residue.chain_id == CHAIN and a.residue.name == "ATP")
prot = sel(lambda a: a.residue.chain_id == CHAIN and a.residue.polymer_type != 0
           and a.element.name != "H")
nlo = sel(lambda a: a.residue.chain_id == CHAIN and N_LOBE[0] <= a.residue.number <= N_LOBE[1])
clo = sel(lambda a: a.residue.chain_id == CHAIN and C_LOBE[0] <= a.residue.number <= C_LOBE[1])
cyc = sel(lambda a: a.residue.chain_id == CYCLIN_CHAIN and a.residue.polymer_type != 0)
hin = sel(lambda a: a.residue.chain_id == CHAIN and HINGE[0] <= a.residue.number <= HINGE[1])

c_atp, c_kin = atp.mean(axis=0), prot.mean(axis=0)
lobe = clo.mean(axis=0) - nlo.mean(axis=0)
# the kinase -> cyclin axis is what sets the roll: it goes across the frame, so the
# two subunits sit side by side rather than one behind the other
axis_kc = cyc.mean(axis=0) - c_kin


def occlusion(z):
    """protein heavy atoms within RAY_RADIUS of the ray ATP -> camera"""
    v = prot - c_atp
    t = v @ z
    perp = np.linalg.norm(v - t[:, None] * z[None, :], axis=1)
    return int(((t > 1.0) & (t < RAY_REACH) & (perp < RAY_RADIUS)).sum())


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


def matrix(_ignored=None):
    """x from the kinase->cyclin axis, y from the ATP->hinge axis, z forced."""
    x = axis_kc / np.linalg.norm(axis_kc)
    y_raw = hin.mean(axis=0) - c_atp
    y = y_raw - (y_raw @ x) * x
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    M = np.column_stack([x, y, z, np.zeros(3)])
    R = M[:, :3]
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-8), "not orthonormal"
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-8), "not right-handed"
    return M


z_naive = c_atp - c_kin
z_naive /= np.linalg.norm(z_naive)
cands = fibonacci_sphere(N_SAMPLE)
scores = np.array([occlusion(z) for z in cands])

# the constrained camera; z is whatever x cross y gives
M_constrained = matrix()
z_best = M_constrained[:, 2].copy()

# ---- one framing sphere that works for all three genes.  Centre it on the
# superposed catalytic domain and size it to the largest kinase+cyclin.
GENES = ["CDK2", "CDK4", "CDK6"]
kd_lo, kd_hi = STRUCTURES["CDK2"]["kd"]
centre = sel(lambda a: a.residue.chain_id == CHAIN
             and kd_lo <= a.residue.number <= kd_hi).mean(axis=0)

radius, per_gene, models_for_zoom = 0.0, {}, {}
for g in GENES:
    cfg = STRUCTURES[g]
    if g == ALIGN_REFERENCE:
        mg = m
    else:
        mg = run(session, f"open {HERE}/structures/{cfg['file']}")[0]
        run(session, f"matchmaker #{mg.id[0]}/{cfg['kinase_chain']} to #{m.id[0]}/{CHAIN}")
    models_for_zoom[g] = mg
    chains = [cfg["kinase_chain"]] + cfg["cyclin_chains"]
    pts = np.array([a.scene_coord for a in mg.atoms
                    if a.residue.chain_id in chains and a.residue.polymer_type != 0])
    r = float(np.linalg.norm(pts - centre, axis=1).max())
    per_gene[g] = r
    radius = max(radius, r)

half = np.radians(FOV) / 2
distance = radius / np.sin(half) * FIT_MARGIN

def full_matrix(z):
    M = matrix(z)
    origin = centre + z * distance        # camera looks along -z, so sit on the +z side
    M[:, 3] = origin
    return M

# ---- ligand-centred zoom cameras, one per gene.  Same ROTATION as the canonical
# camera, so the zoom matches the wide panels and the rest of the figure set; only the
# origin changes, to sit in front of that gene's own ligand.  CDK4 is apo, so it uses
# the superposed 1QMZ ATP -- the same thing its wide panel draws to mark the pocket.
cf = [float(v) for v in CANONICAL_VIEW.split(",")]
can_x = np.array([cf[0], cf[4], cf[8]])
can_y = np.array([cf[1], cf[5], cf[9]])
can_z = np.array([cf[2], cf[6], cf[10]])
can_z = can_z / np.linalg.norm(can_z)
zoom_d = ZOOM_HALF_WIDTH / np.tan(np.radians(FOV) / 2)

zoom_lines = []
for g in GENES:
    cfg = STRUCTURES[g]
    mg = models_for_zoom[g]
    if cfg["ligand"]:
        lig = np.array([a.scene_coord for a in mg.atoms
                        if a.residue.chain_id == cfg["kinase_chain"]
                        and a.residue.name == cfg["ligand"]])
        src = f"{cfg['ligand']} (native)"
    else:
        lig = np.array([a.scene_coord for a in m.atoms
                        if a.residue.chain_id == CHAIN and a.residue.name == "ATP"])
        src = "ATP (borrowed from 1QMZ, this structure is apo)"
    c_lig = lig.mean(axis=0)
    Mz = np.column_stack([can_x, can_y, can_z, c_lig + can_z * zoom_d])
    zoom_lines.append(f"  {g}: {src}, centroid {np.round(c_lig, 2).tolist()}")
    zoom_lines.append(f'  ZOOM_CAMERA_{g} = "'
                      + ",".join(f"{v:.4f}" for v in Mz.flatten()) + '"')

lines = [
    f"ATP centroid {np.round(c_atp, 2).tolist()}   kinase centroid {np.round(c_kin, 2).tolist()}",
    f"catalytic-domain centre (shared frame): {np.round(centre, 2).tolist()}",
    "radius from that centre to the furthest kinase+cyclin atom:",
    *[f"    {g}: {r:6.1f} A" for g, r in per_gene.items()],
    f"shared framing radius {radius:.1f} A -> camera distance {distance:.1f} A "
    f"(fov {FOV} deg, margin {FIT_MARGIN})",
    f"occlusion = protein heavy atoms within {RAY_RADIUS} A of the ray, 1-{RAY_REACH} A out",
    "",
    f"hinge {HINGE[0]}-{HINGE[1]} centroid {np.round(hin.mean(axis=0), 2).tolist()}",
    "",
    f"CHOSEN axis (x = kinase->cyclin, y = ATP->hinge, z = x cross y): "
    f"{np.round(z_best, 4).tolist()}, occlusion {occlusion(z_best)} atoms",
    f"  hinge is {float((hin.mean(axis=0) - c_atp) @ M_constrained[:, 1]):+.1f} A above the ATP in screen y, "
    f"depth offset {float((hin.mean(axis=0) - c_atp) @ M_constrained[:, 2]):+.1f} A (0 = beside, not behind)",
    f"  cyclin is {float((cyc.mean(axis=0) - c_atp) @ M_constrained[:, 0]):+.1f} A to the right",
    "",
    "for context only -- these are NOT used, the layout constraints decide the axis:",
    f"  naive axis (kinase centroid -> ATP): occlusion {occlusion(z_naive)} atoms",
    f"  {N_SAMPLE} sampled directions: occlusion min {scores.min()}  median {int(np.median(scores))}  max {scores.max()}",
    "  the least-occluded axis gives the most open pocket but puts the hinge behind the ligand",
    "",
    'ACTIVE_SITE_ROTATION = "' + ",".join(f"{v:.4f}" for v in matrix(z_best).flatten()) + '"',
    'ACTIVE_SITE_CAMERA = "' + ",".join(f"{v:.4f}" for v in full_matrix(z_best).flatten()) + '"',
    "",
    "",
    f"ligand-centred zoom cameras (half-width {ZOOM_HALF_WIDTH} A, distance {zoom_d:.1f} A,",
    " canonical rotation so the zoom matches the wide panels):",
    *zoom_lines,
    "",
    "for reference, the naive axis would have been:",
    '  "' + ",".join(f"{v:.4f}" for v in matrix(z_naive).flatten()) + '"',
]
text = "\n".join(lines)
print(text)
with open(os.path.join(HERE, "active_site_camera.txt"), "w") as fh:
    fh.write(text + "\n")
