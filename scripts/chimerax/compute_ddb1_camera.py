"""
Derive the shared camera orientation for the CDK12/CDK13 DDB1 panels.

    ChimeraX --nogui --exit --cmd "runscript compute_ddb1_camera.py"

The three subunits sit roughly in a line, DDB1 - CDK - cyclin K.  The figure
wants that line across the frame with the viewer side-on to it:

    camera x (right)         DDB1 centroid -> cyclin K centroid, so DDB1 ends up
                             on the left and cyclin K on the right
    camera z (toward viewer) normal of the plane through the three centroids, so
                             the subunits sit side by side instead of behind
                             one another
    camera y (up)            z x x

Looking *along* the DDB1->CDK axis instead puts DDB1 in front of the kinase and
hides it.

Only the ROTATION is emitted; the translation is left at zero.  The panels then
run `view <spec>` straight after `view matrix camera`, which moves the camera
along its own axis to frame the target without touching the orientation.  That
is much more robust than computing a distance here, which has to agree with
ChimeraX's field-of-view and window aspect ratio to land correctly.

Both mirror images are printed.  SIGN = +1 is the one matching the reference
layout (DDB1 left, cyclin K right).
"""
import numpy as np
from chimerax.core.commands import run
import os

HERE = os.path.dirname(os.path.abspath(__file__))
run(session, "close session")
m = run(session, f"open {HERE}/structures/8bug.cif")[0]

def centroid(chain_id):
    r = [c for c in m.chains if c.chain_id == chain_id][0].existing_residues
    return np.array([a.scene_coord for a in r.atoms]).mean(axis=0)

c_ddb1, c_cdk, c_cyc = centroid("A"), centroid("B"), centroid("C")

x = c_cyc - c_ddb1
x /= np.linalg.norm(x)
n = np.cross(c_cdk - c_ddb1, c_cyc - c_ddb1)
n /= np.linalg.norm(n)

lines = []
for sign in (+1, -1):
    z = n * sign
    y = np.cross(z, x)
    M = np.column_stack([x, y, z, np.zeros(3)])
    assert np.allclose(M[:, :3].T @ M[:, :3], np.eye(3))
    tag = "+1 (DDB1 left -- use this one)" if sign > 0 else "-1 (mirror image)"
    lines.append(f'SIGN {tag}\n  CAMERA_ROTATION = "'
                 + ",".join(f"{v:.4f}" for v in M.flatten()) + '"')

text = "\n".join(lines)
print(text)
with open(os.path.join(HERE, "ddb1_camera.txt"), "w") as fh:
    fh.write(text + "\n")
