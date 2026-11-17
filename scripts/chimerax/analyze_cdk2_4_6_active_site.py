"""
Run inside ChimeraX:
  ChimeraX --nogui --exit --cmd "runscript analyze_cdk2_4_6_active_site.py"

Defines the ATP pocket from the structures rather than by eye, and measures
every CDK2/4/6 hit against it.

CDK2 (1QMZ) has ATP bound and is the superposition reference for the whole repo.
CDK6 (2F2C) has aminopurvalanol in the same cleft.  **CDK4 (2W96) is apo**, so
it has no ligand to frame an active site on.  Superposing all three onto 1QMZ
chain A puts 1QMZ's ATP in a common frame, which gives:

  * one framing object that exists for all three genes (the reference ATP), so
    the three zoom panels come out in directly comparable crops, and
  * a distance-to-ATP for every hit residue, which is what decides whether a hit
    is an active-site hit at all.

Writes cdk2_4_6_active_site_report.txt and cdk2_4_6_hit_distances.csv.
"""
import csv, os, sys
import numpy as np
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
# Optional variant tag.  Set CDK246_TAG to run the whole chain a second time against a
# different hit set without overwriting the first: every file this script reads or writes
# picks up the suffix.  Used for the stringent-FDR variant (CDK246_TAG=_strictFDR).
TAG = os.environ.get("CDK246_TAG", "")

sys.path.insert(0, HERE)
from cdk_config import STRUCTURES, ALIGN_REFERENCE

GENES = ["CDK2", "CDK4", "CDK6"]
POCKET_CUTOFF = 5.0     # A from ATP -> lines the pocket
SITE_CUTOFF = 12.0      # A from ATP -> "in the pocket"
CROP_CUTOFF = 16.0      # A from ATP -> inside the zoomed crop (reporting only)
LABEL_OFFSET = 4.0      # A; radial label offset out of the pocket.  2.5 was not enough
                        # to separate labels on adjacent residues (G165/V64/A167)

THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V','SEP':'S','TPO':'T','PTR':'Y','MSE':'M',
}

hits = list(csv.DictReader(open(os.path.join(HERE, f"cdk2_4_6_hit_residues{TAG}.csv"))))
out_lines, dist_rows, pocket_rows = [], [], []

def emit(s=""):
    print(s)
    out_lines.append(s)

run(session, "close session")

# ---- reference first; it never moves, so its ATP defines the common frame
ref = STRUCTURES[ALIGN_REFERENCE]
refm = run(session, f"open {os.path.join(HERE,'structures',ref['file'])}")[0]
atp = np.array([a.scene_coord for a in refm.atoms
                if a.residue.chain_id == ref["kinase_chain"] and a.residue.name == "ATP"])
emit(f"reference {ref['pdb']} chain {ref['kinase_chain']}: ATP has {len(atp)} atoms")
emit(f"ATP centroid (scene): {np.round(atp.mean(axis=0), 2).tolist()}")
emit(f"pocket cutoff {POCKET_CUTOFF} A | active-site cutoff {SITE_CUTOFF} A")
emit()

models = {ALIGN_REFERENCE: refm}
for gene in GENES:
    if gene == ALIGN_REFERENCE:
        continue
    cfg = STRUCTURES[gene]
    m = run(session, f"open {os.path.join(HERE,'structures',cfg['file'])}")[0]
    models[gene] = m
    res = run(session, f"matchmaker #{m.id[0]}/{cfg['kinase_chain']} "
                       f"to #{refm.id[0]}/{ref['kinase_chain']}")
    emit(f"{gene}: superposed {cfg['pdb']} chain {cfg['kinase_chain']} onto "
         f"{ref['pdb']}/{ref['kinase_chain']}")

emit()
for gene in GENES:
    cfg, m = STRUCTURES[gene], models[gene]
    kc = cfg["kinase_chain"]
    byres = {}
    for r in m.residues:
        if r.chain_id == kc:
            byres.setdefault(r.number, r)

    # pocket: residues lining the ATP site in this structure's own frame
    pocket = []
    for num, r in sorted(byres.items()):
        if r.name not in THREE_TO_ONE:
            continue
        c = np.array([a.scene_coord for a in r.atoms])
        if len(c) and np.linalg.norm(c[:, None, :] - atp[None, :, :], axis=2).min() < POCKET_CUTOFF:
            pocket.append(num)

    emit(f"=== {gene}  ({cfg['pdb']} chain {kc})  ligand={cfg['ligand']}")
    emit(f"    modelled range {min(byres)}-{max(byres)}")
    emit(f"    ATP pocket ({POCKET_CUTOFF} A): {len(pocket)} residues -> "
         + ",".join(str(p) for p in pocket))
    # the pocket is what the zoom panels frame on: conserved, structurally
    # equivalent across the three genes, and larger than ATP alone so nearby
    # hits stay in view
    pocket_rows.append(dict(gene=gene, chain=kc,
                            residues=",".join(str(p) for p in pocket),
                            n_residues=len(pocket)))

    # pocket centroid, for the radial label offsets below
    pc = np.array([a.scene_coord for num in pocket for a in byres[num].atoms]).mean(axis=0)

    gh = [h for h in hits if h["gene"] == gene]
    inside, outside, unmodelled, mismatched = [], [], [], []
    for h in gh:
        pos, wt = int(h["position"]), h["wt"]
        r = byres.get(pos)
        if r is None or r.name not in THREE_TO_ONE:
            unmodelled.append(f"{wt}{pos}")
            dist_rows.append(dict(gene=gene, position=pos, wt=wt, variants=h["variants"],
                                  drug_class=h["drug_class"], modelled=False,
                                  wt_matches="", min_dist_to_ATP="", active_site=False,
                                  in_crop=False))
            continue
        obs = THREE_TO_ONE[r.name]
        ok = (obs == wt)
        if not ok:
            mismatched.append(f"{wt}{pos} (model has {obs})")
        c = np.array([a.scene_coord for a in r.atoms])
        d = float(np.linalg.norm(c[:, None, :] - atp[None, :, :], axis=2).min())
        is_site = d < SITE_CUTOFF
        # push each label radially out of the pocket so neighbouring labels separate
        out_v = c.mean(axis=0) - pc
        n = np.linalg.norm(out_v)
        off = (out_v / n * LABEL_OFFSET) if n > 1e-6 else np.zeros(3)
        (inside if is_site else outside).append((f"{wt}{pos}", d, h["drug_class"]))
        dist_rows.append(dict(gene=gene, position=pos, wt=wt, variants=h["variants"],
                              drug_class=h["drug_class"], modelled=True,
                              wt_matches=ok, min_dist_to_ATP=round(d, 2),
                              active_site=is_site, in_crop=(d < CROP_CUTOFF),
                              label_dx=round(float(off[0]), 2),
                              label_dy=round(float(off[1]), 2),
                              label_dz=round(float(off[2]), 2)))

    emit(f"    hits: {len(gh)} positions | {len(inside)} within {SITE_CUTOFF} A of ATP, "
         f"{len(outside)} distal, {len(unmodelled)} not modelled")
    emit(f"    wt mismatches: {mismatched if mismatched else 'none'}")
    if unmodelled:
        emit(f"    not modelled: {', '.join(unmodelled)}")
    for label, group in (("ACTIVE SITE", inside), ("distal", outside)):
        if group:
            emit(f"    {label}:")
            for nm, d, cl in sorted(group, key=lambda t: t[1]):
                emit(f"        {nm:7s} {d:6.2f} A   {cl}")
    emit()

with open(os.path.join(HERE, f"cdk2_4_6_hit_distances{TAG}.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["gene", "position", "wt", "variants", "drug_class",
                                       "modelled", "wt_matches", "min_dist_to_ATP",
                                       "active_site", "in_crop",
                                       "label_dx", "label_dy", "label_dz"])
    w.writeheader()
    w.writerows(dist_rows)

with open(os.path.join(HERE, f"cdk2_4_6_pocket_residues{TAG}.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["gene", "chain", "n_residues", "residues"])
    w.writeheader()
    w.writerows(pocket_rows)

emit("=== hits inside the zoomed crop (< %.0f A of ATP), which get labeled" % CROP_CUTOFF)
for gene in GENES:
    gr = [r for r in dist_rows if r["gene"] == gene]
    inc = [r for r in gr if r["in_crop"]]
    off = [r for r in gr if r["modelled"] and not r["in_crop"]]
    un = [r for r in gr if not r["modelled"]]
    emit(f"  {gene}: {len(inc)} labeled, {len(off)} off-crop, {len(un)} unmodelled")
    emit("     labeled : " + ", ".join(f'{r["wt"]}{r["position"]}' for r in
                                       sorted(inc, key=lambda r: r["min_dist_to_ATP"])))
    emit("     off-crop: " + (", ".join(f'{r["wt"]}{r["position"]}' for r in
                                        sorted(off, key=lambda r: r["min_dist_to_ATP"])) or "none"))

with open(os.path.join(HERE, f"cdk2_4_6_active_site_report{TAG}.txt"), "w") as fh:
    fh.write("\n".join(out_lines) + "\n")
print("\nwrote cdk2_4_6_active_site_report.txt and cdk2_4_6_hit_distances.csv")
