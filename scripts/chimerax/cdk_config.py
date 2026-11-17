"""
Shared configuration for CDK-family ChimeraX visualizations.

Structure choices, colors and rendering settings follow notebooks/10_PyMOL
(viz1.cxs / viz2_cdk2.cxs / viz3.cxs and align_catalytic_v2.py), extended to
cover CDK4 and CDK6, which were not in the earlier figures.

Fields per gene
---------------
pdb            PDB ID (or "AF-<acc>" for an AlphaFold model)
file           local filename under ./structures/
url            where `fetch_structures.py` gets it
kinase_chain   chain carrying the CDK catalytic domain
cyclin_chains  partner cyclin chain(s); rendered as a surface
other_chains   additional chains in the file (MAT1, DDB1, substrate peptide);
               hidden by default
kd             catalytic (KLIFS) domain residue range, used for the "fade
               everything outside the catalytic core" transparency.  CDK7/8/9/
               12/13/19 are the ranges from align_catalytic_v2.py; CDK2/4/6 were
               derived from the superposition by map_kd_ranges.py (all six known
               ranges map onto CDK2 8-149, which was then mapped out to CDK4/6)
ligand         inhibitor / nucleotide residue name shown as sticks (or None)
display_range  optional; restricts the cartoon to this residue range.  Used for
               the CDK19 AlphaFold model, whose C-terminal tail is disordered
borrowed_cyclin optional; {donor, transparency}.  Superposes the donor CDK onto
               this one and shows the donor's cyclin, drawn semi-transparent
               because it is not part of this structure
"""

STRUCTURES = {
    "CDK7": dict(
        pdb="6XD3", file="6xd3.cif", url="https://files.rcsb.org/download/6XD3.cif",
        kinase_chain="J", cyclin_chains=["I"], other_chains=["H"],
        kd=(16, 159), ligand="V0G", cyclin_name="Cyclin H",
        note="CDK7/CycH/MAT1 + THZ1; chain H = MAT1"),
    "CDK8": dict(
        pdb="5CEI", file="5cei.cif", url="https://files.rcsb.org/download/5CEI.cif",
        kinase_chain="A", cyclin_chains=["B"], other_chains=[],
        kd=(25, 177), ligand="50R", cyclin_name="Cyclin C",
        note="CDK8/CycC + compound 22"),
    "CDK19": dict(
        pdb="AF-Q9BWU1", file="AF-Q9BWU1-F1-model_v4.pdb",
        url="https://alphafold.ebi.ac.uk/files/AF-Q9BWU1-F1-model_v4.pdb",
        kinase_chain="A", cyclin_chains=[], other_chains=[],
        kd=(25, 177), ligand=None, cyclin_name=None,
        display_range=(1, 363),
        # CDK19 is the CDK8 paralog and shares cyclin C, but no CDK19-cyclin C
        # structure exists.  Borrow cyclin C from 5CEI by superposing CDK8 onto
        # CDK19, and render it more transparent to flag that it is inferred.
        borrowed_cyclin=dict(donor="CDK8", transparency=60),
        note="no experimental CDK19 structure; AlphaFold monomer, no cyclin "
             "partner. Only 1-363 is confidently modelled (pLDDT >= 50); the "
             "364-502 tail is disordered spaghetti and is hidden, exactly as "
             "residues 351-502 were hidden in align_catalytic.cxc"),
    "CDK9": dict(
        pdb="8K5R", file="8k5r.cif", url="https://files.rcsb.org/download/8K5R.cif",
        kinase_chain="A", cyclin_chains=["B"], other_chains=[],
        kd=(23, 171), ligand="VQE", cyclin_name="Cyclin T1",
        note="CDK9/CycT1 + KB-0742"),
    "CDK12": dict(
        pdb="8BUG", file="8bug.cif", url="https://files.rcsb.org/download/8BUG.cif",
        kinase_chain="B", cyclin_chains=["C"], other_chains=["A"],
        kd=(731, 881), ligand="RPW", cyclin_name="Cyclin K",
        note="DDB1/CDK12/CycK + CR8 glue; using copy 1 (chains A,B,C) because the RPW\n             ligand is only modelled there; chain A = DDB1"),
    "CDK13": dict(
        pdb="7NXJ", file="7nxj.cif", url="https://files.rcsb.org/download/7NXJ.cif",
        kinase_chain="A", cyclin_chains=["B"], other_chains=[],
        kd=(709, 859), ligand="5I1", cyclin_name="Cyclin K",
        note="CDK13/CycK + THZ1"),
    "CDK2": dict(
        pdb="1QMZ", file="1qmz.cif", url="https://files.rcsb.org/download/1QMZ.cif",
        kinase_chain="A", cyclin_chains=["B"], other_chains=["E"],
        kd=(8, 149), ligand="ATP", cyclin_name="Cyclin A",
        note="pCDK2/CycA + ATP + substrate peptide (chain E); alignment reference"),
    "CDK4": dict(
        pdb="2W96", file="2w96.cif", url="https://files.rcsb.org/download/2W96.cif",
        kinase_chain="B", cyclin_chains=["A"], other_chains=[],
        kd=(10, 160), ligand=None, cyclin_name="Cyclin D1",
        note="CDK4/CycD1, 2.3 A, apo active site"),
    "CDK6": dict(
        pdb="2F2C", file="2f2c.cif", url="https://files.rcsb.org/download/2F2C.cif",
        kinase_chain="B", cyclin_chains=["A"], other_chains=[],
        kd=(18, 167), ligand="AP9", cyclin_name="V-cyclin",
        note="CDK6/KSHV V-cyclin + aminopurvalanol; no human cyclin D-CDK6 "
             "structure exists in the PDB, so the cyclin here is viral"),
}

# structure-based superposition reference (same choice as align_catalytic.cxc)
ALIGN_REFERENCE = "CDK2"

# per-CDK cartoon colors; CDK7-CDK13 and CDK2 unchanged from the earlier sessions
CDK_COLORS = {
    "CDK7":  "#00ffff",
    "CDK8":  "#1f77b4",
    "CDK19": "#6495ed",
    "CDK9":  "#9400d3",
    "CDK12": "#2ca02c",
    "CDK13": "#556b2f",
    "CDK2":  "#daa520",
    "CDK4":  "#e377c2",   # new
    "CDK6":  "#ff7f0e",   # new
}

CYCLIN_COLOR = "#f5deb3"          # wheat, as in viz1/viz2/viz3
CYCLIN_TRANSPARENCY = 0           # percent; surfaces were opaque in viz1-3
OTHER_CHAIN_COLOR = "#d3d3d3"

# every variant in the *_filtered_top_only tables is an enriched (resistance)
# hit, so a single uniform red is used rather than the red/blue pair from viz3.
HIT_COLOR = "#ff0000"
HIT_STYLE = "sphere"

NONCATALYTIC_TRANSPARENCY = 80    # matches cartoon_transparency 0.8 in the PyMOL script

# reproduces the render settings stored in viz1/viz2/viz3
RENDER_PREAMBLE = """\
set bgColor white
lighting soft
lighting multiShadow 64 msMapSize 1024
graphics silhouettes true width 2
set subdivision 3
"""

# Canonical camera, taken verbatim from viz2_cdk2.cxs.  Because every structure
# is superposed onto 1QMZ chain A (which is never moved), this one matrix gives
# all nine panels an identical orientation and scale.
CANONICAL_VIEW = ",".join([
    "-0.50925593", "-0.15998363", "-0.84561436", "-119.48200801",
    "0.34444448", "0.86255213", "-0.37062357", "-17.74379172",
    "0.78868017", "-0.48000945", "-0.38415430", "-66.99018438",
])

# Uniform pull-back applied after the canonical camera.  The viz2 camera was
# framed on CDK2 alone; the larger CDK4/CDK6/CDK19 models clip at that distance,
# so every panel backs off by the same factor and stays directly comparable.
ZOOM_OUT = 0.72

# render window aspect; `save ... width N` follows it
WINDOW_SIZE = (1400, 1400)

# common names for the CCD ligand codes above, for the provenance table
LIGAND_NAMES = {
    "ATP": "ATP", "V0G": "THZ1", "50R": "compound 22", "VQE": "KB-0742",
    "RPW": "CR8", "5I1": "THZ531", "AP9": "aminopurvalanol",
}
