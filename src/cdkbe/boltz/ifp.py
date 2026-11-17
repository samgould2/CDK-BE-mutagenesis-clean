"""Loaders, the interaction fingerprint and the weighted PCA.

Ported from `cdk9_ifp.ipynb` cells 1, 2, 14, 15, 19 and 23. The numbers are the
notebook's: the channel rules, the half-maximal distances and transition
widths, the 0.15 retention floor, the 0.15 weighting offset and the ATP-pocket
residue set are all carried over unchanged.

Two things in the notebook could not be carried over as written:

* cell 14's `RESULTS_DIR` is an absolute path into a collaborator's home
  directory. Paths here resolve from the repository root.
* the structures were copied in gzipped, and `MMCIFParser` wants a readable
  handle, so `parse_cif` opens either spelling.

**Two different ATP-pocket sets exist in the notebook.** Cell 15 builds
`CDK9_ATP_POCKET` out of its `CDK9_REGIONS` map, which has 33 residues; cell 19
then defines its own `ATP_POCKET` with 37, adding 47, 49, 65, 67 and 150. Every
PCA in the notebook -- and so every published panel -- uses cell 19's set, which
is also the one the README documents. `ATP_POCKET` here is that one;
`REGIONS` keeps cell 15's labelled map, which nothing draws but which is what
names the regions.
"""
from __future__ import annotations

import gzip
import json
import re
import shutil
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .. import io as _io
from .. import paths

RESULTS = paths.RAW / "cdk9_boltz" / "results"
#: the fingerprint is slow to build from 108 structures, so it is cached here
CACHE = paths.PROCESSED / "cdk9_boltz"

DRUG_ORDER = ["atuveciclib", "kb0742", "protac", "voruciclib"]
DRUG_LABELS = {"atuveciclib": "Atuveciclib", "kb0742": "KB-0742",
               "protac": "KI-CDK9d-32", "voruciclib": "Voruciclib"}
VARIANTS = ["WT", "L156F", "I157V"]
SEEDS = [635925368, 1821385089, 2036692060]

N_STRUCTURES = 108
N_AFFINITY = 36

#: cell 19's set -- the one every PCA in the notebook restricts to
ATP_POCKET = {
    25, 26, 27, 28, 29, 30, 31, 32, 33,
    46, 47, 48, 49,
    65, 66, 67,
    103, 104, 105, 106, 107, 108, 109,
    147, 148, 149, 150, 151,
    153, 154, 155, 156, 157,
    166, 167, 168, 169,
}

#: cell 15's labelled regions. Their union is 32 residues, NOT `ATP_POCKET`'s
#: 37: nothing here names 47, 49, 65, 67 or 150.
REGIONS = {
    "p_loop": {25, 26, 27, 28, 29, 30, 31, 32, 33},
    "beta3": {46, 48},
    "alphaC": {66},
    "gatekeeper": {103},
    "hinge": {104, 105, 106},
    "hinge_adjacent": {107, 108, 109},
    "catalytic_loop": {147, 148, 149, 151},
    "lower_pocket": {153, 154, 155, 156, 157},
    "dfg_region": {166, 167, 168, 169},
}

#: a feature is kept when its maximum score across all 108 structures reaches this
KEEP_MIN = 0.15
#: added to the mean within-group SD before inverting, so a feature that never
#: moves cannot acquire unbounded weight
WEIGHT_OFFSET = 0.15


# ---------------------------------------------------------------- confidence

def _seed_dirs():
    return sorted(p for p in RESULTS.glob("*/*/seed_*") if p.is_dir())


@lru_cache(maxsize=1)
def metrics_and_affinity():
    """-> (metrics, affinity): one row per structure, one row per seed."""
    per_model, per_seed = [], []
    for seed_dir in _seed_dirs():
        drug = seed_dir.parent.parent.name
        variant = seed_dir.parent.name
        seed = int(seed_dir.name.replace("seed_", ""))
        for cpath in sorted(seed_dir.rglob("confidence_*_model_*.json")):
            conf = json.loads(cpath.read_text())
            per_model.append(dict(
                drug=drug, variant=variant, seed=seed,
                model=int(re.search(r"_model_(\d+)\.json$", cpath.name).group(1)),
                cif_path=str(_cif_for(cpath)), confidence_path=str(cpath), **conf))
        for apath in sorted(seed_dir.rglob("affinity_*.json")):
            per_seed.append(dict(drug=drug, variant=variant, seed=seed,
                                 affinity_path=str(apath),
                                 **json.loads(apath.read_text())))
    m = (pd.DataFrame(per_model)
         .sort_values(["drug", "variant", "seed", "model"]).reset_index(drop=True))
    a = (pd.DataFrame(per_seed)
         .sort_values(["drug", "variant", "seed"]).reset_index(drop=True))
    return m, a


def model_summary() -> pd.DataFrame:
    """Confidence metrics averaged over the structural samples of each seed."""
    m, _ = metrics_and_affinity()
    return m.groupby(["drug", "variant", "seed"], as_index=False).agg(
        confidence_score_mean=("confidence_score", "mean"),
        confidence_score_sd=("confidence_score", "std"),
        iptm_mean=("iptm", "mean"), iptm_sd=("iptm", "std"),
        ligand_iptm_mean=("ligand_iptm", "mean"), ligand_iptm_sd=("ligand_iptm", "std"),
        complex_plddt_mean=("complex_plddt", "mean"),
        complex_plddt_sd=("complex_plddt", "std"),
        complex_iplddt_mean=("complex_iplddt", "mean"),
        complex_iplddt_sd=("complex_iplddt", "std"),
        complex_pde_mean=("complex_pde", "mean"),
        complex_ipde_mean=("complex_ipde", "mean"))


def delta_binding_probability(variant: str = "L156F") -> pd.DataFrame:
    """P(bind, mutant) - P(bind, WT), paired within a seed.

    Pairing by seed is what makes the difference meaningful: the two
    predictions being subtracted used the same sampling identifier.
    """
    _, affinity = metrics_and_affinity()
    wide = (affinity[["drug", "variant", "seed", "affinity_probability_binary"]]
            .pivot(index=["drug", "seed"], columns="variant",
                   values="affinity_probability_binary").reset_index())
    d = wide[["drug", "seed", "WT", variant]].copy()
    d["variant"] = variant
    d["delta_binding_probability"] = d[variant] - d["WT"]
    return d[["drug", "seed", "variant", "WT", variant,
              "delta_binding_probability"]]


def seed_level(variant: str = "WT") -> pd.DataFrame:
    """Mean iPTM against predicted binding probability, one row per seed."""
    m, affinity = metrics_and_affinity()
    s = m.groupby(["drug", "variant", "seed"], as_index=False).agg(
        iptm=("iptm", "mean"))
    s = s.merge(affinity[["drug", "variant", "seed", "affinity_probability_binary"]],
                on=["drug", "variant", "seed"], how="inner", validate="one_to_one")
    s = s.rename(columns={"affinity_probability_binary": "binding_probability"})
    return s[s["variant"] == variant].copy()


# ----------------------------------------------------------------- structures

def _cif_for(confidence_path: Path) -> Path:
    """The structure beside a confidence file, in whichever spelling is there."""
    stem = confidence_path.name.replace("confidence_", "").replace(".json", ".cif")
    return _io.resolve(confidence_path.with_name(stem))


def parse_cif(path) -> dict:
    """Non-hydrogen coordinates, split into protein (chain A) and ligand (B)."""
    from Bio.PDB import MMCIFParser

    path = Path(path)
    parser = MMCIFParser(QUIET=True)
    if path.suffix == ".gz":
        # MMCIFParser needs a seekable plain-text handle
        with gzip.open(path, "rb") as gz, \
                tempfile.NamedTemporaryFile(suffix=".cif") as tmp:
            shutil.copyfileobj(gz, tmp)
            tmp.flush()
            structure = parser.get_structure("model", tmp.name)
    else:
        structure = parser.get_structure("model", str(path))

    model = next(structure.get_models())
    p = dict(xyz=[], res=[], resn=[], name=[], elem=[])
    lig = dict(xyz=[], name=[], elem=[])
    for chain in model:
        for residue in chain:
            resnum, resname = residue.id[1], residue.resname.strip()
            for atom in residue:
                elem = (atom.element or "").strip().upper()
                if elem == "H":
                    continue
                if chain.id == "A":
                    p["xyz"].append(atom.coord.astype(float))
                    p["res"].append(resnum)
                    p["resn"].append(resname)
                    p["name"].append(atom.name.strip())
                    p["elem"].append(elem)
                elif chain.id == "B":
                    lig["xyz"].append(atom.coord.astype(float))
                    lig["name"].append(atom.name.strip())
                    lig["elem"].append(elem)
    return {
        "prot_xyz": np.asarray(p["xyz"], float), "prot_res": np.asarray(p["res"], int),
        "prot_resn": np.asarray(p["resn"]), "prot_name": np.asarray(p["name"]),
        "prot_elem": np.asarray(p["elem"]),
        "lig_xyz": np.asarray(lig["xyz"], float), "lig_name": np.asarray(lig["name"]),
        "lig_elem": np.asarray(lig["elem"]),
    }


@lru_cache(maxsize=1)
def structures() -> tuple:
    """All 108 parsed structures, in (drug, variant, seed, model) order."""
    out = []
    for cpath in sorted(p for d in _seed_dirs()
                        for p in d.rglob("confidence_*_model_*.json")):
        cif = _cif_for(cpath)
        rel = cpath.relative_to(RESULTS)
        s = parse_cif(cif)
        s.update(drug=rel.parts[0], variant=rel.parts[1],
                 seed=int(rel.parts[2].replace("seed_", "")),
                 model=int(re.search(r"_model_(\d+)\.json$", cpath.name).group(1)),
                 cif_path=str(cif))
        out.append(s)
    if len(out) != N_STRUCTURES:
        raise RuntimeError(f"expected {N_STRUCTURES} structures, parsed {len(out)}")
    return tuple(out)


@lru_cache(maxsize=1)
def residue_labels() -> dict:
    """{residue number: "ASP167"}, from the first WT structure."""
    ref = next(s for s in structures() if s["variant"] == "WT")
    return {int(r): f"{n}{int(r)}"
            for r, n in zip(ref["prot_res"], ref["prot_resn"])}


# --------------------------------------------------------- the fingerprint

def soft_distance(d, mid, tau):
    if not np.isfinite(d):
        return 0.0
    return 1.0 / (1.0 + np.exp(np.clip((d - mid) / tau, -50, 50)))


VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
             "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98}

BB = {"N", "CA", "C", "O", "OXT"}
CHARGED_PROTEIN_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"),
                         ("GLU", "OE2"), ("LYS", "NZ"), ("ARG", "NE"),
                         ("ARG", "NH1"), ("ARG", "NH2")}
AROMATIC_PROTEIN_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}

#: channel -> (d_mid, tau, protein-atom rule, ligand-atom rule)
CHANNELS = {
    "carbon_contact": (4.0, 0.6,
                       lambda rn, nm, el: nm not in BB and str(el).upper() == "C",
                       lambda el: str(el).upper() == "C"),
    "polar_bb": (3.3, 0.4,
                 lambda rn, nm, el: nm in {"N", "O", "OXT"},
                 lambda el: str(el).upper() in {"N", "O"}),
    "polar_sc": (3.3, 0.4,
                 lambda rn, nm, el: nm not in BB and str(el).upper() in {"N", "O", "S"},
                 lambda el: str(el).upper() in {"N", "O"}),
    "charged_polar": (4.0, 0.6,
                      lambda rn, nm, el: (rn, nm) in CHARGED_PROTEIN_ATOMS,
                      lambda el: str(el).upper() in {"N", "O"}),
    "aromatic_C": (4.5, 0.6,
                   lambda rn, nm, el: nm in AROMATIC_PROTEIN_ATOMS.get(rn, set()),
                   lambda el: str(el).upper() == "C"),
}


def steric_clash_score(distances, prot_elements, lig_elements,
                       ratio_mid=0.80, tau=0.05):
    """How far the most overlapping atom pair is inside its van der Waals sum.

    This is the single worst pair, not a count of close contacts.
    """
    if distances.size == 0:
        return 0.0
    rp = np.array([VDW_RADII.get(str(e).upper(), 1.70) for e in prot_elements])
    rl = np.array([VDW_RADII.get(str(e).upper(), 1.70) for e in lig_elements])
    ratio = distances / (rp[:, None] + rl[None, :])
    x = np.clip((np.nanmin(ratio) - ratio_mid) / tau, -50, 50)
    return 1.0 / (1.0 + np.exp(x))


def _build_fingerprint() -> pd.DataFrame:
    rows = []
    all_res = sorted({int(r) for s in structures()
                      for r in np.asarray(s["prot_res"])})
    for s in structures():
        P, L = s["prot_xyz"], s["lig_xyz"]
        if len(P) == 0 or len(L) == 0:
            raise ValueError(f"empty chain in {s['cif_path']}")
        d = np.sqrt(((P[:, None, :] - L[None, :, :]) ** 2).sum(axis=-1))
        resn, name, elem = s["prot_resn"], s["prot_name"], s["prot_elem"]
        lig_elem = s["lig_elem"]
        # the per-atom rules depend only on the structure, so they are evaluated
        # once per channel rather than once per residue per channel
        chan_mask = {ch: (np.array([psel(a, b, c)
                                    for a, b, c in zip(resn, name, elem)]),
                          np.array([lsel(e) for e in lig_elem]))
                     for ch, (_, _, psel, lsel) in CHANNELS.items()}
        rec = {k: s[k] for k in ("drug", "variant", "seed", "model")}
        for r in all_res:
            rmask = s["prot_res"] == r
            for ch, (mid, tau, _, _) in CHANNELS.items():
                pm, lm = chan_mask[ch]
                pm = rmask & pm
                dmin = d[np.ix_(pm, lm)].min() if (pm.any() and lm.any()) else np.inf
                rec[f"{r}:{ch}"] = soft_distance(dmin, mid, tau)
            rec[f"{r}:steric_clash"] = steric_clash_score(
                d[np.ix_(rmask, np.ones(len(lig_elem), bool))],
                np.asarray(elem)[rmask], np.asarray(lig_elem))
        rows.append(rec)
    ifp = (pd.DataFrame(rows)
           .set_index(["drug", "variant", "seed", "model"]).sort_index())
    return ifp.loc[:, ifp.max(axis=0) >= KEEP_MIN]


def fingerprint(rebuild: bool = False) -> pd.DataFrame:
    """The retained residue x descriptor matrix, 108 rows.

    Built from the structures and then cached, because parsing 108 CIFs and
    forming every atom-pair distance takes minutes. `rebuild=True` ignores the
    cache; `tests/test_cdk9_boltz.py` rebuilds a slice of it to check the cache
    still matches what the structures give.
    """
    cached = CACHE / "ifp.csv.gz"
    if not rebuild:
        try:
            return _io.read(cached, index_col=[0, 1, 2, 3])
        except FileNotFoundError:
            pass
    ifp = _build_fingerprint()
    CACHE.mkdir(parents=True, exist_ok=True)
    ifp.to_csv(cached)
    return ifp


# ------------------------------------------------------------------- the PCA

def pocket_matrix(variant: str = "WT") -> pd.DataFrame:
    """The fingerprint restricted to one variant and to the ATP pocket."""
    ifp = fingerprint()
    sub = ifp[ifp.index.get_level_values("variant") == variant]
    cols = [c for c in sub.columns if int(c.split(":")[0]) in ATP_POCKET]
    return sub[cols]


#: An SVD gives each component an arbitrary sign: flipping a PC and its loadings
#: together describes the same fit. Which sign LAPACK returns depends on the
#: machine, so these anchors pin the published orientation. Each names a feature
#: and the sign its loading has in the published panel -- PC1's largest-magnitude
#: loading in supplementary 6i, PC2's in 6j.
SIGN_ANCHORS = {0: ("167:charged_polar", -1.0), 1: ("154:polar_sc", -1.0)}


def _orient(vt, columns):
    """Flip each component so its anchor feature keeps the published sign."""
    vt = vt.copy()
    cols = list(columns)
    for pc, (feature, want) in SIGN_ANCHORS.items():
        if pc >= len(vt) or feature not in cols:
            continue
        loading = vt[pc][cols.index(feature)]
        if loading != 0 and np.sign(loading) != want:
            vt[pc] = -vt[pc]
    return vt


def weighted_pca(X: pd.DataFrame, group_level="drug", orient: bool = True):
    """-> (scores, explained %, loadings, weights).

    Contacts that fluctuate between predictions of the same ligand are
    downweighted: each feature's within-group SD is averaged over groups, and
    the feature is scaled by 1 / (0.15 + that). Without it the leading PCs
    describe sampling noise in the flexible loops rather than binding mode.

    `orient` applies `SIGN_ANCHORS`, so the panels come out the way they were
    published rather than mirrored on whichever axis this machine's LAPACK
    happened to flip. Pass False for the raw decomposition.
    """
    sd = X.groupby(level=group_level).transform("std").mean(axis=0)
    w = 1.0 / (WEIGHT_OFFSET + sd.values)
    Xw = (X.values - X.values.mean(axis=0)) * np.sqrt(w)[None, :]
    _, sv, vt = np.linalg.svd(Xw, full_matrices=False)
    if orient:
        vt = _orient(vt, X.columns)
    scores = Xw @ vt[:2].T
    var = (sv ** 2 / (sv ** 2).sum())[:2] * 100
    return scores, var, vt, w


def pc_loadings(pc: int = 0, variant: str = "WT") -> pd.Series:
    """PC loadings for the WT pocket PCA, labelled `ASP167:charged_polar`.

    The sign of a PCA axis is arbitrary: positive and negative describe
    opposite ends of this fit, not favourable and unfavourable interactions.
    """
    X = pocket_matrix(variant)
    _, _, vt, _ = weighted_pca(X)
    lab = residue_labels()
    return pd.Series(
        vt[pc],
        index=[f"{lab[int(c.split(':')[0])]}:{c.split(':')[1]}" for c in X.columns],
        name=f"PC{pc + 1}_loading")
