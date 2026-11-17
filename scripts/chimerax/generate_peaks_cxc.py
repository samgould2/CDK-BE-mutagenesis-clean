"""
Project the resistance *peaks* of the CDK alignment onto a single structure.

generate_consensus_cxc.py draws the individual alignment columns that recur
across CDKs.  This script instead draws the peaks that
notebooks/09_analysis_and_viz/canonical_resistance_figures.py detects in the
smoothed enrichment profile -- the same peaks marked with a red x in
`MSA_resistance_peaks_kinase_domain.pdf` -- so the structural figure and the
profile figure show the same objects.

A peak is a region, not a column: its extent is the half-prominence width from
`scipy.signal.peak_widths`.  Every alignment column inside that window that
carries an enriched call is mapped onto the reference kinase, and all of those
residues are coloured as one peak, so each cluster of spheres on the structure
is one bump in the profile.

Spheres are coloured by how many CDKs contribute enriched calls to the peak;
peaks whose apex sits on a numbered KLIFS position use the red ramp, those in a
KLIFS gap the blue-grey ramp, matching generate_consensus_cxc.py.

    python3 generate_peaks_cxc.py --set core7
    python3 generate_peaks_cxc.py --set all9 --reference CDK2
"""
import argparse, collections, os, sys

import numpy as np
import pandas as pd
from Bio import AlignIO
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
from cdk_config import (STRUCTURES, CDK_COLORS, CYCLIN_COLOR, OTHER_CHAIN_COLOR,
                        NONCATALYTIC_TRANSPARENCY, RENDER_PREAMBLE, CANONICAL_VIEW,
                        WINDOW_SIZE)

KLIFS_DIR = os.path.join(ROOT, "data", "reference", "KLIFS")
ALN_FILE = os.path.join(KLIFS_DIR, "clustalo-I20250722-154001-0730-44590675-p1m.aln-clustal_num")

GENES_ALL = ["CDK7", "CDK8", "CDK19", "CDK9", "CDK12", "CDK13", "CDK2", "CDK4", "CDK6"]
EXCLUDE = {"core7": ["CDK8", "CDK19"], "all9": []}

# same ramps as the consensus panel, indexed by number of contributing CDKs
POCKET_RAMP = {1: "#fee5d9", 2: "#fcae91", 3: "#fb6a4a", 4: "#cb181d", 5: "#67000d"}
GAP_RAMP = {1: "#eff3ff", 2: "#c6dbef", 3: "#6baed6", 4: "#2171b5", 5: "#08306b"}

# identical to canonical_resistance_figures.py, so the peaks are the same ones
SIGMA, HEIGHT, DISTANCE, PROMINENCE = 5, 0.2, 10, 0.3


def ramp(d, n):
    return d[min(max(n, 1), 5)]


def alignment_maps():
    aln = AlignIO.read(ALN_FILE, "clustal")
    res2col, col2res, aln_char = {}, {}, {}
    for rec in aln:
        d, inv, n = {}, {}, 0
        for c, ch in enumerate(str(rec.seq)):
            if ch != "-":
                n += 1
                d[n], inv[c] = c, n
        res2col[rec.id], col2res[rec.id] = d, inv
        aln_char[rec.id] = str(rec.seq)
    return aln.get_alignment_length(), res2col, col2res, aln_char


def load_peaks(aset, reference):
    """One dict per detected peak, with the reference kinase's residues inside it."""
    calls = pd.read_csv(os.path.join(
        ROOT, "data", "processed", "lfc_fdr_tables", "resistance_patterns", aset,
        "resistance_calls_canonical.csv"))
    aln_len, res2col, col2res, aln_char = alignment_maps()

    y = np.zeros(aln_len)
    np.add.at(y, calls.aln_col.values, 1)
    ys = gaussian_filter1d(y, sigma=SIGMA)
    idx, _ = find_peaks(ys, height=HEIGHT, distance=DISTANCE, prominence=PROMINENCE)
    _, _, left, right = peak_widths(ys, idx, rel_height=0.5)

    # KLIFS label per alignment column, via each gene's own pocket residues
    KL = pd.read_csv(os.path.join(KLIFS_DIR, "KLIFS_compiled.csv"))
    genes = [g for g in GENES_ALL if g not in EXCLUDE[aset]]
    col_klifs = {}
    for _, row in KL.iterrows():
        for g in genes:
            if pd.notna(row[g]):
                c = res2col[g].get(int(row[g]))
                if c is not None:
                    col_klifs.setdefault(c, row["KLIFS_idx"])

    by_col = calls.groupby("aln_col")
    peaks = []
    for rank, (apex, lo, hi) in enumerate(zip(idx, left, right), start=1):
        lo, hi = int(np.floor(lo)), int(np.ceil(hi))
        cols = [c for c in range(lo, hi + 1) if c in by_col.groups]
        sub = calls[calls.aln_col.isin(cols)]
        resis = sorted({col2res[reference][c] for c in cols
                        if c in col2res[reference]})
        labels = [col_klifs[c] for c in cols if c in col_klifs]
        apex_label = col_klifs.get(int(apex))
        # the apex itself may be a gap in the reference, so fall back to the
        # nearest column inside the peak that the reference does have
        near = sorted((abs(c - apex), c) for c in cols if c in col2res[reference])
        apex_resi = col2res[reference][near[0][1]] if near else None
        apex_aa = aln_char[reference][near[0][1]] if near else ""
        peaks.append(dict(
            rank=rank, apex=int(apex), lo=lo, hi=hi, height=float(ys[apex]),
            cols=cols, resis=resis, apex_resi=apex_resi, apex_aa=apex_aa,
            n_tests=len(sub),
            n_cdks=sub.Gene.nunique(), cdks=",".join(sorted(sub.Gene.unique())),
            in_klifs=apex_label is not None,
            label=apex_label or (labels[0] if labels else "C-lobe"),
            variants=",".join(sorted({f"{g}-{s}" for g, s in zip(sub.Gene, sub.SAV)}))))
    return peaks, calls


def build(peaks, reference, labeled):
    cfg = STRUCTURES[reference]
    kc, cyc, other = cfg["kinase_chain"], cfg["cyclin_chains"], cfg["other_chains"]
    lo, hi = cfg["kd"]
    L = []
    L.append("# " + "=" * 74)
    L.append(f"# Resistance peaks of the CDK alignment projected onto {reference} "
             f"({cfg['pdb']}).")
    L.append("# Generated by generate_peaks_cxc.py -- edit that script, not this file.")
    L.append(f"# {len(peaks)} peaks; each sphere cluster is one bump in the smoothed")
    L.append("# enrichment profile (MSA_resistance_peaks_kinase_domain.pdf).")
    L.append("# Sphere colour = number of CDKs contributing; reds = the peak apex is a")
    L.append("# numbered KLIFS position, blues = the apex falls in a KLIFS gap.")
    L.append("# " + "=" * 74)
    L.append("close session")
    L.append(f"cd {HERE}")
    L.append("")
    L.append(f"open structures/{cfg['file']} name {reference}")
    L.append("hide #1 target acsb")
    L.append(f"show #1/{kc} cartoon")
    L.append(f"color #1/{kc} {CDK_COLORS[reference]} target c")
    L.append(f"transparency #1/{kc} & ~:{lo}-{hi} {NONCATALYTIC_TRANSPARENCY} cartoons")
    if cyc:
        cs = "#1/" + ",".join(cyc)
        L.append(f"# {cfg['cyclin_name']} as a see-through surface")
        L.append(f"show {cs} atoms")
        L.append(f"surface {cs}")
        L.append(f"hide {cs} atoms")
        L.append(f"color {cs} {CYCLIN_COLOR} target s")
        L.append(f"transparency {cs} 55 surfaces")
    if cfg["ligand"]:
        ls = f"#1/{kc} & :{cfg['ligand']}"
        L.append(f"show {ls} target ab")
        L.append(f"style {ls} stick")
        L.append(f"color {ls} #909090 target ab")
        L.append(f"color {ls} byhetero")
        if reference == "CDK2":
            L.append(f"show #1/{kc} & :MG atoms")
            L.append(f"style #1/{kc} & :MG sphere")
            L.append(f"size #1/{kc} & :MG atomRadius 0.8")
    if other:
        L.append("color #1/" + ",".join(other) + f" {OTHER_CHAIN_COLOR}")

    L.append("\n# ------------------------------------------------- resistance peaks")
    for p in peaks:
        if not p["resis"]:
            L.append(f"# peak {p['rank']} (col {p['apex']}, {p['label']}): no "
                     f"{reference} residue in the window -- gap in the reference")
            continue
        col = ramp(POCKET_RAMP if p["in_klifs"] else GAP_RAMP, p["n_cdks"])
        sel = f"#1/{kc}:" + ",".join(str(r) for r in p["resis"])
        L.append(f"# peak {p['rank']}: col {p['apex']} ({p['label']}), height "
                 f"{p['height']:.2f}, {p['n_cdks']} CDKs [{p['cdks']}], "
                 f"{len(p['resis'])} {reference} residues")
        # the peak is a region, so its whole extent is painted onto the ribbon;
        # spheres mark only the apex, otherwise the clusters merge into one blob
        # and the fold behind them stops being readable
        L.append(f"color {sel} {col} target c")
        L.append(f"transparency {sel} 0 cartoons")
        apex_sel = f"#1/{kc}:{p['apex_resi']}"
        L.append(f"show {apex_sel} atoms")
        L.append(f"cartoon {apex_sel} suppressBackboneDisplay false")
        L.append(f"style {apex_sel} sphere")
        L.append(f"color {apex_sel} {col} target ac")
        L.append(f"name peak{p['rank']} {sel}")
    allres = sorted({r for p in peaks for r in p["resis"]})
    if allres:
        L.append(f"name peak_sites #1/{kc}:" + ",".join(str(r) for r in allres))

    L.append("\n# ---------------------------------------------------- look and camera")
    L.extend(RENDER_PREAMBLE.strip().splitlines())
    L.append("show #* models")
    L.append("label delete")
    L.append("hide #* target p")
    L.append("hide @@element=H target ab")
    L.append(f"view matrix camera {CANONICAL_VIEW}")
    L.append("zoom 0.88")

    if labeled:
        L.append("\n# --------------------------------------------------- peak labels")
        L.append("# one label per peak, on its apex residue; the rank disambiguates")
        L.append("# the several peaks that share the generic 'C-lobe' name")
        for p in peaks:
            if not p["resis"]:
                continue
            L.append(f'label #1/{kc}:{p["apex_resi"]} residues '
                     f'text "P{p["rank"]} {p["label"]}" '
                     f'height 2.4 color black bgColor none')
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", default="core7", choices=["core7", "all9"])
    ap.add_argument("--reference", default="CDK2")
    ap.add_argument("--tag", help="output directory suffix (default: peaks)")
    a = ap.parse_args()

    peaks, calls = load_peaks(a.set, a.reference)
    out = os.path.join(HERE, f"cxc_{a.set}_{a.tag or 'peaks'}")
    os.makedirs(out, exist_ok=True)
    open(os.path.join(out, "peaks.cxc"), "w").write(build(peaks, a.reference, False))
    open(os.path.join(out, "peaks_labeled.cxc"), "w").write(
        build(peaks, a.reference, True))

    tag = os.path.basename(out)
    fig = os.path.join(ROOT, "figures", "chimerax", tag)
    os.makedirs(fig, exist_ok=True)
    rl = ["# Render the peak panels; needs a GPU context",
          f"cd {HERE}", f"windowsize {WINDOW_SIZE[0]} {WINDOW_SIZE[1]}"]
    for name in ("peaks", "peaks_labeled"):
        rl.append(f"\nopen {tag}/{name}.cxc")
        rl.append(f"save {fig}/{name}.png width 3000 supersample 3 "
                  f"transparentBackground true")
    rl.append("\nclose session")
    open(os.path.join(out, "render_all.cxc"), "w").write("\n".join(rl) + "\n")

    tdir = os.path.join(ROOT, "data", "processed", "lfc_fdr_tables",
                        "resistance_patterns", a.set)
    # provenance table alongside the figure
    pd.DataFrame([{k: v for k, v in p.items() if k != "cols"} for p in peaks]).assign(
        resis=lambda d: d.resis.apply(lambda r: ",".join(map(str, r)))).to_csv(
        os.path.join(tdir, "MSA_peak_sites.csv"), index=False)

    # second copy in the column layout label_consensus.py reads, so the peak panel
    # can borrow the same leader-line labeller as the consensus panel
    pd.DataFrame([dict(
        aln_col=p["apex"], n_cdks=p["n_cdks"], tier="recurrent",
        in_KLIFS=p["in_klifs"], KLIFS_idx=p["label"] if p["in_klifs"] else "",
        segment="" if p["in_klifs"] else p["label"],
        region="C-lobe post-pocket" if p["label"] == "C-lobe" else "",
        cdks=p["cdks"], variants=p["variants"],
        residues=f'{a.reference}:{p["apex_aa"]}{p["apex_resi"]}')
        for p in peaks if p["resis"]]).to_csv(
        os.path.join(tdir, "MSA_peak_sites_labels.csv"), index=False)

    print(f"{len(peaks)} peaks projected onto {a.reference} "
          f"({STRUCTURES[a.reference]['pdb']})")
    for p in peaks:
        print(f"  peak {p['rank']:>2}  col {p['apex']:>5}  win {p['lo']}-{p['hi']:<5} "
              f"h {p['height']:4.2f}  {p['n_cdks']} CDKs  "
              f"{len(p['resis']):>2} {a.reference} res  {p['label']:<11} {p['cdks']}")
    missing = [p["rank"] for p in peaks if not p["resis"]]
    if missing:
        print(f"  peaks with no {a.reference} residue (gap in the reference): {missing}")
    print(f"wrote {out}/peaks.cxc, peaks_labeled.cxc, render_all.cxc")


if __name__ == "__main__":
    main()
