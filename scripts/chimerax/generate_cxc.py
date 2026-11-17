"""
Generate ChimeraX command files for the top CDK resistance variants.

Writes, into ./cxc/ :
    00_all_cdks_aligned.cxc   all nine CDKs superposed on CDK2 (1QMZ chain A)
    CDK<n>.cxc                one single-kinase panel per gene

Every script uses the same structures, the same colors, the same lighting and
the same camera matrix, so panels can be placed side by side without any
per-figure fiddling.  Conventions carried over from notebooks/10_PyMOL:

    * CDK catalytic domain  : cartoon, per-gene color
    * outside the catalytic domain : same color, 80% transparent
    * cyclin partner        : molecular surface, wheat
    * inhibitor / nucleotide: sticks, element-colored on gray carbons
    * resistance hits       : red CPK spheres
    * white background, soft lighting, silhouettes

Usage:  python3 generate_cxc.py
        then, in ChimeraX:  open cxc/00_all_cdks_aligned.cxc
"""
import argparse, csv, os, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cdk_config import (STRUCTURES, ALIGN_REFERENCE, CDK_COLORS, CYCLIN_COLOR,
                        OTHER_CHAIN_COLOR, HIT_COLOR, HIT_STYLE,
                        NONCATALYTIC_TRANSPARENCY, RENDER_PREAMBLE, CANONICAL_VIEW,
                        ZOOM_OUT, WINDOW_SIZE, LIGAND_NAMES)

OUT = os.path.join(HERE, "cxc")
HITS_CSV = os.path.join(HERE, "hit_residues.csv")
# reference first so it is always model #1 and never gets moved by matchmaker
ORDER = [ALIGN_REFERENCE] + [g for g in STRUCTURES if g != ALIGN_REFERENCE]


def load_hits():
    """gene -> {'all': [...positions...], 'unmodelled': [...]} plus labels."""
    rows = list(csv.DictReader(open(HITS_CSV)))
    unmodelled = set()
    if os.path.exists(os.path.join(HERE, "unmodelled_positions.txt")):
        for line in open(os.path.join(HERE, "unmodelled_positions.txt")):
            line = line.strip()
            if line and not line.startswith("#"):
                g, p = line.split()
                unmodelled.add((g, int(p)))
    hits = collections.defaultdict(list)
    for r in rows:
        hits[r["gene"]].append(r)
    return hits, unmodelled


def spec(model, chains):
    return f"#{model}/" + ",".join(chains)


def model_block(mid, gene, hits, unmodelled, lines):
    """Emit the load + representation commands for one CDK."""
    cfg = STRUCTURES[gene]
    kc, cyc, other = cfg["kinase_chain"], cfg["cyclin_chains"], cfg["other_chains"]
    used = [kc] + cyc + other
    lo, hi = cfg["kd"]
    color = CDK_COLORS[gene]

    lines.append(f"\n# ---------------------------------------------------------------- {gene}")
    lines.append(f"# {cfg['pdb']}: {cfg['note']}")
    lines.append(f"open structures/{cfg['file']} name {gene}")
    lines.append(f"hide #{mid} target acsb")
    lines.append(f"hide #{mid} models")   # re-shown at the end, keeps intermediate steps cheap

    # ---- kinase cartoon
    dr = cfg.get("display_range")
    if dr:
        lines.append(f"# only residues {dr[0]}-{dr[1]} are confidently modelled")
        lines.append(f"show #{mid}/{kc}:{dr[0]}-{dr[1]} cartoon")
    else:
        lines.append(f"show #{mid}/{kc} cartoon")
    lines.append(f"color #{mid}/{kc} {color} target c")
    # everything outside the catalytic core fades back (cartoon_transparency 0.8
    # in the old PyMOL script)
    lines.append(f"transparency #{mid}/{kc} & ~:{lo}-{hi} {NONCATALYTIC_TRANSPARENCY} cartoons")

    # ---- cyclin surface
    if cyc:
        cs = spec(mid, cyc)
        lines.append(f"# {cfg['cyclin_name']} as surface")
        lines.append(f"show {cs} atoms")
        lines.append(f"surface {cs}")
        lines.append(f"hide {cs} atoms")
        lines.append(f"color {cs} {CYCLIN_COLOR} target s")
    else:
        lines.append(f"# no cyclin partner available for {gene}")

    # ---- ligand
    if cfg["ligand"]:
        ls = f"#{mid}/{kc} & :{cfg['ligand']}"
        # target `ab` matters: the blanket `hide ... target acsb` above turns
        # bond display off, and `show ... atoms` alone would leave the ligand
        # drawn as unbonded dots
        lines.append(f"show {ls} target ab")
        lines.append(f"style {ls} stick")
        lines.append(f"color {ls} #909090 target ab")
        lines.append(f"color {ls} byhetero")
        if gene == "CDK2":   # keep the catalytic Mg that comes with the ATP
            lines.append(f"show #{mid}/{kc} & :MG atoms")
            lines.append(f"style #{mid}/{kc} & :MG sphere")
            lines.append(f"size #{mid}/{kc} & :MG atomRadius 0.8")

    # ---- resistance hits
    rows = sorted(hits.get(gene, []), key=lambda r: int(r["position"]))
    def shown(r):
        pos = int(r["position"])
        if (gene, pos) in unmodelled:
            return False
        return not dr or dr[0] <= pos <= dr[1]

    present = [r for r in rows if shown(r)]
    absent = [r for r in rows if not shown(r)]
    if present:
        resis = ",".join(str(p) for p in sorted({int(r["position"]) for r in present}))
        hs = f"#{mid}/{kc}:{resis}"
        lines.append(f"# {len(present)} resistance positions: "
                     + ", ".join(sorted({r["variants_at_position"] for r in present}))[:900])
        lines.append(f"show {hs} atoms")
        # ChimeraX hides backbone atoms under a cartoon by default, which makes
        # glycine hits (no side chain past CA) disappear entirely -- only the
        # ribbon gets recolored.  Turn that off for the hit residues.
        lines.append(f"cartoon {hs} suppressBackboneDisplay false")
        lines.append(f"style {hs} {HIT_STYLE}")
        lines.append(f"color {hs} {HIT_COLOR} target ac")
        lines.append(f"name {gene}_hits {hs}")
    if absent:
        lines.append("# not shown - disordered / not modelled in this structure: "
                     + ", ".join(f"{r['variant']}" for r in absent))

    # ---- other chains stay available but hidden
    if other:
        lines.append(f"# hidden partners in this file: {', '.join(other)}")
        lines.append(f"color {spec(mid, other)} {OTHER_CHAIN_COLOR}")
    unused = [c for c in _chains_in(gene) if c not in used]
    if unused:
        lines.append(f"# extra copies in the asymmetric unit, left hidden: {', '.join(unused)}")

    lines.append(f"name {gene}_kinase #{mid}/{kc}")
    if cyc:
        lines.append(f"name {gene}_cyclin {spec(mid, cyc)}")


# chains present in each file, filled in by validate_structures.py's output
_CHAINS = {
    "CDK7": list("HIJ"), "CDK8": list("AB"), "CDK19": ["A"], "CDK9": list("AB"),
    "CDK12": list("ABCDEFGHI"), "CDK13": list("ABCD"), "CDK2": list("ABCDEF"),
    "CDK4": list("AB"), "CDK6": list("AB"),
}


def _chains_in(gene):
    return _CHAINS[gene]


def borrowed_block(mid, gene, host_mid, lines):
    """Superpose a donor CDK onto `gene` and show the donor's cyclin."""
    cfg = STRUCTURES[gene]
    b = cfg["borrowed_cyclin"]
    donor = STRUCTURES[b["donor"]]
    lines.append(f"\n# ------------------------------- borrowed cyclin for {gene}")
    lines.append(f"# {gene} has no structure with a cyclin bound.  {b['donor']} "
                 f"({donor['pdb']}) is superposed onto it and its")
    lines.append(f"# {donor['cyclin_name']} shown at {b['transparency']}% transparency, "
                 "to make clear the cyclin is inferred rather than observed.")
    lines.append(f"open structures/{donor['file']} name {gene}_borrowed_{donor['cyclin_name'].replace(' ', '')}")
    lines.append(f"hide #{mid} target acsb")
    lines.append(f"matchmaker #{mid}/{donor['kinase_chain']} to #{host_mid}/{cfg['kinase_chain']}")
    cs = spec(mid, donor["cyclin_chains"])
    lines.append(f"show {cs} atoms")
    lines.append(f"surface {cs}")
    lines.append(f"hide {cs} atoms")
    lines.append(f"color {cs} {CYCLIN_COLOR} target s")
    lines.append(f"transparency {cs} {b['transparency']} target s")
    lines.append(f"name {gene}_borrowed_cyclin {cs}")


def header(title, extra=""):
    return [f"# {'=' * 74}", f"# {title}", "# Generated by generate_cxc.py -- edit that script, not this file.",
            f"# {extra}" if extra else "# ", f"# {'=' * 74}", "close session", "cd " + HERE, ""]


def footer(lines, focus=None):
    lines.append("\n# ---------------------------------------------------- look and camera")
    lines.extend(RENDER_PREAMBLE.strip().splitlines())
    lines.append("show #* models")
    lines.append("# mmCIF chain breaks otherwise draw dashed pseudobonds labelled")
    lines.append("# \"N residues\", which land on top of the kinase in the render")
    lines.append("label delete")
    lines.append("hide #* target p")
    lines.append("# 8BUG (CDK12) ships riding hydrogens; the crystal structures do")
    lines.append("# not, so drop them everywhere to keep the panels consistent")
    lines.append("hide @@element=H target ab")
    lines.append(f"view matrix camera {CANONICAL_VIEW}")
    lines.append("# uniform pull-back so no panel clips at the frame edge")
    lines.append(f"zoom {ZOOM_OUT}")
    lines.append("# `view` alone re-frames on everything shown; the matrix above")
    lines.append("# reproduces the orientation used in viz2_cdk2.cxs.")
    if focus:
        lines.append(f"# zoom on the catalytic cleft:  view {focus}")
    lines.append("")
    lines.append("# save figures/<name>.png width 3000 supersample 3 transparentBackground true")
    lines.append("# save <name>.cxs")


def write_structure_table():
    """Write the structure/chain provenance table next to the rendered figures."""
    fig = os.path.normpath(os.path.join(HERE, "..", "..", "figures", "chimerax"))
    os.makedirs(fig, exist_ok=True)

    rows = []
    for gene in ORDER:
        c = STRUCTURES[gene]
        cyc = c.get("borrowed_cyclin")
        if cyc:
            donor = STRUCTURES[cyc["donor"]]
            cyclin = (f"{donor['cyclin_name']} (borrowed from {donor['pdb']} chain "
                      f"{donor['cyclin_chains'][0]}, {cyc['transparency']}% transparent)")
            cyclin_chain = f"{donor['pdb']}:{donor['cyclin_chains'][0]}"
        elif c["cyclin_chains"]:
            cyclin = c["cyclin_name"]
            cyclin_chain = ",".join(c["cyclin_chains"])
        else:
            cyclin, cyclin_chain = "none", ""
        rows.append(dict(
            gene=gene, pdb=c["pdb"], cdk_chain=c["kinase_chain"],
            cyclin=cyclin, cyclin_chain=cyclin_chain,
            ligand=(f"{LIGAND_NAMES.get(c['ligand'], c['ligand'])} ({c['ligand']})"
                    if c["ligand"] else "none"),
            hidden_chains=",".join(c["other_chains"]),
            catalytic_domain=f"{c['kd'][0]}-{c['kd'][1]}",
            displayed_range=(f"{c['display_range'][0]}-{c['display_range'][1]}"
                             if c.get("display_range") else "full chain"),
            cartoon_color=CDK_COLORS[gene],
            note=" ".join(c["note"].split()),
        ))

    fields = ["gene", "pdb", "cdk_chain", "cyclin", "cyclin_chain", "ligand",
              "hidden_chains", "catalytic_domain", "displayed_range",
              "cartoon_color", "note"]
    with open(os.path.join(fig, "structure_table.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    md = ["# Structures used for the CDK resistance-variant figures",
          "",
          "Generated by `notebooks/11_ChimeraX/generate_cxc.py` from `cdk_config.py`.",
          "One structure per gene, preferring an experimental structure with the",
          "cyclin partner bound.",
          "",
          "| gene | PDB | CDK chain | cyclin (surface) | ligand | catalytic domain | notes |",
          "|------|-----|-----------|------------------|--------|------------------|-------|"]
    for r in rows:
        md.append(f"| {r['gene']} | **{r['pdb']}** | {r['cdk_chain']} | {r['cyclin']} "
                  f"| {r['ligand']} | {r['catalytic_domain']} | {r['note']} |")
    md += ["",
           "Chains present in the file but hidden in the figures (MAT1, DDB1, the CDK2",
           "substrate peptide, extra copies in the asymmetric unit) are listed in",
           "`structure_table.csv`.",
           "",
           "All nine are superposed onto CDK2 (1QMZ chain A) with `matchmaker`;",
           "pruned RMSDs are 0.81-1.12 A.",
           ""]
    open(os.path.join(fig, "structure_table.md"), "w").write("\n".join(md))
    print(f"wrote structure_table.csv / .md to {fig}")


def main():
    global OUT, HITS_CSV, ORDER
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hits", help="hit table in hit_residues.csv format "
                                   "(default: hit_residues.csv next to this script)")
    ap.add_argument("--out", help="output directory for the .cxc files (default: ./cxc)")
    ap.add_argument("--genes", help="comma-separated subset of CDKs to emit "
                                    "(default: all nine)")
    a = ap.parse_args()
    if a.hits:
        HITS_CSV = os.path.abspath(a.hits)
    if a.out:
        OUT = os.path.abspath(a.out)
    if a.genes:
        keep = [g.strip() for g in a.genes.split(",")]
        unknown = [g for g in keep if g not in STRUCTURES]
        if unknown:
            sys.exit(f"unknown gene(s): {unknown}")
        # the alignment reference is always loaded, even if it is not in `keep`
        ORDER = [ALIGN_REFERENCE] + [g for g in keep if g != ALIGN_REFERENCE]
    print(f"hits:   {HITS_CSV}")
    print(f"output: {OUT}")
    os.makedirs(OUT, exist_ok=True)
    hits, unmodelled = load_hits()

    # ---------------------------------------------------------- combined figure
    lines = header("All CDK-family resistance variants, superposed on CDK2 (1QMZ chain A)",
                   "Model numbering: " + ", ".join(f"#{i+1}={g}" for i, g in enumerate(ORDER)))
    for i, gene in enumerate(ORDER):
        model_block(i + 1, gene, hits, unmodelled, lines)
    lines.append("\n# ------------------------------------------------ structural alignment")
    lines.append(f"# every kinase domain is superposed onto {ALIGN_REFERENCE} "
                 f"({STRUCTURES[ALIGN_REFERENCE]['pdb']} chain "
                 f"{STRUCTURES[ALIGN_REFERENCE]['kinase_chain']}); the reference never moves,")
    lines.append("# so cyclins and ligands are carried into a common frame too.")
    ref_spec = f"#1/{STRUCTURES[ALIGN_REFERENCE]['kinase_chain']}"
    for i, gene in enumerate(ORDER):
        if gene == ALIGN_REFERENCE:
            continue
        lines.append(f"matchmaker #{i+1}/{STRUCTURES[gene]['kinase_chain']} to {ref_spec}")
    next_id = len(ORDER) + 1
    for i, gene in enumerate(ORDER):
        if "borrowed_cyclin" in STRUCTURES[gene]:
            borrowed_block(next_id, gene, i + 1, lines)
            next_id += 1
    footer(lines, focus="CDK2_hits")
    combined = "\n".join(lines) + "\n"
    open(os.path.join(OUT, "00_all_cdks_aligned.cxc"), "w").write(combined)

    # same figure with the cyclin surfaces off -- with nine kinases superposed the
    # surfaces stack up and hide the far side (cf. top_hits_MSA_no_surface.png)
    nosurf = combined.replace("# ---------------------------------------------------- look and camera",
                              "hide #* target s\n"
                              "# ---------------------------------------------------- look and camera")
    nosurf = nosurf.replace("All CDK-family resistance variants, superposed on CDK2 (1QMZ chain A)",
                            "All CDK-family resistance variants, superposed, cyclin surfaces hidden")
    open(os.path.join(OUT, "01_all_cdks_no_surface.cxc"), "w").write(nosurf)

    # -------------------------------------------------------- single-CDK panels
    for gene in ORDER:
        cfg = STRUCTURES[gene]
        lines = header(f"{gene} resistance variants ({cfg['pdb']})", cfg["note"].replace("\n", " "))
        if gene == ALIGN_REFERENCE:
            model_block(1, gene, hits, unmodelled, lines)
        else:
            # the reference is loaded but hidden, purely to put this panel in the
            # same frame as every other panel
            lines.append("# reference loaded hidden so this panel shares the common orientation")
            ref = STRUCTURES[ALIGN_REFERENCE]
            lines.append(f"open structures/{ref['file']} name {ALIGN_REFERENCE}_ref")
            lines.append("hide #1 target acsb")
            lines.append("hide #1 models")
            model_block(2, gene, hits, unmodelled, lines)
            lines.append(f"\nmatchmaker #2/{cfg['kinase_chain']} to #1/{ref['kinase_chain']}")
        host = 1 if gene == ALIGN_REFERENCE else 2
        if "borrowed_cyclin" in cfg:
            borrowed_block(host + 1, gene, host, lines)
        footer(lines, focus=f"{gene}_hits")
        # keep the reference hidden in single panels
        shown = "#*" if gene == ALIGN_REFERENCE else (
            "#2,3" if "borrowed_cyclin" in cfg else "#2")
        out = "\n".join(lines).replace("show #* models", f"show {shown} models")
        open(os.path.join(OUT, f"{gene}.cxc"), "w").write(out + "\n")

    write_structure_table()

    # ------------------------------------------------------------- batch render
    tag = os.path.basename(OUT)                     # "cxc" for the default run
    fig = os.path.normpath(os.path.join(HERE, "..", "..", "figures", "chimerax"))
    if tag != "cxc":                                # keep alternate hit sets out of the
        fig = os.path.join(fig, tag)                # original figures/chimerax/
    os.makedirs(fig, exist_ok=True)
    rl = header(f"Render every panel to {os.path.relpath(fig, os.path.join(HERE, '..', '..'))} "
                "at publication size",
                "ChimeraX needs a real GPU context for `save`, so run this from the "
                f"ChimeraX GUI\n#   (or `ChimeraX --exit --cmd \"open {tag}/render_all.cxc\"`), "
                "not --nogui.")
    rl = rl[:-2]   # drop the trailing `close session` / blank from the header
    rl.append(f"cd {HERE}")
    rl.append("# fixed aspect ratio so every panel comes out the same shape")
    rl.append(f"windowsize {WINDOW_SIZE[0]} {WINDOW_SIZE[1]}")
    for name in ["00_all_cdks_aligned", "01_all_cdks_no_surface"] + ORDER:
        rl.append(f"\nopen {tag}/{name}.cxc")
        rl.append(f"save {fig}/{name}.png width 3000 supersample 3 "
                  "transparentBackground true")
    rl.append("\nclose session")
    open(os.path.join(OUT, "render_all.cxc"), "w").write("\n".join(rl) + "\n")

    print(f"wrote {len(ORDER)+3} scripts to {OUT}")
    for f in sorted(os.listdir(OUT)):
        print("   ", f)


if __name__ == "__main__":
    main()
