"""
Parse the top-hit LFC tables into a tidy, per-residue table for structural
visualization.

Input : data/processed/lfc_fdr_tables/{ABE,CBE}_master_LFC_table_filtered_top_only.csv
Output: hit_residues.csv  (one row per gene / residue / editor)

HGVSp entries can encode several substitutions from one guide
(e.g. "M816I_D817N", "D185N_W187*").  Each substitution becomes its own row,
inheriting that guide's max_LFC / min_FDR, so a multi-edit guide highlights
every residue it changes.
"""
import csv, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
TBL = os.path.join(HERE, "..", "..", "data", "processed", "lfc_fdr_tables")
SUBST = re.compile(r"^([A-Z])(\d+)([A-Z]|\*)$")


def parse_hgvsp(hgvsp):
    """'M816I_D817N' -> [('M',816,'I'), ('D',817,'N')]"""
    out = []
    for part in hgvsp.split("_"):
        m = SUBST.match(part.strip())
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def main():
    rows = []
    for editor in ("CBE", "ABE"):
        path = os.path.join(TBL, f"{editor}_master_LFC_table_filtered_top_only.csv")
        with open(path) as fh:
            for r in csv.DictReader(fh):
                for wt, pos, mut in parse_hgvsp(r["HGVSp"]):
                    rows.append(dict(
                        gene=r["gene"], editor=editor, position=pos, wt=wt, mut=mut,
                        variant=f"{wt}{pos}{mut}",
                        kind="nonsense" if mut == "*" else "missense",
                        max_LFC=float(r["max_LFC"]), min_FDR=float(r["min_FDR"]),
                        gRNA_id=r["gRNA_id"], HGVSp=r["HGVSp"], screen=r["Screen_ID"],
                    ))

    # collapse to one row per gene/position/editor, keeping the strongest guide
    best = {}
    for r in rows:
        key = (r["gene"], r["position"], r["editor"])
        if key not in best or r["max_LFC"] > best[key]["max_LFC"]:
            best[key] = r
    # note every variant seen at that position/editor
    seen = collections.defaultdict(set)
    for r in rows:
        seen[(r["gene"], r["position"], r["editor"])].add(r["variant"])
    for key, r in best.items():
        r["variants_at_position"] = ";".join(sorted(seen[key]))

    out = sorted(best.values(), key=lambda r: (r["gene"], r["position"], r["editor"]))
    fields = ["gene", "editor", "position", "wt", "mut", "variant",
              "variants_at_position", "kind", "max_LFC", "min_FDR",
              "gRNA_id", "HGVSp", "screen"]
    dest = os.path.join(HERE, "hit_residues.csv")
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out)

    print(f"wrote {dest}: {len(out)} gene/position/editor rows "
          f"from {len(rows)} substitutions")
    per = collections.Counter((r["gene"], r["editor"]) for r in out)
    for gene in sorted({g for g, _ in per}):
        print(f"  {gene:6s} CBE={per[(gene,'CBE')]:2d}  ABE={per[(gene,'ABE')]:2d}")


if __name__ == "__main__":
    main()
