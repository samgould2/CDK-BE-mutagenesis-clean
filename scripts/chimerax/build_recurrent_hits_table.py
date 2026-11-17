"""
Build a hit_residues.csv-format table from the cross-CDK resistance analysis.

Companion to build_hits_table.py, which reads the *_filtered_top_only tables.
This one reads the per-call table written by
notebooks/09_analysis_and_viz/CDK_resistance_patterns.ipynb, so the structural
panels can show exactly the variant set that analysis calls resistant --
including the option to show only the positions that recur across CDKs.

Input : data/processed/lfc_fdr_tables/resistance_patterns/<set>/resistance_calls_canonical.csv
Output: hit_residues_<set>.csv  (one row per gene / residue / editor)

    python3 build_recurrent_hits_table.py --set core7
    python3 build_recurrent_hits_table.py --set all9 --min-cdks 2

Feed the result to generate_cxc.py:

    python3 generate_cxc.py --hits hit_residues_core7.csv --out cxc_core7 \
        --genes CDK7,CDK9,CDK12,CDK13,CDK2,CDK4,CDK6
"""
import argparse, collections, csv, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SUBST = re.compile(r"^([A-Z])(\d+)([A-Z]|\*)$")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", default="core7",
                    help="analysis set directory under resistance_patterns/ (default: core7)")
    ap.add_argument("--calls", help="explicit path to resistance_calls_canonical.csv")
    ap.add_argument("--min-cdks", type=int, default=0,
                    help="keep only residues whose alignment column carries resistance in at "
                         "least this many CDKs (0 = keep every resistance residue)")
    ap.add_argument("-o", "--out", help="output csv (default: hit_residues_<set>.csv)")
    a = ap.parse_args()

    calls_path = a.calls or os.path.join(
        ROOT, "data", "processed", "lfc_fdr_tables", "resistance_patterns", a.set,
        "resistance_calls_canonical.csv")
    if not os.path.exists(calls_path):
        sys.exit(f"no calls table at {calls_path}\n"
                 f"run CDK_resistance_patterns.ipynb first")
    rows = list(csv.DictReader(open(calls_path)))
    print(f"read {len(rows)} enriched calls from {os.path.relpath(calls_path, ROOT)}")

    # recurrence is counted per alignment column, which is what makes a site
    # comparable across CDKs; a residue is kept if its column clears the cutoff
    by_col = collections.defaultdict(set)
    for r in rows:
        by_col[r["aln_col"]].add(r["Gene"])
    if a.min_cdks:
        rows = [r for r in rows if len(by_col[r["aln_col"]]) >= a.min_cdks]
        print(f"  {len(rows)} calls remain at >= {a.min_cdks} CDKs per alignment column")

    out = {}
    seen = collections.defaultdict(set)
    for r in rows:
        m = SUBST.match(r["SAV"])
        if not m:
            continue
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        key = (r["Gene"], pos, r["Editor"])
        seen[key].add(r["SAV"])
        lfc, fdr = float(r["LFC"]), float(r["FDR"])
        prev = out.get(key)
        if prev is None or lfc > prev["max_LFC"]:
            out[key] = dict(gene=r["Gene"], editor=r["Editor"], position=pos, wt=wt, mut=mut,
                            variant=r["SAV"],
                            kind="nonsense" if mut == "*" else "missense",
                            max_LFC=lfc, min_FDR=fdr, gRNA_id=r["gRNA_id"],
                            HGVSp=r["SAV"], screen=r["Screen"])
        else:
            prev["min_FDR"] = min(prev["min_FDR"], fdr)
    for key, r in out.items():
        r["variants_at_position"] = ";".join(sorted(seen[key]))

    dest = a.out or os.path.join(HERE, f"hit_residues_{a.set}.csv")
    fields = ["gene", "editor", "position", "wt", "mut", "variant",
              "variants_at_position", "kind", "max_LFC", "min_FDR",
              "gRNA_id", "HGVSp", "screen"]
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(out.values(), key=lambda r: (r["gene"], r["position"], r["editor"])))
    print(f"wrote {dest}: {len(out)} gene/position/editor rows")
    per = collections.Counter((r["gene"], r["editor"]) for r in out.values())
    for gene in sorted({g for g, _ in per}):
        print(f"  {gene:6s} CBE={per[(gene, 'CBE')]:3d}  ABE={per[(gene, 'ABE')]:3d}")


if __name__ == "__main__":
    main()
