"""
Per-residue hit table for the CDK2/4/6 active-site panels.

Input : data/processed/lfc_fdr_tables/CDK2_4_6_hits/CDK2_4_6_resistance_variants.csv
        (written by notebooks/09_analysis_and_viz/logo_plots_CDK2_4_6_canonical.ipynb)
Output: cdk2_4_6_hit_residues.csv  -- one row per gene / residue

Same shape as cdk12_13_ddb1_hit_residues.csv, which serves the CDK12/13 panels.

Which hits: with TOP_N = None

# Optional variant tag + stricter FDR.  Set CDK246_TAG / CDK246_FDR to build a second hit
# set without overwriting the first.  The source table is already filtered to FDR < 0.1 &
# LFC >= 1 in >=1 condition; CDK246_FDR re-applies the same rule at a tighter FDR.
TAG = os.environ.get("CDK246_TAG", "")
FDR_CUTOFF = os.environ.get("CDK246_FDR")
FDR_CUTOFF = float(FDR_CUTOFF) if FDR_CUTOFF else None
CONDITIONS = ["RIBO", "ABEMA", "ATIRMO", "PALBO", "TAGTO", "INX"]
LFC_MIN = 1.0, **every** resistance hit for that gene -- i.e. every
guide in the notebook's table, which is already filtered to FDR < 0.1 & LFC >= 1 in
at least one of the six conditions, base_RAW >= 100 and sensor WT < 90%.  Set TOP_N
to an integer to keep only the top N guides by `active_median_LFC` instead.

A compound guide (e.g. CDK4 I51V_T53A) contributes a row per substituted residue, so
the position count exceeds the guide count.

Also recorded per position:
  drug_class  CDK4/6i only | CDK2i only | both  -- which inhibitor class(es) the
              strongest guide at that position is resistant to.  The panels color
              by this, so an active-site hit that is CDK2i-specific is visually
              distinct from one that is pan-class.
"""
import csv, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "data", "processed", "lfc_fdr_tables",
                   "CDK2_4_6_hits", "CDK2_4_6_resistance_variants.csv")
OUT = os.path.join(HERE, "cdk2_4_6_hit_residues%s.csv" % os.environ.get("CDK246_TAG",""))

GENES = ["CDK2", "CDK4", "CDK6"]
# None = every hit in the screen for that gene.  Set an integer to keep only the top N
# guides by active_median_LFC (it was 10); `g[:None]` slices the whole list, so the same
# code path serves both.
TOP_N = None

# Optional variant tag + stricter FDR.  Set CDK246_TAG / CDK246_FDR to build a second hit
# set without overwriting the first.  The source table is already filtered to FDR < 0.1 &
# LFC >= 1 in >=1 condition; CDK246_FDR re-applies the same rule at a tighter FDR.
TAG = os.environ.get("CDK246_TAG", "")
FDR_CUTOFF = os.environ.get("CDK246_FDR")
FDR_CUTOFF = float(FDR_CUTOFF) if FDR_CUTOFF else None
CONDITIONS = ["RIBO", "ABEMA", "ATIRMO", "PALBO", "TAGTO", "INX"]
LFC_MIN = 1.0
SUBST = re.compile(r"^([A-Z])(\d+)([A-Z]|\*)$")


def parse_hgvsp(h):
    """'I51V_T53A' -> [('I',51,'V'), ('T',53,'A')]"""
    out = []
    for part in h.split("_"):
        m = SUBST.match(part.strip())
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def main():
    allrows = list(csv.DictReader(open(SRC)))
    rows = []
    if FDR_CUTOFF is not None:
        keep = [r for r in allrows
                if any(float(r[f"FDR_{c}"]) < FDR_CUTOFF and float(r[f"LFC_{c}"]) >= LFC_MIN
                       for c in CONDITIONS)]
        print(f"FDR < {FDR_CUTOFF} & LFC >= {LFC_MIN} in >=1 condition: "
              f"{len(keep)} of {len(allrows)} guides kept")
        allrows = keep

    for gene in GENES:
        g = [r for r in allrows if r["Gene"] == gene]
        g.sort(key=lambda r: -float(r["active_median_LFC"]))
        for rank, r in enumerate(g[:TOP_N], start=1):
            for wt, pos, mut in parse_hgvsp(r["Top_HGVSp"]):
                rows.append(dict(
                    gene=gene, editor=r["Editor"], position=pos, wt=wt, mut=mut,
                    variant=f"{wt}{pos}{mut}",
                    kind="nonsense" if mut == "*" else "missense",
                    rank=rank,
                    active_median_LFC=float(r["active_median_LFC"]),
                    max_LFC=float(r["max_LFC"]),
                    n_conditions=int(r["n_conditions"]),
                    drug_class=r["drug_class"],
                    gRNA_id=r["gRNA_id"], HGVSp=r["Top_HGVSp"]))

    # one row per gene/position: keep the strongest guide, note every variant seen
    best, seen, classes = {}, collections.defaultdict(set), collections.defaultdict(set)
    for r in rows:
        key = (r["gene"], r["position"])
        if key not in best or r["active_median_LFC"] > best[key]["active_median_LFC"]:
            best[key] = r
        seen[key].add(r["variant"])
        classes[key].add(r["drug_class"])
    for key, r in best.items():
        r["variants"] = ";".join(sorted(seen[key]))
        # a position carrying guides of both classes is reported as 'both'
        r["drug_class"] = ("both" if len(classes[key]) > 1 or "both" in classes[key]
                           else classes[key].pop())

    out = sorted(best.values(), key=lambda r: (GENES.index(r["gene"]), r["position"]))
    fields = ["gene", "position", "wt", "mut", "variant", "variants", "kind", "rank",
              "active_median_LFC", "max_LFC", "n_conditions", "drug_class",
              "editor", "gRNA_id", "HGVSp"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    print(f"wrote {OUT}: {len(out)} gene/position rows from top {TOP_N} guides per gene")
    for gene in GENES:
        sub = [r for r in out if r["gene"] == gene]
        cls = collections.Counter(r["drug_class"] for r in sub)
        print(f"  {gene}: {len(sub):2d} positions  {dict(cls)}")
        print("     " + ", ".join(f'{r["wt"]}{r["position"]}({r["editor"]})' for r in sub))


if __name__ == "__main__":
    main()
