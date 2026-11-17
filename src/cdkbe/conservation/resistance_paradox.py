"""
Why resistance alleles sit at residues evolution never changes.

The apparent paradox: the strongest resistance variants are at positions that are
invariant across ~1000 species, where a mutation "should" break the protein. The
resolution is that conservation predicts *functional impact*, not its direction,
and the two outcomes land on different invariant residues of the same pocket.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss
from Bio.Align import substitution_matrices

from . import screen_conservation_integration as sci
from .. import paths

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
plt.rc("font", family="Helvetica")

B62 = substitution_matrices.load("BLOSUM62")
CATALYTIC = {17: "VAIK Lys", 24: "αC-Glu", 68: "HRD-His", 70: "HRD-Asp",
             75: "cat-loop Asn", 81: "DFG-Asp", 82: "DFG-Phe", 83: "DFG-Gly"}


def main():
    d = pd.read_csv(f"{CONS}/screen_conservation_merged.csv")
    d = d[d.consequence == "missense"].copy()
    d["blosum"] = [B62[r, a] if (r, a) in B62 else np.nan
                   for r, a in zip(d.ref_aa, d.alt_aa)]
    h = pd.read_excel(f"{CONS}/TABLE_conservation_confidence_tiers.xlsx",
                      "Highest_confidence_variants")
    res = h[h.outcome == "resistance"]

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.2),
                             gridspec_kw=dict(wspace=0.42, left=0.05, right=0.985,
                                              top=0.80, bottom=0.20))
    a0, a1, a2, a3 = axes

    # --- a: invariant sites give BOTH outcomes
    order = ["Tolerant in evolution", "Intermediate", "CDK-specific, fixed",
             "Universally immutable"]
    sub = d[d.site_class.isin(order)]
    x = np.arange(len(order))
    for off, col, key, lab in [(-0.2, "#2166AC", "fitness_hit", "fitness"),
                               (0.2, "#B2182B", "resistance_hit", "resistance")]:
        r = sub.groupby("site_class")[key].mean().reindex(order) * 100
        a0.bar(x + off, r, 0.4, color=col, label=lab)
    n = sub.groupby("site_class").size().reindex(order)
    a0.set_xticks(x)
    a0.set_xticklabels([f"{o}\n(n={n[o]})" for o in order], fontsize=7.5,
                       rotation=20, ha="right")
    a0.set_ylabel("Hit rate (%)", fontsize=10)
    a0.legend(fontsize=8.5, frameon=False)
    a0.set_title("1. Invariant sites give BOTH outcomes\n"
                 "conservation predicts impact, not direction",
                 fontsize=10.5, weight="bold")
    sns.despine(ax=a0)
    a0.text(-0.20, 1.06, "a", transform=a0.transAxes, fontsize=15, weight="bold")

    # --- b: resistance alleles are not depleted
    rv = d.assign(variant=d.gene + "_" + d.HGVSp)
    rv = rv[rv.variant.isin(res.gene + "_" + res.HGVSp)]
    g = rv.groupby(["gene", "HGVSp"])["min_fitness_LFC"].median().sort_values()
    colors = ["#B2182B" if v > -1 else "#777777" for v in g]
    a1.barh(range(len(g)), g.values, color=colors)
    a1.axvline(-1, ls="--", color="black", lw=0.9)
    a1.text(-1.06, len(g) - 0.4, "fitness-hit\nthreshold", fontsize=7,
            ha="right", va="top")
    a1.set_yticks(range(len(g)))
    a1.set_yticklabels([f"{a} {b}" for a, b in g.index], fontsize=8.5)
    a1.set_xlabel("Fitness LFC (most depleted contrast)", fontsize=9.5)
    a1.set_title("2. Those alleles keep the protein working\n"
                 "7/8 are not depleted without drug", fontsize=10.5, weight="bold")
    sns.despine(ax=a1)
    a1.text(-0.20, 1.06, "b", transform=a1.transAxes, fontsize=15, weight="bold")

    # --- c: different invariant residues of the same pocket
    p = d[d.KLIFS_position.notna()].copy()
    regs = ["I", "II", "III", "g.l", "aC", "b.l", "V", "GK", "hinge",
            "aD", "aE", "VI", "c.l", "VII", "xDFG"]
    tab = []
    for key, lab in [("fitness_hit", "fitness"), ("resistance_hit", "resistance")]:
        v = p[p[key]].region.value_counts()
        tab.append((v.reindex(regs).fillna(0) / max(v.sum(), 1) * 100).rename(lab))
    tab = pd.concat(tab, axis=1)
    xx = np.arange(len(regs))
    a2.bar(xx - 0.2, tab["fitness"], 0.4, color="#2166AC", label="fitness")
    a2.bar(xx + 0.2, tab["resistance"], 0.4, color="#B2182B", label="resistance")
    a2.set_xticks(xx); a2.set_xticklabels(regs, fontsize=7.5, rotation=60)
    a2.set_ylabel("% of pocket hits in region", fontsize=9.5)
    a2.legend(fontsize=8.5, frameon=False)
    a2.set_title("3. But at DIFFERENT invariant residues\n"
                 "resistance: b1/b7 pocket lining;  fitness: catalytic loop",
                 fontsize=10.5, weight="bold")
    sns.despine(ax=a2)
    a2.text(-0.20, 1.06, "c", transform=a2.transAxes, fontsize=15, weight="bold")

    # --- d: it is NOT explained by milder substitutions
    pl = pd.concat([
        pd.DataFrame({"b": d.loc[d.resistance_hit, "blosum"], "g": "resistance\nhits"}),
        pd.DataFrame({"b": d.loc[d.fitness_hit, "blosum"], "g": "fitness\nhits"}),
        pd.DataFrame({"b": d.loc[~d.resistance_hit & ~d.fitness_hit, "blosum"],
                      "g": "non-hits"})]).dropna()
    o2 = ["non-hits", "fitness\nhits", "resistance\nhits"]
    sns.boxplot(data=pl, x="g", y="b", ax=a3, order=o2, hue="g", hue_order=o2,
                palette=["#CCCCCC", "#2166AC", "#B2182B"], legend=False, fliersize=1)
    pv = ss.mannwhitneyu(d.loc[d.resistance_hit, "blosum"].dropna(),
                         d.loc[d.fitness_hit, "blosum"].dropna())[1]
    a3.set_xlabel(""); a3.set_ylabel("BLOSUM62 substitution score", fontsize=9.5)
    a3.set_xticklabels(o2, fontsize=8.5)
    a3.set_title(f"4. NOT because the swaps are milder\n"
                 f"resistance vs fitness p = {pv:.2f} (n.s.)",
                 fontsize=10.5, weight="bold")
    sns.despine(ax=a3)
    a3.text(-0.20, 1.06, "d", transform=a3.transAxes, fontsize=15, weight="bold")

    fig.suptitle("Why the most conserved residues yield resistance alleles",
                 fontsize=13.5, weight="bold", y=0.955)
    fig.savefig(f"{FIG}/resistance_conservation_paradox.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/resistance_conservation_paradox.png", dpi=280, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/resistance_conservation_paradox.{{pdf,png}}")
    return d, h


if __name__ == "__main__":
    main()
