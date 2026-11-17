"""Figures + tables for deep ortholog conservation of the nine screened CDKs."""
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from adjustText import adjust_text

from . import cdk_ortholog_conservation as co
from . import kinome_conservation as kc
from .. import paths

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
plt.rc("font", family="Helvetica")


CDKS = list(co.CDKS)

CLADE_ORDER = ["Mammalia", "Aves", "Amphibia", "Actinopterygii", "Chordata",
               "Arthropoda", "Nematoda", "Other Metazoa", "Fungi",
               "Viridiplantae", "Other Eukaryota", "Unclassified"]
CLADE_COLORS = dict(zip(CLADE_ORDER, sns.color_palette("Spectral", len(CLADE_ORDER))))


def load_all():
    cons = pd.read_csv(f"{CONS}/ortholog_conservation_all_CDKs.csv")
    klifs = pd.read_csv(os.path.join(KLIFS, "KLIFS_compiled.csv"))
    klifs["KLIFS_position"] = klifs["KLIFS_idx"].str.split(":").str[1].astype(int)
    kin = pd.read_csv(f"{CONS}/TABLE_kinome_KLIFS_conservation.csv")

    # long form: one row per (CDK, KLIFS position)
    rows = []
    for _, r in klifs.iterrows():
        for c in CDKS:
            if pd.notna(r[c]):
                rows.append({"gene": c, "KLIFS_position": int(r["KLIFS_position"]),
                             "human_position": int(r[c]), "klifs_aa": r[f"{c}_seq"]})
    kmap = pd.DataFrame(rows)
    merged = (kmap.merge(cons, on=["gene", "human_position"], how="left")
                  .merge(kin[["KLIFS_position", "KLIFS_label", "region", "consensus_aa",
                              "conservation", "mutability_class"]]
                         .rename(columns={"conservation": "kinome_conservation",
                                          "consensus_aa": "kinome_consensus"}),
                         on="KLIFS_position", how="left"))
    return cons, merged, kin


def build_tables(cons, merged, kin):
    merged["ortholog_conservation"] = merged["conservation"]

    def call(r):
        o, k = r["ortholog_conservation"], r["kinome_conservation"]
        if pd.isna(o) or pd.isna(k):
            return "unmapped"
        if o >= 0.75 and k >= 0.50:
            return "Universally immutable (catalytic core)"
        if o >= 0.75 and k < 0.50:
            return "CDK-specific, evolutionarily fixed"
        if o < 0.45:
            return "Tolerant in evolution"
        return "Intermediate"
    merged["site_class"] = merged.apply(call, axis=1)

    out = merged[["gene", "KLIFS_position", "KLIFS_label", "region", "human_position",
                  "human_residue", "klifs_aa", "kinome_consensus", "consensus_residue",
                  "ortholog_conservation", "pct_identical_to_human", "kinome_conservation",
                  "mutability_class", "site_class", "n_seqs_with_residue",
                  "n_clades_covered", "coverage"]].round(4)
    out.to_csv(f"{CONS}/TABLE_CDK_ortholog_conservation_KLIFS.csv", index=False)

    comp = []
    for g in CDKS:
        d = pd.read_csv(f"{ORTH}/{g}_orthologs.csv")
        v = d.tax_bin.value_counts()
        v.name = g
        comp.append(v)
    comp = pd.concat(comp, axis=1).fillna(0).astype(int)
    comp = comp.reindex([c for c in CLADE_ORDER if c in comp.index])

    with pd.ExcelWriter(f"{CONS}/TABLE_CDK_ortholog_conservation.xlsx",
                        engine="openpyxl") as w:
        cons.round(4).to_excel(w, "Per_residue_all_CDKs", index=False)
        out.to_excel(w, "KLIFS_pocket_per_CDK", index=False)
        (out[out.site_class == "Universally immutable (catalytic core)"]
         .to_excel(w, "Cannot_be_mutated", index=False))
        (out[out.site_class == "CDK-specific, evolutionarily fixed"]
         .to_excel(w, "CDK_specific_fixed", index=False))
        (out[out.site_class == "Tolerant in evolution"]
         .to_excel(w, "Can_be_mutated", index=False))
        comp.to_excel(w, "Homolog_set_composition")
    print(f"wrote {CONS}/TABLE_CDK_ortholog_conservation.xlsx")
    return out, comp


def figure3(merged, comp, out):
    fig = plt.figure(figsize=(17, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], width_ratios=[1.15, 1],
                          hspace=0.30, wspace=0.22,
                          left=0.065, right=0.975, top=0.94, bottom=0.07)

    # ---- a: how deep the search went
    axA = fig.add_subplot(gs[0, 0])
    bottom = np.zeros(len(CDKS))
    for clade in comp.index:
        vals = comp.loc[clade, CDKS].values.astype(float)
        axA.bar(CDKS, vals, bottom=bottom, label=clade,
                color=CLADE_COLORS.get(clade, "#999999"), edgecolor="white", lw=0.5)
        bottom += vals
    axA.set_ylabel("Homologous sequences (one per organism)", fontsize=10)
    axA.set_title("Depth of the homolog search", fontsize=12, weight="bold")
    axA.set_ylim(0, bottom.max() * 1.62)
    axA.legend(fontsize=7, ncol=3, loc="upper left", framealpha=0.92,
               borderpad=0.5, columnspacing=1.0, handlelength=1.1)
    axA.tick_params(axis="x", rotation=45, labelsize=9)
    for i, g in enumerate(CDKS):
        axA.text(i, bottom[i] + 6, f"{int(bottom[i])}", ha="center", fontsize=8)
    sns.despine(ax=axA)
    axA.text(-0.09, 1.04, "a", transform=axA.transAxes, fontsize=15, weight="bold")

    # ---- b: ortholog conservation at the 85 KLIFS positions, per CDK
    axB = fig.add_subplot(gs[1, :])
    # one KLIFS position is a gap in every CDK, so reindex to the full 85
    piv = (merged.pivot_table(index="gene", columns="KLIFS_position",
                              values="ortholog_conservation")
                 .reindex(index=CDKS, columns=range(1, 86)))
    sns.heatmap(piv, ax=axB, cmap="RdYlBu_r", vmin=0.2, vmax=1.0,
                cbar_kws=dict(label="Ortholog conservation\n(cross-species, per CDK)",
                              shrink=0.8, pad=0.010))
    axB.set_xticks(np.arange(85) + 0.5)
    kin_tab = pd.read_csv(f"{CONS}/TABLE_kinome_KLIFS_conservation.csv")
    kinome_cons = kin_tab.sort_values("KLIFS_position")["consensus_aa"].values
    axB.set_xticklabels(kinome_cons, fontsize=6.5, family="monospace", rotation=0)
    axB.tick_params(axis="x", length=0)
    axB.set_yticklabels(CDKS, rotation=0, fontsize=10)
    axB.set_ylabel(""); axB.set_xlabel("KLIFS pocket position (kinome consensus residue)",
                                       fontsize=10)
    axB.set_title("Cross-species conservation of each screened CDK, mapped onto the KLIFS pocket",
                  fontsize=12, weight="bold", pad=8)
    axB.text(-0.065, 1.03, "c", transform=axB.transAxes, fontsize=15, weight="bold")

    # ---- c: the key scatter — kinome vs ortholog conservation
    axC = fig.add_subplot(gs[0, 1])
    palette = {"Universally immutable (catalytic core)": "#B2182B",
               "CDK-specific, evolutionarily fixed": "#2166AC",
               "Intermediate": "#BBBBBB", "Tolerant in evolution": "#4D9221"}
    for cls, col in palette.items():
        s = out[out.site_class == cls]
        axC.scatter(s["kinome_conservation"], s["ortholog_conservation"], s=16,
                    c=col, label=f"{cls} (n={len(s)})", alpha=0.75,
                    edgecolors="white", linewidths=0.3)
    axC.axhline(0.75, ls="--", lw=0.8, color="grey")
    axC.axvline(0.50, ls="--", lw=0.8, color="grey")
    axC.set_xlabel("Conservation across the human kinome\n(can another kinase use a different residue?)",
                   fontsize=9.5)
    axC.set_ylabel("Conservation across species within one CDK\n(has evolution ever changed it?)",
                   fontsize=9.5)
    axC.set_title("Which residues can and cannot be mutated", fontsize=12, weight="bold")
    axC.legend(fontsize=7, loc="lower left", framealpha=0.94)
    axC.set_ylim(top=1.055)
    # label the most kinome-divergent of the CDK-specific fixed sites
    lab = (out[(out.site_class == "CDK-specific, evolutionarily fixed")
               & (out.ortholog_conservation >= 0.97)]
           .nsmallest(7, "kinome_conservation"))
    texts = [axC.text(r.kinome_conservation, r.ortholog_conservation,
                      f"{r.gene} {r.human_residue}{r.human_position}", fontsize=6.4)
             for r in lab.itertuples()]
    adjust_text(texts, ax=axC,
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
    sns.despine(ax=axC)
    axC.text(-0.12, 1.04, "b", transform=axC.transAxes, fontsize=15, weight="bold")

    fig.savefig(f"{FIG}/CDK_ortholog_conservation.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/CDK_ortholog_conservation.png", dpi=280, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/CDK_ortholog_conservation.{{pdf,png}}")


def figure4(cons, merged):
    """Per-residue conservation tracks along each CDK."""
    fig, axes = plt.subplots(len(CDKS), 1, figsize=(15, 1.45 * len(CDKS)))
    for ax, g in zip(axes, CDKS):
        d = cons[cons.gene == g]
        dstart, dend = co.kinase_domain(g)
        ax.axvspan(dstart, dend, color="#FDE0DD", zorder=0)
        ax.fill_between(d.human_position, d.conservation, color="#3B5B92", lw=0)
        pk = merged[(merged.gene == g) & merged.human_position.notna()]
        ax.scatter(pk.human_position, [1.06] * len(pk), s=5, marker="v",
                   color="#B2182B", clip_on=False)
        ax.set_ylim(0, 1.12); ax.set_xlim(0, d.human_position.max())
        ax.set_ylabel(g, fontsize=10, rotation=0, ha="right", va="center")
        ax.set_yticks([0, 0.5, 1.0]); ax.tick_params(labelsize=7)
        n = int(d.n_seqs_total.iloc[0])
        ax.text(0.995, 0.86, f"n = {n} species", transform=ax.transAxes,
                ha="right", fontsize=7.5, color="#555555")
        sns.despine(ax=ax)
    axes[-1].set_xlabel("Residue position (human numbering)", fontsize=10)
    axes[0].set_title("Per-residue cross-species conservation of the screened CDKs\n"
                      "pink band = Pfam kinase domain;  red triangles = KLIFS pocket positions",
                      fontsize=12.5, weight="bold", pad=12)
    fig.tight_layout()
    fig.savefig(f"{FIG}/CDK_ortholog_conservation_tracks.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/CDK_ortholog_conservation_tracks.png", dpi=260, bbox_inches="tight")
    print(f"wrote {paths.rel(FIG)}/CDK_ortholog_conservation_tracks.{{pdf,png}}")


if __name__ == "__main__":
    cons, merged, kin = load_all()
    out, comp = build_tables(cons, merged, kin)
    figure3(merged, comp, out)
    figure4(cons, merged)
