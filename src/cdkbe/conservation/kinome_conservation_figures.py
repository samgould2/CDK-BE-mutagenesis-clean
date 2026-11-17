"""Figures for kinome-wide KLIFS pocket conservation."""
import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import logomaker

from . import kinome_conservation as kc

from ._paths import (ALN, BLAST, CONS, EDITING, FIG, KLIFS, LFC, ORTH,  # noqa: F401
                      OUT, ROOT)
plt.rc("font", family="Helvetica")


REGION_COLORS = {
    "I": "#8DA0CB", "g.l": "#4F6FB5", "II": "#8DA0CB", "III": "#8DA0CB",
    "aC": "#E78AC3", "b.l": "#B3B3B3", "IV": "#8DA0CB", "V": "#8DA0CB",
    "GK": "#E41A1C", "hinge": "#FC8D62", "linker": "#FDCDAC", "aD": "#66C2A5",
    "aE": "#66C2A5", "VI": "#66C2A5", "c.l": "#1B7837", "VII": "#66C2A5",
    "VIII": "#66C2A5", "xDFG": "#984EA3", "a.l": "#CAB2D6",
}
GROUP_ORDER = ["AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", "Other", "Atypical"]

# Landmarks worth calling out on the axis
LANDMARKS = {4: "G-loop", 6: "G-loop", 11: "β2 Val", 17: "VAIK Lys", 24: "αC-Glu",
             45: "Gatekeeper", 46: "Hinge1", 48: "Hinge3", 68: "HRD-His",
             70: "HRD-Asp\n(catalytic)", 75: "Cat-loop Asn", 81: "DFG-Asp",
             82: "DFG-Phe", 83: "DFG-Gly"}


def _xaxis(ax, master, letters=True):
    ax.set_xlim(0.5, 85.5)
    ax.set_xticks(range(1, 86))
    if letters:
        ax.set_xticklabels(master["consensus_aa"], fontsize=6.5, family="monospace")
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis="x", length=0)


def region_bar(ax, master):
    """Thin ribbon showing KLIFS structural regions."""
    start = 1
    for name, n in kc.REGIONS:
        ax.add_patch(mpatches.Rectangle((start - 0.5, 0), n, 1,
                                        color=REGION_COLORS[name], ec="white", lw=0.6))
        if n >= 3:
            ax.text(start + n / 2 - 0.5, 0.5, name, ha="center", va="center",
                    fontsize=6.5, color="white", weight="bold")
        start += n
    ax.set_xlim(0.5, 85.5); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def figure1(master, per_group, cdk_tab, mat):
    fig = plt.figure(figsize=(17, 11.5))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.5, 0.22, 2.0, 2.6, 1.35],
                          hspace=0.12, left=0.075, right=0.965, top=0.955, bottom=0.06)

    # ---- A: sequence logo of the human kinome pocket
    freqs = json.load(open(os.path.join(OUT, "klifs_position_aa_frequencies.json")))
    fm = pd.DataFrame({int(k): v for k, v in freqs.items()}).T.fillna(0.0)
    fm = fm.reindex(range(1, 86)).fillna(0.0)
    info_mat = logomaker.transform_matrix(fm, from_type="probability", to_type="information")
    axA = fig.add_subplot(gs[0])
    logomaker.Logo(info_mat, ax=axA, color_scheme="chemistry", show_spines=False)
    axA.set_ylabel("Information\n(bits)", fontsize=10)
    axA.set_xlim(0.5, 85.5); axA.set_xticks([]); axA.set_ylim(0, 4.4)
    axA.tick_params(labelsize=8)
    axA.set_title("Evolutionary constraint across the human kinome at the 85 KLIFS pocket positions "
                  f"(n = {len(mat)} kinases)", fontsize=13, weight="bold", pad=10)
    axA.text(-0.055, 1.02, "a", transform=axA.transAxes, fontsize=15, weight="bold")

    # ---- region ribbon
    axR = fig.add_subplot(gs[1]); region_bar(axR, master)

    # ---- B: kinome conservation bars
    axB = fig.add_subplot(gs[2])
    colors = [REGION_COLORS[r] for r in master["region"]]
    axB.bar(master["KLIFS_position"], master["conservation"], color=colors,
            edgecolor="black", lw=0.4, width=0.82)
    for y, lab, c in [(0.70, "immutable", "#B2182B"), (0.50, "strongly constrained", "#EF8A62"),
                      (0.30, "moderately constrained", "#999999"), (0.15, "tolerant", "#BBBBBB")]:
        axB.axhline(y, ls="--", lw=0.8, color=c, zorder=0)
        axB.text(85.8, y, lab, fontsize=7, va="center", color=c)
    axB.set_ylabel("Conservation\n(1 − H/log₂20)", fontsize=10)
    axB.set_ylim(0, 1.0); _xaxis(axB, master, letters=False)
    axB.tick_params(labelsize=8)
    for pos, lab in LANDMARKS.items():
        v = master.loc[master.KLIFS_position == pos, "conservation"].iloc[0]
        axB.annotate(lab, (pos, v + 0.02), fontsize=6.8, rotation=90,
                     ha="center", va="bottom", color="black")
    axB.set_ylim(0, 1.22)
    axB.text(-0.055, 1.0, "b", transform=axB.transAxes, fontsize=15, weight="bold")

    # ---- C: conservation by kinase group
    axC = fig.add_subplot(gs[3])
    wide = per_group.pivot(index="group", columns="KLIFS_position", values="conservation")
    order = [g for g in GROUP_ORDER if g in wide.index]
    wide = wide.loc[order]
    wide.loc["CDK family"] = cdk_tab.set_index("KLIFS_position")["conservation"]
    n_by = per_group.groupby("group")["n_seqs"].first().to_dict()
    n_by["CDK family"] = int(cdk_tab["n_seqs"].iloc[0])
    sns.heatmap(wide, ax=axC, cmap="RdYlBu_r", vmin=0, vmax=1, cbar_kws=dict(
        label="Conservation", shrink=0.75, pad=0.012, aspect=14))
    axC.set_yticklabels([f"{g}  (n={n_by[g]})" for g in wide.index],
                        rotation=0, fontsize=9)
    axC.set_xlabel(""); axC.set_ylabel("")
    axC.set_xticks(np.arange(85) + 0.5)
    axC.set_xticklabels(master["consensus_aa"], fontsize=6.5, family="monospace", rotation=0)
    axC.tick_params(axis="x", length=0)
    axC.axhline(len(order), color="black", lw=1.8)
    axC.text(-0.055, 1.01, "c", transform=axC.transAxes, fontsize=15, weight="bold")

    # ---- D: numbering strip + screened-CDK residue examples
    axD = fig.add_subplot(gs[4]); axD.axis("off")
    for p in range(1, 86):
        axD.text(p, 0.90, str(p), fontsize=6.2, ha="center", va="center", rotation=90)
    axD.text(0.0, 0.90, "KLIFS #", fontsize=8, ha="right", va="center",
             transform=axD.get_yaxis_transform())
    for j, cdk in enumerate(["CDK7", "CDK9", "CDK12", "CDK2"]):
        y = 0.66 - j * 0.165
        axD.text(0.0, y, cdk, fontsize=8, ha="right", va="center",
                 transform=axD.get_yaxis_transform())
        for p in range(1, 86):
            r = master.loc[master.KLIFS_position == p, f"{cdk}_residue"].iloc[0]
            if r:
                axD.text(p, y, r[1:], fontsize=5.4, ha="center", va="center", rotation=90,
                         color="#444444")
    axD.set_xlim(0.5, 85.5); axD.set_ylim(0, 1)
    axD.text(-0.055, 1.0, "d", transform=axD.transAxes, fontsize=15, weight="bold")

    fig.savefig(os.path.join(FIG, "kinome_KLIFS_conservation.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "kinome_KLIFS_conservation.png"), dpi=280, bbox_inches="tight")
    print("wrote figures/kinome_KLIFS_conservation.{pdf,png}")
    return fig


def figure2(master, spec):
    """Group-specific constraint, corrected against a size-matched null."""
    order = [g for g in GROUP_ORDER if g in set(spec["group"])] + ["CDK family"]
    zc = spec.pivot(index="group", columns="KLIFS_position", values="z_conservation").loc[order]
    zk = spec.pivot(index="group", columns="KLIFS_position", values="z_KL").loc[order]
    n_by = spec.groupby("group")["n_kinases"].first().to_dict()
    ylabs = [f"{g}  (n={n_by[g]})" for g in order]

    fig, axes = plt.subplots(4, 1, figsize=(17, 12.5),
                             gridspec_kw=dict(height_ratios=[2.4, 0.22, 2.4, 2.3], hspace=0.40),
                             constrained_layout=False)
    ax0, axR, ax1, ax2 = axes

    sns.heatmap(zc, ax=ax0, cmap="PuOr_r", center=0, vmin=-8, vmax=8,
                cbar_kws=dict(label="z (excess conservation)", shrink=0.85, pad=0.012))
    ax0.set_xticks([]); ax0.set_xlabel(""); ax0.set_ylabel("")
    ax0.set_yticklabels(ylabs, rotation=0, fontsize=9)
    ax0.axhline(len(order) - 1, color="black", lw=1.8)
    ax0.set_title("Group-specific constraint at KLIFS pocket positions, corrected for subset size\n"
                  "(purple = locked down within the group beyond what a random size-matched set of kinases shows)",
                  fontsize=12.5, weight="bold", pad=12)
    ax0.text(-0.075, 1.04, "a", transform=ax0.transAxes, fontsize=15, weight="bold")

    region_bar(axR, master)

    sns.heatmap(zk, ax=ax1, cmap="Reds", vmin=0, vmax=25,
                cbar_kws=dict(label="z (KL divergence vs kinome)", shrink=0.85, pad=0.012))
    ax1.set_xticks(np.arange(85) + 0.5)
    ax1.set_xticklabels(master["consensus_aa"], fontsize=6.5, family="monospace", rotation=0)
    ax1.tick_params(axis="x", length=0)
    ax1.set_yticklabels(ylabs, rotation=0, fontsize=9)
    ax1.set_xlabel("")
    ax1.set_ylabel("")
    ax1.axhline(len(order) - 1, color="black", lw=1.8)
    ax1.set_title("Specificity-determining positions — the group prefers a different residue than the kinome does",
                  fontsize=11, weight="bold", pad=8)
    ax1.text(-0.075, 1.06, "b", transform=ax1.transAxes, fontsize=15, weight="bold")

    # ---- CDK family detail
    cdk = spec[spec.group == "CDK family"].set_index("KLIFS_position").sort_index()
    cols = ["#7B3294" if z > 4 else ("#E66101" if z < -2 else "#BBBBBB")
            for z in cdk["z_conservation"]]
    ax2.bar(cdk.index, cdk["z_conservation"], color=cols, edgecolor="black", lw=0.4, width=0.82)
    ax2.axhline(0, color="black", lw=0.8)
    ax2.axhline(4, ls="--", lw=0.8, color="#7B3294")
    ax2.set_ylabel("z (excess conservation)\nCDK family vs size-matched kinome", fontsize=9.5)
    ax2.set_xlim(0.5, 85.5)
    ax2.set_xticks(range(1, 86))
    ax2.set_xticklabels([f"{master.loc[master.KLIFS_position==p,'consensus_CDKfamily'].iloc[0]}{p}"
                         for p in range(1, 86)], fontsize=5.6, rotation=90, family="monospace")
    ax2.set_xlabel("KLIFS position (CDK-family consensus residue)", fontsize=10)
    top = cdk.nlargest(10, "z_KL").index
    for p in top:
        r7 = master.loc[master.KLIFS_position == p, "CDK7_residue"].iloc[0]
        r9 = master.loc[master.KLIFS_position == p, "CDK9_residue"].iloc[0]
        ax2.annotate(f"{r7} / {r9}", (p, cdk.loc[p, "z_conservation"]), fontsize=6,
                     ha="center", va="bottom", rotation=90, color="#7B3294", weight="bold")
    ax2.set_ylim(top=cdk["z_conservation"].max() * 1.7)
    ax2.text(-0.075, 1.02, "c", transform=ax2.transAxes, fontsize=15, weight="bold")
    ax2.set_title("CDK family detail — purple = held fixed more tightly than kinases in general;\n"
                  "labels give the CDK7 / CDK9 residue at the 10 strongest specificity-determining sites",
                  fontsize=10.5, pad=8)
    sns.despine(ax=ax2)

    fig.subplots_adjust(left=0.115, right=0.965, top=0.935, bottom=0.075)
    fig.savefig(os.path.join(FIG, "kinome_KLIFS_group_specific_constraint.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "kinome_KLIFS_group_specific_constraint.png"), dpi=280, bbox_inches="tight")
    print("wrote figures/kinome_KLIFS_group_specific_constraint.{pdf,png}")
    return fig


if __name__ == "__main__":
    master, per_group, cdk_tab, mat, info = kc.main()
    groups = {g: list(i) for g, i in info.groupby("group").groups.items() if len(i) >= 10}
    groups["CDK family"] = list(info.index[info["family"].str.startswith("CDK", na=False)])
    spec = kc.specificity_table(mat, info, groups, n_perm=1000)
    spec.to_csv(os.path.join(OUT, "klifs_group_specificity.csv"), index=False)
    figure1(master, per_group, cdk_tab, mat)
    figure2(master, spec)
