"""Per-hit bar chart beside a sequence logo of what that guide actually edits.

This layout is shared by figures 2b, 2i, 3b and 4p. In the original it was
written out once per figure, in CDK7_SY5609_hits.ipynb, logo_plots.ipynb,
logo_plots_CDK12_13_canonical.ipynb and clinical.ipynb, with the logo builder
(`logo_df_maker`) duplicated alongside it each time.

Left: LFC per replicate per condition for each hit guide, bars with seaborn's
default 95% bootstrap CI across the three replicates.
Right: one logo per guide, over the edited window plus `context` residues
either side. Each column sums to 100 — the observed edits at their measured
%Reads, and the reference residue taking whatever is left — so glyph height is
the fraction of sensor reads carrying that residue. Edited positions are red.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec

from . import editing, io, paths
from .plotting import EDITOR_COLORS

CONTEXT = 1      # residues shown either side of the edited window
MIN_PERC = 2     # a variant below this share of sensor reads is not drawn

#: label colour for a value `label_colors` does not name -- supp 9b picks out
#: CDK7 in red and leaves the other three genes black
DEFAULT_LABEL_COLOR = "black"


def protein(gene: str) -> str:
    cdks = io.read(paths.REFERENCE / "library_generation" / "cdks_info_updated.csv")
    return cdks.loc[cdks["Gene"] == gene, "Protein"].values[0]


def logo_frame(savs: pd.DataFrame, gRNA_id: str, context: int = CONTEXT,
               min_perc: float = MIN_PERC):
    """-> (logo frame indexed by codon, first codon, last codon, variants, sequence).

    The frame has one row per codon and one column per amino acid seen there.
    """
    gene = savs.loc[savs["gRNA_id"] == gRNA_id, "Gene"].values[0]
    seq = protein(gene)

    d = savs[(savs["gRNA_id"] == gRNA_id) & (savs["HGVSp"] != "WT")
             & (savs["%Reads"] >= min_perc)].copy()
    if d.empty:
        raise ValueError(f"{gRNA_id}: no variant clears {min_perc}% of sensor reads")
    d["codon"] = [editing.first_number(h) for h in d["HGVSp"]]
    d["ref"] = [h[0] for h in d["HGVSp"]]
    d["alt"] = [h[-1] for h in d["HGVSp"]]

    lo, hi = int(d["codon"].min()), int(d["codon"].max())
    # the protein sequence must agree with what the variant calls say is there
    for c in (lo, hi):
        assert seq[c - 1] == d.loc[d["codon"] == c, "ref"].values[0], \
            f"{gRNA_id}: codon {c} is {seq[c - 1]} in the sequence but " \
            f"{d.loc[d['codon'] == c, 'ref'].values[0]} in the variant call"

    start, end = max(lo - context, 1), min(hi + context, len(seq))
    window = seq[start - 1:end]
    residues = sorted(set(window) | set(d["ref"]) | set(d["alt"]))

    frame = pd.DataFrame({"codon": list(range(start, end + 1)),
                          **{r: np.zeros(len(window)) for r in residues}})
    for _, v in d.iterrows():
        frame.loc[frame["codon"] == v["codon"], v["alt"]] += v["%Reads"]
    # whatever is left of the 100% belongs to the reference residue
    for i, codon in enumerate(range(start, end + 1)):
        remaining = 100 - frame.loc[frame["codon"] == codon, residues].values.sum()
        frame.loc[frame["codon"] == codon, window[i]] += remaining
    return frame.set_index("codon"), lo, hi, d, window


def logo_group(gs, row0: int, sub: pd.DataFrame, plot_df: pd.DataFrame, savs: dict,
               palette, hue_order, label_colors: dict, color_by: str = "Editor",
               xlim=None, xlabel: str = "LFC (vs. DMSO)", sizes: dict | None = None,
               context: int = CONTEXT, min_perc: float = MIN_PERC):
    """Draw one group of hits into rows row0.. of an existing gridspec.

    Column 0 of those rows is one shared bar axis; column 1 gets a logo per
    guide. Figures that stack several groups (figure 3b has two) call this once
    per group; `hits_with_logos` is the single-group case.

    `label_colors` maps the value of `color_by` to a colour -- the published
    panels colour the variant labels by editor (2b, 2i, 4p) or by gene (3b).
    """
    import logomaker

    sizes = sizes or dict(logo=9, ylabel=12, xtick=12, xlabel=12, tick_len=3.5, tick_w=0.8)
    sub = sub.reset_index(drop=True)
    order = list(sub["key"])

    axb = plt.subplot(gs[row0:row0 + len(sub), 0])
    rows = plot_df[plot_df["key"].isin(order)]
    sns.barplot(data=rows, x="LFC", y="key", hue="Condition", order=order,
                hue_order=hue_order, palette=palette, ax=axb, linewidth=1,
                edgecolor="black", legend=False)
    sns.stripplot(data=rows, x="LFC", y="key", hue="Condition", order=order,
                  hue_order=hue_order, palette=palette, ax=axb, dodge=True,
                  linewidth=1, edgecolor="black", size=4, legend=False)

    for i, row in sub.iterrows():
        frame, _, _, d, _ = logo_frame(savs[row["Editor"]], row["gRNA_id"],
                                       context, min_perc)
        ax1 = plt.subplot(gs[row0 + i, 1])
        logo = logomaker.Logo(frame, vpad=0.02, ax=ax1)
        logo.style_xticks(anchor=0, spacing=1, rotation=0, fontsize=sizes["logo"])
        logo.ax.set_yticks([0, 50, 100])
        logo.ax.set_yticklabels([])
        logo.style_spines(spines=["right", "top"], visible=False)
        logo.ax.tick_params(axis="both", which="major", labelsize=sizes["logo"])
        for _, v in d[d["%Reads"] > 0].iterrows():
            logo.style_single_glyph(c=v["alt"], p=v["codon"], color="tab:red")

    axb.set_yticklabels(list(sub["label"]))
    for i, lbl in enumerate(axb.get_yticklabels()):
        lbl.set_color(label_colors.get(sub.loc[i, color_by], DEFAULT_LABEL_COLOR))
    if xlim is not None:
        axb.set_xlim(xlim)
    axb.tick_params(axis="y", which="major", labelsize=sizes["ylabel"])
    axb.tick_params(axis="x", which="major", labelsize=sizes["xtick"],
                    length=sizes["tick_len"], width=sizes["tick_w"])
    axb.set_ylabel("")
    axb.set_xlabel(xlabel, fontsize=sizes["xlabel"])
    for i in range(len(sub) - 1):
        axb.axhline(i + 0.5, linestyle="dashed", color="black", linewidth=0.5)
    axb.axvline(0, linestyle="dashed", color="black", linewidth=0.5)
    axb.spines[["right", "top"]].set_visible(False)
    return axb


def hits_with_logos(hits: pd.DataFrame, plot_df: pd.DataFrame, savs: dict,
                    palette, hue_order, title: str = "", xlabel: str = "LFC (vs. DMSO)",
                    context: int = CONTEXT, min_perc: float = MIN_PERC,
                    width: float = 9.0, row_h: float = 1.0, legend: bool = True,
                    sizes: dict | None = None, label_colors: dict | None = None,
                    color_by: str = "Editor", legend_loc: str = "axes",
                    legend_fontsize: float = 9):
    """One group of hits as a standalone figure.

    `hits` needs key, gRNA_id, Editor and label; `plot_df` needs key, LFC, Condition.

    `sizes` overrides the font and tick sizes -- the supplementary panels are
    printed larger than the main-figure ones and use the notebooks' `_large`
    style. `label_colors`/`color_by` choose what the variant labels are
    coloured by: editor for 2b/2i/4n and supp 8c, gene for 3b and supp 9b.
    """
    hits = hits.reset_index(drop=True)
    n = len(hits)
    fig = plt.figure(figsize=(width, row_h * n + 1.2))
    gs = gridspec.GridSpec(n, 2, width_ratios=[1.5, 1], wspace=0.15, hspace=0.55)
    axb = logo_group(gs, 0, hits, plot_df, savs, palette, hue_order,
                     label_colors=label_colors or EDITOR_COLORS, color_by=color_by,
                     xlabel=xlabel, context=context, min_perc=min_perc, sizes=sizes)
    if legend:
        # the group drawer suppresses legends; the standalone panels show one.
        # `legend_loc="figure"` puts it across the top instead of inside the bar
        # axis, which is what the supplementary panels do -- they are tall
        # enough that an in-axis legend lands on top of a row.
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="black")
                   for c in palette]
        if legend_loc == "figure":
            fig.legend(handles, list(hue_order), fontsize=legend_fontsize,
                       ncol=3, frameon=False, loc="upper center",
                       bbox_to_anchor=(0.5, 1 - 0.60 / fig.get_figheight()))
        else:
            axb.legend(handles, list(hue_order), fontsize=legend_fontsize,
                       loc="lower right", frameon=False)
    if title:
        fig.suptitle(title, fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
    else:
        # no title, so no reserved strip -- otherwise the long variant labels clip
        fig.tight_layout()
    return fig
