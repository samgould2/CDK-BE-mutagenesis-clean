"""Figure 4m and supplementary 13e-13g: the T47D CDK variant ± palbociclib RNA-seq.

One experiment, four panels. Six T47D lines -- the two non-targeting controls
and four base-edited CDK variants -- grown in DMSO or palbociclib, three
replicates each, 36 libraries in all (nf-core/rnaseq, STAR + Salmon):

    4m    cell-cycle marker gene heatmap, DMSO above palbociclib
    S13e  PCA of the 36 libraries
    S13f  GSEA dotplot, five Hallmark sets by four variants (source only)
    S13g  volcano grid, each variant against its matched NT control

This is **not** the T47D RNA-seq in `rnaseq_t47d.py`. That one is figure 5:
four CDK4 genotypes, vehicle only, twelve libraries. The two were run months
apart against different GENCODE releases -- 78,932 genes here against 78,458
there -- which is why figure 4m could not be drawn from the matrix this
repository already had; see `docs/DISCREPANCIES.md` item 32, which this module
closes.

They do share three sample **names**, and not the libraries behind them:
`T47D_CBE_NT_veh_1`, `_2` and `_3` appear in both matrices. Within either
experiment two CBE NT replicates correlate at 0.994 or better on log2 TPM;
the same name across the two correlates at 0.93. Joining the matrices on
sample name would silently pair unrelated libraries, so nothing here does --
`tests/test_t47d_cdk_palbo.py` pins the distinction.

The delivered package is `data/processed/rnaseq/T47D_CDK_palbo_RNAseq/`: the
Salmon gene-level TPM matrix, the eight DESeq2 contrasts, the MSigDB Hallmark
sets the highlight and the GSEA use, the three notebooks that drew the panels,
and the two z-score tables the marker-gene notebook wrote. The notebooks are
kept unchanged as the record of what was run -- they read literal `.csv` paths
relative to their own folder, so they do not run here (item 17); this module is
their plotting code moved into the package, reading through `cdkbe.io`.

Three things are worth stating because they change what comes out:

* **The marker-gene panel z-scores within each treatment arm**, not across the
  whole matrix. Each panel therefore shows how the genotypes compare under that
  treatment, and the DMSO -> palbociclib shift itself is deliberately not shown.
  Both panels share one colour scale.
* **Two CDK6 lines are relabelled on every panel.** `CDK6_G53S` prints as
  `CDK6 G53S_M54I` and `CDK6_R46W` as `CDK6 R44C_R46W`, because each line
  carries a second edit the sample name does not record. The delivered PCA
  notebook predates that relabelling and prints the short names; the published
  panel does not, and `DISPLAY_NAMES` here is what the published panel says.
* **The volcano's y axis is ticked every 20.** The delivered notebook leaves
  the locator to matplotlib, which ticks every 10; the published panel ticks
  every 20. Same data, same limits -- see `docs/DISCREPANCIES.md` item 35.

S13f is the one panel of the four that is not redrawn. Its dots encode NES by
colour and FDR q-value by area, and the published values came out of a GSEA run
whose permutation seed was not delivered, so a redraw here would print numbers
that are close to the published ones rather than the published ones. What this
module ships instead is the input that run consumed: the pre-ranked gene list
per variant and the five Hallmark sets the panel draws. `docs/COVERAGE.md`
marks it **source only**.
"""
from __future__ import annotations

import inspect
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import io, paths

BASE = paths.PROCESSED / "rnaseq" / "T47D_CDK_palbo_RNAseq"
EXPR_FILE = BASE / "salmon.merged.gene_tpm.csv"
DESEQ2_DIR = BASE / "deseq2_results"
HALLMARK_CACHE = BASE / "hallmark_gene_sets_cache.json"

#: genotype order on every panel, the two non-targeting controls first
GENOTYPES = ["ABE_NT", "CBE_NT", "CDK2_S261G", "CDK4_I51V",
             "CDK6_G53S", "CDK6_R46W"]
CONTROLS = ["ABE_NT", "CBE_NT"]
#: variant -> the non-targeting arm it is compared against
CONTROL_FOR = {"CDK2_S261G": "ABE_NT", "CDK4_I51V": "ABE_NT",
               "CDK6_G53S": "CBE_NT", "CDK6_R46W": "CBE_NT"}
#: file suffix -> the label the panels print. `veh` is DMSO.
TREATMENT_LABELS = {"veh": "DMSO", "palbo": "Palbociclib"}
ARMS = ["veh", "palbo"]

#: the labels the published panels print, which are not the sample names: each
#: CDK6 line carries a second edit that the sample name does not record
DISPLAY_NAMES = {"CDK6_G53S": "CDK6 G53S_M54I", "CDK6_R46W": "CDK6 R44C_R46W"}

SAMPLE_RE = r"^T47D_(?P<genotype>.+)_(?P<treatment>palbo|veh)_\d+$"


def pretty_label(g: str) -> str:
    """`CDK2_S261G` -> `CDK2 S261G`, `ABE_NT` -> `ABE NT`, with the relabels."""
    if g in DISPLAY_NAMES:
        return DISPLAY_NAMES[g]
    return re.sub(r"\b(CDK\d+)_(?=[A-Z]\d+[A-Z])", r"\1 ", str(g)).replace("_NT", " NT")


# ----------------------------------------------------------------- reading --

def expression():
    """-> (log2(TPM+1) genes x samples, sample metadata, gene_id -> symbol).

    Samples come back in panel order: genotype as `GENOTYPES` lists them, DMSO
    before palbociclib, replicate within that.
    """
    expr = io.read(EXPR_FILE, index_col=0)
    gene_names = expr.pop("gene_name")
    parsed = expr.columns.str.extract(SAMPLE_RE)
    metadata = pd.DataFrame({"genotype": parsed["genotype"].values,
                             "treatment": parsed["treatment"].values},
                            index=expr.columns)
    order = (metadata.assign(g=metadata["genotype"].map({g: i for i, g in enumerate(GENOTYPES)}),
                             t=metadata["treatment"].map({"veh": 0, "palbo": 1}))
             .sort_values(["g", "t"]).index.tolist())
    return np.log2(expr[order] + 1), metadata.loc[order], gene_names


def deseq2(variant: str, arm: str) -> pd.DataFrame:
    """One DESeq2 contrast, variant against its matched NT arm.

    Rows with no `padj` -- DESeq2's independent filtering and outlier calls --
    are dropped, as the delivered notebook does, so a count taken from this
    frame is the count the published panel prints.
    """
    ctrl = CONTROL_FOR[variant]
    d = io.read(DESEQ2_DIR / f"{variant}_vs_{ctrl}_{arm}.csv")
    d = d.dropna(subset=["padj", "log2FoldChange"]).copy()
    d["y"] = -np.log10(d["padj"])
    return d


def hallmark() -> dict:
    """The MSigDB Hallmark sets delivered with the analysis -> {name: [genes]}."""
    return json.loads(io.resolve(HALLMARK_CACHE).read_text())


# ------------------------------------------------------ 4m: marker heatmap --

#: the panel's genes, top to bottom
MARKER_GENES = ["MKI67", "PCNA",                                    # proliferation
                "CCNE1", "CCNE2", "MCM2", "MCM7", "RRM2", "TYMS",   # G1/S
                "CCNA2",                                            # S/G2
                "CCNB1", "CCNB2", "CDK1", "TOP2A", "AURKA", "PLK1",  # G2/M
                "CDKN1A", "CDKN2A"]                                 # p21, p16

FS_TITLE, FS_GENE, FS_XTICK, FS_CBAR, FS_CBTICK = 15, 9, 10, 12, 11
#: inches per gene row / per sample column, and the colour-bar height
ROW_H, COL_W, CBAR_H = 0.25, 0.3, 0.14


def _marker_label(g: str) -> str:
    """x-axis label: controls on one line, variants with the gene above it."""
    if g in CONTROLS:
        return g.replace("_", " ")
    gene, variant = str(g).split("_", 1)
    return f"{gene}\n{DISPLAY_NAMES.get(g, gene + ' ' + variant).split(' ', 1)[1]}"


def _zscore_rows(m: pd.DataFrame) -> pd.DataFrame:
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)


def marker_expression():
    """-> (log2 TPM of the marker genes, samples x genes, sample metadata).

    Where a symbol maps to more than one Ensembl id, the more highly expressed
    id wins -- the delivered notebook's rule, and it matters for CDKN2A.
    """
    log_expr, metadata, gene_names = expression()
    sym_to_id = {}
    for gid, gname in gene_names.items():
        if gname in MARKER_GENES and (
                gname not in sym_to_id
                or log_expr.loc[gid].mean() > log_expr.loc[sym_to_id[gname]].mean()):
            sym_to_id[gname] = gid
    present = [g for g in MARKER_GENES if g in sym_to_id]
    marker = log_expr.loc[[sym_to_id[g] for g in present]]
    marker.index = present
    return marker, metadata


def marker_zscores() -> dict:
    """-> {"veh": genes x samples, "palbo": ...}, z-scored within each arm."""
    marker, metadata = marker_expression()
    out = {}
    for arm in ARMS:
        cols = [s for s in marker.columns if metadata.loc[s, "treatment"] == arm]
        out[arm] = _zscore_rows(marker[cols])
    return out


def _genotype_blocks(genotypes):
    """[(first_col, last_col, genotype), ...] for runs of the same genotype."""
    blocks, start = [], 0
    for i in range(1, len(genotypes) + 1):
        if i == len(genotypes) or genotypes[i] != genotypes[start]:
            blocks.append((start, i - 1, genotypes[start]))
            start = i
    return blocks


def panel_4m(outdir=None, **_):
    """The stacked marker-gene heatmap: DMSO above palbociclib, one colour bar."""
    _, metadata = marker_expression()
    arm_z = marker_zscores()
    vmax = np.ceil(max(np.nanmax(np.abs(z.values)) for z in arm_z.values()) * 2) / 2

    n_rows = len(arm_z["veh"])
    n_cols = arm_z["veh"].shape[1]
    heat_h = ROW_H * n_rows
    heights = [heat_h, 0.95, heat_h, 0.75, CBAR_H]
    fig = plt.figure(figsize=(COL_W * n_cols + 1.2, sum(heights) + 0.4))
    gs = fig.add_gridspec(5, 3, height_ratios=heights,
                          width_ratios=[0.2, 0.6, 0.2], hspace=0)
    axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[2, :])]
    cax = fig.add_subplot(gs[4, 1])

    im = None
    for ax, arm in zip(axes, ARMS):
        z = arm_z[arm]
        im = ax.imshow(z.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       aspect="auto")
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(z.index, fontsize=FS_GENE, style="italic")
        ax.tick_params(axis="y", length=0)

        geno = metadata.loc[z.columns, "genotype"].values
        blocks = _genotype_blocks(geno)
        ax.set_xticks([(a + b) / 2 for a, b, _ in blocks])
        ax.set_xticklabels([_marker_label(g) for _, _, g in blocks],
                           fontsize=FS_XTICK, rotation=0, ha="center",
                           va="top", linespacing=1.2)
        ax.tick_params(axis="x", length=0, pad=3)
        for _, b, _ in blocks[:-1]:
            ax.axvline(b + 0.5, color="black", lw=0.3)
        n_ctrl = int(np.isin(geno, CONTROLS).sum())
        if 0 < n_ctrl < len(geno):
            ax.axvline(n_ctrl - 0.5, color="black", lw=1.0)   # controls | variants
        for sp in ax.spines.values():
            sp.set_edgecolor("black")
            sp.set_linewidth(0.3)
        ax.set_title(TREATMENT_LABELS[arm], fontsize=FS_TITLE, pad=6)

    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("z-score (log2 TPM, within arm)", fontsize=FS_CBAR)
    cb.ax.tick_params(labelsize=FS_CBTICK, width=0.3)
    cb.outline.set_linewidth(0.3)
    _save(fig, outdir, _STEM["4m"])
    return fig, pd.concat({TREATMENT_LABELS[a]: arm_z[a] for a in ARMS},
                          names=["treatment"])


# ----------------------------------------------------------- S13e: the PCA --

N_TOP = 500                 # most-variable genes the PCA is computed over
PCA_COLORS = {
    "ABE_NT": "#8C8C86", "CBE_NT": "#55554F",   # the darker grey tells them apart
    "CDK2_S261G": "#5B84B1", "CDK4_I51V": "#9A82B0",
    "CDK6_G53S": "#D9702E", "CDK6_R46W": "#7FB49B",
}
TREATMENT_MARKERS = {"veh": "o", "palbo": "^"}

#: An SVD fixes each component only up to its sign, and which sign LAPACK
#: returns depends on the build -- this machine gets PC2 the other way up from
#: the one the published panel was drawn on, which would mirror the plot and
#: relabel the axis from -15..10 to -10..15. Each component is therefore
#: oriented against what drives it: PC1 (52%) is the CDK6 R44C_R46W line
#: standing apart from the other five, PC2 (25%) is DMSO against palbociclib.
#: Both margins are large -- +32 and +10 score units -- so neither rule is
#: close to a coin toss. Written as component index -> (metadata column, the
#: value whose samples sit on the positive side).
PC_ORIENTATION = {0: ("genotype", "CDK6_R46W"), 1: ("treatment", "veh")}


def pca():
    """-> (PC1/PC2 per sample with its genotype and arm, % variance per PC).

    Top-`N_TOP` most-variable genes of log2(TPM+1), mean-centred, SVD -- the
    delivered notebook's own PCA rather than DESeq2's -- with each component
    oriented by `PC_ORIENTATION`.
    """
    log_expr, metadata, _ = expression()
    X = log_expr.T
    top = X.var(axis=0).sort_values(ascending=False).head(N_TOP).index
    Xc = X[top] - X[top].mean(axis=0)
    U, S, _ = np.linalg.svd(Xc.values, full_matrices=False)
    scores = U * S
    for i, (column, value) in PC_ORIENTATION.items():
        anchor = (metadata[column].reindex(X.index) == value).to_numpy()
        if scores[anchor, i].mean() < 0:
            scores[:, i] = -scores[:, i]
    var_pct = (S ** 2) / np.sum(S ** 2) * 100
    out = (pd.DataFrame({"sample": X.index, "PC1": scores[:, 0],
                         "PC2": scores[:, 1]})
           .merge(metadata.reset_index(names="sample"), on="sample"))
    return out, var_pct


def panel_S13e(outdir=None, **_):
    """The PCA of all 36 libraries: colour is genotype, marker is treatment."""
    from matplotlib.lines import Line2D

    pca_df, var_pct = pca()
    fig, ax = plt.subplots(figsize=(5, 5))
    for (genotype, treatment), grp in pca_df.groupby(["genotype", "treatment"]):
        ax.scatter(grp["PC1"], grp["PC2"],
                   color=PCA_COLORS.get(genotype, "gray"),
                   marker=TREATMENT_MARKERS.get(treatment, "o"),
                   s=70, zorder=3, linewidths=0.5, edgecolors="black")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(f"PC1: {var_pct[0]:.0f}% variance", fontsize=11)
    ax.set_ylabel(f"PC2: {var_pct[1]:.0f}% variance", fontsize=11)

    # one column, a bold DMSO block above a bold Palbociclib block, genotypes
    # under each header in PCA_COLORS order
    shown = [g for g in PCA_COLORS if g in set(pca_df["genotype"])]
    handles, labels, header_idx = [], [], []
    for arm in ARMS:
        if handles:
            handles.append(Line2D([], [], linestyle="none"))
            labels.append("")
        header_idx.append(len(labels))
        handles.append(Line2D([], [], linestyle="none"))
        labels.append(TREATMENT_LABELS[arm])
        for g in shown:
            handles.append(Line2D([], [], linestyle="none",
                                  marker=TREATMENT_MARKERS[arm], markersize=7,
                                  markerfacecolor=PCA_COLORS[g],
                                  markeredgecolor="black", markeredgewidth=0.5))
            labels.append(pretty_label(g))
    leg = ax.legend(handles, labels, fontsize=8, frameon=False, ncol=1,
                    loc="center left", bbox_to_anchor=(1.0, 0.5),
                    handletextpad=0.4)
    for i in header_idx:
        leg.get_texts()[i].set_fontweight("bold")
    fig.tight_layout()
    _save(fig, outdir, _STEM["S13e"])
    return fig, pca_df


# ------------------------------------------------------- S13g: the volcanos --

LFC_THRESH, PADJ_THRESH = 1.0, 0.05
X_LIM = 6          # fixed log2FC axis range
N_LABELS = 8       # most-significant proliferation genes labelled per panel
Y_TICK_STEP = 20   # the published panel's locator; the notebook leaves it to matplotlib

UP_COLOR, DOWN_COLOR, NS_COLOR, THRESH_COLOR = "#d6604d", "#4393c3", "#d0d0d0", "#999999"
FS_ROW, FS_VTITLE, FS_AXIS, FS_TICK, FS_GENELAB, FS_LEGEND = 20, 15, 14, 12, 10, 11

#: the highlight is one group: the union of these two Hallmark sets
PROLIF_SETS = ["E2F Targets", "G2-M Checkpoint"]
PROLIF_LABEL = "Cell proliferation\n(Hallmark E2F + G2-M)"
PROLIF_COLOR = "#7B3294"


def proliferation_genes() -> set:
    hm = hallmark()
    return set().union(*(hm[s] for s in PROLIF_SETS))


def volcano_counts() -> pd.DataFrame:
    """The DE counts each volcano prints, one row per variant and arm."""
    rows = []
    for arm in ARMS:
        for variant in CONTROL_FOR:
            d = deseq2(variant, arm)
            sig = d["padj"] < PADJ_THRESH
            rows.append({
                "variant": variant, "label": pretty_label(variant),
                "vs_control": CONTROL_FOR[variant],
                "treatment": TREATMENT_LABELS[arm],
                "genes_tested": len(d),
                "up": int((sig & (d["log2FoldChange"] > LFC_THRESH)).sum()),
                "down": int((sig & (d["log2FoldChange"] < -LFC_THRESH)).sum()),
            })
    return pd.DataFrame(rows)


def labelled_genes() -> pd.DataFrame:
    """The proliferation genes each volcano names, and why they were picked."""
    prolif = proliferation_genes()
    rows = []
    for arm in ARMS:
        for variant in CONTROL_FOR:
            d = deseq2(variant, arm)
            cc = d[d["gene_name"].isin(prolif)]
            keep = cc[(cc["padj"] < PADJ_THRESH)
                      & (cc["log2FoldChange"].abs() > LFC_THRESH)
                      & (cc["log2FoldChange"].abs() <= X_LIM)].nsmallest(N_LABELS, "padj")
            for _, r in keep.iterrows():
                rows.append({"variant": variant, "label": pretty_label(variant),
                             "treatment": TREATMENT_LABELS[arm],
                             "gene_name": r["gene_name"],
                             "log2FoldChange": r["log2FoldChange"],
                             "padj": r["padj"]})
    return pd.DataFrame(rows)


def panel_S13g(outdir=None, **_):
    """The 2x4 volcano grid: DMSO on top, palbociclib below."""
    try:
        from adjustText import adjust_text
    except ImportError:                                   # labels will overlap
        adjust_text = None
    # `expand=` is adjustText 1.0 and later; 0.8 spells the same thing as two
    # separate keywords and raises on the tuple. pyproject asks for >=0.8, so
    # the call has to work either way -- it only moves labels, and no test
    # depends on where they land.
    spread = ({"expand": (1.2, 1.4)}
              if adjust_text is not None
              and "expand" in inspect.signature(adjust_text).parameters else {})

    prolif = proliferation_genes()
    results = {(v, a): deseq2(v, a) for a in ARMS for v in CONTROL_FOR}
    y_max = max(r.loc[np.isfinite(r["y"]), "y"].max() for r in results.values())
    y_lim = np.ceil(y_max) + 2

    nrows, ncols = len(ARMS), len(CONTROL_FOR)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows),
                             squeeze=False)
    for r, arm in enumerate(ARMS):
        for c, variant in enumerate(CONTROL_FOR):
            ax, rc = axes[r, c], results[(variant, arm)]
            sig = rc["padj"] < PADJ_THRESH
            up = sig & (rc["log2FoldChange"] > LFC_THRESH)
            down = sig & (rc["log2FoldChange"] < -LFC_THRESH)
            ns = ~(up | down)

            ax.scatter(rc.loc[ns, "log2FoldChange"], rc.loc[ns, "y"], s=4,
                       color=NS_COLOR, alpha=0.5, linewidths=0, rasterized=True)
            up_pts = ax.scatter(rc.loc[up, "log2FoldChange"], rc.loc[up, "y"],
                                s=12, color=UP_COLOR, alpha=0.85, linewidths=0,
                                rasterized=True)
            down_pts = ax.scatter(rc.loc[down, "log2FoldChange"], rc.loc[down, "y"],
                                  s=12, color=DOWN_COLOR, alpha=0.85,
                                  linewidths=0, rasterized=True)

            cc = rc[rc["gene_name"].isin(prolif)]
            ax.scatter(cc["log2FoldChange"], cc["y"], s=24, color=PROLIF_COLOR,
                       edgecolors="black", linewidths=0.3, zorder=5,
                       label=(PROLIF_LABEL if (r, c) == (0, 0) else None))

            keep = cc[(cc["padj"] < PADJ_THRESH)
                      & (cc["log2FoldChange"].abs() > LFC_THRESH)
                      & (cc["log2FoldChange"].abs() <= X_LIM)].nsmallest(N_LABELS, "padj")
            texts = [ax.text(row["log2FoldChange"], row["y"], row["gene_name"],
                             fontsize=FS_GENELAB, style="italic")
                     for _, row in keep.iterrows()]

            ax.set_xlim(-X_LIM, X_LIM)
            ax.set_ylim(0, y_lim)
            ax.set_yticks(np.arange(0, y_lim + 1, Y_TICK_STEP))
            for x in (LFC_THRESH, -LFC_THRESH):
                ax.axvline(x, color=THRESH_COLOR, lw=1.0, linestyle="--", alpha=0.8)
            ax.axhline(-np.log10(PADJ_THRESH), color=THRESH_COLOR, lw=1.0,
                       linestyle="--", alpha=0.8)
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(labelsize=FS_TICK)
            ax.set_xlabel(r"log$_2$ Fold Change", fontsize=FS_AXIS)
            ax.set_ylabel(r"-log$_{10}$(padj)", fontsize=FS_AXIS)
            ax.set_title(f"NT vs. {pretty_label(variant)}", fontsize=FS_VTITLE)

            counts = ax.legend(handles=[up_pts, down_pts],
                               labels=[f"{int(up.sum())} DE UP",
                                       f"{int(down.sum())} DE DOWN"],
                               loc="upper right", fontsize=FS_LEGEND,
                               frameon=True, framealpha=1, edgecolor="#cccccc",
                               fancybox=True, borderpad=0.6, handletextpad=0.5)
            ax.add_artist(counts)
            if adjust_text is not None and texts:
                adjust_text(texts, ax=ax,
                            arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
                            **spread)

    fig.tight_layout(rect=(0.04, 0, 1, 0.94), h_pad=3)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper right",
               bbox_to_anchor=(1.0, 1.0), fontsize=FS_LEGEND, frameon=False,
               markerscale=1.5)
    for r, arm in enumerate(ARMS):
        bb = axes[r, 0].get_position()
        fig.text(0.01, (bb.y0 + bb.y1) / 2, TREATMENT_LABELS[arm], rotation=90,
                 ha="left", va="center", fontsize=FS_ROW)
    _save(fig, outdir, _STEM["S13g"])
    return fig, volcano_counts()


# ------------------------------------------------- S13f: the GSEA dot plot --

#: the five Hallmark sets the published dot plot draws, top to bottom
GSEA_SETS = ["E2F Targets", "G2-M Checkpoint", "Mitotic Spindle",
             "p53 Pathway", "Apoptosis"]
#: the dot plot is of the palbociclib arm
GSEA_ARM = "palbo"


def gsea_ranking(variant: str) -> pd.DataFrame:
    """The whole pre-ranked gene list for one variant, highest statistic first.

    A pre-ranked GSEA ranks on DESeq2's Wald statistic, and the enrichment
    score of a set depends on where its genes sit in the *whole* list, so this
    is the input a rerun needs. It is ~19,000 rows per variant.
    """
    d = deseq2(variant, GSEA_ARM).sort_values("stat", ascending=False)
    out = d[["gene_name", "stat", "log2FoldChange", "padj"]].reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out


def gsea_ranks() -> pd.DataFrame:
    """Where each gene set's genes sit in each variant's ranking.

    The dot plot's own numbers -- NES and FDR q-value -- are not in this
    repository, so what the source table carries is the statistic they were
    computed from: every gene of the five Hallmark sets the panel draws, its
    Wald statistic, and its rank in that variant's full list of `n_genes`.
    The full ranking itself is `gsea_ranking`, and the DESeq2 tables it reads
    ship whole; this sheet is the part of it the panel is about.
    """
    hm = hallmark()
    members = set().union(*(hm[s] for s in GSEA_SETS))
    out = []
    for variant in CONTROL_FOR:
        ranked = gsea_ranking(variant)
        t = ranked[ranked["gene_name"].isin(members)].copy()
        t.insert(0, "variant", variant)
        t.insert(1, "label", pretty_label(variant))
        t.insert(2, "vs_control", CONTROL_FOR[variant])
        t.insert(3, "treatment", TREATMENT_LABELS[GSEA_ARM])
        t["n_genes"] = len(ranked)
        for s in GSEA_SETS:
            t[f"in_{s}"] = t["gene_name"].isin(hm[s])
        out.append(t)
    return pd.concat(out, ignore_index=True)


def gsea_sets() -> pd.DataFrame:
    """The five Hallmark sets, one row per gene."""
    hm = hallmark()
    return pd.DataFrame([{"gene_set": s, "gene_name": g}
                         for s in GSEA_SETS for g in sorted(hm[s])])


# ---------------------------------------------------------------- registry --

_STEM = {
    "4m": "marker_gene_panel_DMSO_palbo_stacked",
    "S13e": "PCA",
    "S13g": "volcano_variant_vs_control",
}


def _save(fig, outdir, stem):
    if outdir is not None and stem:
        fig.savefig(str(Path(outdir) / f"{stem}.pdf"), bbox_inches="tight",
                    transparent=True)


def source(name: str) -> dict:
    """-> the frames a reader needs to redraw this panel, keyed by sheet name."""
    if name == "4m":
        z = marker_zscores()
        marker, metadata = marker_expression()
        return {"zscore_DMSO": z["veh"].reset_index(names="gene"),
                "zscore_Palbociclib": z["palbo"].reset_index(names="gene"),
                "log2_TPM": marker.reset_index(names="gene"),
                "samples": metadata.reset_index(names="sample")}
    if name == "S13e":
        pca_df, var_pct = pca()
        pca_df = pca_df.copy()
        pca_df.insert(1, "label", [pretty_label(g) for g in pca_df["genotype"]])
        variance = pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(len(var_pct))],
                                 "percent_variance": var_pct})
        return {"pca": pca_df, "variance_explained": variance}
    if name == "S13f":
        return {"ranks": gsea_ranks(), "gene_sets": gsea_sets()}
    if name == "S13g":
        return {"counts": volcano_counts(), "labelled_genes": labelled_genes()}
    raise KeyError(name)


MAIN_PANELS = {"4m": lambda: source("4m")}
MAIN_FIGURES = {"4m": panel_4m}
PANELS = {k: (lambda _p=k: source(_p)) for k in ("S13e", "S13f", "S13g")}
FIGURES = {"S13e": panel_S13e, "S13g": panel_S13g}

#: `scripts/render_figures.py` reads this to name the file a panel produces
PARTS = {k: [v] for k, v in _STEM.items()}
