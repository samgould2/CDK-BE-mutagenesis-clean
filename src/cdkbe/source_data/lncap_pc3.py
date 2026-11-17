"""The LNCaP and PC3 CDK12 panels: supplementary 11c, 11D-11H and 11K.

These seven panels were produced outside this repository and came as complete
packages: the values that were plotted, the script that plotted them and the
published image, one folder per panel, all under
`data/raw/collaborator_panels/`. This module redraws each panel from those
values, so the figure is reproducible here rather than only archived here.

    11c  qPCR of four DNA-repair genes over a CDK12-inhibition time course
    11D  PCA of the LNCaP RNA-seq, from nf-core/rnaseq's own DESeq2 QC step
    11E  the four LNCaP DESeq2 contrasts as volcanoes
    11F  leading-edge heatmap of HALLMARK_ANDROGEN_RESPONSE in the NT arm
    11G  preranked GSEA of GO:0006302, double-strand break repair
    11H  amino-acid logo plots of the CDK12-edited PC3 pools
    11K  the three PC3 intrinsic (untreated) contrasts as volcanoes

**Nothing here recomputes a statistic.** The DESeq2 fits, the GSEA runs and the
qPCR ddCq are the delivered authors' own, and their scripts ship beside the
data. What this module does is read the delivered values back and draw them
with the same parameters those scripts used, which is what makes the panels
reproducible from this repository alone. Two corrections the authors record in
their own READMEs are carried with the data rather than re-litigated here: a
`lfc_shrink` reference-level bug in the PC3 contrasts, and a reverse-strand
CIGAR bug in the logo-plot pileup. `docs/DISCREPANCIES.md` item 25 has both.

The panel letters are not guesswork either. The delivered folders are named
`RNAseq_panel_d` through `RNAseq_panel_g` after the panels they make, and the
deposition README puts the qPCR heatmap before them.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import paths

BASE = paths.RAW / "collaborator_panels"
LNCAP = BASE / "LNCaP"
PC3 = BASE / "PC3"

#: panel -> the delivered folder it is drawn from
FOLDER = {
    "S11c": LNCAP / "qPCR_heatmap",
    "S11D": LNCAP / "RNAseq_panel_d",
    "S11E": LNCAP / "RNAseq_panel_e",
    "S11F": LNCAP / "RNAseq_panel_f",
    "S11G": LNCAP / "RNAseq_panel_g",
    "S11H": PC3 / "Cdk12_logo_plots",
    "S11K": PC3 / "PC3_panel_k",
}


def _sections(path: Path) -> list:
    """Split a delivered `*_source_data.txt` into its `## Section` blocks."""
    blocks, current = [], []
    for line in path.read_text().splitlines():
        if line.startswith("## "):
            if current:
                blocks.append(current)
            current = []
            continue
        if line.startswith("#") or not line.strip():
            continue
        current.append(line)
    if current:
        blocks.append(current)
    return blocks


def _table(lines: list) -> pd.DataFrame:
    head = lines[0].split("\t")
    rows = [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]
    d = pd.DataFrame(rows)
    for c in d.columns:
        # "NA" marks the reference row of a comparison, not a missing number,
        # so it has to become NaN before the column will convert at all --
        # otherwise the whole column stays text and every p-value compares as
        # a string.
        col = d[c].replace("NA", np.nan)
        try:
            d[c] = pd.to_numeric(col)
        except (TypeError, ValueError):
            d[c] = col
    return d


# -------------------------------------------------- 11c: qPCR time course ---

QPCR_GENES = ["ATR", "BRCA1", "BRCA2", "ERCC4"]
QPCR_TIMEPOINTS = [0, 4, 8, 12]


def qpcr_replicates() -> pd.DataFrame:
    """Section 1 of the delivered file: dCq, ddCq and fold change per well."""
    return _table(_sections(
        FOLDER["S11c"] / "LNCaP_CDK12_qPCR_DNArepair_heatmap_source_data.txt")[0])


def qpcr_summary() -> pd.DataFrame:
    """Section 2: the mean, SEM and Bonferroni-corrected p the panel prints."""
    return _table(_sections(
        FOLDER["S11c"] / "LNCaP_CDK12_qPCR_DNArepair_heatmap_source_data.txt")[1])


def supp_11c(outdir=None, **_):
    """-> (figure, the summary table). Genes down the rows, hours across."""
    import seaborn as sns

    s = qpcr_summary().set_index(["Gene", "Timepoint_hr"])
    fc = np.array([[s.loc[(g, t), "Mean_FoldChange"] for t in QPCR_TIMEPOINTS]
                   for g in QPCR_GENES])
    sem = np.array([[s.loc[(g, t), "SEM"] for t in QPCR_TIMEPOINTS]
                    for g in QPCR_GENES])
    fig, ax = plt.subplots(figsize=(2.8, 3.2))
    sns.heatmap(fc, ax=ax, cmap="RdBu_r", center=1.0, vmin=0.0, vmax=2.0,
                annot=False, linewidths=0.8, linecolor="black", cbar=True,
                cbar_kws={"shrink": 0.7, "pad": 0.02},
                xticklabels=[f"{t} hr" for t in QPCR_TIMEPOINTS],
                yticklabels=QPCR_GENES)
    for i, gene in enumerate(QPCR_GENES):
        for j, t in enumerate(QPCR_TIMEPOINTS):
            star = s.loc[(gene, t), "Significance"]
            text = f"{fc[i, j]:.2f}\n±{sem[i, j]:.2f}"
            if isinstance(star, str) and star not in ("ns", "NA"):
                text += f"\n{star}"
            # the delivered script picks the text colour from the value, not
            # from the cell's luminance as seaborn's own `annot` would
            colour = "white" if (fc[i, j] < 0.55 or fc[i, j] > 1.6) else "black"
            ax.text(j + 0.5, i + 0.5, text, ha="center", va="center",
                    fontsize=6.5, color=colour)
    ax.set_xticklabels([f"{t} hr" for t in QPCR_TIMEPOINTS], fontsize=8,
                       rotation=90)
    ax.set_yticklabels(QPCR_GENES, fontsize=8, rotation=90, ha="center",
                       va="bottom")
    ax.tick_params(axis="both", length=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    cbar = ax.collections[0].colorbar
    cbar.set_label("Fold Change", fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.4)
    _save(fig, outdir, "LNCaP_CDK12_qPCR_DNArepair_heatmap")
    return fig, qpcr_summary()


# ------------------------------------------------------------- 11D: PCA ----

PCA_COLOURS = {"NT_Control": "#2166ac", "NT_HQ461": "#92c5de",
               "I733V_Control": "#d6604d", "I733V_HQ461": "#f4a582"}
PCA_LABELS = {"NT_Control": "NT Control", "NT_HQ461": "NT + HQ461",
              "I733V_Control": "I733V Control", "I733V_HQ461": "I733V + HQ461"}


def pca_values() -> pd.DataFrame:
    """PC1/PC2 per sample, as nf-core/rnaseq's DESeq2 QC step wrote them."""
    d = pd.read_csv(FOLDER["S11D"] / "deseq2.pca.vals.txt", sep="\t")
    d = d.rename(columns={d.columns[0]: "sample"})
    d["sample"] = d["sample"].astype(str)
    d["condition"] = [next((c for c in PCA_COLOURS if c in s), "other")
                      for s in d["sample"]]
    return d


def supp_11D(outdir=None, **_):
    """-> (figure, the PC1/PC2 table)."""
    d = pca_values()
    pc1, pc2 = [c for c in d.columns if c.startswith("PC")][:2]
    fig, ax = plt.subplots(figsize=(5, 4.5))
    for cond, grp in d.groupby("condition"):
        ax.scatter(grp[pc1], grp[pc2], label=PCA_LABELS.get(cond, cond),
                   color=PCA_COLOURS.get(cond, "gray"), s=80, zorder=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(pc1, fontsize=11)
    ax.set_ylabel(pc2, fontsize=11)
    ax.set_title("PCA — LNCaP CDK12 RNA-seq", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    _save(fig, outdir, "PCA")
    return fig, d


# ------------------------------------------------------- 11E, 11K: volcanoes --

LFC_THRESH, PADJ_THRESH = 1.0, 0.05
#: pydeseq2 can return padj of 0; the PC3 script floors it so -log10 is finite
PADJ_FLOOR = 1e-300

#: panel -> [(csv stem, panel title), ...] in the order the published figure
#: draws them, taken from the `comparisons`/`COMPS` list in each script
VOLCANOES = {
    "S11E": [
        ("NT_HQ461_vs_NT_Control", "NT: HQ461 vs Control"),
        ("I733V_HQ461_vs_I733V_Control", "Cdk12 I733V: HQ461 vs Control"),
        ("NT_HQ461_vs_I733V_HQ461", "HQ461: NT vs Cdk12 I733V"),
        ("NT_Control_vs_I733V_Control", "Control: NT vs Cdk12 I733V"),
    ],
    "S11K": [
        ("NT_vs_I733V_Control", "Control: NT vs I733V"),
        ("NT_vs_G822E_Control", "Control: NT vs G822E"),
        ("NT_vs_G731K_Control", "Control: NT vs G731K"),
    ],
}
VOLCANO_FIGSIZE = {"S11E": (20, 5), "S11K": (16, 5)}
VOLCANO_FONT = {"S11E": (12, 13), "S11K": (10, 12)}


def deseq2(panel: str, stem: str) -> pd.DataFrame:
    """One delivered DESeq2 contrast."""
    return pd.read_csv(FOLDER[panel] / f"{stem}.csv")


def volcano_counts(panel: str) -> pd.DataFrame:
    """The DE call each volcano's legend prints: padj < 0.05 and |LFC| > 1."""
    rows = []
    for stem, title in VOLCANOES[panel]:
        d = deseq2(panel, stem).dropna(subset=["padj", "log2FoldChange"])
        up = ((d["padj"] < PADJ_THRESH) & (d["log2FoldChange"] > LFC_THRESH)).sum()
        down = ((d["padj"] < PADJ_THRESH)
                & (d["log2FoldChange"] < -LFC_THRESH)).sum()
        rows.append({"panel": panel, "comparison": stem, "title": title,
                     "tested": len(d), "up": int(up), "down": int(down)})
    return pd.DataFrame(rows)


def _volcano(panel: str, outdir=None, **_):
    specs = VOLCANOES[panel]
    frames = {stem: deseq2(panel, stem).dropna(
        subset=["padj", "log2FoldChange"]).copy() for stem, _ in specs}
    if panel == "S11K":
        for d in frames.values():
            d["padj"] = d["padj"].clip(lower=PADJ_FLOOR)
    all_lfc = np.concatenate([d["log2FoldChange"].to_numpy() for d in frames.values()])
    all_y = np.concatenate([-np.log10(d["padj"].to_numpy()) for d in frames.values()])
    x_lim = np.ceil(max(abs(np.nanmin(all_lfc)), abs(np.nanmax(all_lfc)))) + 0.5
    y_lim = np.ceil(np.nanmax(all_y)) + 2
    tick_fs, label_fs = VOLCANO_FONT[panel]

    fig, axes = plt.subplots(1, len(specs), figsize=VOLCANO_FIGSIZE[panel])
    for ax, (stem, title) in zip(np.atleast_1d(axes), specs):
        d = frames[stem]
        y = -np.log10(d["padj"])
        up = (d["padj"] < PADJ_THRESH) & (d["log2FoldChange"] > LFC_THRESH)
        down = (d["padj"] < PADJ_THRESH) & (d["log2FoldChange"] < -LFC_THRESH)
        ns = ~(up | down)
        ax.scatter(d.loc[ns, "log2FoldChange"], y[ns], s=4, color="#AAAAAA",
                   alpha=0.4, linewidths=0, rasterized=True)
        tag = "DE UP" if panel == "S11E" else "UP"
        ax.scatter(d.loc[up, "log2FoldChange"], y[up], s=25, color="tab:red",
                   alpha=0.7, linewidths=0.3, edgecolors="black",
                   rasterized=True, label=f"{int(up.sum())} {tag}")
        tag = "DE DOWN" if panel == "S11E" else "DOWN"
        ax.scatter(d.loc[down, "log2FoldChange"], y[down], s=25,
                   color="tab:blue", alpha=0.7, linewidths=0.3,
                   edgecolors="black", rasterized=True,
                   label=f"{int(down.sum())} {tag}")
        ax.set_xlim(-x_lim, x_lim)
        ax.set_ylim(0, y_lim)
        ax.legend(fontsize=tick_fs, frameon=True, framealpha=0.9,
                  edgecolor="#CCCCCC")
        for v in (LFC_THRESH, -LFC_THRESH):
            ax.axvline(v, color="black", lw=0.5, ls="--", alpha=0.4)
        ax.axhline(-np.log10(PADJ_THRESH), color="black", lw=0.5, ls="--",
                   alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", which="major", labelsize=tick_fs)
        ax.set_xlabel("log$_2$ Fold Change", fontsize=label_fs)
        ax.set_ylabel("-log$_{10}$(padj)", fontsize=label_fs)
        ax.set_title(title, fontsize=label_fs)
    fig.tight_layout()
    _save(fig, outdir, _STEM[panel])
    return fig, volcano_counts(panel)


# ------------------------------------------- 11F: androgen leading edge -----

def androgen_table() -> pd.DataFrame:
    """Leading-edge genes with their per-sample z-scores, log2FC and padj."""
    path = (FOLDER["S11F"]
            / "LNCaP_NT_Hallmark_androgen_leadingedge_heatmap_source_data.txt")
    return _table(_sections(path)[0])


def supp_11F(outdir=None, **_):
    """-> (figure, the leading-edge table)."""
    import seaborn as sns

    d = androgen_table()
    genes = d["gene"].tolist()
    z = d[[c for c in d.columns if c.endswith("_zscore")]].to_numpy(float)
    fig, ax = plt.subplots(figsize=(2.6, 0.17 * len(genes) + 0.9))
    sns.heatmap(z, ax=ax, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                linewidths=0.8, linecolor="black", cbar=True,
                cbar_kws={"shrink": 0.4, "pad": 0.04}, yticklabels=genes)
    ax.set_xticks([1.5, 4.5])
    ax.set_xticklabels(["Control", "HQ461"], fontsize=8, rotation=0)
    ax.set_yticklabels(genes, fontsize=8, rotation=0, style="italic")
    ax.tick_params(length=0)
    ax.set_title("Hallmark Androgen Response\n(leading edge, NT)", fontsize=8)
    cbar = ax.collections[0].colorbar
    cbar.set_label("Row z-score", fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.4)
    _save(fig, outdir, "LNCaP_NT_Hallmark_androgen_leadingedge_heatmap")
    return fig, d


# -------------------------------------------------------- 11G: GSEA run ----

def gsea_table() -> pd.DataFrame:
    """The ranked list with set membership and the running enrichment score."""
    path = (FOLDER["S11G"]
            / "LNCaP_GSEA_GO_DSB_repair_HQ461_NT_vs_I733V_source_data.txt")
    return _table(_sections(path)[0])


def gsea_statistics() -> dict:
    """NES and the nominal p, read from the delivered file's own header."""
    path = (FOLDER["S11G"]
            / "LNCaP_GSEA_GO_DSB_repair_HQ461_NT_vs_I733V_source_data.txt")
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if "NES =" in line:
            for part in line.lstrip("# ").split(";"):
                k, _, v = part.partition("=")
                out[k.strip()] = float(v)
    return out


def supp_11G(outdir=None, **_):
    """-> (figure, the ranked list)."""
    d = gsea_table()
    stats = gsea_statistics()
    x = np.arange(len(d))
    es = d["running_ES"].to_numpy(float)
    rnk = d["wald_stat"].to_numpy(float)
    hits = x[d["in_gene_set"].to_numpy(int) == 1]

    fig, axes = plt.subplots(3, 1, figsize=(3.0, 2.8), sharex=True,
                             gridspec_kw={"height_ratios": [3, 0.5, 1.3],
                                          "hspace": 0.05})
    ax = axes[0]
    ax.plot(x, es, color="#2166ac", lw=1.5)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylabel("Enrichment score", fontsize=7)
    ax.set_title("Double-strand break repair\nHQ461: NT vs Cdk12 I733V",
                 fontsize=8)
    ax.text(0.03, 0.05,
            f"NES = {stats['NES']:.2f}\np = {stats['nominal p']:.3f}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=7)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    ax.vlines(hits, 0, 1, color="black", lw=0.6)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)

    ax = axes[2]
    ax.fill_between(x, rnk, 0, color="grey", lw=0)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xlim(0, len(x))
    ax.set_ylabel("Wald stat", fontsize=7)
    ax.set_xlabel("Gene rank", fontsize=7)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.01, 0.95, "higher in NT", transform=ax.transAxes, fontsize=6.5,
            va="top")
    ax.text(0.99, 0.05, "higher in I733V", transform=ax.transAxes,
            fontsize=6.5, ha="right")
    _save(fig, outdir, "LNCaP_GSEA_GO_DSB_repair_HQ461_NT_vs_I733V")
    return fig, d


# --------------------------------------------------- 11H: CDK12 aa logos ---

#: the wild-type residue at every position plotted, from the pileup script.
#: Anything else above 1% of reads is drawn in red.
LOGO_WT = {731: "G", 732: "I", 733: "I",
           816: "M", 817: "D", 818: "H", 819: "D", 820: "L", 821: "M", 822: "G"}
LOGO_ORDER = ["NT control (amp1)", "G731K gRNA", "I733V gRNA",
              "NT control (amp2)", "G822E/K gRNA"]
LOGO_TITLES = {"NT control (amp1)": "NT control",
               "NT control (amp2)": "NT control",
               "G731K gRNA": "G731K gRNA", "I733V gRNA": "I733V gRNA",
               "G822E/K gRNA": "G822E/K gRNA"}
LOGO_MIN_PCT = 1.0


def logo_table() -> pd.DataFrame:
    """Per-sample, per-residue amino-acid percentages from the read pileup."""
    return pd.read_csv(FOLDER["S11H"] / "output" / "aa_frequency_table.csv")


def supp_11H(outdir=None, **_):
    """-> (figure, the frequency table)."""
    import logomaker

    d = logo_table()
    fig, axes = plt.subplots(len(LOGO_ORDER), 1,
                             figsize=(6, 2.3 * len(LOGO_ORDER)))
    for ax, sample in zip(axes, LOGO_ORDER):
        sub = (d[d["sample"] == sample].drop(columns=["sample"])
               .set_index("residue").sort_index())
        logo = logomaker.Logo(sub, vpad=0.02, ax=ax)
        logo.style_xticks(anchor=0, spacing=1, rotation=0, fontsize=11)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels([0, 25, 50, 75, 100], fontsize=10)
        ax.set_ylabel("% reads", fontsize=11)
        ax.set_xlabel("CDK12 residue", fontsize=11)
        ax.set_title(LOGO_TITLES[sample], fontsize=12, loc="left")
        for residue in sub.index:
            wt = LOGO_WT.get(residue)
            for aa in sub.columns:
                if aa != wt and sub.loc[residue, aa] >= LOGO_MIN_PCT:
                    logo.style_single_glyph(c=aa, p=residue, color="tab:red")
        logo.style_spines(spines=["right", "top"], visible=False)
    fig.tight_layout()
    _save(fig, outdir, "CDK12_aa_logo_plots")
    return fig, d


# --------------------------------------------------------------- registry ---

_STEM = {
    "S11c": "LNCaP_CDK12_qPCR_DNArepair_heatmap",
    "S11D": "PCA",
    "S11E": "volcano_all_comparisons",
    "S11F": "LNCaP_NT_Hallmark_androgen_leadingedge_heatmap",
    "S11G": "LNCaP_GSEA_GO_DSB_repair_HQ461_NT_vs_I733V",
    "S11H": "CDK12_aa_logo_plots",
    "S11K": "PC3_panel_k_volcano_intrinsic",
}


def _save(fig, outdir, stem):
    if outdir is not None:
        fig.savefig(str(Path(outdir) / f"{stem}.pdf"), bbox_inches="tight",
                    transparent=True)


def supp_11E(outdir=None, **kw):
    """The four LNCaP contrasts, drawn side by side as the published panel is."""
    return _volcano("S11E", outdir=outdir, **kw)


def supp_11K(outdir=None, **kw):
    """The three PC3 intrinsic contrasts, side by side."""
    return _volcano("S11K", outdir=outdir, **kw)


PANELS = {
    "S11c": qpcr_summary,
    "S11D": pca_values,
    "S11E": lambda: volcano_counts("S11E"),
    "S11F": androgen_table,
    "S11G": gsea_table,
    "S11H": logo_table,
    "S11K": lambda: volcano_counts("S11K"),
}

FIGURES = {
    "S11c": supp_11c, "S11D": supp_11D, "S11E": supp_11E, "S11F": supp_11F,
    "S11G": supp_11G, "S11H": supp_11H, "S11K": supp_11K,
}

#: `scripts/render_figures.py` reads this to name the file a panel produces
PARTS = {k: [v] for k, v in _STEM.items()}
