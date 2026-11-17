"""Draw the lollipop panels (figures 2a, 2h, 3a).

Each wrapper holds the arguments its published panel was drawn with. The hit
set each produces is the `labelled_hit` column of the matching source-data
table, so the figure and the workbook cannot drift apart.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .. import counts as _counts
from .. import editing, heatmap, lfc, lollipop, paths
from .screens import _TABLES, _s3_hits, figure_4b as _hit_counts

#: pool each subpool's guides live in, for the gRNA denominator in the title
SUBPOOL_POOL = {"subpool1": "F1-R1", "subpool2": "F2-R2", "subpool3": "F3-R3"}

#: The denominator printed in every lollipop title was a hard-coded literal in
#: `lollipop_plots_canonical.ipynb` (cells 8-10), not a computed number.
#: subpool1 and subpool2 match this repository's targeting-guide counts exactly
#: (2,439 and 4,533), but subpool3's 1,655 is the *whole* F3-R3 pool -- it
#: includes the 79 non-targeting and 71 intron controls, where the other two
#: exclude them. The targeting count for subpool3 is 1,505.
#:
#: These literals are used so the panels reproduce the published titles. Set
#: `n_library=None` to get the computed, internally consistent value instead.
#: See docs/DISCREPANCIES.md item 12.
PUBLISHED_N_LIBRARY = {"subpool1": 2439, "subpool2": 4533, "subpool3": 1655}

#: panel -> (subpool, LFC table key, CBE condition, ABE condition, title, LFC cutoff)
PANEL_SPEC = {
    "2a": ("subpool1", "subpool1_DMSO", "KI-CDK9d-32_1000nM", "KI-CDK9d-32_1000nM",
           "KI-CDK9d-32 1,000 nM", 2),
    "2h": ("subpool1", "SY-5609_DMSO", "CBE_SY-5609_10nM", "ABE_SY-5609_10nM",
           "SY-5609 10 nM", 2),
    "3a": ("subpool2", "subpool2_DMSO", "CBE_HQ461", "ABE_HQ461", "HQ461", 2),
}


def draw(panel: str, outdir=None, **kw):
    """-> (figure, {editor: filtered frame})."""
    subpool, table, cond_cbe, cond_abe, title, lfc_cutoff = PANEL_SPEC[panel]
    abe_p, cbe_p, abe_s, cbe_s = _TABLES[table]
    edits = editing.edit_dict(sorted({abe_s, cbe_s}))
    m_cbe, m_abe, targ_cbe, targ_abe = lfc.load_merged(abe_p, cbe_p, abe_s, cbe_s, edits)
    n_lib = kw.pop("n_library", PUBLISHED_N_LIBRARY[subpool])
    if n_lib is None:                      # computed rather than as published
        lib = _counts.library()
        n_lib = int(((lib["Pool"] == SUBPOOL_POOL[subpool])
                     & (lib["classification"] == "targeting")).sum())
    fig, kept = lollipop.lollipop(
        subpool,
        targ={"CBE": targ_cbe, "ABE": targ_abe},
        merged={"CBE": m_cbe, "ABE": m_abe},
        conditions={"CBE": cond_cbe, "ABE": cond_abe},
        title=title, lfc_cutoff=lfc_cutoff, n_library=n_lib, **kw)
    if outdir is not None:
        fig.savefig(str(outdir / f"Figure_{panel}_lollipop.pdf"), dpi=300, transparent=True)
    return fig, kept


def _wrap(panel: str, spec: dict):
    """A named, documented drawer for one spec-table panel.

    The registries used to hold bare lambdas, which left `FIGURES[k].__doc__`
    empty and every generated listing of the figures blank for these panels.
    The description is composed from the spec row, so it states what the panel
    actually draws rather than repeating a title by hand.
    """
    subpool, table, cond_cbe, cond_abe, title, cutoff = spec[panel]
    f = lambda **kw: draw(panel, **kw)          # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    f.__doc__ = (f"{subpool} lollipop, {title}. CBE {cond_cbe} / ABE {cond_abe} "
                 f"against {table}, labelled at |LFC| > {cutoff}.")
    return f


FIGURES = {k: _wrap(k, PANEL_SPEC) for k in PANEL_SPEC}


# ---------------------------------------------------------------- figure 4a / 4b

#: the published per-gene heatmaps are the narrow, with-T0, large-text variant
HEATMAP_KW = dict(per_gene=True, width=6.6, include_T0=True, suffix="_large")


def draw_4a(outdir=None, genes=heatmap.GENES, **kw):
    """One clustered heatmap per gene. -> {gene: (figure, rows)}."""
    hits = _s3_hits()
    out = {}
    for g in genes:
        sub = hits[hits["Gene"] == g].reset_index(drop=True)
        base = "FDR < 0.1 & LFC >= 1 in >=1 of 6 conditions"
        fig = heatmap.cross_drug_heatmap(
            sub, f"{g}  |  {base}  |  n = {len(sub)}", **{**HEATMAP_KW, **kw})
        if outdir is not None:
            fig.savefig(str(outdir / f"Figure_4a_{g}_cross_drug_heatmap.pdf"),
                        dpi=300, transparent=True, bbox_inches="tight")
        out[g] = (fig, sub)
    return out


def draw_4b(outdir=None, **kw):
    """Resistance hits per drug by gene."""
    counts_df = _hit_counts()
    fig = heatmap.hits_per_drug(counts_df, n_hits=len(_s3_hits()), suffix="_large", **kw)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_4b_resistance_hits_per_drug.pdf"),
                    dpi=300, transparent=True, bbox_inches="tight")
    return fig, counts_df


FIGURES["4a"] = draw_4a
FIGURES["4b"] = draw_4b


# ---------------------------------------------------------------- figure 2i

#: CDK7 SY-5609 hits, from CDK7_SY5609_hits.ipynb
SY5609_CONDITIONS = {"SY-5609_10nM": "SY-5609 10 nM", "SY-5609_100nM": "SY-5609 100 nM"}
SY5609_PALETTE = ["lightcyan", "cyan"]


def _sy5609_hits(gene="CDK7", fdr=0.1, lfc_cutoff=2, min_base=100, wt_cutoff=90):
    """The hit table figure 2i draws, rebuilt from the LFC and editing tables."""
    from .. import editing as _ed
    rows = []
    for editor in ("ABE", "CBE"):
        d = lfc.io.read(paths.LFC / "barcode_counts" / "DMSO_base" /
                        f"SY_5609_{editor}_barcode_DMSO_base.csv")
        d = d[(d["Gene"] == gene) & (d["classification"] == "targeting")].copy()
        d["Editor"] = editor
        h = lfc.io.read(paths.EDITING / f"{editor}_subpool1_HGVSp_canonical.csv")
        wt = h[h["HGVSp"] == "WT"].set_index("gRNA_id")["%Reads"]
        top = h[h["HGVSp"] != "WT"].drop_duplicates("gRNA_id").set_index("gRNA_id")
        # a guide with no canonical sensor reads has no editing call, so cannot be labelled
        d = d[d["gRNA_id"].isin(h["gRNA_id"])]
        d["WT_perc"] = d["gRNA_id"].map(wt).fillna(0)
        d["Top_HGVSp"] = d["gRNA_id"].map(top["HGVSp"])
        d["Top_perc"] = d["gRNA_id"].map(top["%Reads"])
        for cond in SY5609_CONDITIONS:
            d[f"LFC_{cond}"] = d[f"LFC_median_{editor}_{cond}"]
            d[f"FDR_{cond}"] = d[f"FDR_{editor}_{cond}"]
        rows.append(d)
    data = pd.concat(rows, ignore_index=True)
    data = data[(data["base_RAW"] >= min_base) & (data["WT_perc"] < wt_cutoff)]

    enriched = np.column_stack([(data[f"FDR_{c}"] < fdr) & (data[f"LFC_{c}"] >= lfc_cutoff)
                                for c in SY5609_CONDITIONS]).any(axis=1)
    hits = data[enriched].copy()
    hits["max_LFC"] = hits[[f"LFC_{c}" for c in SY5609_CONDITIONS]].max(axis=1)
    hits = hits.sort_values("max_LFC", ascending=False).reset_index(drop=True)
    # three separate guides all call L18F, so the row key carries the gRNA id too
    hits["key"] = hits["gRNA_id"] + "_" + hits["Editor"]
    hits["label"] = hits["Top_HGVSp"]
    return hits


def _sy5609_replicates(hits):
    rows = []
    for _, v in hits.iterrows():
        for cond, name in SY5609_CONDITIONS.items():
            for rep in (1, 2, 3):
                rows.append(dict(gRNA_id=v["gRNA_id"], key=v["key"], Editor=v["Editor"],
                                 Condition=name, Replicate=rep,
                                 LFC=v[f'{v["Editor"]}_{cond}_REP{rep}']))
    return pd.DataFrame(rows)


def draw_2i(outdir=None, **kw):
    """CDK7 SY-5609 resistance hits: LFC bars beside each guide's edit logo."""
    from .. import counts as _c, logo_panel
    hits = _sy5609_hits()
    plot_df = _sy5609_replicates(hits)
    library = _c.library()
    savs = {}
    for editor in ("ABE", "CBE"):
        s = lfc.io.read(paths.EDITING / "SAVs_canonical" / f"{editor}_subpool1_canonical.csv")
        savs[editor] = s.merge(library[["gRNA_id", "Gene"]], on="gRNA_id")
    fig = logo_panel.hits_with_logos(
        hits, plot_df, savs, palette=SY5609_PALETTE,
        hue_order=list(SY5609_CONDITIONS.values()),
        title=("CDK7 SY-5609 resistance hits | FDR < 0.1 & LFC >= 2 in either dose\n"
               f"{len(hits)} gRNAs, base mean count >= 100, sensor editing > 10%"),
        **kw)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_2i_CDK7_SY5609_hits.pdf"), transparent=True)
    return fig, hits


FIGURES["2i"] = draw_2i


# ---------------------------------------------------------------- figure 3b

CDK12_13_CONDITIONS = {"HQ461": "HQ461", "BSJ-4-116": "BSJ-4-116",
                       "CDK12-IN-2": "CDK12-IN-2"}
CDK12_13_PALETTE = ["darkolivegreen", "tab:green", "lightgreen"]
CDK12_13_GENE_COLORS = {"CDK12": "black", "CDK13": "tab:blue"}


def _cdk12_13_replicates(hits):
    """Per-replicate LFC for the CDK12/13 hits, from the upstream LFC tables.

    The saved hit table keeps medians and FDRs but not the per-replicate
    columns the bars' confidence intervals come from.
    """
    reps = {}
    for editor in ("ABE", "CBE"):
        d = lfc.io.read(paths.LFC / "barcode_counts" / "DMSO_base" /
                        f"CDK12_13_{editor}_barcode_DMSO_base.csv")
        reps[editor] = d.set_index("gRNA_id")
    rows = []
    for _, v in hits.iterrows():
        src = reps[v["Editor"]]
        if v["gRNA_id"] not in src.index:
            continue
        r = src.loc[v["gRNA_id"]]
        for cond, name in CDK12_13_CONDITIONS.items():
            for rep in (1, 2, 3):
                col = f'{v["Editor"]}_{cond}_REP{rep}'
                if col in src.columns:
                    rows.append(dict(gRNA_id=v["gRNA_id"], key=v["key"],
                                     Editor=v["Editor"], Condition=name,
                                     Replicate=rep, LFC=r[col]))
    return pd.DataFrame(rows)


def draw_3b(outdir=None, top_hits=12, p_cutoff=0.01, suffix="_large", width=10.0):
    """CDK12/13: the top hits and the HQ461-specific ones, as two stacked groups."""
    from matplotlib import gridspec
    from matplotlib.patches import Patch
    import matplotlib.pyplot as plt
    from .. import counts as _c, heatmap, logo_panel

    t = lfc.io.read(paths.LFC / "CDK12_13_hits" / "CDK12_13_resistance_variants.csv")
    t["key"] = t["gRNA_id"] + "_" + t["Editor"]
    # the gene goes on its own line above the variant, as published
    t["label"] = t["Gene"] + "\n" + t["Top_HGVSp"]
    hq = (t["p_HQ461"] < p_cutoff) & (t["LFC_HQ461"] > 0) & (t["n_conditions"] == 1)
    specific = t[hq].reset_index(drop=True)
    top = t[~t["key"].isin(specific["key"])].head(top_hits).reset_index(drop=True)

    plot_df = _cdk12_13_replicates(pd.concat([top, specific]))
    pad = 0.05 * (plot_df["LFC"].max() - plot_df["LFC"].min())
    xlim = (plot_df["LFC"].min() - pad, plot_df["LFC"].max() + pad)

    library = _c.library()
    savs = {}
    for editor in ("ABE", "CBE"):
        s = lfc.io.read(paths.EDITING / "SAVs_canonical" / f"{editor}_subpool2_canonical.csv")
        savs[editor] = s.merge(library[["gRNA_id", "Gene"]], on="gRNA_id")

    st = heatmap.STYLES[suffix]
    panels = [(f"Top {top_hits} resistance hits (excluding HQ461-specific)", top),
              ("HQ461-specific resistance hits (resistant in HQ461 only)", specific)]
    # a spacer row between the two groups, so their bars do not touch
    GAP, TITLE_GAP = 1.0, 0.9
    kinds, starts = [], {}
    for pi, (_, sub) in enumerate(panels):
        kinds += ["gap"] * (1 if pi else 0)
        starts[pi] = len(kinds)
        kinds += ["row"] * len(sub)
    heights = [GAP if k == "gap" else 1.0 for k in kinds]

    fig = plt.figure(figsize=(width, sum(heights) + TITLE_GAP * len(panels) + 1.2))
    gs = gridspec.GridSpec(len(kinds), 2, width_ratios=[1.5, 1],
                           height_ratios=heights, wspace=0.15, hspace=0.55)
    for pi, (title, sub) in enumerate(panels):
        axb = logo_panel.logo_group(
            gs, starts[pi], sub, plot_df, savs, CDK12_13_PALETTE,
            list(CDK12_13_CONDITIONS.values()),
            label_colors=CDK12_13_GENE_COLORS, color_by="Gene", xlim=xlim, sizes=st)
        axb.set_title(f"{title}  (n = {len(sub)})", fontsize=st["panel_title"],
                      loc="left", pad=10)

    H = fig.get_figheight()
    handles = [Patch(facecolor=c, edgecolor="black", label=n)
               for c, n in zip(CDK12_13_PALETTE, CDK12_13_CONDITIONS.values())]
    fig.legend(handles=handles, fontsize=st["legend"], ncol=3, frameon=False,
               loc="upper center", bbox_to_anchor=(0.5, 1 - 0.52 / H))
    fig.suptitle(
        f"CDK12/13 resistance variants | top {top_hits} broader hits and the "
        f"HQ461-specific hits\nresistant = empirical p < {p_cutoff} & LFC > 0  |  "
        "CDK13 labels in blue", fontsize=st["title"], y=1 - 0.08 / H)
    fig.subplots_adjust(left=0.20, right=0.99, top=1 - 1.35 / H, bottom=0.80 / H)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_3b_CDK12_13_hits.pdf"), transparent=True)
    return fig, pd.concat([top, specific], ignore_index=True)


FIGURES["3b"] = draw_3b


# ---------------------------------------------------------------- figure 2b

CDK9_CONDITIONS = {"KI-CDK9d-32_100nM": "KI-CDK9d-32 100nM",
                   "KI-CDK9d-32_1000nM": "KI-CDK9d-32 1,000nM",
                   "KI-CDK9d-32N_1250nM": "KI-CDK9d-32N 1,250nM",
                   "KI-CDK9d-32N_5000nM": "KI-CDK9d-32N 5,000nM",
                   "KB-0742_1500nM": "KB-0742 1,500nM"}
CDK9_PALETTE = ["indigo", "darkviolet", "mediumorchid", "thistle", "magenta"]
#: the dose the hits are ranked by
CDK9_RANK_BY = "KI-CDK9d-32N_5000nM"


def _cdk9_top_hits(n=10, fdr=0.1, gene="CDK9", min_perc=2, canonical=True):
    """The CDK9 variants figure 2b shows, ranked by the strongest dose.

    `canonical` picks which editing tables supply the labels and logos.
    The default is the canonical tables, which keep only the editor's own edit
    and are the standard everywhere else in this repository.

    **The published panel used the non-canonical ones.** `logo_plots.ipynb`
    predates canonical filtering and reads `SAVs/` and `*_HGVSp`, which give
    slightly wider logo windows -- a bystander C>T on an ABE guide still counts
    there. Pass `canonical=False` to reproduce the published panel exactly.
    See docs/DISCREPANCIES.md.
    """
    rows = []
    for editor in ("ABE", "CBE"):
        d = lfc.io.read(paths.LFC / "barcode_counts" / "DMSO_base" /
                        f"{editor}_subpool1_barcode_DMSO_base.csv")
        d["Editor"] = editor
        rows.append(d)
    data = pd.concat(rows)
    sel = data[(data[f"FDR_{CDK9_RANK_BY}"] < fdr)
               & (data["classification"] == "targeting")]
    sel = sel.sort_values(f"LFC_median_{CDK9_RANK_BY}", ascending=False)
    sel = sel[sel["Gene"] == gene].head(n).reset_index(drop=True)
    sel["key"] = sel["gRNA_id"] + "_" + sel["Editor"]

    # the label is the guide's top variant call, from the canonical HGVSp table
    labels = []
    for _, v in sel.iterrows():
        suffix = "_HGVSp_canonical" if canonical else "_HGVSp"
        h = lfc.io.read(paths.EDITING / f'{v["Editor"]}_subpool1{suffix}.csv')
        s = h[(h["gRNA_id"] == v["gRNA_id"]) & (h["HGVSp"] != "WT")
              & (h["%Reads"] >= min_perc)]
        labels.append(s["HGVSp"].iloc[0] if len(s) else "")
    sel["label"] = labels
    return sel


def _cdk9_replicates(hits):
    rows = []
    for _, v in hits.iterrows():
        for cond, name in CDK9_CONDITIONS.items():
            for rep in (1, 2, 3):
                col = f"{cond}_REP{rep}"
                if col in hits.columns:
                    rows.append(dict(gRNA_id=v["gRNA_id"], key=v["key"],
                                     Editor=v["Editor"], Condition=name,
                                     Replicate=rep, LFC=v[col]))
    return pd.DataFrame(rows)


def draw_2b(outdir=None, n=10, canonical=True, **kw):
    """CDK9 top resistance variants across the five CDK9-inhibitor arms.

    See `_cdk9_top_hits` for what `canonical` changes. The default uses the
    canonical editing tables, like the rest of the repository; the published
    panel used the non-canonical ones.
    """
    from .. import counts as _c, logo_panel
    hits = _cdk9_top_hits(n=n, canonical=canonical)
    plot_df = _cdk9_replicates(hits)
    library = _c.library()
    savs = {}
    for editor in ("ABE", "CBE"):
        sav_path = (paths.EDITING / "SAVs_canonical" / f"{editor}_subpool1_canonical.csv"
                    if canonical else paths.EDITING / "SAVs" / f"{editor}_subpool1.csv")
        s = lfc.io.read(sav_path)
        savs[editor] = s.merge(library[["gRNA_id", "Gene"]], on="gRNA_id")
    fig = logo_panel.hits_with_logos(
        hits, plot_df, savs, palette=CDK9_PALETTE,
        hue_order=list(CDK9_CONDITIONS.values()), xlabel="LFC",
        width=7.0, row_h=0.8, legend=False, **kw)
    if outdir is not None:
        # the figure is 7in but the published PDF measures ~8in: saved tight, so
        # the long compound labels (D104N_C106Y) are not clipped
        fig.savefig(str(outdir / "Figure_2b_CDK9_top_variants.pdf"),
                    transparent=True, bbox_inches="tight")
    return fig, hits


FIGURES["2b"] = draw_2b


# ---------------------------------------------------------------- figure 4p

CDK46_CONDITIONS = {"RIBO": "Ribociclib", "ABEMA": "Abemaciclib", "PALBO": "Palbociclib"}
CDK46_PALETTE = ["lightcoral", "indianred", "maroon"]
#: the CDK6 residues reported as clinical resistance variants
CLINICAL_CODONS = {"CBE": {"168", "46"}, "ABE": {"168", "46", "47"}}
CLINICAL_MIN_EDIT = 10      # % of sensor reads for a guide to be included at all
#: the six drug arms of the CDK2/4/6 screen. The pre-selection below asks for
#: significance in at least one of them -- not only in the three CDK4/6
#: inhibitors the panel draws, which would drop R46W (its significant arm is
#: one of the other three).
CLINICAL_CONDITIONS = ["PALBO", "RIBO", "ABEMA", "ATIRMO", "INX", "TAGTO"]
#: A guide reaches the panel only if it is a resistance hit in at least one of
#: those arms. 0.1 is this repository's threshold for a resistance hit
#: throughout (it is written on the face of supplementary 8c and 9b), and it is
#: what reproduces the published panel: of the 14 guides editing a clinically
#: reported codon, exactly five clear it, one per variant, and they are the
#: five the manuscript shows, in the order it shows them. Raising it to 0.05
#: leaves none. See `docs/DISCREPANCIES.md` item 5b.
CLINICAL_MAX_FDR = 0.1


def _clinical_guides(canonical=True):
    """CDK6 guides editing the clinically reported residues, per editor."""
    from .. import counts as _c
    library = _c.library()
    out = {}
    for editor, codons in CLINICAL_CODONS.items():
        path = (paths.EDITING / "SAVs_canonical" / f"{editor}_subpool3_canonical.csv"
                if canonical else paths.EDITING / "SAVs" / f"{editor}_subpool3.csv")
        s = lfc.io.read(path).merge(library[["gRNA_id", "Gene"]], on="gRNA_id")
        sel = s[(s["Codon"].astype(str).isin(codons)) & (s["Gene"] == "CDK6")
                & (s["%Reads"] >= CLINICAL_MIN_EDIT)]
        out[editor] = (sel.sort_values("%Reads", ascending=False)["gRNA_id"]
                       .astype(str) + f"_{editor}").tolist()
    return out


def _clinical_hits(n=15, canonical=True, min_perc=1, max_fdr=CLINICAL_MAX_FDR):
    """The guides the clinical-variants panel draws, in the order it draws them.

    Two filters, then a sort. A guide has to edit a clinically reported CDK6
    codon at `CLINICAL_MIN_EDIT`% of sensor reads, and it has to be a
    resistance hit -- FDR below `max_fdr` in at least one of the screen's six
    drug arms. What is left is ordered by palbociclib LFC.

    `n` is a cap carried over from the notebook, which took the top 15 by LFC
    and had no significance filter. With the filter in place it does not bind:
    five guides clear it.
    """
    guides = _clinical_guides(canonical)
    wanted = set(guides["CBE"]) | set(guides["ABE"])
    rows = []
    for editor in ("ABE", "CBE"):
        d = lfc.io.read(paths.LFC / "barcode_counts" / "DMSO_base" /
                        f"CDK2_4_6_{editor}_barcode_DMSO_base.csv")
        d = d.rename(columns={f"LFC_median_{editor}_PALBO": "LFC_median_PALBO"})
        # the FDR columns carry the editor in their name; drop it so one
        # concatenated frame can be filtered on a single set of columns
        d = d.rename(columns={f"FDR_{editor}_{c}": f"FDR_{c}"
                              for c in CLINICAL_CONDITIONS})
        d["Editor"] = editor
        rows.append(d)
    data = pd.concat(rows)
    data["key"] = data["gRNA_id"] + "_" + data["Editor"]
    sel = data[data["key"].isin(wanted)].copy()
    fdr = [f"FDR_{c}" for c in CLINICAL_CONDITIONS]
    sel["best_FDR"] = sel[fdr].min(axis=1)
    sel = sel[sel["best_FDR"] < max_fdr]
    sel = sel.sort_values("LFC_median_PALBO", ascending=False)
    sel = sel.head(n).reset_index(drop=True)

    suffix = "_HGVSp_canonical" if canonical else "_HGVSp"
    labels = []
    for _, v in sel.iterrows():
        h = lfc.io.read(paths.EDITING / f'{v["Editor"]}_CDK2_4_6{suffix}.csv')
        s = h[(h["gRNA_id"] == v["gRNA_id"]) & (h["HGVSp"] != "WT")
              & (h["%Reads"] >= min_perc)]
        labels.append(f'CDK6\n{s["HGVSp"].iloc[0]}' if len(s) else "CDK6")
    sel["label"] = labels
    return sel


def _clinical_replicates(hits):
    rows = []
    for _, v in hits.iterrows():
        for cond, name in CDK46_CONDITIONS.items():
            for rep in (1, 2, 3):
                col = f'{v["Editor"]}_{cond}_REP{rep}'
                if col in hits.columns:
                    rows.append(dict(gRNA_id=v["gRNA_id"], key=v["key"],
                                     Editor=v["Editor"], Condition=name,
                                     Replicate=rep, LFC=v[col]))
    return pd.DataFrame(rows)


def draw_4p(outdir=None, n=15, canonical=True, min_perc=1, **kw):
    """CDK6 clinical resistance variants across the three CDK4/6 inhibitors."""
    from .. import counts as _c, logo_panel
    hits = _clinical_hits(n=n, canonical=canonical, min_perc=min_perc)
    plot_df = _clinical_replicates(hits)
    library = _c.library()
    savs = {}
    for editor in ("ABE", "CBE"):
        path = (paths.EDITING / "SAVs_canonical" / f"{editor}_subpool3_canonical.csv"
                if canonical else paths.EDITING / "SAVs" / f"{editor}_subpool3.csv")
        savs[editor] = lfc.io.read(path).merge(library[["gRNA_id", "Gene"]], on="gRNA_id")
    fig = logo_panel.hits_with_logos(
        hits, plot_df, savs, palette=CDK46_PALETTE,
        hue_order=list(CDK46_CONDITIONS.values()), xlabel="LFC",
        width=7.0, row_h=0.6, legend=False, min_perc=min_perc, **kw)
    if outdir is not None:
        fig.savefig(str(outdir / "Figure_4p_CDK6_clinical_variants.pdf"),
                    transparent=True, bbox_inches="tight")
    return fig, hits


FIGURES["4p"] = draw_4p
