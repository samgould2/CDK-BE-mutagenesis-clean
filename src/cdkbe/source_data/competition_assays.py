"""Flow competition assays: figure 2f and supplementary 7c, 7e, 10f, 10g, 11b, 11i.

A competition assay mixes edited (fluorescent) cells with unedited ones and
reads the edited fraction by flow at two timepoints. The number each panel
shows is that fraction at the end of treatment divided by the same fraction at
T0 -- a fold change above 1 means the variant was enriched under drug, which is
the resistance phenotype the screens predicted.

`data/raw/competition_assays/` holds the delivered fold-change tables, three
replicates per condition, and `README.txt` there is the collaborator's own
file-by-file guide to them. Nothing in this module re-derives a fold change:
the replicate values are the rawest form of this experiment that was delivered
(the FlowJo exports they come from stay in the lab's Dropbox, and the
normalisation between the two is not written down anywhere -- see
`docs/DISCREPANCIES.md` item 23). Everything else on the panel -- the mean, the
SEM, the significance stars and the colour -- is computed here.

The five heatmap panels are reproduced against the published SVGs, which are
vendored under `tests/reference/competition_assays/`. Every cell's mean, SEM,
star and fill colour is pinned to those files, so the panel cannot drift from
what the manuscript shows without a test failing.

Two panels are not heatmaps. Supplementary 7c and 7e are the mouse B-ALL
competition assays, and this repository has their three replicates but not the
published figure, so the numbers are exact while the panel style is this
repository's own convention rather than a reproduction of the original.
`docs/DISCREPANCIES.md` item 24 says so.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns

from .. import paths

ASSAYS = paths.RAW / "competition_assays"

#: the colour scale every published heatmap in this family uses. Read back out
#: of the vendored SVGs: the colourbar runs 0 to 5.5 with evenly spaced ticks
#: (so a plain linear norm, not a two-slope one), and the colours it draws are
#: RdBu_r truncated to 0.389-1.0, which is exactly what seaborn's `center=1`
#: does to that colormap at this vmin/vmax.
CMAP, VMIN, VMAX, CENTER = "RdBu_r", 0.0, 5.5, 1.0
CBAR_LABEL = "Fold Change\n(Relative to T0)"

#: seaborn picks the annotation colour from the cell's luminance; these are the
#: two values it uses, and the published SVGs carry both.
DARK_TEXT, LIGHT_TEXT = "#262626", "#ffffff"


def _stars(p: float) -> str:
    """The convention the delivered PC3 table states in its own header."""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ---------------------------------------------------------------- reading ---

def matrix_file(name: str) -> pd.DataFrame:
    """One delivered fold-change table, as long form.

    The files are a small matrix: a header naming the treatment once per
    replicate, then one row per variant. Quoting is inconsistent between them
    ("Control " with a trailing space, bare `HQ461`), so both the variant and
    the treatment are stripped of quotes and whitespace here.

    -> columns variant, treatment, replicate, fold_change
    """
    path = ASSAYS / name
    lines = [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.strip()]
    header = [c.strip().strip('"').strip() for c in lines[0].split("\t")][1:]
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        variant = parts[0].strip().strip('"').strip()
        seen: dict = {}
        for treatment, value in zip(header, parts[1:]):
            seen[treatment] = seen.get(treatment, 0) + 1
            rows.append({"variant": variant, "treatment": treatment,
                         "replicate": seen[treatment], "fold_change": float(value)})
    return pd.DataFrame(rows)


def pc3_table() -> pd.DataFrame:
    """The PC3 G822E file, which is delivered in long form already.

    Unlike the matrix files this one carries its own summary and p-values in a
    second section; `pc3_delivered_summary` reads that, and the tests check it
    against the summary computed here.
    """
    path = ASSAYS / "PC3" / "PC3_CDK12_G822E_CompAssay_heatmap_source_data.txt"
    lines = path.read_text().splitlines()
    start = lines.index("Genotype\tTreatment\tReplicate\tFoldChange_vs_T0")
    rows = []
    for ln in lines[start + 1:]:
        if not ln.strip() or ln.startswith("#") or ln.startswith("##"):
            break
        g, t, r, v = ln.split("\t")
        rows.append({"variant": g.strip(), "treatment": t.strip(),
                     "replicate": int(r), "fold_change": float(v)})
    return pd.DataFrame(rows)


def pc3_delivered_summary() -> pd.DataFrame:
    """Section 2 of the PC3 file: the summary its own author computed."""
    path = ASSAYS / "PC3" / "PC3_CDK12_G822E_CompAssay_heatmap_source_data.txt"
    lines = path.read_text().splitlines()
    start = lines.index(
        "Genotype\tTreatment\tn\tMean\tSD\tSEM\tp_vs_NT\tSignificance")
    rows = []
    for ln in lines[start + 1:]:
        if not ln.strip():
            break
        g, t, n, mean, sd, sem, p, sig = ln.split("\t")
        rows.append({"variant": g.strip(), "treatment": t.strip(), "n": int(n),
                     "mean": float(mean), "sd": float(sd), "sem": float(sem),
                     "p_vs_NT": None if p == "NA" else float(p),
                     "significance": None if sig == "NA" else sig})
    return pd.DataFrame(rows)


def ball_file(name: str) -> pd.DataFrame:
    """A mouse B-ALL summary: one column per arm, one row per replicate."""
    lines = [ln for ln in (ASSAYS / "B-ALL" / name).read_text().splitlines()
             if ln.strip()]
    header = [c.strip().strip('"').strip() for c in lines[0].split("\t")]
    rows = []
    for i, ln in enumerate(lines[1:], start=1):
        for treatment, value in zip(header, ln.split("\t")):
            rows.append({"treatment": treatment, "replicate": i,
                         "fold_change": float(value)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- panels ---

#: Each heatmap panel: where the replicates come from, which rows to show and
#: in what order, how to label them, which row the stars are measured against,
#: and where the thick rules go. `rows` is explicit rather than taken from the
#: file because two panels show fewer rows than their file holds -- 10g's file
#: carries an ABE and a CBE non-targeting arm and the published panel shows
#: only the ABE one, as a row labelled "NT".
PANEL_SPEC = {
    "2f": {
        "files": ["A549/Cdk9 FC New Variant.txt"],
        "variants": ["ABE NT", "CBE NT", "L156F", "I157V"],
        "labels": ["ABE NT", "CBE NT", "L156F", "I157V"],
        "treatments": ["Control", "KI-Cdk9d-32"],
        "reference": "ABE NT",
        # both non-targeting arms are controls: the published panel annotates
        # neither, and the test is against the ABE arm (the CBE arm gives
        # I157V a single star where the panel shows two).
        "controls": ["ABE NT", "CBE NT"],
        "transposed": True,          # variants across, treatments down
        "stars_on": ["KI-Cdk9d-32"],  # the published panel marks only this row
        "rules": [],
        "tick_fontsize": 10, "annot_fontsize": 8, "figsize": (5.0, 2.4),
        "rotate_xticks": 0,
    },
    "S10f": {
        "files": ["A549/Cdk12 FC.txt", "A549/Cdk13 FC.txt"],
        "variants": ["NT", "G731K", "I733V", "I711V"],
        "labels": ["NT", "CDK12$^{G731K}$", "CDK12$^{I733V}$", "CDK13$^{I711V}$"],
        "treatments": ["Control", "HQ461"],
        "reference": "NT",
        "transposed": False,
        "stars_on": ["Control", "HQ461"],
        "rules": [1, 3],
        "tick_fontsize": 8, "annot_fontsize": 7, "figsize": (2.4, 4.25),
        "rotate_xticks": 90,
    },
    "S10g": {
        "files": ["A549/Cdk12 FC ZLY.txt", "A549/Cdk13 FC ZLY.txt"],
        "variants": ["ABE NT", "G731K", "I733V", "I711V"],
        "labels": ["NT", "CDK12$^{G731K}$", "CDK12$^{I733V}$", "CDK13$^{I711V}$"],
        "treatments": ["Control", "ZLY-025"],
        "reference": "ABE NT",
        "transposed": False,
        "stars_on": ["Control", "ZLY-025"],
        "rules": [1, 3],
        "tick_fontsize": 8, "annot_fontsize": 7, "figsize": (2.4, 4.25),
        "rotate_xticks": 90,
    },
    "S11b": {
        "files": ["LNCaP/LNCaP Cdk12 FC.txt"],
        "variants": ["NT", "I733V"],
        "labels": ["NT", "CDK12$^{I733V}$"],
        "treatments": ["Control", "HQ461"],
        "reference": "NT",
        "transposed": False,
        "stars_on": ["Control", "HQ461"],
        "rules": [1],
        "tick_fontsize": 8, "annot_fontsize": 7, "figsize": (2.4, 4.25),
        "rotate_xticks": 90,
    },
    "S11I": {
        "files": [],                  # long form; see `pc3_table`
        "variants": ["NT", "CDK12 G822E"],
        "labels": ["NT", "CDK12$^{G822E}$"],
        "treatments": ["Control", "HQ461", "ZLY-025"],
        "reference": "NT",
        "transposed": False,
        "stars_on": ["Control", "HQ461", "ZLY-025"],
        "rules": [1],
        "tick_fontsize": 8, "annot_fontsize": 7, "figsize": (2.4, 4.25),
        "rotate_xticks": 90,
    },
}

#: the two mouse B-ALL panels: (file, the treated arm, the colour)
BALL_SPEC = {
    "S7c": ("mCdk9 Summary Plot.txt", "KI-Cdk9d-32", "tab:purple", "CDK9 L156F"),
    "S7e": ("mCdk7 Comp Assay Summary.txt", "SY-5609", "cyan", "CDK7 L18F"),
}


def replicates(panel: str) -> pd.DataFrame:
    """Every replicate behind a heatmap panel, in the panel's own row order."""
    spec = PANEL_SPEC[panel]
    if panel == "S11I":
        d = pc3_table()
    else:
        parts = [matrix_file(f) for f in spec["files"]]
        d = pd.concat(parts, ignore_index=True).drop_duplicates(
            subset=["variant", "treatment", "replicate"])
    d = d[d["variant"].isin(spec["variants"])
          & d["treatment"].isin(spec["treatments"])]
    order = {v: i for i, v in enumerate(spec["variants"])}
    torder = {t: i for i, t in enumerate(spec["treatments"])}
    d = d.assign(_v=d["variant"].map(order), _t=d["treatment"].map(torder))
    return (d.sort_values(["_v", "_t", "replicate"])
             .drop(columns=["_v", "_t"]).reset_index(drop=True))


def source(panel: str) -> pd.DataFrame:
    """The panel's source-data sheet: one row per cell, with the statistics.

    Mean and SEM are over the three replicates; the p-value is a two-sided
    unpaired t-test of the variant against the panel's non-targeting row within
    the same treatment, which is the test the delivered PC3 table names in its
    header and the one that reproduces every star on all five panels.
    """
    spec = PANEL_SPEC[panel]
    d = replicates(panel)
    rows = []
    for variant in spec["variants"]:
        for treatment in spec["treatments"]:
            x = d[(d["variant"] == variant)
                  & (d["treatment"] == treatment)]["fold_change"].to_numpy()
            ref = d[(d["variant"] == spec["reference"])
                    & (d["treatment"] == treatment)]["fold_change"].to_numpy()
            p = (np.nan if variant == spec["reference"]
                 else scipy.stats.ttest_ind(x, ref).pvalue)
            rows.append({
                "panel": panel, "variant": variant, "treatment": treatment,
                "n": len(x), "mean": x.mean(),
                "sd": x.std(ddof=1), "sem": x.std(ddof=1) / np.sqrt(len(x)),
                "p_vs_reference": p,
                "significance": "" if np.isnan(p) else _stars(p),
            })
    return pd.DataFrame(rows)


def _grid(panel: str):
    """-> (values, annotations, row labels, column labels) laid out as drawn."""
    spec = PANEL_SPEC[panel]
    s = source(panel).set_index(["variant", "treatment"])
    values, annot = [], []
    for variant in spec["variants"]:
        vrow, arow = [], []
        for treatment in spec["treatments"]:
            r = s.loc[(variant, treatment)]
            vrow.append(r["mean"])
            text = f"{r['mean']:.2f}\n±{r['sem']:.2f}"
            if (r["significance"] and treatment in spec["stars_on"]
                    and variant not in spec.get("controls", [spec["reference"]])):
                text += f"\n{r['significance']}"
            arow.append(text)
        values.append(vrow)
        annot.append(arow)
    values = np.asarray(values)
    annot = np.asarray(annot, dtype=object)
    rows, cols = spec["labels"], list(spec["treatments"])
    if spec["transposed"]:
        values, annot = values.T, annot.T
        rows, cols = cols, rows
    return values, annot, rows, cols


def heatmap(panel: str, outdir=None, **_):
    """-> (figure, the panel's source-data sheet)."""
    spec = PANEL_SPEC[panel]
    values, annot, rows, cols = _grid(panel)
    fig, ax = plt.subplots(figsize=spec["figsize"])
    sns.heatmap(values, cmap=CMAP, vmin=VMIN, vmax=VMAX, center=CENTER,
                annot=annot, fmt="", linewidths=0.8, linecolor="k",
                annot_kws={"fontsize": spec["annot_fontsize"]},
                xticklabels=cols, yticklabels=rows, ax=ax,
                cbar_kws={"label": CBAR_LABEL})
    ax.set_xticklabels(cols, rotation=spec["rotate_xticks"],
                       fontsize=spec["tick_fontsize"],
                       ha="center" if spec["rotate_xticks"] else "center")
    ax.set_yticklabels(rows, rotation=0 if spec["transposed"] else 90,
                       va="center", fontsize=spec["tick_fontsize"])
    for y in spec["rules"]:
        ax.hlines(y, *ax.get_xlim(), color="k", lw=2.5)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks(np.arange(int(VMIN), int(VMAX) + 1))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(CBAR_LABEL, fontsize=10)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(Path(outdir) / f"{_STEM[panel]}.pdf"), transparent=True,
                    bbox_inches="tight")
    return fig, source(panel)


# ------------------------------------------------------- mouse B-ALL pair ---

def ball_source(panel: str) -> pd.DataFrame:
    """Supplementary 7c/7e: the three replicates, with the arm comparison."""
    name, treated, _, genotype = BALL_SPEC[panel]
    d = ball_file(name)
    d.insert(0, "genotype", genotype)
    control = d[d["treatment"] == "Control"]["fold_change"].to_numpy()
    drug = d[d["treatment"] == treated]["fold_change"].to_numpy()
    t, p = scipy.stats.ttest_ind(drug, control)
    d.attrs["ttest"] = (t, p)
    return d


def ball_panel(panel: str, outdir=None, **_):
    """-> (figure, the replicates).

    No artwork was delivered for these two, so this is the bar-and-points style
    the other in vivo panels in this repository use rather than a reproduction
    of a published panel.
    """
    name, treated, colour, genotype = BALL_SPEC[panel]
    d = ball_source(panel)
    order = ["Control", treated]
    fig, ax = plt.subplots(figsize=(3.2, 4.2))
    sns.barplot(data=d, x="treatment", y="fold_change", order=order,
                palette=["lightgrey", colour], errorbar="se", capsize=0.2,
                err_kws={"linewidth": 1.5}, ax=ax, hue="treatment",
                legend=False)
    sns.stripplot(data=d, x="treatment", y="fold_change", order=order,
                  color="k", size=5, ax=ax)
    t, p = d.attrs["ttest"]
    ax.set_title(f"{genotype}\nt-test p = {p:.3g}", fontsize=11)
    ax.set_ylabel("Fold change (relative to T0)", fontsize=12)
    ax.set_xlabel("")
    ax.axhline(1, color="0.4", ls="--", lw=1)
    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="both", which="major", labelsize=11)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(Path(outdir) / f"{_STEM[panel]}.pdf"), transparent=True)
    return fig, d


# --------------------------------------------------------------- registry ---

#: the file name each panel writes, kept close to the delivered artwork's own
_STEM = {
    "2f": "CDK9_L156F_I157V_KI_CompAssay_heatmap",
    "S7c": "mCDK9_L156F_protac_competition",
    "S7e": "mCDK7_L18F_SY5609_competition",
    "S10f": "CDK12_CDK13_HQ461_CompAssay_heatmap",
    "S10g": "CDK12_CDK13_ZLY025_CompAssay_heatmap",
    "S11b": "LNCaP_CDK12_I733V_HQ461_CompAssay_heatmap",
    "S11I": "PC3_CDK12_G822E_CompAssay_heatmap",
}


def _wrap(panel: str, fn):
    def drawer(outdir=None, **kw):
        return fn(panel, outdir=outdir, **kw)
    drawer.__name__ = drawer.__qualname__ = f"panel_{panel.replace('.', '_')}"
    drawer.__doc__ = f"Competition-assay panel {panel}."
    return drawer


MAIN_PANELS = {"2f": lambda: source("2f")}
MAIN_FIGURES = {"2f": _wrap("2f", heatmap)}

PANELS = {k: (lambda _p=k: source(_p)) for k in
          ("S10f", "S10g", "S11b", "S11I")}
PANELS.update({k: (lambda _p=k: ball_source(_p)) for k in ("S7c", "S7e")})

FIGURES = {k: _wrap(k, heatmap) for k in ("S10f", "S10g", "S11b", "S11I")}
FIGURES.update({k: _wrap(k, ball_panel) for k in ("S7c", "S7e")})

#: `scripts/render_figures.py` reads this to name the file a panel produces
PARTS = {k: [v] for k, v in _STEM.items()}
