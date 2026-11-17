"""Resistance-variant panels for CDK8/19 and the compound-mutant screen (supp 8c, 9b).

Ported from `logo_plots_CDK8_19_canonical.ipynb` cell 16 and
`logo_plots_KB_compound_mut.ipynb` cell 19. Both are the `hits_with_logos`
layout the main figures use (2b, 2i, 4p), at the notebooks' `_large` style --
the published supplementary panels are printed bigger than the main ones.

The hit tables are shipped, so nothing re-derives them: the notebooks write
`{OUT}_resistance_variants.csv` beside the figure and those tables are in
`data/processed/lfc_fdr_tables/`. The per-replicate LFCs the bars show are read
back out of the screen's own LFC table.

Two differences from each other worth keeping straight:

* 8c colours the variant labels **by editor** (ABE purple, CBE blue), since
  both editors contribute hits;
* 9b is a CBE-only screen, so it colours them **by gene** instead, with CDK7
  picked out in red -- that is the point of the panel.
"""
from __future__ import annotations

import pandas as pd

from .. import counts as _counts
from .. import io, logo_panel, paths

#: the notebooks' '_large' style: the supplementary panels are printed bigger
LARGE = dict(ylabel=17, xtick=15, xlabel=18, logo=11, tick_len=8.0, tick_w=1.8)

SPEC = {
    "S8c": dict(
        out="CDK8_19",
        hits_table=paths.LFC / "CDK8_19_hits" / "CDK8_19_resistance_variants.csv",
        lfc_tables={"ABE": paths.LFC / "barcode_counts/DMSO_base/ABE_subpool1_barcode_DMSO_base.csv",
                    "CBE": paths.LFC / "barcode_counts/DMSO_base/CBE_subpool1_barcode_DMSO_base.csv"},
        sav_source=("shipped", "subpool1"),
        conditions={"SEL120_4000nM": "SEL120 4,000 nM",
                    "Senexin B_2000nM": "Senexin B 2,000 nM",
                    "Senexin B_15000nM": "Senexin B 15,000 nM"},
        palette=["tab:blue", "lightsteelblue", "darkslateblue"],
        color_by="Editor", label_colors=logo_panel.EDITOR_COLORS,
        title=("CDK8/19 resistance variants | all {n} guides with FDR < 0.1 & LFC >= 2 "
               "in >=1 of the 3 conditions\n"
               "base mean count >= 100, sensor editing > 10%, ordered by median LFC "
               "across resistant conditions  |  ABE labels purple, CBE blue")),
    "S9b": dict(
        out="KB_compound_mut",
        hits_table=paths.LFC / "KB_compound_mut_hits" / "KB_compound_mut_resistance_variants.csv",
        lfc_tables={"CBE": paths.LFC / "proto_counts/DMSO_base/Compound_mutant_proto_DMSO_base.csv"},
        sav_source=("rebuild", "KB_compound_mut"),
        conditions={"KB_2000": "KB-0742 2,000 nM", "KB_4000": "KB-0742 4,000 nM"},
        palette=["lightsteelblue", "darkslateblue"],
        color_by="Gene", label_colors={"CDK7": "tab:red"},
        title=("KB-0742 resistance variants, compound mutant screen (CBE)  |  all {n} "
               "guides with FDR < 0.1 at >=1 dose & LFC >= 2 at >=1 dose\n"
               "base mean count >= 100, sensor editing > 10%, ordered by median LFC "
               "across resistant doses  |  CDK7 labels red")),
}


def hits(key: str) -> pd.DataFrame:
    """The shipped hit table, with the row key and label the panel uses."""
    spec = SPEC[key]
    h = io.read(spec["hits_table"]).copy()
    h["key"] = h["gRNA_id"] + "_" + h["Editor"]
    h["label"] = h["Gene"] + "\n" + h["Top_HGVSp"]
    return h


def replicates(key: str) -> pd.DataFrame:
    """Per-replicate LFCs for each hit, read back out of the screen's LFC table."""
    spec = SPEC[key]
    tables = {e: io.read(p).set_index("gRNA_id") for e, p in spec["lfc_tables"].items()}
    rows = []
    for v in hits(key).itertuples():
        t = tables[v.Editor]
        for cond, name in spec["conditions"].items():
            for rep in (1, 2, 3):
                rows.append(dict(gRNA_id=v.gRNA_id, key=v.key, Gene=v.Gene,
                                 Editor=v.Editor, Condition=name, Replicate=rep,
                                 LFC=t.loc[v.gRNA_id, f"{cond}_REP{rep}"]))
    out = pd.DataFrame(rows)
    assert not out["LFC"].isna().any(), f"{key}: a hit is missing a replicate"
    return out


def _savs(key: str) -> dict:
    """Single amino-acid variants per editor, in the convention its panel used.

    8c reads the shipped `SAVs_canonical` tables, as the main-figure panels do.
    9b rebuilds them with the editor-own-edits filter its notebook applies --
    the two conventions give different tables; see docs/DISCREPANCIES.md 14.
    """
    from .. import editing as _editing

    spec = SPEC[key]
    how, stem = spec["sav_source"]
    library = _counts.library()[["gRNA_id", "Gene"]]
    out = {}
    for editor in spec["lfc_tables"]:
        if how == "shipped":
            s = io.read(paths.EDITING / "SAVs_canonical" / f"{editor}_{stem}_canonical.csv")
        else:
            s = _editing.savs_from_full_canonical(f"{editor}_{stem}", own_edits_only=True)
        out[editor] = s.merge(library, on="gRNA_id")
    return out


def draw(key: str, outdir=None, **kw):
    """-> (figure, the hit table)."""
    spec = SPEC[key]
    h = hits(key)
    fig = logo_panel.hits_with_logos(
        h, replicates(key), _savs(key), palette=spec["palette"],
        hue_order=list(spec["conditions"].values()),
        title=spec["title"].format(n=len(h)),
        sizes=LARGE, label_colors=spec["label_colors"], color_by=spec["color_by"],
        legend_loc="figure", legend_fontsize=14, **kw)
    if outdir is not None:
        fig.savefig(str(outdir / f"Supp_{key[1:]}_{spec['out']}_resistance_variants.pdf"),
                    transparent=True)
    return fig, h


def _wrap(panel):
    sp = SPEC[panel]
    f = lambda **kw: draw(panel, **kw)          # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    # the spec's `title` is the published panel title, with its own line break
    # and the {n} placeholder the hit count fills in; the first clause of it is
    # the description.
    f.__doc__ = (sp["title"].split("|")[0].strip()
                 + ": LFC bars beside each guide's edit logo, over "
                 + ", ".join(sp["conditions"].values()) + ".")
    return f


FIGURES = {k: _wrap(k) for k in SPEC}
PANELS = {k: (lambda _k=k: hits(_k)) for k in SPEC}
