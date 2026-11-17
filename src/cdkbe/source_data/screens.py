"""Screen panels: lollipops and the cross-drug resistance summaries.

lollipop_plots_canonical.ipynb cell 15
    -> figures 2a, 2h, 3a   (the plotted points are `targ_*` after `guide_filter`)
logo_plots_CDK12_13_canonical.ipynb cells 19-20
    -> figure 3b   (top 12 hits plus the HQ461-specific ones)
logo_plots_CDK2_4_6_canonical.ipynb cells 29, 33
    -> figures 4a, 4b
"""
from __future__ import annotations
import pandas as pd

from .. import editing, io, lfc, paths

DMSO = paths.LFC / "barcode_counts" / "DMSO_base"

#: (abe table, cbe table, abe editing sample, cbe editing sample)
_TABLES = {
    "subpool1_DMSO": (DMSO / "ABE_subpool1_barcode_DMSO_base.csv",
                      DMSO / "CBE_subpool1_barcode_DMSO_base.csv",
                      "ABE_subpool1", "CBE_subpool1"),
    "SY-5609_DMSO":  (DMSO / "SY_5609_ABE_barcode_DMSO_base.csv",
                      DMSO / "SY_5609_CBE_barcode_DMSO_base.csv",
                      "ABE_subpool1", "CBE_subpool1"),
    "subpool2_DMSO": (DMSO / "CDK12_13_ABE_barcode_DMSO_base.csv",
                      DMSO / "CDK12_13_CBE_barcode_DMSO_base.csv",
                      "ABE_CDK12_13", "CBE_CDK12_13"),
    "subpool3_DMSO": (DMSO / "CDK2_4_6_ABE_barcode_DMSO_base.csv",
                      DMSO / "CDK2_4_6_CBE_barcode_DMSO_base.csv",
                      "ABE_CDK2_4_6", "CBE_CDK2_4_6"),
}

_OUT = ["Editor", "gRNA_id", "Gene", "HGVSp", "Codon", "LFC_median", "FDR",
        "WT_perc", "pct_reads_top_variant", "base_RAW", "labelled_hit"]


def _lollipop(table: str, genes, cond_cbe: str, cond_abe: str,
              lfc_cutoff: float, fdr_cutoff: float = 0.1) -> pd.DataFrame:
    abe_p, cbe_p, abe_s, cbe_s = _TABLES[table]
    edits = editing.edit_dict(sorted({abe_s, cbe_s}))
    _, _, targ_cbe, targ_abe = lfc.load_merged(abe_p, cbe_p, abe_s, cbe_s, edits)
    out = []
    for editor, targ, cond in (("CBE", targ_cbe, cond_cbe), ("ABE", targ_abe, cond_abe)):
        s = lfc.guide_filter(targ, genes)
        d = pd.DataFrame({
            "Editor": editor, "gRNA_id": s["gRNA_id"], "Gene": s["Gene"],
            "HGVSp": s["HGVSp"], "Codon": s["Codon"],
            "LFC_median": s[f"LFC_median_{cond}"], "FDR": s[f"FDR_{cond}"],
            "WT_perc": s["WT_perc"], "pct_reads_top_variant": s["%Reads"],
            "base_RAW": s["base_RAW"]})
        d["labelled_hit"] = (d["FDR"] <= fdr_cutoff) & (d["LFC_median"] >= lfc_cutoff)
        out.append(d)
    return pd.concat(out)[_OUT]


def figure_2a() -> pd.DataFrame:
    """Subpool 1 lollipop, KI-CDK9d-32 1,000 nM."""
    return _lollipop("subpool1_DMSO", ["CDK7", "CDK8", "CDK19", "CDK9"],
                     "KI-CDK9d-32_1000nM", "KI-CDK9d-32_1000nM", lfc_cutoff=2)


def figure_2h() -> pd.DataFrame:
    """Subpool 1 lollipop, SY-5609 10 nM."""
    return _lollipop("SY-5609_DMSO", ["CDK7", "CDK8", "CDK19", "CDK9"],
                     "CBE_SY-5609_10nM", "ABE_SY-5609_10nM", lfc_cutoff=2)


def figure_3a() -> pd.DataFrame:
    """Subpool 2 lollipop, HQ461."""
    return _lollipop("subpool2_DMSO", ["CDK12", "CDK13"],
                     "CBE_HQ461", "ABE_HQ461", lfc_cutoff=2)


def figure_3b(top_hits: int = 12, p_cutoff: float = 0.01) -> pd.DataFrame:
    """CDK12/13: the top hits plus the HQ461-specific ones.

    `hits_all` is already ordered in the saved table, so "top 12" is its first
    twelve rows once the HQ461-specific guides are held out.
    """
    t = io.read(paths.LFC / "CDK12_13_hits" / "CDK12_13_resistance_variants.csv")
    hq = (t["p_HQ461"] < p_cutoff) & (t["LFC_HQ461"] > 0) & (t["n_conditions"] == 1)
    specific = t[hq].reset_index(drop=True)
    key = "key" if "key" in t.columns else "gRNA_id"
    top = t[~t[key].isin(specific[key])].head(top_hits).reset_index(drop=True)
    top = top.assign(panel_group=f"Top {top_hits} resistance hits (excluding HQ461-specific)")
    specific = specific.assign(panel_group="HQ461-specific resistance hits")
    return pd.concat([top, specific], ignore_index=True)


CONDITIONS_S3 = {"RIBO": "Ribociclib", "ABEMA": "Abemaciclib", "ATIRMO": "Atirmociclib",
                 "PALBO": "Palbociclib", "TAGTO": "Tagtociclib", "INX": "INX-315"}
GENES_S3 = ["CDK2", "CDK4", "CDK6"]


def _s3_hits() -> pd.DataFrame:
    return io.read(paths.LFC / "CDK2_4_6_hits" / "CDK2_4_6_resistance_variants.csv")


def figure_4a() -> pd.DataFrame:
    """Per-gene cross-drug heatmap rows: every hit, with LFC/FDR/p per drug and T0."""
    t = _s3_hits()
    cols = ["gRNA_id", "Gene", "Editor", "Top_HGVSp", "Top_perc", "WT_perc", "base_RAW",
            "active_median_LFC", "median_LFC", "max_LFC", "max_cond", "n_conditions"]
    cols += [f"{s}_{c}" for c in CONDITIONS_S3 for s in ("LFC", "FDR", "p")]
    cols += ["LFC_T0", "p_T0", "FDR_T0"]
    return t[[c for c in cols if c in t.columns]]


def figure_4b(p_cutoff: float = 0.01) -> pd.DataFrame:
    """Resistance hits per drug by gene (a hit counts once per drug it resists in)."""
    t = _s3_hits()
    rows = []
    for cond, name in CONDITIONS_S3.items():
        active = t[(t[f"p_{cond}"] < p_cutoff) & (t[f"LFC_{cond}"] > 0)]
        for g in GENES_S3:
            rows.append(dict(Drug=name, Condition=cond, Gene=g,
                             n_hits=int((active["Gene"] == g).sum())))
    return pd.DataFrame(rows)


def figure_2i() -> dict:
    """CDK7 SY-5609 resistance hits: the bar heights, and the per-replicate LFCs
    behind their confidence intervals."""
    from . import screen_figures as _sf
    hits = _sf._sy5609_hits()
    keep = ["gRNA_id", "Editor", "Gene", "Top_HGVSp", "Top_perc", "WT_perc", "base_RAW"]
    keep += [f"{stat}_{c}" for c in _sf.SY5609_CONDITIONS for stat in ("LFC", "FDR")]
    return {"hits": hits[[c for c in keep if c in hits.columns]],
            "per_replicate_LFC": _sf._sy5609_replicates(hits)}


def figure_2b() -> dict:
    """CDK9 top resistance variants: the guides shown and their per-replicate LFCs."""
    from . import screen_figures as _sf
    hits = _sf._cdk9_top_hits()
    keep = ["gRNA_id", "Editor", "Gene", "label", "classification",
            *[f"LFC_median_{c}" for c in _sf.CDK9_CONDITIONS],
            *[f"FDR_{c}" for c in _sf.CDK9_CONDITIONS]]
    return {"hits": hits[[c for c in keep if c in hits.columns]],
            "per_replicate_LFC": _sf._cdk9_replicates(hits)}


def figure_4p() -> dict:
    """CDK6 clinical resistance variants across the three CDK4/6 inhibitors."""
    from . import screen_figures as _sf
    hits = _sf._clinical_hits()
    keep = ["gRNA_id", "Editor", "Gene", "label", "classification", "LFC_median_PALBO"]
    return {"hits": hits[[c for c in keep if c in hits.columns]],
            "per_replicate_LFC": _sf._clinical_replicates(hits)}
