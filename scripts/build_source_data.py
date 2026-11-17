#!/usr/bin/env python3
"""Write the manuscript source-data workbooks: one .xlsx per figure, one sheet per panel.

    python3 scripts/build_source_data.py [--out source_data_tables]

Each sheet is the frame the panel plots, reached by replaying that figure's own
filters against data/processed/. An INDEX sheet lists every panel with its row
and column counts and the original notebook cell it came from.
"""
from __future__ import annotations
import argparse, collections, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import pandas as pd                                            # noqa: E402
from cdkbe.source_data import (DRAWN, PANELS, SUPP_DRAWN,       # noqa: E402
                               SUPP_PANELS)

DESCRIPTIONS = {
    "1e": ("Codon and variant coverage vs sensor-editing cutoff, per CDK",
           "editing_analysis_v2.ipynb cell 16"),
    "1f": ("Resistance variants by residue, per CDK (max LFC over that CDK's drugs)",
           "fig1_analysis.ipynb cells 4-13"),
    "2a": ("Lollipop, KI-CDK9d-32 1,000 nM, subpool 1 (CDK7/8/19/9); median LFC vs DMSO",
           "lollipop_plots_canonical.ipynb cell 15"),
    "2h": ("Lollipop, SY-5609 10 nM, subpool 1; median LFC vs DMSO",
           "lollipop_plots_canonical.ipynb cell 15"),
    "2m": ("Xenograft tumour volume fold change (CDK9 I157V Day 5, CDK7 L18F Day 7)",
           "invivo_xenograft_CDK9_CDK7_combined.ipynb cell 5"),
    "2n": ("Xenograft Ki67 (CDK9 I157V, CDK7 L18F)",
           "invivo_xenograft_CDK9_CDK7_combined.ipynb cell 7"),
    "3a": ("Lollipop, HQ461, subpool 2 (CDK12/13); median LFC vs DMSO",
           "lollipop_plots_canonical.ipynb cell 15"),
    "3b": ("CDK12/13 resistance hits: top 12 plus HQ461-specific",
           "logo_plots_CDK12_13_canonical.ipynb cells 19-20"),
    "3f": ("PC3 DEG counts per variant and drug, incl. DMSO-vs-NT baseline",
           "reanalysis/scripts/03_figures.py"),
    "3g": ("PC3 drug response retained (% of NT signature shift)",
           "reanalysis/scripts/03_figures.py"),
    "3h": ("PC3 elongation profiles across gene-length deciles",
           "reanalysis/scripts/06_elongation_combined.py L92"),
    "4a": ("CDK2/4/6 cross-drug resistance heatmap; LFC/FDR/p per drug plus T0",
           "logo_plots_CDK2_4_6_canonical.ipynb cell 29"),
    "4b": ("Resistance hits per drug by gene",
           "logo_plots_CDK2_4_6_canonical.ipynb cell 33"),
    "4n": ("Xenograft tumour volume fold change, CDK6 R46W, palbociclib Day 7",
           "invivo_xenograft_CDK6_palbociclib.ipynb cell 20"),
    "3d": ("CDK12/13 variant validation MTT: LC50 fold shift vs matched NT control, "
           "extra-sum-of-squares F test, Holm adjusted",
           "FACS_dose_response_2024_8_22.ipynb cells 136, 142-143"),
    "4n": ("Xenograft tumour volume fold change, CDK6 R46W, palbociclib Day 7",
           "invivo_xenograft_CDK6_palbociclib.ipynb cell 20"),
    "4o": ("Xenograft Ki67, CDK6 R46W, palbociclib",
           "invivo_xenograft_CDK6_palbociclib.ipynb cell 22"),
}

SUPP_DESCRIPTIONS = {
    "S1a": ("Cumulative population doublings, subpool 1 screen", "CPDS.ipynb cell 2"),
    "S1b": ("CPDs, SY-5609 retest screen", "CPDS.ipynb cell 4"),
    "S1c": ("CPDs, CDK12/13 screen", "CPDS.ipynb cell 6"),
    "S1d": ("CPDs, CDK2/4/6 screen", "CPDS.ipynb cell 7"),
    "S3a": ("CBE sensor editing by protospacer position, per CDK",
            "editing_by_CDK_protospacer_location.ipynb cell 14"),
    "S3b": ("ABE sensor editing by protospacer position, per CDK",
            "editing_by_CDK_protospacer_location.ipynb cell 14"),
    "S4a": ("LFC correlation across subpool-1 comparisons, per editor",
            "LFC_correlation_clustermaps.ipynb cell 9"),
    "S4b": ("LFC correlation, subpool 2", "LFC_correlation_clustermaps.ipynb cell 9"),
    "S4c": ("LFC correlation, subpool 3", "LFC_correlation_clustermaps.ipynb cell 9"),
    "S7f": ("B-ALL leukemia burden, CDK9 L156F vs KI-CDK9d-32", "invivo.ipynb cell 2"),
    "S7g": ("The same, normalised to each genotype's vehicle arm", "invivo.ipynb cell 3"),
    "S7h": ("B-ALL leukemia burden, CDK7 L18F vs SY-5609", "invivo.ipynb cell 7"),
    "S7i": ("The same, normalised to each genotype's vehicle arm", "invivo.ipynb cell 8"),
    "S9d": ("Allele-specific MTT, SY-5609", "FACS_dose_response_2024_8_22.ipynb cell 131"),
    "S9e": ("Allele-specific MTT, KI-CDK9d-32", "FACS_dose_response_2024_8_22.ipynb cell 132"),
    "S9f": ("Allele-specific MTT, KB-0742", "FACS_dose_response_2024_8_22.ipynb cell 133"),
    "S9g": ("KB-0742 viability at the top four doses", "FACS_dose_response_2024_8_22.ipynb cell 134"),
    "S15b": ("OncoKB resistance mutations inside vs outside the KLIFS pocket, per kinase",
             "KLIFS_analysis.ipynb cells 17-21"),
    "S15c": ("KLIFS pocket resistance: CDK screen by position, and the OncoKB kinases "
             "on the same 85 positions", "CDK_resistance_patterns.ipynb cells 11-14"),
    "S16b": ("KLIFS position x CDK early-essentiality detail",
             "klifs_detail_figures.py build()"),
    "S16c": ("Depletion significance pooled to KLIFS subdomains",
             "essentiality_peaks.py within_pocket_test()"),
    "S16d": ("Depletion significance per KLIFS pocket column",
             "essentiality_peaks.py within_pocket_test()"),
    "S16e": ("Per-residue early essentiality (the montage's underlying scores)",
             "essentiality_figures.py fig_depleted_montage()"),
    "S2a": ("Dose response, subpool 1 compounds (A549); no background subtraction",
            "FACS_dose_response_2024_8_22.ipynb cells 157-160"),
    "S2b": ("Dose response, CDK12/13 compounds (A549)",
            "FACS_dose_response_2024_8_22.ipynb cells 161-164"),
    "S2c": ("Dose response, CDK4/6 inhibitors (NT CBE)",
            "FACS_dose_response_2024_8_22.ipynb cells 165-168"),
    "S5e": ("CDK9 L156F RNA-seq: DE gene counts per contrast",
            "RNAseq_analysis.ipynb cells 2-4"),
    "S8c": ("CDK8/19 resistance variants with their sensor logos",
            "logo_plots_CDK8_19_canonical.ipynb cell 16"),
    "S9b": ("Compound-mutant KB-0742 resistance variants with their sensor logos",
            "logo_plots_KB_compound_mut.ipynb cell 19"),
}
# `PANEL_COMPOSITE` maps a panel to (composite stem, its letter within it)
for _k, (_c, _letter) in __import__("cdkbe.source_data.beatseq_qc",
                                    fromlist=["PANEL_COMPOSITE"]).PANEL_COMPOSITE.items():
    _where = f"{_c} ({_letter})" if _letter else _c
    SUPP_DESCRIPTIONS[_k] = (
        f"BEAT-seq {_where.replace('_', ' ')}: the frame the delivered "
        "composite was drawn from",
        f"CDK-BEAT-seq plotting/figures/{_c}")
for _p, (_sub, _tab, _cbe, _abe, _title, _lfc) in __import__(
        "cdkbe.source_data.supp_figures", fromlist=["SUPP_SPEC"]).SUPP_SPEC.items():
    SUPP_DESCRIPTIONS[_p] = (f"Lollipop, {_title}, {_sub}; median LFC vs DMSO",
                             "lollipop_plots_canonical.ipynb cell 15")


def described(panel, descriptions, drawers):
    """-> (what the sheet is, where it came from).

    `DESCRIPTIONS` above is hand-written and does not cover every panel. Where
    it has nothing, the drawer's own docstring says what the panel is and its
    module says where the code lives, which is better than a blank cell and
    cannot go stale against the code the way a second copy would.
    """
    desc, origin = descriptions.get(panel, ("", ""))
    drawer = drawers.get(panel)
    if not desc and drawer is not None and drawer.__doc__:
        desc = drawer.__doc__.strip().splitlines()[0].rstrip(".")
    if not origin and drawer is not None:
        module = getattr(drawer, "__module__", "")
        if module.startswith("cdkbe"):
            origin = f"{module}.{getattr(drawer, '__name__', '')}"
    return desc, origin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "source_data_tables"))
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    import re

    def as_frame(obj):
        """A panel may hand back a Series; a sheet is always a frame."""
        if isinstance(obj, pd.Series):
            return obj.to_frame(name=obj.name or "value")
        return obj

    def tidy(df):
        """Drop the spreadsheet noise a frame carries in from its source.

        A column that is both unnamed -- pandas calls it `Unnamed: 16` when
        the original .xlsx had a blank header -- and entirely empty is a
        trailing cell range somebody once formatted, not data. Named columns
        are kept whatever they hold, and so is an empty column somebody
        named, because that absence may be the point.
        """
        junk = [c for c in df.columns
                if str(c).startswith("Unnamed:") and df[c].isna().all()]
        return df.drop(columns=junk) if junk else df

    def keeps_index(df):
        """Whether the row labels carry meaning and so have to be written.

        `to_excel(index=False)` was dropping the index of every frame, which
        for a matrix -- a heatmap's genes, a correlation's conditions, a
        contingency table's categories -- threw away the row labels and left a
        block of numbers nobody could redraw from. An anonymous integer index
        is still dropped: it is only position.
        """
        return bool(df.index.name) or not pd.api.types.is_integer_dtype(df.index)

    def write_book(xlsx, panels, registry, descriptions, drawers):
        index = []
        with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
            for panel in panels:
                res = registry[panel]()
                frames = res if isinstance(res, dict) else {"": res}
                for suffix, df in frames.items():
                    sheet = f"{panel}_{suffix}"[:31] if suffix else panel
                    df = tidy(as_frame(df))
                    labelled = keeps_index(df)
                    if labelled and not df.index.name:
                        # the sheet's first column is these labels, and a
                        # column with no header is not self-describing
                        df = df.rename_axis("label")
                    df.to_excel(w, sheet_name=sheet, index=labelled)
                    desc, origin = described(panel, descriptions, drawers)
                    index.append(dict(sheet=sheet, panel=panel, rows=len(df),
                                      columns=len(df.columns) + int(labelled),
                                      row_labels=df.index.name or
                                      ("row labels" if labelled else ""),
                                      description=desc, from_original=origin))
            pd.DataFrame(index).to_excel(w, sheet_name="INDEX", index=False)
            wb = w.book
            wb.move_sheet("INDEX", offset=-len(wb.sheetnames) + 1)
        print(f"{xlsx.name}: {len(index)} sheets")

    by_supp = collections.defaultdict(list)
    for panel in SUPP_PANELS:
        by_supp[int(re.match(r"S(\d+)", panel).group(1))].append(panel)
    for fignum, panels in sorted(by_supp.items()):
        write_book(out / f"Supplementary_Figure_{fignum}.xlsx",
                   sorted(panels), SUPP_PANELS, SUPP_DESCRIPTIONS, SUPP_DRAWN)

    by_fig = collections.defaultdict(list)
    for panel in sorted(PANELS):
        by_fig[panel[0]].append(panel)

    for fignum, panels in sorted(by_fig.items()):
        xlsx = out / f"Figure_{fignum}.xlsx"
        write_book(xlsx, panels, PANELS, DESCRIPTIONS, DRAWN)


if __name__ == "__main__":
    main()
