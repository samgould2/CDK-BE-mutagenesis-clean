"""PC3 RNA-seq panels (figure 3f-h).

the reanalysis scripts, under data/processed/rnaseq/.../reanalysis/scripts:
  03_figures.py       writes Fig2b_DEG_counts_with_baseline.csv beside the figure
                      and Fig3_pct_response_per_replicate.csv  -> 3f, 3g
  06_elongation_combined.py line 92 writes Fig10_replicate_decile_profiles.csv -> 3h

These scripts save the plotted frame next to the PDF, so there is nothing to
re-derive; the functions below just locate and load it.

Note on 3f: the panel carries the DMSO-vs-NT baseline sub-panel, whose counts
exist only in the `Fig2b` (with-baseline) table, not in `Fig2`.
"""
from __future__ import annotations
import pandas as pd

from .. import io, paths

PC3 = paths.PROCESSED / "rnaseq" / "PC3_RNAseq" / "reanalysis" / "figures"


def figure_3f() -> pd.DataFrame:
    """DEG counts per variant and drug, including the DMSO-vs-NT baseline."""
    return io.read(PC3 / "Fig2b_DEG_counts_with_baseline.csv")


def figure_3g() -> dict:
    """Drug response retained. Two sheets: per-replicate percentages and per-sample score."""
    return {
        "drug_response_retained": io.read(PC3 / "Fig3_pct_response_per_replicate.csv"),
        "signature_score_per_sample": io.read(PC3 / "Fig3_signature_score_per_sample.csv"),
    }


def figure_3h() -> pd.DataFrame:
    """Elongation profiles across gene-length deciles (mean +/- SD over 3 replicates)."""
    return io.read(PC3 / "Fig10_replicate_decile_profiles.csv")
