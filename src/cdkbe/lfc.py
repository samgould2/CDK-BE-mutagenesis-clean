"""LFC/FDR tables joined to the sensor-editing call for each gRNA."""
from __future__ import annotations
import pandas as pd

from . import io, paths

#: guide-level filters shared by the lollipop and logo figures
MIN_BASE_COUNTS = 100     # raw counts in the baseline arm
WT_CUTOFF_MAIN = 90       # % WT reads above which a guide has no usable edit
FDR_CUTOFF = 0.1


def load_merged(abe_path, cbe_path, abe_edit: str, cbe_edit: str, edits: dict):
    """Merge an ABE/CBE LFC pair with their editing tables.

    Returns (m_cbe, m_abe, targ_cbe, targ_abe). The `targ_*` frames keep only
    targeting guides that actually have sensor information -- a guide with no
    editing call cannot be placed on a residue axis.
    """
    abe = io.read(abe_path)
    cbe = io.read(cbe_path)
    m_abe = pd.merge(abe, edits[abe_edit], on="gRNA_id", how="outer")
    targ_abe = m_abe[(m_abe["classification"] == "targeting") & (~m_abe["HGVSp"].isna())]
    m_cbe = pd.merge(cbe, edits[cbe_edit], on="gRNA_id", how="outer")
    targ_cbe = m_cbe[(m_cbe["classification"] == "targeting") & (~m_cbe["HGVSp"].isna())]
    return m_cbe, m_abe, targ_cbe, targ_abe


def guide_filter(targ: pd.DataFrame, genes, min_base_counts=MIN_BASE_COUNTS,
                 wt_cutoff=WT_CUTOFF_MAIN) -> pd.DataFrame:
    """The filter every residue-level figure applies before plotting."""
    return targ[(targ["base_RAW"] >= min_base_counts)
                & (targ["Gene"].isin(genes))
                & (targ["WT_perc"] < wt_cutoff)]


def hits(targ: pd.DataFrame, condition: str, lfc_cutoff: float,
         fdr_cutoff: float = FDR_CUTOFF) -> pd.DataFrame:
    """Guides called as resistant in one condition."""
    lfc, fdr = f"LFC_median_{condition}", f"FDR_{condition}"
    return targ[(targ[fdr] <= fdr_cutoff) & (targ[lfc] >= lfc_cutoff)]


def table_path(kind: str, baseline: str, name: str):
    """`lfc.table_path('barcode_counts', 'DMSO_base', 'ABE_subpool1_barcode_DMSO_base')`"""
    return paths.LFC / kind / baseline / f"{name}.csv"
