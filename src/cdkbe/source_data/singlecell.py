"""Figure 5: the CROP-seq single-cell readout.

These panels were produced by a separate analysis (CDK-BEAT-seq). Its derived
tables are copied into `data/processed/beatseq/` so this repository is
self-contained; the delivered single-cell objects behind them (~2.4 GB of
.h5ad) are not shipped and nothing here reads them.

Panel sources, from that analysis's `fig_panels.py`:

    5a  library_composition.csv            panel_i, the library donuts
    5b  umap_cells_A549.csv.gz             panel_a, raw embedding
    5c  umap_cells_T47D.csv.gz             panel_a, raw embedding
    5d  nmd_readout_classified.csv         panel_c1s
    5e  editing_target_transcript_phaseadj.csv   panel_c3
    5j  severity.csv                       panel_d12
    5k  severity_vs_resistance.csv         panel_e

5f-5i come from the T47D CDK4 RNA-seq instead and live in `rnaseq_t47d`.
"""
from __future__ import annotations
import pandas as pd

from .. import io, paths

BEATSEQ = paths.PROCESSED / "beatseq"

#: the UMAP panels plot the raw embedding, not the cell-cycle-regressed one
EMBEDDING = "raw"
UMAP_COLUMNS = ["cell_barcode", "cell_line", "arm", "guide_target", "guide_identity",
                "phase", "S_score", "G2M_score", "umap1_raw", "umap2_raw"]


def figure_5a() -> pd.DataFrame:
    """CROP-seq library composition, by editor and guide category."""
    return io.read(BEATSEQ / "library_composition.csv")


def _umap(cell_line: str) -> pd.DataFrame:
    d = io.read(BEATSEQ / f"umap_cells_{cell_line}.csv.gz", index_col=0)
    d = d.reset_index().rename(columns={"index": "cell_barcode"})
    return d[[c for c in UMAP_COLUMNS if c in d.columns]]


def figure_5b() -> pd.DataFrame:
    """Per-cell UMAP, A549: coordinates, CDK targeted, cell-cycle phase."""
    return _umap("A549")


def figure_5c() -> pd.DataFrame:
    """Per-cell UMAP, T47D."""
    return _umap("T47D")


def figure_5d() -> pd.DataFrame:
    """Editing readout by variant class."""
    return io.read(BEATSEQ / "nmd_readout_classified.csv")


def figure_5e() -> pd.DataFrame:
    """Target-transcript editing against the cell cycle, phase-adjusted."""
    return io.read(BEATSEQ / "editing_target_transcript_phaseadj.csv")


def figure_5j() -> pd.DataFrame:
    """Transcriptional severity per gRNA, by CDK."""
    return io.read(BEATSEQ / "severity.csv")


def figure_5k() -> pd.DataFrame:
    """Severity against the maximum drug LFC (the resistance axis)."""
    d = io.read(BEATSEQ / "severity_vs_resistance.csv")
    return d.dropna(subset=["max_drug_lfc", "severity"])
