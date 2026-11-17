"""Early-essentiality panels (supplementary 16b-16e).

The analysis itself lives in `cdkbe.essentiality`, ported from the sixteen
scripts that were in `notebooks/09_analysis_and_viz/`. This module only picks
out the four panels the supplementary figure uses and hands back the frame
behind each.

None of them re-runs the permutation null: `essentiality_peaks.main()` writes
`klifs_within_pocket_tests.csv` and `klifs_subdomain_tests.csv`, both of which
ship in `data/processed/lfc_fdr_tables/T0_plasmid_essentiality/`, so the panels
read those. Rebuilding them means 2,000 circular shifts per gene and is a
`python3 -m cdkbe.essentiality.essentiality_peaks` job, not a figure step.

16e is a montage of ChimeraX renders. They cannot be regenerated without
ChimeraX, so they ship under `data/reference/structure_renders/essentiality`
and are treated as input.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from .. import io, paths

TABLES = paths.LFC / "T0_plasmid_essentiality"


def within_pocket_tests() -> pd.DataFrame:
    """Per-KLIFS-column depletion test, ranked (the `wp` frame)."""
    return io.read(TABLES / "klifs_within_pocket_tests.csv")


def subdomain_tests() -> pd.DataFrame:
    """The same test pooled to KLIFS subdomains (the `wsd` frame)."""
    return io.read(TABLES / "klifs_subdomain_tests.csv")


def klifs_position_table() -> pd.DataFrame:
    """The 9 x 85 KLIFS position table behind 16b."""
    from ..essentiality import klifs_detail_figures as kdf
    t, _ = kdf.build()          # build() returns (table, None)
    return t


def supp_16b(outdir=None, **kw):
    """KLIFS position detail: nine CDKs by eighty-five pocket positions."""
    from ..essentiality import klifs_detail_figures as kdf
    t = klifs_position_table()
    fig = kdf.fig_heatmap(t)
    if outdir is not None:
        fig.savefig(str(outdir / "Supp_16b_klifs_position_detail.pdf"),
                    bbox_inches="tight")
    return fig, t


def supp_16cd(outdir=None, **kw):
    """Significance of the KLIFS calls: subdomains, then the pocket columns.

    16c and 16d are the two halves of one published panel
    (`peak_significance_wide_flipped.pdf`), drawn with `orientation='v'` and
    thirty columns, so they are produced together.
    """
    from ..essentiality import essentiality_peaks as ep
    wp, wsd = within_pocket_tests(), subdomain_tests()
    fig = ep.fig_klifs_significance(wp, wsd,
                                    fname="peak_significance_wide_flipped.pdf",
                                    figsize=(17.0, 6.6), orientation="v",
                                    n_cols=30)
    if outdir is not None:
        fig.savefig(str(outdir / "Supp_16cd_peak_significance.pdf"),
                    bbox_inches="tight")
    return fig, {"subdomain_tests": wsd, "pocket_columns": wp}


def supp_16e(outdir=None, **kw):
    """Depleted residues on all nine CDK structures, five across."""
    from ..essentiality import essentiality_figures as ef
    fig = ef.fig_depleted_montage(cols=5, figsize=(16.5, 7.6),
                                  fname="depleted_residues_montage_wide.pdf")
    if outdir is not None:
        fig.savefig(str(outdir / "Supp_16e_depleted_montage.pdf"),
                    bbox_inches="tight")
    return fig, io.read(TABLES / "residue_essentiality.csv")


FIGURES = {"S16b": supp_16b, "S16c": supp_16cd, "S16d": supp_16cd, "S16e": supp_16e}
PANELS = {
    "S16b": klifs_position_table,
    "S16c": subdomain_tests,
    "S16d": within_pocket_tests,
    "S16e": lambda: io.read(TABLES / "residue_essentiality.csv"),
}
