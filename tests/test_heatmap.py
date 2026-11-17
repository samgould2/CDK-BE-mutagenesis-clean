"""The cross-drug heatmap must keep the decisions its original comments record.

Two are easy to undo by accident and both change the figure:

1. the row linkage is computed on the six drug columns only, *before* the T0
   column is appended, so a `_withT0` figure keeps exactly the row order of its
   plain counterpart and baseline fitness never drives the drug clustering;
2. correlation distance is used unless some row has zero variance, where it is
   undefined, in which case it falls back to euclidean.
"""
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from cdkbe import heatmap, plotting
from cdkbe.source_data import screens
from cdkbe.source_data import screen_figures as sf


@pytest.fixture(scope="module", autouse=True)
def style():
    plotting.use_style()


@pytest.fixture(scope="module")
def hits():
    return screens._s3_hits()


def test_t0_column_does_not_change_row_order(hits):
    """A with-T0 heatmap must order rows exactly like the plain one."""
    sub = hits[hits["Gene"] == "CDK4"].reset_index(drop=True)
    orders = []
    for include in (False, True):
        fig = heatmap.cross_drug_heatmap(sub, "t", per_gene=True, width=6.6,
                                         include_T0=include, suffix="_large")
        labels = [t.get_text() for t in fig.axes[2].get_yticklabels()]
        orders.append(labels)
        matplotlib.pyplot.close(fig)
    assert orders[0] == orders[1]


def test_hit_counts_match_source_data(hits):
    """4b's bars are the counts in its source-data sheet."""
    fig, counts_df = sf.draw_4b()
    expected = screens.figure_4b()
    pd.testing.assert_frame_equal(
        counts_df.sort_values(["Drug", "Gene"]).reset_index(drop=True),
        expected.sort_values(["Drug", "Gene"]).reset_index(drop=True))
    matplotlib.pyplot.close(fig)


@pytest.mark.parametrize("gene", heatmap.GENES)
def test_heatmap_rows_are_that_genes_hits(gene, hits):
    res = sf.draw_4a(genes=[gene])
    fig, sub = res[gene]
    assert len(sub) == int((hits["Gene"] == gene).sum())
    assert (sub["Gene"] == gene).all()
    matplotlib.pyplot.close(fig)


def test_flat_row_falls_back_to_euclidean():
    """A row with no variance has no correlation distance; it must still draw."""
    cols = [f"LFC_{c}" for c in heatmap.CONDITIONS]
    n = 4
    d = pd.DataFrame({c: np.linspace(0, 1, n) for c in cols})
    d.loc[0, cols] = 0.5                       # flat row
    for c in heatmap.CONDITIONS:
        d[f"p_{c}"] = 0.5
    d["Gene"] = "CDK4"
    d["Editor"] = "ABE"
    d["Top_HGVSp"] = [f"A{i}T" for i in range(n)]
    fig = heatmap.cross_drug_heatmap(d, "flat", per_gene=True, width=6.0)
    assert fig.axes
    matplotlib.pyplot.close(fig)
