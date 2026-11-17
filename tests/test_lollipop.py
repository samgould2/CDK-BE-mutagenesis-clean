"""The lollipop figures must label exactly the guides the source data calls hits.

Three ~150-line plotters that differed only in gene list, figure size, panel
widths and x ticks were merged into one. The risk in that merge is that the
figure and the source-data workbook drift apart, so these tests derive the hit
set from the drawn frames and compare it with the `labelled_hit` column.
"""
import matplotlib
matplotlib.use("Agg")

import pytest

from cdkbe import lollipop, plotting
from cdkbe.source_data import PANELS
from cdkbe.source_data import screen_figures as sf

PANEL_AXES = {"2a": 5, "2h": 5, "3a": 3}     # genes + the controls panel


@pytest.fixture(scope="module", autouse=True)
def style():
    plotting.use_style()


@pytest.fixture(scope="module")
def drawn():
    return {p: sf.draw(p) for p in PANEL_AXES}


@pytest.mark.parametrize("panel", sorted(PANEL_AXES))
def test_panel_count(drawn, panel):
    fig, _ = drawn[panel]
    assert len(fig.get_axes()) == PANEL_AXES[panel]


@pytest.mark.parametrize("panel", sorted(PANEL_AXES))
def test_labelled_hits_match_source_data(drawn, panel):
    _, kept = drawn[panel]
    _, _, cond_cbe, cond_abe, _, cutoff = sf.PANEL_SPEC[panel]
    on_figure = set()
    for editor, cond in (("CBE", cond_cbe), ("ABE", cond_abe)):
        g = kept[editor]
        hit = g[(g[f"FDR_{cond}"] <= 0.1) & (g[f"LFC_median_{cond}"] >= cutoff)]
        on_figure |= set(hit["gRNA_id"])
    sd = PANELS[panel]()
    in_table = set(sd.loc[sd["labelled_hit"], "gRNA_id"])
    assert on_figure == in_table, (
        f"{panel}: {len(on_figure ^ in_table)} guides differ between "
        "the figure and its source-data sheet")


def test_find_ranges():
    assert lollipop.find_ranges([1, 2, 3, 7, 8, 20]) == [(1, 3), (7, 8), (20, 20)]
    assert lollipop.find_ranges([]) == []
    assert lollipop.find_ranges([5, 5, 6]) == [(5, 6)]


def test_control_classes_order():
    """The controls panel x order is Intron, NT, WT-heavy -- as published."""
    import pandas as pd
    m = pd.DataFrame({"classification": ["intron", "targeting", "targeting"],
                      "Gene": ["CDK9", "NT", "CDK9"],
                      "WT_perc": [0.0, 0.0, 99.0],
                      "LFC_median_x": [0.1, 0.2, 0.3]})
    out = lollipop.control_classes(m, "LFC_median_x", 98, "CBE")
    assert list(dict.fromkeys(out["Guide Class"]))[:2] == ["Intron", "NT"]
