"""The in vivo panels must reproduce the p-values printed on the published figures.

These four panels are the ones where the figure itself carries the numbers, so
they can be checked against the manuscript directly rather than against another
table. Every value below was read off the published panel.
"""
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import pytest

from cdkbe import plotting
from cdkbe.source_data import invivo_figures as ivf

#: p-values as printed on each published panel, to two significant figures
PUBLISHED = {
    "2m": [0.18, 0.045, 0.45, 0.68, 0.36, 0.016, 0.004, 0.30],
    "2n": [0.92, 0.048, 0.016, 0.75, 0.012, 0.84, 0.008, 0.73],
    # the CDK6 xenograft pair; these were figure 4o and 4p before the
    # November relabelling (docs/DISCREPANCIES.md item 28)
    "4n": [0.44, 0.34, 0.025, 0.29],
    "4o": [0.12, 0.69, 0.018, 0.17],
}


@pytest.fixture(scope="module", autouse=True)
def style():
    plotting.use_style()


@pytest.mark.parametrize("panel", sorted(PUBLISHED))
def test_pvalues_match_published_figure(panel):
    results = ivf.FIGURES[panel]()
    stats = pd.concat([s for _, s in results])
    got = [float(f"{p:.2g}") for p in stats["p_shown"]]
    want = PUBLISHED[panel]
    assert len(got) == len(want), f"{panel}: {len(got)} comparisons, expected {len(want)}"
    for g, w in zip(got, want):
        assert abs(g - w) <= max(0.005, 0.05 * w), f"{panel}: got {g}, figure says {w}"


@pytest.mark.parametrize("panel", sorted(PUBLISHED))
def test_panel_draws(panel):
    for fig, _ in ivf.FIGURES[panel]():
        assert fig.get_axes(), f"{panel}: figure has no axes"
