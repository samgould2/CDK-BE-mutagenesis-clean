"""The bar-plus-logo panel, shared by figures 2b, 2i, 3b and 4n.

Each logo column must sum to 100: the observed edits at their measured %Reads,
and the reference residue taking whatever is left. If that invariant breaks the
glyph heights stop meaning "fraction of sensor reads" and the panel is quietly
wrong rather than visibly broken.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from cdkbe import counts, io, logo_panel, paths, plotting
from cdkbe.source_data import screen_figures as sf


@pytest.fixture(scope="module", autouse=True)
def style():
    plotting.use_style()


@pytest.fixture(scope="module")
def savs():
    library = counts.library()
    out = {}
    for editor in ("ABE", "CBE"):
        s = io.read(paths.EDITING / "SAVs_canonical" / f"{editor}_subpool1_canonical.csv")
        out[editor] = s.merge(library[["gRNA_id", "Gene"]], on="gRNA_id")
    return out


def test_logo_columns_sum_to_100(savs):
    hits = sf._sy5609_hits()
    for _, row in hits.iterrows():
        frame, lo, hi, d, seq = logo_panel.logo_frame(savs[row["Editor"]], row["gRNA_id"])
        sums = frame.sum(axis=1).to_numpy()
        assert np.allclose(sums, 100, atol=1e-6), \
            f"{row['gRNA_id']}: logo columns sum to {sums.min():.1f}-{sums.max():.1f}"


def test_logo_window_covers_the_edits_plus_context(savs):
    hits = sf._sy5609_hits()
    row = hits.iloc[0]
    frame, lo, hi, d, seq = logo_panel.logo_frame(savs[row["Editor"]], row["gRNA_id"],
                                                  context=1)
    assert frame.index.min() <= lo and frame.index.max() >= hi
    assert len(seq) == len(frame)


def test_sy5609_hits_are_the_published_nine():
    hits = sf._sy5609_hits()
    assert len(hits) == 9
    assert list(hits["Top_HGVSp"])[:3] == ["L18F", "L18F", "L18F"], \
        "three separate guides all call L18F, which is why the row key carries the gRNA id"
    assert hits["key"].is_unique


def test_panel_draws():
    fig, hits = sf.draw_2i()
    # one shared bar axis plus one logo per hit
    assert len(fig.get_axes()) == len(hits) + 1
    plt.close(fig)
