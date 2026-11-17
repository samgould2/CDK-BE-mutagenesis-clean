"""Reading data back out of a published figure.

`cdkbe.svg_values` exists for the panels this repository has as a published
figure but not as a table. Inverting a plot is the kind of thing that looks
right and is wrong, so the module is pinned against a case where the answer is
already known: the T47D cell-cycle panels, whose published SVGs *and*
underlying FlowJo export are both here. If the inversion cannot reproduce
those exactly, nothing it says elsewhere is worth anything.

The second half checks the two tables recovered this way -- figure 4g and
supplementary 13a -- against the figures they came from, so the files on disk
cannot drift away from their source.
"""
from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cdkbe import paths, svg_values as sv  # noqa: E402
from cdkbe.source_data import t47d_cell_cycle as cc  # noqa: E402
from cdkbe.source_data import selection_allele_frequency as sa  # noqa: E402

ARTWORK = paths.RAW / "collaborator_panels" / "T47D_MCF7"
ROOT = Path(__file__).resolve().parents[1]

#: the published panels whose values this repository independently knows
KNOWN = [("T47D_cell_cycle_stacked_Sphase_ABE.svg", cc.ABE_LINES),
         ("T47D_cell_cycle_stacked_Sphase_CBE.svg", cc.CBE_LINES)]


@pytest.mark.parametrize("stem,lines", KNOWN)
def test_the_inversion_reproduces_values_this_repository_already_knows(stem, lines):
    """Bar heights off the artwork == gate frequencies off the FlowJo export."""
    got = [value for _, value in sv.bars(sv.load(ARTWORK / stem))]
    means = cc.means().set_index(["CellLine", "Treatment"])
    want = [float(means.loc[(line, treatment), "S"])
            for treatment in cc.TREATMENT_ORDER for line in lines]
    # matplotlib writes coordinates to six decimal places, which on this
    # axis is a little under 2e-7 percentage points -- that, not the method,
    # is the whole of the difference
    assert got == pytest.approx(want, abs=1e-6)


@pytest.mark.parametrize("stem,lines", KNOWN)
def test_the_axes_background_and_spines_are_not_mistaken_for_bars(stem, lines):
    assert len(sv.bars(sv.load(ARTWORK / stem))) == 2 * len(lines)


def test_a_non_linear_axis_is_refused_rather_than_inverted():
    """A log axis with plain numeric labels is the dangerous case: every tick
    parses as a number, and reading it as linear would be silently wrong."""
    from matplotlib.ticker import ScalarFormatter

    fig, ax = plt.subplots()
    ax.plot([1, 10, 100], [1, 10, 100])
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    buf = io.StringIO()
    fig.savefig(buf, format="svg")
    plt.close(fig)
    root = sv.load(io.StringIO(buf.getvalue()))
    with pytest.raises(ValueError, match="not linear"):
        sv.scale(root, "y")


def test_a_round_trip_through_matplotlib_comes_back_unchanged():
    """Draw known numbers, read them back, get the same numbers."""
    want = [3.5, 41.25, 98.125]
    fig, ax = plt.subplots()
    ax.bar([0, 1, 2], want, edgecolor="black")
    buf = io.StringIO()
    fig.savefig(buf, format="svg")
    plt.close(fig)
    got = [v for _, v in sv.bars(sv.load(io.StringIO(buf.getvalue())))]
    assert got == pytest.approx(want, abs=1e-6)


# ------------------------------------------------- figure 4g and supp 13a

@pytest.mark.parametrize("panel", ["4g", "S13a"])
def test_the_shipped_table_still_matches_the_artwork(panel):
    artwork = sa.SPEC[panel][0]
    assert artwork.exists()
    table = sa.source(panel).set_index("Allele")
    name = {"4g": "cdk_variant_T47D_plot_t0_tf.svg",
            "S13a": "cdk_variant_MCF7_plot_t0_tf.svg"}[panel]
    for label, points in sv.labelled_series(ARTWORK / name):
        pre, post = points[0][1], points[1][1]
        assert table.loc[label, "Pre-Selection"] == pytest.approx(pre, abs=1e-6)
        assert table.loc[label, "Post-Selection"] == pytest.approx(post, abs=1e-6)


@pytest.mark.parametrize("panel", ["4g", "S13a"])
def test_the_panel_draws_the_eight_variants_in_the_published_order(panel):
    fig, table = sa.panel(panel)
    assert list(table["Allele"]) == [
        "CDK2 S261G", "CDK4 I51V", "CDK6 M54T", "CDK6 L56P",
        "CDK6 G53S", "CDK6 M54I", "CDK6 R46W", "CDK6 R44C"]
    assert len(fig.axes) == 1
    plt.close(fig)


@pytest.mark.parametrize("panel", ["4g", "S13a"])
def test_the_redraw_plots_the_values_the_table_holds(panel):
    """Draw it, read the drawing back, and compare with the table."""
    fig, table = sa.panel(panel)
    buf = io.StringIO()
    fig.savefig(buf, format="svg")
    plt.close(fig)
    got = {label: points for label, points
           in [(lab, pts) for lab, pts in
               sv.labelled_series(io.StringIO(buf.getvalue()))]}
    for _, row in table.iterrows():
        (_, pre), (_, post) = got[row["Allele"]]
        assert pre == pytest.approx(row["Pre-Selection"], abs=1e-6)
        assert post == pytest.approx(row["Post-Selection"], abs=1e-6)


def test_the_t47d_values_are_whole_percentages():
    """As figure 4f's hand-quantified table is, which is what it should be."""
    table = sa.source("4g")
    for column in ("Pre-Selection", "Post-Selection"):
        assert (table[column] % 1 == 0).all()


def test_the_one_fractional_value_is_the_mcf7_cdk2_variant_after_selection():
    """Recorded so it reads as measured rather than as a slip."""
    table = sa.source("S13a").set_index("Allele")
    fractional = table[(table % 1 != 0).any(axis=1)]
    assert list(fractional.index) == ["CDK2 S261G"]
    assert fractional.loc["CDK2 S261G", "Post-Selection"] == pytest.approx(4.264871)
    assert fractional.loc["CDK2 S261G", "Pre-Selection"] == 0


def test_the_recovery_script_regenerates_what_is_on_disk():
    out = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "recover_selection_tables.py"),
         "--check"], capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "DIFFERS" not in out.stdout
