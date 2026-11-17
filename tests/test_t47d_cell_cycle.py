"""Figure 4k, 4l and supplementary 13c, 13d.

One flow experiment, four panels. The delivered package carried its FlowJo
export, the notebook that drew the panels and the published SVGs, so these
tests compare the redraw against the manuscript's own artwork character for
character and fill for fill -- the SVGs are vendored under
`reference/t47d_cell_cycle/`.

The thing most worth pinning is the Bonferroni family. The delivered notebook
corrects across **every comparison drawn on one panel**, two variants by two
treatments, so the multiplier is four. Correcting within each treatment
instead would halve it and turn CDK2 S261G's DMSO bracket from `ns` into a
star.
"""
from __future__ import annotations

import io
import re
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cdkbe.source_data import t47d_cell_cycle as cc  # noqa: E402

REFERENCE = Path(__file__).resolve().parent / "reference" / "t47d_cell_cycle"
PANELS = ["4k", "4l", "S13c", "S13d"]


def published(panel: str) -> str:
    return (REFERENCE / f"published_{cc._STEM[panel]}.svg").read_text()


def redrawn(panel: str) -> str:
    # Start from matplotlib's own defaults rather than whatever style an
    # earlier test left behind: `plotting.use_style()` saves transparently,
    # which drops the figure and axes background patches and so changes the
    # list of fills this file compares.
    rc = {k: v for k, v in matplotlib.rcParamsDefault.items()
          if k not in {"backend", "interactive"}}
    rc["svg.fonttype"] = "none"
    with matplotlib.rc_context(rc):
        fig, _ = cc.panel(panel)
        buf = io.StringIO()
        fig.savefig(buf, format="svg", bbox_inches="tight")
        plt.close(fig)
    return buf.getvalue()


def texts(svg: str) -> list:
    return re.findall(r"<text[^>]*>([^<]*)</text>", svg)


def fills(svg: str) -> list:
    return re.findall(r'style="fill: #([0-9a-f]{6})', svg)


@pytest.mark.parametrize("panel", PANELS)
def test_every_label_and_star_matches_the_published_panel(panel):
    assert texts(redrawn(panel)) == texts(published(panel))


@pytest.mark.parametrize("panel", PANELS)
def test_every_fill_matches_the_published_panel(panel):
    assert fills(redrawn(panel)) == fills(published(panel))


def test_the_export_is_six_lines_by_three_treatments_by_three_wells():
    d = cc.replicates()
    assert len(d) == 54
    assert set(d["CellLine"]) == set(cc.ABE_LINES) | set(cc.CBE_LINES)
    assert set(d["Treatment"]) == {"DMSO", "Palbo", "ss"}
    assert d.groupby(["CellLine", "Treatment"]).size().unique().tolist() == [3]


def test_the_four_gates_do_not_sum_to_a_hundred():
    """They are siblings on one parent, which is why the stack rescales."""
    total = cc.replicates()[cc.PHASE_ORDER].sum(axis=1)
    assert (total < 100).all()
    assert total.min() > 90


def test_the_bonferroni_family_is_the_whole_panel():
    for lines, control in ((cc.ABE_LINES, "ABE_NT"), (cc.CBE_LINES, "CBE_NT")):
        sig = cc.significance(lines, control)
        assert len(sig) == 4
        assert set(sig["m"]) == {4}
        assert sig["p_adj"].tolist() == pytest.approx(
            (sig["p_value"] * 4).clip(upper=1.0).tolist())


def test_the_published_brackets_reproduce():
    """Exactly what the two %S-phase panels print, in their drawn order."""
    abe = cc.significance(cc.ABE_LINES, "ABE_NT").set_index(
        ["CellLine", "Treatment"])["sig"]
    assert abe[("CDK2_S261G", "DMSO")] == "ns"
    assert abe[("CDK4_I51V", "DMSO")] == "**"
    assert abe[("CDK2_S261G", "Palbo")] == "**"
    assert abe[("CDK4_I51V", "Palbo")] == "**"
    cbe = cc.significance(cc.CBE_LINES, "CBE_NT").set_index(
        ["CellLine", "Treatment"])["sig"]
    assert cbe[("CDK6_G53S", "DMSO")] == "ns"
    assert cbe[("CDK6_R46W", "DMSO")] == "ns"
    assert cbe[("CDK6_G53S", "Palbo")] == "*"
    assert cbe[("CDK6_R46W", "Palbo")] == "**"


def test_correcting_within_each_treatment_would_change_a_call():
    """Why the family matters, pinned so it cannot be loosened by accident."""
    sig = cc.significance(cc.ABE_LINES, "ABE_NT").set_index(
        ["CellLine", "Treatment"])
    raw = sig.loc[("CDK2_S261G", "DMSO"), "p_value"]
    assert raw * 2 < 0.05 < raw * 4


def test_the_cdk6_lines_are_relabelled_because_they_carry_a_second_edit():
    assert cc._label("CDK6_G53S") == "CDK6 G53S_M54I"
    assert cc._label("CDK6_R46W") == "CDK6 R44C_R46W"
    assert cc._label("CDK2_S261G") == "CDK2 S261G"


def test_serum_starvation_is_measured_but_not_drawn():
    assert "ss" in set(cc.replicates()["Treatment"])
    assert cc.TREATMENT_ORDER == ["DMSO", "Palbo"]
    for panel in PANELS:
        assert "serum free" not in texts(published(panel))


def test_the_delivered_package_is_in_the_repository():
    assert cc.DATA_CSV.is_file()
    for name in ("cell_cycle_stacked_barplot.ipynb", "figure_style.py"):
        assert (cc.FOLDER / name).is_file(), name


def test_every_panel_is_registered_once():
    assert set(cc.MAIN_FIGURES) == {"4k", "4l"}
    assert set(cc.FIGURES) == {"S13c", "S13d"}
    assert set(cc._STEM) == set(cc.GROUPS) == {"4k", "4l", "S13c", "S13d"}
