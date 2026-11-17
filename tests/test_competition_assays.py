"""Figure 2f and supplementary 7c, 7e, 10f, 10g, 11b, 11i.

Five of these seven are heatmaps whose published SVG was delivered with the
data, and those SVGs are vendored under `reference/competition_assays/`. A
matplotlib SVG written with `svg.fonttype: none` carries every cell's fill as a
hex string and every annotation as real text, so the comparison here is
cell-for-cell and character-for-character against the manuscript's own artwork
-- not against a previous run of this code.

Four of the five come back identical. The fifth, supplementary 11b, differs in
its four fill colours and nothing else: the published LNCaP panel was centred
near 1.38 rather than 1.0, which this file pins so the difference cannot be
mistaken for drift. `docs/DISCREPANCIES.md` item 23 asks which is wanted.

This repository has no published figure for supplementary 7c and 7e, so they
are checked
against their numbers only.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cdkbe.source_data import competition_assays as ca  # noqa: E402

REFERENCE = Path(__file__).resolve().parent / "reference" / "competition_assays"
HEATMAPS = ["2f", "S10f", "S10g", "S11b", "S11I"]
#: the one panel whose published colours this code does not reproduce
DIFFERENT_COLOURS = {"S11b"}


def published(panel: str) -> str:
    return (REFERENCE / f"published_{panel}.svg").read_text()


def redrawn(panel: str) -> str:
    import io

    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig, _ = ca.heatmap(panel)
        buf = io.StringIO()
        fig.savefig(buf, format="svg", bbox_inches="tight")
        plt.close(fig)
    return buf.getvalue()


def cells(svg: str) -> list:
    block = re.search(r'id="QuadMesh_1">(.*?)\n   </g>', svg, re.S)
    return re.findall(r"fill: #([0-9a-f]{6}); stroke", block.group(1))


def texts(svg: str) -> list:
    return re.findall(r"<text[^>]*>([^<]*)</text>", svg)


@pytest.mark.parametrize("panel", HEATMAPS)
def test_every_label_and_number_matches_the_published_panel(panel):
    assert texts(redrawn(panel)) == texts(published(panel))


@pytest.mark.parametrize("panel", [p for p in HEATMAPS
                                   if p not in DIFFERENT_COLOURS])
def test_every_cell_colour_matches_the_published_panel(panel):
    assert cells(redrawn(panel)) == cells(published(panel))


def test_the_lncap_panel_is_the_one_with_a_different_colour_scale():
    """Pinned, not fixed: the numbers agree and only the four fills differ."""
    assert texts(redrawn("S11b")) == texts(published("S11b"))
    assert cells(redrawn("S11b")) != cells(published("S11b"))
    assert cells(published("S11b")) == ["f6f7f7", "e3edf3", "f9efe9", "ef9979"]


def test_the_family_shares_one_colour_scale():
    """0 to 5.5 centred at 1, which is what the published colourbars show."""
    assert (ca.VMIN, ca.VMAX, ca.CENTER, ca.CMAP) == (0.0, 5.5, 1.0, "RdBu_r")


@pytest.mark.parametrize("panel", HEATMAPS)
def test_three_replicates_in_every_cell(panel):
    d = ca.source(panel)
    assert set(d["n"]) == {3}


def test_the_pc3_panel_agrees_with_the_summary_its_author_computed():
    """The only delivered file that carries its own statistics -- so use them."""
    mine = ca.source("S11I").set_index(["variant", "treatment"])
    theirs = ca.pc3_delivered_summary().set_index(["variant", "treatment"])
    for key, row in theirs.iterrows():
        assert mine.loc[key, "mean"] == pytest.approx(row["mean"], abs=5e-5)
        assert mine.loc[key, "sd"] == pytest.approx(row["sd"], abs=5e-5)
        assert mine.loc[key, "sem"] == pytest.approx(row["sem"], abs=5e-5)
        if row["significance"] is not None:
            assert mine.loc[key, "significance"] == row["significance"]
            assert mine.loc[key, "p_vs_reference"] == pytest.approx(
                row["p_vs_NT"], rel=1e-2)


def test_the_two_cdk13_files_repeat_the_cdk12_files_control_arm():
    """Which is why supplementary 10f and 10g draw one NT row, not two."""
    for a, b, nt in (("A549/Cdk12 FC.txt", "A549/Cdk13 FC.txt", "NT"),
                     ("A549/Cdk12 FC ZLY.txt", "A549/Cdk13 FC ZLY.txt",
                      "ABE NT")):
        x = ca.matrix_file(a)
        y = ca.matrix_file(b)
        xs = x[x["variant"] == nt].reset_index(drop=True)
        ys = y[y["variant"] == nt].reset_index(drop=True)
        assert xs["fold_change"].tolist() == ys["fold_change"].tolist()


def test_the_zly_panel_drops_the_cbe_arm_the_file_carries():
    """The published panel's NT row is the ABE arm; the file holds both."""
    delivered = set(ca.matrix_file("A549/Cdk12 FC ZLY.txt")["variant"])
    assert {"ABE NT", "CBE NT"} <= delivered
    assert ca.PANEL_SPEC["S10g"]["variants"][0] == "ABE NT"
    assert ca.PANEL_SPEC["S10g"]["labels"][0] == "NT"
    assert "CBE NT" not in ca.PANEL_SPEC["S10g"]["variants"]


def test_figure_2f_annotates_neither_non_targeting_arm():
    """Both are controls; the published panel marks only L156F and I157V."""
    annot = ca._grid("2f")[1]
    stars = [t for row in annot for t in row if "\n" in t and t.count("\n") == 2]
    assert len(stars) == 2
    assert all(t.rsplit("\n", 1)[1] in {"***", "**"} for t in stars)


def test_figure_2f_is_measured_against_the_abe_arm():
    """Against the CBE arm I157V would carry one star, and the panel has two."""
    import scipy.stats

    d = ca.replicates("2f")

    def arm(variant, treatment):
        return d[(d["variant"] == variant)
                 & (d["treatment"] == treatment)]["fold_change"]

    vs_abe = scipy.stats.ttest_ind(arm("I157V", "KI-Cdk9d-32"),
                                   arm("ABE NT", "KI-Cdk9d-32")).pvalue
    vs_cbe = scipy.stats.ttest_ind(arm("I157V", "KI-Cdk9d-32"),
                                   arm("CBE NT", "KI-Cdk9d-32")).pvalue
    assert ca._stars(vs_abe) == "**"
    assert ca._stars(vs_cbe) == "*"


@pytest.mark.parametrize("panel", ["S7c", "S7e"])
def test_the_mouse_panels_have_three_replicates_per_arm(panel):
    d = ca.ball_source(panel)
    assert d.groupby("treatment").size().tolist() == [3, 3]
    t, p = d.attrs["ttest"]
    assert p < 0.01          # both arms separate; the panel is about that


@pytest.mark.parametrize("panel", ["S7c", "S7e"])
def test_the_mouse_panels_draw(panel):
    fig, _ = ca.ball_panel(panel)
    assert fig.axes
    plt.close(fig)


def test_every_panel_is_registered_once():
    both = set(ca.MAIN_FIGURES) | set(ca.FIGURES)
    assert both == set(ca._STEM)
    assert not (set(ca.MAIN_FIGURES) & set(ca.FIGURES))
    assert all(not k.startswith("S") for k in ca.MAIN_FIGURES)
    assert all(k.startswith("S") for k in ca.FIGURES)


def test_the_delivered_files_are_in_the_repository():
    """No panel here reaches outside the repo for its data."""
    for spec in ca.PANEL_SPEC.values():
        for name in spec["files"]:
            assert (ca.ASSAYS / name).is_file(), name
    assert (ca.ASSAYS / "README.txt").is_file()


def test_means_and_sems_are_over_the_replicates_not_read_from_a_file():
    """A spot check that the panel numbers really are computed here."""
    d = ca.replicates("S11b")
    x = d[(d["variant"] == "I733V")
          & (d["treatment"] == "HQ461")]["fold_change"].to_numpy()
    row = ca.source("S11b").set_index(["variant", "treatment"]).loc[
        ("I733V", "HQ461")]
    assert row["mean"] == pytest.approx(x.mean())
    assert row["sem"] == pytest.approx(x.std(ddof=1) / np.sqrt(3))
    assert f"{row['mean']:.2f}" == "3.17"
