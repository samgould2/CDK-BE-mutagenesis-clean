"""Figures 2d, 2k and 4f: the three main panels recovered last.

Each was traced to the cell that drew it by its own `fig.savefig` name, which is
the file `docs/source_data_mapping_main.csv` lists for that panel -- so the
provenance is the author's own code, not a resemblance:

    2d  FACS_dose_response_2024_8_22.ipynb cells 126 and 127
        -> figures/L156F_KI-CDK9d-32_wide.pdf, figures/L156F_KB-0742_wide.pdf
    2k  the same notebook, cell 130      -> figures/L18F_SY_wide.pdf
    4f  analysis.ipynb cell 137          -> figures/CDK4_validation_v3.pdf

2d and 2k turn out to be the same plate and the same fit as supplementary
9d-9f, with two of the four lines and a shorter canvas. 4f is different in kind:
its numbers were literals in the notebook, so the test here is that the
transcription still matches them.

`tests/reference/main_panels/` holds the text of each published PDF, extracted
once and vendored, so the comparison keeps working without the original tree.
The published panels carry only axis labels, so text identity is a check on the
axes and scales -- the fits are pinned separately, by value.
"""
from __future__ import annotations

import collections
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from cdkbe.source_data import allele_frequency as af  # noqa: E402
from cdkbe.source_data import allele_specific_mtt as mtt  # noqa: E402

REF = Path(__file__).resolve().parent / "reference" / "main_panels"

needs_pdftotext = pytest.mark.skipif(
    shutil.which("pdftotext") is None, reason="poppler's pdftotext not installed")


def tokens_of(pdf) -> collections.Counter:
    out = subprocess.run(["pdftotext", "-q", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    return collections.Counter(out.split())


def reference(stem: str) -> collections.Counter:
    return collections.Counter((REF / f"{stem}.txt").read_text().split())


# ------------------------------------------------------------------ 2d / 2k

def test_main_curves_use_the_supplementary_plate():
    """2d/2k and supplementary 9d-9f are one experiment, so one plate."""
    assert mtt.PLATE == "2025_8_10.csv"
    for panel, specs in mtt.MAIN_CURVES.items():
        for drug, _, cts, _, _ in specs:
            assert set(cts) <= set(mtt.CELL_TYPES), (panel, drug)
            assert len(cts) == 2, f"{panel}/{drug} draws two lines, not {len(cts)}"


#: LC50 in nM per (panel, drug, cell line), from the shipped plate and the
#: notebook's cell-125 fit. Pinned so a change in the fit or the plate shows up.
LC50 = {
    ("2d", "KI-CDK9d-32", "NT"): 65.2,
    ("2d", "KI-CDK9d-32", "L156F"): 203.3,
    ("2d", "KB-0742", "NT"): 1414.4,
    ("2d", "KB-0742", "L156F"): 1141.7,
    ("2k", "SY-5609", "NT"): 6.8,
    ("2k", "SY-5609", "L18F"): 45.5,
}


def test_main_panels_are_not_in_the_supplementary_registry():
    """2d/2k share a module with supplementary 9d-9f; they must not leak."""
    assert set(mtt.MAIN_FIGURES) == {"2d", "2k"}
    assert not {k for k in mtt.FIGURES if not k.startswith("S")}
    assert not {k for k in mtt.PANELS if not k.startswith("S")}


@pytest.mark.parametrize("panel", sorted(mtt.MAIN_CURVES))
def test_lc50s_reproduce(panel):
    out = mtt.MAIN_FIGURES[panel]()
    for drug, (fig, params) in out.items():
        plt.close(fig)
        for _, r in params.iterrows():
            want = LC50[(panel, drug, r["Cell Line"])]
            assert round(r["LC50_nM"], 1) == pytest.approx(want, abs=0.15), \
                f"{panel}/{drug}/{r['Cell Line']}: {r['LC50_nM']:.1f} != {want}"


def test_the_resistance_direction_is_the_published_one():
    """CDK9 L156F resists KI-CDK9d-32 and CDK7 L18F resists SY-5609.

    KB-0742 is the stated exception: L156F cells are slightly MORE sensitive,
    which is the point figure 2d makes, so a sign flip there is a real failure
    rather than a rounding difference.
    """
    assert LC50[("2d", "KI-CDK9d-32", "L156F")] > LC50[("2d", "KI-CDK9d-32", "NT")]
    assert LC50[("2k", "SY-5609", "L18F")] > LC50[("2k", "SY-5609", "NT")]
    assert LC50[("2d", "KB-0742", "L156F")] < LC50[("2d", "KB-0742", "NT")]


@needs_pdftotext
@pytest.mark.parametrize("panel", sorted(mtt.MAIN_CURVES))
def test_curves_match_the_published_pdfs(panel, tmp_path):
    for drug, xlim, cts, ylim, stem in mtt.MAIN_CURVES[panel]:
        fig, _ = mtt.curve_panel(drug, xlim, cell_types=cts, ylim=ylim,
                                 figsize=(5, 3))
        out = tmp_path / f"{stem}.pdf"
        fig.savefig(out, transparent=True)
        plt.close(fig)
        assert tokens_of(out) == reference(stem), f"{panel}/{drug} vs {stem}"


def test_the_supplementary_panels_are_unchanged():
    """The four-line panels must not have moved when 2d/2k were added."""
    fig, params = mtt.curve_panel(*mtt.CURVES["S9d"])
    plt.close(fig)
    assert sorted(params["Cell Line"]) == sorted(mtt.CELL_TYPES)
    assert fig.get_size_inches().tolist() == [5.0, 6.0]


# ----------------------------------------------------------------------- 4f

def test_4f_table_is_the_notebooks_literals():
    d = af.source()
    assert list(d.columns) == ["Allele", "Pre-Selection", "Post-Selection"]
    assert len(d) == 11
    # the first and last of the three lists in cell 137, as a spot check that
    # the transcription kept its row order
    assert d.iloc[0].tolist() == ["CDK4 E56K", 54, 72]
    assert d.iloc[-1].tolist() == ["CDK6 E21K", 74, 77]
    assert set(d["Allele"].str.split().str[0]) == {"CDK2", "CDK4", "CDK6"}


def test_4f_every_variant_is_enriched_by_selection():
    d = af.source()
    assert (d["Post-Selection"] >= d["Pre-Selection"]).all()


@needs_pdftotext
def test_4f_matches_the_published_pdf(tmp_path):
    fig, _ = af.draw()
    out = tmp_path / "CDK4_validation_v3.pdf"
    fig.savefig(out, transparent=True)
    plt.close(fig)
    assert tokens_of(out) == reference("CDK4_validation_v3")
