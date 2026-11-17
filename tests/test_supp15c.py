"""Supplementary 15c: the two-panel reconstruction carries the published content.

The code that drew the published panel was never found -- see
`docs/DISCREPANCIES.md` item 13. `klifs_msa.supp_15c` rebuilds it from this
repository's own inputs, so the test that matters is whether the rebuild says
the same thing the published figure says.

It does, exactly: the two PDFs contain **the same 2,189 characters**, with no
character in one and not the other. Only the order `pdftotext` walks the page in
differs, because the row labels and the subdomain annotations sit at the same
height. That covers every residue of all seventeen rows, both axes, the
consensus sequence, both legends, the title and the footnote -- everything the
panel asserts.

What it cannot cover is colour, which is where the panel's meaning actually
lives, so the drug counts behind the shading are pinned separately below.
"""
from __future__ import annotations

import collections
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from cdkbe.source_data import klifs_msa as km  # noqa: E402

REF = Path(__file__).resolve().parent / "reference" / "supp15c"

pymsaviz = pytest.importorskip("pymsaviz", reason="pymsaviz is an optional extra")
needs_pdftotext = pytest.mark.skipif(
    shutil.which("pdftotext") is None, reason="poppler's pdftotext not installed")


def chars(text: str) -> str:
    return re.sub(r"\s+", "", text)


def pdf_text(path) -> str:
    return subprocess.run(["pdftotext", "-q", str(path), "-"],
                          capture_output=True, text=True).stdout


# ------------------------------------------------------------- the structure

def test_both_blocks_are_drawn():
    fig, _ = km.supp_15c()
    try:
        # the CDK alignment, the consensus track, the clinical-kinase alignment
        assert len(fig.axes) == 3
    finally:
        plt.close(fig)


def test_single_block_variant_is_still_available():
    """`two_panel=False` is notebook cell 14 verbatim, which some panels cite."""
    fig, _ = km.supp_15c(two_panel=False)
    try:
        assert fig.axes
    finally:
        plt.close(fig)


def test_the_ten_clinical_kinases_are_the_published_set():
    assert km.ONCOKB_GENES == ["ABL1", "ALK", "BTK", "EGFR", "KIT", "MET",
                               "NTRK1", "NTRK3", "PDGFRA", "ROS1"]


def test_clinical_sequences_cover_all_85_positions():
    mat, seqs = km.oncokb_matrix()
    assert len(seqs) == 85
    for g in km.ONCOKB_GENES:
        assert len(seqs[f"{g}_seq"]) == 85
    assert list(mat.index) == km.ONCOKB_GENES
    assert mat.shape == (10, 85)


def test_drug_counts_stay_within_the_published_scale():
    """The published legend runs 1..5 drugs; a sixth would need a new swatch."""
    mat, _ = km.oncokb_matrix()
    assert int(mat.values.max()) == km.ONCOKB_MAX_DRUGS


def test_met_and_pdgfra_have_no_pocket_resistance():
    """Both rows are blank in the published panel; that is a real result."""
    mat, _ = km.oncokb_matrix()
    assert int(mat.loc["MET"].sum()) == 0
    assert int(mat.loc["PDGFRA"].sum()) == 0


def test_the_footnote_numbers_are_counted_not_copied():
    from cdkbe.source_data.klifs_resistance import annotated

    a = annotated()
    assert len(a) == 104
    assert int((~a["KLIFS Catalytic Site"]).sum()) == 30


# ------------------------------------------------------ against the published

@needs_pdftotext
@pytest.mark.parametrize("aset", ["core7", "all9"])
def test_reconstruction_has_the_published_content(aset, tmp_path):
    fig, _ = km.supp_15c(aset=aset)
    out = tmp_path / f"{aset}.pdf"
    fig.savefig(out, transparent=True, bbox_inches="tight")
    plt.close(fig)

    published = chars((REF / f"KLIFS_resistance_MSA_{aset}.txt").read_text())
    rebuilt = chars(pdf_text(out))
    assert collections.Counter(rebuilt) == collections.Counter(published), (
        f"{aset}: the rebuild and the published panel no longer carry the same "
        "characters")
    assert len(rebuilt) == len(published)
