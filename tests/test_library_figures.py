"""Figure 1's panels.

1f carries a discrepancy worth pinning: the published PDF puts all nine panels
on one left edge, while the notebook cell that names that file offsets them to
align the kinase domains. The figure predates the change, so both layouts are
supported and both are checked here.
"""
import re
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from cdkbe import plotting
from cdkbe.source_data import library_figures as lf


@pytest.fixture(scope="module", autouse=True)
def style():
    plotting.use_style()


def _label_spread(fig, tmp_path):
    p = tmp_path / "f.pdf"
    fig.savefig(str(p))
    out = subprocess.run(["pdftotext", "-bbox", str(p), "-"],
                         capture_output=True, text=True).stdout
    xs = [float(m.group(1))
          for m in re.finditer(r'xMax="([\d.]+)"[^>]*>CDK\d+</word>', out)]
    return max(xs) - min(xs)


def test_published_layout_shares_one_left_edge(tmp_path):
    fig, _ = lf.figure_1f(align_domains=False)
    assert _label_spread(fig, tmp_path) < 1.0
    plt.close(fig)


def test_aligned_layout_offsets_the_panels(tmp_path):
    fig, _ = lf.figure_1f(align_domains=True)
    assert _label_spread(fig, tmp_path) > 100
    plt.close(fig)


def test_1f_draws_a_panel_per_gene():
    fig, d = lf.figure_1f()
    assert len(fig.get_axes()) == len(lf.GENES)
    assert set(d["gene"]) <= set(lf.GENES)
    plt.close(fig)


def test_1e_three_panels_and_a_ceiling_per_gene():
    fig, d = lf.figure_1e(coverage=lf.lib.figure_1e(cutoffs=[0, 50, 100]))
    assert len(fig.get_axes()) == 3
    plt.close(fig)
    designed = lf.designed_coverage()
    assert set(designed.index) == set(lf.GENES)
    # the designed library cannot cover more than every codon
    assert (designed[["combined", "ABE", "CBE"]] <= 100).all().all()


def test_measured_coverage_exceeds_the_designed_library():
    """The dashed line is the designed library, not a ceiling.

    Measured coverage sits ABOVE it at low cutoffs, by 8-19% depending on the
    gene, because bystander editing reaches codons no guide was designed for.
    That gap is the point of the panel, so a regression that clamped the
    measured curve to the design would look plausible and be wrong.
    """
    d = lf.lib.figure_1e(cutoffs=[0]).set_index("Gene")
    designed = lf.designed_coverage()
    for g in lf.GENES:
        excess = d.loc[g, "total_codons"] - designed.loc[g, "combined"]
        assert excess > 0, f"{g}: measured coverage did not exceed the design"
        assert excess < 40, f"{g}: implausible excess of {excess:.0f} points"
