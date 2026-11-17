"""Figure 1's panels, neither of which had a saved intermediate.

1e's `plotting_df` was built in a loop over 101 cutoffs x 9 genes into a
pre-allocated wide frame and never written to disk, so it is a full
re-derivation; these values come from replaying that loop.
"""
import pytest

from cdkbe.source_data import library as lb

#: total_codons at a few cutoffs, from the original notebook's logic
EXPECTED = {
    0:  {"CDK12": 98.3, "CDK13": 97.8, "CDK19": 96.4, "CDK2": 99.3, "CDK4": 99.7,
         "CDK6": 99.4, "CDK7": 97.7, "CDK8": 99.4, "CDK9": 100.0},
    10: {"CDK12": 90.3, "CDK13": 89.7, "CDK19": 86.1, "CDK2": 94.0, "CDK4": 98.0,
         "CDK6": 91.7, "CDK7": 83.5, "CDK8": 90.3, "CDK9": 97.6},
}


@pytest.fixture(scope="module")
def coverage():
    return lb.figure_1e(cutoffs=sorted(EXPECTED))


@pytest.mark.parametrize("cutoff", sorted(EXPECTED))
def test_codon_coverage(coverage, cutoff):
    d = coverage[coverage["editing_cutoff"] == cutoff].set_index("Gene")
    for gene, want in EXPECTED[cutoff].items():
        assert round(float(d.loc[gene, "total_codons"]), 1) == want, gene


def test_coverage_falls_as_the_cutoff_rises():
    d = lb.figure_1e(cutoffs=[0, 25, 50]).pivot(
        index="editing_cutoff", columns="Gene", values="total_codons")
    assert (d.diff().dropna() <= 0).all().all()


def test_resistance_variants():
    d = lb.figure_1f()
    assert len(d) == 97, "97 resistance variants across the nine CDKs"
    assert set(d["Editor"]) == {"ABE", "CBE"}
    assert (d["Codon"] != "WT").all()
    assert d["max_LFC"].min() > lb.LFC_CUTOFF
