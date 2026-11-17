"""The PC3 panels must keep printing the numbers the published figures print.

Every value below was read off the panel PDFs the re-analysis scripts saved
(`reanalysis/figures/Fig2b_DEG_counts_with_baseline.pdf`,
`Fig3_response_retained.pdf`, `Fig10_elongation_combined.pdf`), not off the
tables that feed them -- so this fails if a table is swapped for a different
run, as well as if the drawing code drifts.
"""
import matplotlib
matplotlib.use("Agg")

import pytest

from cdkbe.source_data import rnaseq, rnaseq_figures as rf


#: (drug, genotype) -> (up, down), as printed on the bars of 3f
DEG_LABELS = {
    ("HQ461", "NT"): (3096, 3826), ("HQ461", "I733V"): (295, 606),
    ("HQ461", "G822E"): (517, 641), ("HQ461", "G731K"): (44, 115),
    ("ZLY025", "NT"): (2603, 4163), ("ZLY025", "I733V"): (2068, 1895),
    ("ZLY025", "G822E"): (1964, 1873), ("ZLY025", "G731K"): (47, 134),
    ("Control_vs_NT", "I733V"): (21, 62), ("Control_vs_NT", "G822E"): (244, 408),
    ("Control_vs_NT", "G731K"): (230, 360),
}

#: (drug, genotype) -> the rounded percentage printed above each bar of 3g
PCT_LABELS = {
    ("HQ461", "NT"): 100, ("HQ461", "I733V"): 15,
    ("HQ461", "G822E"): 25, ("HQ461", "G731K"): 10,
    ("ZLY025", "NT"): 100, ("ZLY025", "I733V"): 61,
    ("ZLY025", "G822E"): 49, ("ZLY025", "G731K"): 5,
}


def test_deg_counts_match_the_panel():
    d = rnaseq.figure_3f()
    got = {(r.drug, r.genotype): (int(r.up), int(r.down)) for r in d.itertuples()}
    assert got == DEG_LABELS


def test_core_signature_matches_the_panel_title():
    """The 3g title reads "2,098 up / 3,068 down genes"."""
    up, dn = rf.core_nt_signature()
    assert (len(up), len(dn)) == (2098, 3068)


@pytest.mark.parametrize("key,expected", sorted(PCT_LABELS.items()))
def test_pct_response_matches_the_panel(key, expected):
    drug, geno = key
    pct = rnaseq.figure_3g()["drug_response_retained"]
    v = pct[(pct.genotype == geno) & (pct.drug == drug)].pct_response
    assert len(v) == 3, f"{geno}/{drug}: {len(v)} replicates, expected 3"
    assert round(v.mean()) == expected


def test_decile_profiles_are_complete():
    """10 deciles x (4 genotypes x 2 drugs + 3 variants at baseline) = 110 rows."""
    d = rnaseq.figure_3h()
    assert len(d) == 110
    assert sorted(d.decile.unique()) == list(range(10))
    assert set(d.group) == {"drug", "baseline"}
    assert (d.n_replicates == 3).all()
    # the panel's y-limits are set in the drawing code; nothing may fall outside
    drug, base = d[d.group == "drug"], d[d.group == "baseline"]
    assert drug.mean_log2ratio.between(-0.45, 0.85).all()
    assert base.mean_log2ratio.between(-0.12, 0.16).all()


@pytest.mark.parametrize("panel,n_axes", [("3f", 3), ("3g", 1), ("3h", 3)])
def test_panels_draw(panel, n_axes):
    fig, obj = rf.FIGURES[panel]()
    assert len(fig.axes) == n_axes
    assert len(obj) > 0
