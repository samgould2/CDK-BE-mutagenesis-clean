"""Figure 4p: the CDK6 clinical resistance variants.

The published panel shows five guides and the notebook's own selection -- top
15 by palbociclib LFC among guides editing a clinically reported codon --
produces fourteen. The step that was missing is a significance filter, and
these tests pin that it reproduces the panel exactly: the five guides, one per
variant, in the order the figure draws them.

(This was figure 4n before the November relabelling; `docs/DISCREPANCIES.md`
items 5b and 28.)
"""
from __future__ import annotations

import pytest

from cdkbe.source_data import screen_figures as sf

#: as printed on the published panel, top to bottom
PUBLISHED = ["V47A", "R46W", "R46W_Q48*", "R168H", "A167T"]


def variants(hits):
    return [v.split("\n")[1] for v in hits["label"]]


def test_the_panel_is_the_published_five_in_the_published_order():
    assert variants(sf._clinical_hits()) == PUBLISHED


def test_one_guide_per_variant():
    """Why the panel has no duplicates: the filter leaves only one of each."""
    hits = sf._clinical_hits()
    assert len(hits) == len(set(hits["key"])) == 5
    assert len(set(variants(hits))) == 5


def test_without_the_filter_the_notebook_rule_gives_fourteen():
    """The state this repository was in, and what it could not explain."""
    unfiltered = sf._clinical_hits(max_fdr=1.1)
    assert len(unfiltered) == 14
    for extra in ("Q48*", "R46Q", "R168C"):
        assert extra in variants(unfiltered)
        assert extra not in PUBLISHED


def test_the_threshold_has_to_be_point_one():
    """At 0.05 nothing clears; the smallest FDR among the candidates is 0.053."""
    assert sf.CLINICAL_MAX_FDR == 0.1
    assert len(sf._clinical_hits(max_fdr=0.05)) == 0
    best = sf._clinical_hits(max_fdr=1.1)["best_FDR"]
    assert best.min() > 0.05
    # and the set is not sensitive to exactly where the line falls: the sixth
    # best is R168C at 0.137
    assert len(sf._clinical_hits(max_fdr=0.13)) == 5


def test_all_six_drug_arms_count_not_only_the_three_the_panel_draws():
    """R46W's significant arm is not a CDK4/6 inhibitor."""
    assert set(sf.CLINICAL_CONDITIONS) == {"PALBO", "RIBO", "ABEMA",
                                           "ATIRMO", "INX", "TAGTO"}
    assert set(sf.CDK46_CONDITIONS) < set(sf.CLINICAL_CONDITIONS)
    hits = sf._clinical_hits().set_index("key")
    cdk46 = [f"FDR_{c}" for c in sf.CDK46_CONDITIONS]
    r46w = hits.loc["gRNA_CDK6_targ_8422_CBE"]
    assert r46w[cdk46].min() > 0.5
    assert r46w["best_FDR"] < 0.1


def test_the_rows_are_ordered_by_palbociclib_lfc():
    lfc = sf._clinical_hits()["LFC_median_PALBO"].tolist()
    assert lfc == sorted(lfc, reverse=True)
