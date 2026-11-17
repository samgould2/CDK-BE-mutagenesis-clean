"""The editing rebuild must agree with the shipped tables on the variant call.

The per-sample files in `data/raw/screen_editing/` are a redo CRISPResso run;
the shipped tables came from an earlier run that is not in the repository, so
read counts differ by a couple of percent. What must not differ is which
variant each guide is called as -- that is what every downstream figure uses.
"""
import numpy as np
import pytest

from cdkbe import counts, editing, editing_pipeline as ep

SCREEN = "CBE_subpool1"


def _call_from(full_canonical, editor):
    ref, alt = editing.EDITS[editor]
    keep = full_canonical["DNA Change"].str.findall(editing._DNA_CHANGE).apply(
        lambda e: all((r == ref) and (a == alt) for _, r, a in e))
    h = (full_canonical[keep].groupby(["gRNA_id", "HGVSp"], as_index=False)["#Reads"].sum()
         .sort_values(["gRNA_id", "#Reads"], ascending=[True, False]).reset_index(drop=True))
    h["%Reads"] = h["#Reads"] / h.groupby("gRNA_id")["#Reads"].transform("sum") * 100
    return editing.top_hgvsp(h)


@pytest.fixture(scope="module")
def pair():
    key, library = counts.sample_key(), counts.library()
    shipped = editing.top_hgvsp(editing.hgvsp_from_full_canonical(SCREEN))
    rebuilt = _call_from(ep.build(SCREEN, key, library)["full_canonical"], SCREEN[:3])
    return shipped, rebuilt


def test_same_guides(pair):
    shipped, rebuilt = pair
    assert set(shipped["gRNA_id"]) == set(rebuilt["gRNA_id"])


def test_same_variant_call(pair):
    shipped, rebuilt = pair
    m = shipped[["gRNA_id", "HGVSp"]].merge(rebuilt[["gRNA_id", "HGVSp"]],
                                            on="gRNA_id", suffixes=("_s", "_r"))
    agree = (m["HGVSp_s"] == m["HGVSp_r"]).mean()
    assert agree == 1.0, f"only {agree:.1%} of guides get the same variant call"


def test_wt_fraction_tracks(pair):
    shipped, rebuilt = pair
    m = shipped[["gRNA_id", "WT_perc"]].merge(rebuilt[["gRNA_id", "WT_perc"]],
                                              on="gRNA_id", suffixes=("_s", "_r"))
    assert np.corrcoef(m["WT_perc_s"], m["WT_perc_r"])[0, 1] > 0.99


def test_canonical_is_edit_not_window():
    """The canonical filter is Canonical_edit alone; Canonical_window is not applied."""
    t = ep.build(SCREEN)["full_canonical"]
    assert t["Canonical_edit"].all()
    assert not t["Canonical_window"].all(), \
        "canonical table should still contain out-of-window canonical edits"
