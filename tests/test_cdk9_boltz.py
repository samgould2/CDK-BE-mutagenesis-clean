"""Supplementary 6e-6j: the CDK9 Boltz-2 interaction-fingerprint panels.

The strongest checks here are the two counts the analysis README states
independently of any figure -- 132 retained fingerprint features and 58
retained WT ATP-pocket features. Both come out of the port exactly, from the
108 copied structures, which is what says the fingerprint was carried over
faithfully rather than approximately.

The published panels are PNG/SVG with no extractable numbers, so the panels
themselves are pinned by the values they plot: the explained variances, the
feature sets at the ends of each PC, and the paired binding-probability
differences.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

pytest.importorskip("Bio", reason="biopython parses the predicted structures")

from cdkbe.boltz import ifp  # noqa: E402
from cdkbe.source_data import cdk9_boltz as B  # noqa: E402


# ------------------------------------------------------------- the dataset

def test_the_delivered_design_is_intact():
    m, a = ifp.metrics_and_affinity()
    assert len(m) == ifp.N_STRUCTURES == 108
    assert len(a) == ifp.N_AFFINITY == 36
    assert sorted(m["drug"].unique()) == sorted(ifp.DRUG_ORDER)
    assert sorted(m["variant"].unique()) == sorted(ifp.VARIANTS)
    assert sorted(m["seed"].unique()) == sorted(ifp.SEEDS)
    # three seeds x three structural samples for every ligand-variant pair
    assert (m.groupby(["drug", "variant"]).size() == 9).all()
    # affinity is recorded per seed, not per structural sample
    assert (a.groupby(["drug", "variant"]).size() == 3).all()


def test_the_two_atp_pocket_sets_are_not_the_same():
    """The notebook has two, and the PCA uses the larger one. Easy to mix up.

    Cell 15's labelled regions cover 32 residues; cell 19's flat set covers 37.
    The five extra are the ones no region names.
    """
    union = set().union(*ifp.REGIONS.values())
    assert len(union) == 32
    assert len(ifp.ATP_POCKET) == 37
    assert union < ifp.ATP_POCKET
    assert ifp.ATP_POCKET - union == {47, 49, 65, 67, 150}


# --------------------------------------------------------- the fingerprint

def test_feature_counts_match_the_analysis_readme():
    """132 retained features, 58 of them in the WT pocket -- the README's numbers."""
    F = ifp.fingerprint()
    assert F.shape == (108, 132)
    assert ifp.pocket_matrix("WT").shape == (36, 58)


def test_every_retained_feature_clears_the_floor():
    F = ifp.fingerprint()
    assert (F.max(axis=0) >= ifp.KEEP_MIN).all()


def test_scores_are_bounded_probabilities():
    F = ifp.fingerprint()
    assert F.values.min() >= 0.0
    assert F.values.max() <= 1.0


def test_the_cache_matches_a_rebuild_from_the_structures():
    """The cached matrix is a convenience; the structures are the source."""
    cached = ifp.fingerprint()
    fresh = ifp._build_fingerprint()
    assert list(cached.columns) == list(fresh.columns)
    assert np.allclose(cached.values, fresh.values, atol=1e-9)


# ----------------------------------------------------------------- the PCA

def test_explained_variance_matches_the_published_axes():
    _, var, _, _ = ifp.weighted_pca(ifp.pocket_matrix("WT"))
    assert round(var[0]) == 34
    assert round(var[1]) == 15


def test_the_axis_signs_are_anchored_not_left_to_lapack():
    """An unoriented fit may come back mirrored; the panels must not."""
    X = ifp.pocket_matrix("WT")
    for pc, (feature, want) in ifp.SIGN_ANCHORS.items():
        _, _, vt, _ = ifp.weighted_pca(X)
        assert np.sign(vt[pc][list(X.columns).index(feature)]) == want
    # and flipping is all that orientation does
    _, _, raw, _ = ifp.weighted_pca(X, orient=False)
    _, _, fixed, _ = ifp.weighted_pca(X)
    for pc in ifp.SIGN_ANCHORS:
        assert np.allclose(np.abs(raw[pc]), np.abs(fixed[pc]))


#: the ten features supplementary 6i plots, bottom to top
PC1_ENDS = ["ASP167:charged_polar", "ASP167:polar_sc", "PHE30:aromatic_C",
            "PHE30:carbon_contact", "LYS48:charged_polar", "ASP109:carbon_contact",
            "GLN27:carbon_contact", "PHE105:aromatic_C", "CYS106:polar_bb",
            "PHE105:carbon_contact"]

#: the ten features supplementary 6j plots, bottom to top
PC2_TOP = ["ASN154:polar_sc", "LYS151:charged_polar", "ASN154:steric_clash",
           "ASP109:carbon_contact", "ILE25:steric_clash", "ASP167:polar_bb",
           "LYS48:charged_polar", "LYS48:polar_sc", "CYS106:polar_bb",
           "ASP104:polar_bb"]


def test_6i_is_the_published_pc1_features():
    assert list(B.source_6i().index) == PC1_ENDS


def test_6j_is_the_published_pc2_features():
    assert list(B.source_6j().index) == PC2_TOP


def test_6i_and_6j_are_different_selections():
    """6i takes the two ends of PC1; 6j takes the ten largest on PC2."""
    v = B.source_6i()
    assert (v.values[:5] < 0).all() and (v.values[5:] > 0).all()
    pc2 = B.source_6j()
    full = ifp.pc_loadings(1)
    assert set(pc2.index) == set(full.abs().nlargest(10).index)


# ------------------------------------------------------- binding probability

def test_l156f_lowers_predicted_binding_for_every_ligand():
    """The result supplementary 6f reports; a sign flip here would be serious."""
    d = ifp.delta_binding_probability("L156F")
    assert len(d) == 12                      # 4 ligands x 3 matched seeds
    means = d.groupby("drug")["delta_binding_probability"].mean()
    assert (means < 0).all()
    # atuveciclib and the PROTAC move most, voruciclib least
    assert means.idxmin() in {"atuveciclib", "protac"}
    assert means.idxmax() == "voruciclib"


def test_delta_is_paired_within_a_seed():
    d = ifp.delta_binding_probability("L156F")
    assert np.allclose(d["delta_binding_probability"], d["L156F"] - d["WT"])
    assert sorted(d["seed"].unique()) == sorted(ifp.SEEDS)


def test_seed_level_joins_one_to_one():
    s = ifp.seed_level("WT")
    assert len(s) == 12
    assert s["iptm"].between(0, 1).all()
    assert s["binding_probability"].between(0, 1).all()


# --------------------------------------------------------------- the panels

@pytest.mark.parametrize("panel", ["S6e", "S6f", "S6h", "S6i", "S6j"])
def test_panel_draws(panel):
    fig, data = B.FIGURES[panel]()
    try:
        assert fig.get_axes()
        assert len(data) > 0
    finally:
        plt.close(fig)


def test_6h_axes_span_every_pose_not_only_the_drawn_ones():
    """The invisible non-model-0 scatter is what sets the published limits."""
    fig, _ = B.draw_6h()
    try:
        ax = fig.axes[0]
        assert ax.get_xlim()[0] < -2.0 and ax.get_xlim()[1] > 2.0
        assert ax.get_ylim()[0] < -2.0 and ax.get_ylim()[1] > 2.0
    finally:
        plt.close(fig)
