"""Supplementary panels reproduce the numbers their published versions print.

Each value below was read off the published PDF, not off the table that feeds
the panel, so these fail if a table is swapped as well as if the code drifts.
"""
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from cdkbe import dose_response as dr
from cdkbe import lfc_correlation as lc
from cdkbe.source_data import SUPP_DRAWN, SUPP_PANELS
from cdkbe.source_data import allele_specific_mtt as am
from cdkbe.source_data import cpd, editing_position as ep, leukemia_burden as lb
from cdkbe.source_data import supp_figures as sf
from cdkbe.source_data import validation_mtt as vm


# ------------------------------------------------------------------ lollipops

#: the denominator each published lollipop title prints. subpool3's is the
#: whole pool rather than the targeting count -- see docs/DISCREPANCIES.md 12.
PUBLISHED_N = {"subpool1": "2,439", "subpool2": "4,533", "subpool3": "1,655"}


@pytest.mark.parametrize("panel", sorted(sf.SUPP_SPEC))
def test_lollipop_title_matches_published(panel):
    subpool, _, _, _, title, _ = sf.SUPP_SPEC[panel]
    fig, kept = sf.FIGURES[panel]()
    text = fig._suptitle.get_text()
    plt.close(fig)
    assert title in text, f"{panel}: title does not name the drug"
    assert PUBLISHED_N[subpool] in text, f"{panel}: wrong gRNA denominator"
    assert set(kept) == {"CBE", "ABE"}


def test_subpool3_denominator_is_the_published_inconsistent_one():
    """Pinned deliberately: the published figure quotes the whole F3-R3 pool."""
    from cdkbe.source_data.screen_figures import PUBLISHED_N_LIBRARY
    assert PUBLISHED_N_LIBRARY["subpool3"] == 1655


# ------------------------------------------------------- LFC correlation (S4)

@pytest.mark.parametrize("editor", lc.EDITORS)
@pytest.mark.parametrize("subpool,n_cols,n_guides",
                         [("subpool1", 14, 2730), ("subpool2", 5, 4868),
                          ("subpool3", 8, 1655)])
def test_correlation_shape(editor, subpool, n_cols, n_guides):
    corr, wide = lc.correlation(editor, subpool)
    assert corr.shape == (n_cols, n_cols)
    assert len(wide) == n_guides
    # the matrix is a correlation matrix: unit diagonal, symmetric
    assert (corr.to_numpy().diagonal().round(10) == 1).all()
    assert (corr.round(10).to_numpy() == corr.round(10).to_numpy().T).all()


def test_subpool1_carries_the_sy5609_retest():
    corr, _ = lc.correlation("ABE", "subpool1")
    assert "DMSO vs. T0 (SY-5609 retest)" in corr.columns
    assert "DMSO vs. T0" in corr.columns


# ------------------------------------------------------------------ CPDs (S1)

@pytest.mark.parametrize("panel,n_rows,n_levels",
                         [("S1a", 432, 9), ("S1b", 180, 3),
                          ("S1c", 192, 4), ("S1d", 336, 7)])
def test_cpd_shape(panel, n_rows, n_levels):
    d = cpd.source(panel)
    hue = cpd.SPEC[panel][2]
    assert len(d) == n_rows
    assert d[hue].nunique() == n_levels
    assert set(d["Editor"]) == {"CBE", "ABE"}


def test_subpool1_drops_the_rerun_arms():
    """SY-5609 was rerun as its own screen, so its arms leave panel 1a."""
    d = cpd.source("S1a")
    assert not d["Condition"].str.startswith("SY-5609").any()
    assert "SEL120 HIGH" not in set(d["Condition"])


# -------------------------------------------------------- editing by position

def test_editing_position_table():
    d = ep.table()
    assert len(d) == 83_175
    assert set(d["Editor"]) == {"ABE", "CBE"}
    assert sorted(d["Gene"].unique()) == sorted(ep.CDK_ORDER)
    # every row is an editable base of a 20 nt protospacer
    assert d["Location"].between(1, 20).all()
    # %Reads is a percentage of that guide's own total
    assert d["%Reads"].between(0, 100).all()


def test_editing_window_is_the_library_window():
    assert ep.window() == (4, 8)


# ------------------------------------------------- leukemia burden (S7f-S7i)

@pytest.mark.parametrize("experiment,n_mice", [("CDK9", 19), ("CDK7", 20)])
def test_burden_experiment(experiment, n_mice):
    d = lb.raw(experiment)
    assert len(d) == n_mice
    mutant = lb.EXPERIMENTS[experiment][2]
    assert set(d["Genotype"]) == {"WT", mutant}


def test_normalisation_puts_each_genotype_on_its_own_vehicle():
    """Every treated mouse is a percentage of its *own* genotype's vehicle."""
    for experiment in lb.EXPERIMENTS:
        d = lb.normalised(experiment)
        assert (d[lb.NORMED] > 0).all()
        assert len(d) < len(lb.raw(experiment))     # vehicle arms are the divisor


# ----------------------------------------------- allele-specific MTT (S9d-g)

#: LC50 ordering the panels show: the resistance allele shifts its own drug
def test_resistance_alleles_shift_their_own_drug():
    d = am.plate()
    sy = dr.curve_table(d, am.CELL_TYPES, "SY-5609").set_index("Cell Line")
    ki = dr.curve_table(d, am.CELL_TYPES, "KI-CDK9d-32").set_index("Cell Line")
    # CDK7 L18F resists the CDK7 inhibitor; CDK9 L156F resists the CDK9 degrader
    assert sy.loc["L18F", "LC50_nM"] > sy.loc["NT", "LC50_nM"]
    assert ki.loc["L156F", "LC50_nM"] > ki.loc["NT", "LC50_nM"]
    # and the double mutant resists both
    assert sy.loc["L156F/L18F", "LC50_nM"] > sy.loc["NT", "LC50_nM"]
    assert ki.loc["L156F/L18F", "LC50_nM"] > ki.loc["NT", "LC50_nM"]


# ------------------------------------------- CDK12/13 validation MTT (fig 3d)

#: (drug, variant) -> (fold shift as printed, significance as printed)
PUBLISHED_3D = {
    ("BSJ-4-116", "CDK12 I733V (3)"): (4.0, "****"),
    ("BSJ-4-116", "CDK12 G731K (4)"): (2.4, "****"),
    ("BSJ-4-116", "CDK12 M816I_D817N (5)"): (5.5, "****"),
    ("BSJ-4-116", "CDK12 G822K (6)"): (0.6, "**"),
    ("BSJ-4-116", "CDK12 G822E (7)"): (0.8, "ns"),
    ("BSJ-4-116", "CDK13 I711V (8)"): (3.4, "****"),
    ("CDK12-IN-2", "CDK12 I733V (3)"): (2.5, "***"),
    ("CDK12-IN-2", "CDK12 G731K (4)"): (2.9, "**"),
    ("CDK12-IN-2", "CDK12 M816I_D817N (5)"): (3.5, "****"),
    ("CDK12-IN-2", "CDK12 G822K (6)"): (1.2, "ns"),
    ("CDK12-IN-2", "CDK12 G822E (7)"): (1.6, "ns"),
    ("CDK12-IN-2", "CDK13 I711V (8)"): (2.1, "****"),
    ("HQ461", "CDK12 I733V (3)"): (2.2, "**"),
    ("HQ461", "CDK12 G731K (4)"): (3.4, "**"),
    ("HQ461", "CDK12 M816I_D817N (5)"): (2.8, "****"),
    ("HQ461", "CDK12 G822K (6)"): (2.8, "****"),
    ("HQ461", "CDK12 G822E (7)"): (2.5, "****"),
    ("HQ461", "CDK13 I711V (8)"): (1.7, "***"),
}


def test_figure_3d_statistics_match_the_published_panel():
    st = vm.stats().set_index(["Drug", "Variant"])
    bad = []
    for key, (fold, sig) in PUBLISHED_3D.items():
        row = st.loc[key]
        if round(row["fold_shift"], 1) != fold or row["sig"] != sig:
            bad.append(f"{key}: got {row['fold_shift']:.1f}x {row['sig']}, "
                       f"published {fold}x {sig}")
    assert not bad, "\n".join(bad)


def test_each_variant_is_tested_against_its_own_editor_control():
    for variant, control in vm.CONTROL_OF.items():
        editor = "ABE" if control.endswith("ABE") else "CBE"
        assert control == f"NT {editor}"
    # the ABE variants are the two made with ABE
    abe = [v for v, c in vm.CONTROL_OF.items() if c == "NT ABE"]
    assert set(abe) == {"CDK12 I733V (3)", "CDK13 I711V (8)"}


# ---------------------------------------------------------------- the registry

def test_every_drawn_supplementary_panel_has_source_data():
    assert set(SUPP_DRAWN) <= set(SUPP_PANELS)


@pytest.mark.parametrize("panel", sorted(SUPP_PANELS))
def test_supplementary_source_data_is_non_empty(panel):
    res = SUPP_PANELS[panel]()
    frames = res.values() if isinstance(res, dict) else [res]
    for df in frames:
        assert len(df) > 0, f"{panel}: empty frame"


# ------------------------------------------------ KLIFS resistance (supp 15)

def test_oncokb_pocket_mapping_matches_the_published_footnote():
    """The published 15c reads: 30 of 104 OncoKB mutations fall outside."""
    from cdkbe.source_data import klifs_resistance as kr
    a = kr.annotated()
    assert len(a) == 104
    assert int((~a["KLIFS Catalytic Site"]).sum()) == 30


def test_oncokb_per_gene_counts_match_the_published_bars():
    """Read off observed_resistance_oncokb_KLIFS_inclusion.pdf."""
    from cdkbe.source_data import klifs_resistance as kr
    published = {"ABL1": (14, 0), "ALK": (4, 1), "BTK": (6, 1), "EGFR": (9, 6),
                 "KIT": (3, 5), "MET": (0, 2), "NTRK1": (1, 0), "NTRK3": (3, 0),
                 "PDGFRA": (0, 1), "ROS1": (2, 0)}
    pivot = kr.counts()
    got = {g: (int(pivot.loc[g, True]), int(pivot.loc[g, False])) for g in pivot.index}
    assert got == published


def test_klifs_resistance_matrices_match_the_shipped_summary():
    """58 sites over 39 positions; 13 recurrent in >=2 CDKs, 5 in >=3."""
    from cdkbe import io, paths
    from cdkbe.source_data import klifs_msa as km
    mat, assayed, n_cdk, kl = km.matrices()
    assert (mat > 0).values.sum() == 58
    assert int((n_cdk > 0).sum()) == 39
    assert int((n_cdk >= 2).sum()) == 13
    assert int((n_cdk >= 3).sum()) == 5
    # the shipped per-position summary has one row per hit position
    shipped = io.read(paths.LFC / "resistance_patterns" / "core7"
                      / "KLIFS_resistance_by_position.csv")
    assert len(shipped) == int((n_cdk > 0).sum())
    assert set(shipped["pos"]) == set(n_cdk.index[n_cdk > 0])


# ------------------------------------------------- early essentiality (16)

def test_essentiality_modules_all_import():
    """The sixteen ported scripts must resolve their paths inside the repo."""
    import importlib
    for m in ("t0_plasmid_essentiality", "essentiality_peaks", "essentiality_figures",
              "klifs_detail_figures", "correction_explained", "correction_ladder",
              "dfg_robustness", "gc_by_contrast", "gc_domain_check", "threshold_sweep",
              "assay_validation", "dbnsfp_recheck", "build_kinase_domain_msa",
              "canonical_resistance_figures", "cdk_cyclin_klifs_cartoon",
              "cdk_cyclin_klifs_ribbon"):
        importlib.import_module(f"cdkbe.essentiality.{m}")


def test_klifs_significance_tables():
    from cdkbe.source_data import essentiality_panels as ep
    wsd, wp = ep.subdomain_tests(), ep.within_pocket_tests()
    # the published panel titles: "2 clear FDR < 0.05" and "the 74 pocket columns"
    assert len(wp) == 74
    assert int((wsd.FDR < 0.05).sum()) == 2


def test_essentiality_renders_are_present():
    """16e is a montage of renders that cannot be regenerated without ChimeraX."""
    from cdkbe.essentiality import essentiality_figures as ef
    import os
    missing = [g for g in ef.GENES
               if not os.path.exists(f"{ef.CXDIR}/{g}_depleted.png")]
    assert not missing, f"missing renders: {missing}"


# ------------------------------------------------------- BEAT-seq QC (14)

# The bundle shipped outputs without plotting code, so every one of these
# panels is pinned by a number the published PDF prints on its own face.

def test_14a_cells_before_and_after_qc():
    """93 / 95 / 94 / 93 % passing, and 138,363 cells in the title."""
    from cdkbe.source_data import beatseq_qc as bq
    d = bq.source_14a().set_index("arm")
    assert [round(v) for v in d["pct_passing_qc"]] == [93, 95, 94, 93]
    assert int(d.loc["A549/ABE", "total_pre_qc"]) == 67049
    assert int(bq.overall()["total_post_qc"]) == 138363


def test_14b_guide_assignment_shares():
    """73 / 70 / 86 / 88 % singlets, over 101,805 singlets in all."""
    from cdkbe.source_data import beatseq_qc as bq
    d = bq.source_14b()
    assert [round(v) for v in d["singlet_pct"]] == [73, 70, 86, 88]
    assert int(bq.overall()["singlet"]) == 101805
    # the three shares are a partition of the post-QC cells
    total = d[["singlet_pct", "multiplet_pct", "unassigned_pct"]].sum(axis=1)
    assert (total.round(2) == 100.0).all()


def test_14c_library_coverage_is_the_sum_of_both_guide_classes():
    """91/92, 71/72, 91/92, 71/72 -- targeting plus non-targeting."""
    from cdkbe.source_data import beatseq_qc as bq
    d = bq.source_14c()
    got = [(int(r["recovered"]), int(r["designed"])) for _, r in d.iterrows()]
    assert got == [(91, 92), (71, 72), (91, 92), (71, 72)]


def test_14d_and_14e_depth_matches_the_panel_labels():
    """Both panels print median/mean and n for each arm."""
    from cdkbe.source_data import beatseq_qc as bq
    published = {"A549/ABE": (509, 500, 91), "A549/CBE": (537, 549, 71),
                 "T47D/ABE": (138, 136, 91), "T47D/CBE": (71, 70, 71)}
    d = bq.depth_stats()
    got = {a: (int(r["median"]), int(r["mean"]), int(r["n_guides"]))
           for a, r in d.iterrows()}
    assert got == published


def test_14e_bins_reproduce_the_published_bar_heights():
    """The bin count is recovered, not chosen: 28 and 30 both miss."""
    import numpy as np
    from cdkbe.source_data import beatseq_qc as bq
    published = {
        "A549/ABE": [1, 1, 1, 1, 9, 3, 13, 14, 21, 19, 7],
        "A549/CBE": [2, 4, 17, 23, 14, 8, 3],
        "T47D/ABE": [1, 1, 1, 1, 3, 3, 4, 3, 16, 23, 14, 12, 9],
        "T47D/CBE": [1, 1, 1, 2, 2, 6, 14, 24, 13, 6, 1]}
    c = bq.cells_per_guide()
    edges = bq.bin_edges()
    for arm, want in published.items():
        h, _ = np.histogram(c.loc[c["arm"] == arm, "n_cells"], bins=edges)
        assert [int(x) for x in h if x] == want, arm


def test_14f_active_rate_by_designed_consequence():
    """4/232 missense, 2/34 nonsense, 2/22 start_loss -- 288 guide-arms."""
    from cdkbe.source_data import beatseq_qc as bq
    d = bq.source_14f().set_index("variant_class")
    got = {k: (int(r["active"]), int(r["n_guide_arms"])) for k, r in d.iterrows()}
    assert got == {"missense": (4, 232), "nonsense": (2, 34),
                   "start_loss": (2, 22)}
    assert d["n_guide_arms"].sum() == 288


def test_14g_de_summary_matches_the_panel_fractions():
    """de_summary.pdf prints 2/82, 0/62, 3/82, 3/62 across the four arms."""
    from cdkbe.source_data import beatseq_qc as bq
    d = bq.source_14g().set_index("arm")
    published = {"A549/ABE": (2, 82), "A549/CBE": (0, 62),
                 "T47D/ABE": (3, 82), "T47D/CBE": (3, 62)}
    got = {a: (int(r["active"]), int(r["n_guide_arms"])) for a, r in d.iterrows()}
    assert got == published


def test_14h_rho_recomputes_to_the_delivered_values():
    """The title's +0.06, -0.03, -0.05, +0.05 come out of the guide table."""
    from cdkbe.source_data import beatseq_qc as bq
    fresh = bq.rho_cells_vs_de()
    delivered = bq.de_summary_stats().set_index("arm")["rho_cells_vs_de"]
    assert (fresh.round(3) == delivered.round(3)).all()
    assert len(bq.source_14h()) == 288


def test_14i_pairs_a_guide_with_itself_across_cell_lines():
    """144 pairs, Spearman 0.26, and the seven points the panel labels."""
    from cdkbe.source_data import beatseq_qc as bq
    w = bq.source_14i()
    assert len(w) == 144
    assert round(bq.rho_between_lines(), 2) == 0.26
    labelled = set(w.loc[w["active_either"], "display_label"])
    assert labelled == {"CDK9_L289P", "CDK4_W198*", "CDK8_R29*",
                        "CDK8_M1V_D2G", "CDK13_T66M_A67V", "CDK8_S343P",
                        "CDK6_H73R"}


def test_14j_is_the_eight_active_guide_arms_longest_first():
    from cdkbe.source_data import beatseq_qc as bq
    d = bq.source_14j()
    assert len(d) == 8
    assert list(d["num_de_genes"]) == sorted(d["num_de_genes"], reverse=True)
    top = d.iloc[0]
    assert top["display_label"] == "CDK9_L289P" and top["num_de_genes"] == 50
    assert int(top["n_cells"]) == 122 and round(top["dominant_edit_pct"]) == 66


def test_14k_alone_is_not_redrawn():
    """The gene filter behind the published clustermap is not in this repo."""
    from cdkbe.source_data import SUPP_DRAWN, SUPP_PANELS
    from cdkbe.source_data import beatseq_qc as bq
    assert set(bq.PANEL_COMPOSITE) <= set(SUPP_PANELS)
    assert set(bq.PANEL_COMPOSITE) - set(SUPP_DRAWN) == {"S14k"}


def test_the_clustermap_gene_set_really_is_unrecoverable():
    """Pinned so the claim in beatseq_qc's docstring stays honest.

    The published panel says 2,607 shared genes. The delivered pseudobulk
    matrix has 4,779 columns and no missing values, so no "shared across arms"
    rule cuts it down -- something not in the bundle chose those genes.
    """
    from cdkbe import io, paths
    m = io.read(paths.PROCESSED / "beatseq" / "outputs" / "05_primary_analysis"
                / "outputs" / "tables" / "de" / "pseudobulk_control_norm.csv")
    m = m.set_index(m.columns[0])            # the guide-arm id
    assert m.shape == (324, 4779)
    assert not m.isna().any().any()


@pytest.mark.parametrize("panel", [f"S14{k}" for k in "abcdefghij"])
def test_14_panels_draw(panel):
    import matplotlib.pyplot as plt
    from cdkbe.source_data import beatseq_qc as bq
    fig, data = bq.FIGURES[panel]()
    try:
        assert fig.get_axes()
        assert len(data) > 0
    finally:
        plt.close(fig)


@pytest.mark.parametrize("name", ["experiment_summary", "de_summary"])
def test_the_four_up_composites_redraw(name):
    """The delivered PDFs are 1x4 plates; these put the same four back."""
    import matplotlib.pyplot as plt
    from cdkbe.source_data import beatseq_qc as bq
    fig = bq.COMPOSITE_FIGURES[name]()
    fig = fig[0] if isinstance(fig, tuple) else fig
    try:
        assert len(fig.get_axes()) == 4
        assert [ax.get_title(loc="left")[0] for ax in fig.axes] == list("abcd")
    finally:
        plt.close(fig)


def test_14_panel_letters_come_from_the_composites_own_labels():
    """experiment_summary.pdf is lettered a-d; de_summary.pdf is too."""
    from cdkbe.source_data import beatseq_qc as bq
    assert [bq.PANEL_COMPOSITE[f"S14{k}"] for k in "abcd"] == [
        ("experiment_summary", c) for c in "abcd"]
    assert [bq.PANEL_COMPOSITE[f"S14{k}"] for k in "ghij"] == [
        ("de_summary", c) for c in "abcd"]


# ------------------------------------------------- dose-response grids (2a-c)

#: LC50 in nM as printed on each published DR-curve panel
PUBLISHED_LC50 = {
    "S2a": {"KI-CDK9d-32": 100, "KI-CDK9d-32N": 1384, "KB-0742": 1285,
            "SY-5609": 321, "SEL120": 4893, "Senexin B": 852},
    "S2b": {"BSJ-4-116": 147, "CDK12-IN-2": 2382, "HQ461": 1720},
}


@pytest.mark.parametrize("panel", sorted(PUBLISHED_LC50))
def test_dose_response_lc50s_match_published(panel):
    from cdkbe.source_data import screen_dose_response as sdr
    got = sdr.parameters(panel).set_index("Drug")["LC50_nM"].round().astype(int)
    assert got.to_dict() == PUBLISHED_LC50[panel]


def test_subpool3_reproduces_only_three_of_six():
    """Pinned deliberately -- see docs/DISCREPANCIES.md item 15.

    Ribociclib, Atirmociclib and INX-315 reproduce exactly; the other three
    do not, and no plate in the repository gives the published values.
    """
    from cdkbe.source_data import screen_dose_response as sdr
    published = {"Palbociclib": 6298, "Ribociclib": 2754, "Abemaciclib": 3354,
                 "Atirmociclib": 5542, "INX-315": 1687, "Tagtociclib": 3653}
    got = sdr.parameters("S2c").set_index("Drug")["LC50_nM"].round().astype(int)
    agree = {d for d, v in published.items() if abs(got[d] - v) <= 1}
    assert agree == {"Ribociclib", "Atirmociclib", "INX-315"}


def test_each_dr_panel_uses_its_own_background_convention():
    from cdkbe.source_data import screen_dose_response as sdr
    assert sdr.SPEC["S2a"]["background"] == "none"
    assert sdr.SPEC["S2b"]["background"] == "cell_line_no_cells"
    assert sdr.SPEC["S2c"]["background"] == "cell_line_no_cells"


# ------------------------------------------------------ RNA-seq volcanoes (5e)

def test_rnaseq_de_counts_match_the_panel_legends():
    """RNA_seq_comp.png prints 8/13, 3/3, 820/547, 481/234, 109/204."""
    from cdkbe.source_data import rnaseq_volcanoes as rv
    published = {"L156F_24hours_vs_L156F_no_drug": (8, 13),
                 "L156F_48hours_vs_L156F_no_drug": (3, 3),
                 "NT_24hours_vs_NT_no_drug": (820, 547),
                 "NT_48hours_vs_NT_no_drug": (481, 234),
                 "L156F_no_drug_vs_NT_no_drug": (109, 204)}
    c = rv.counts().set_index("contrast")
    got = {k: (int(r["up"]), int(r["down"])) for k, r in c.iterrows()}
    assert got == published


# --------------------------------------------- resistance variant panels (8c, 9b)

@pytest.mark.parametrize("panel,n_hits,genes", [
    ("S8c", 18, {"CDK8", "CDK19"}),
    ("S9b", 13, {"CDK7", "CDK8", "CDK9", "CDK19"}),
])
def test_resistance_variant_hits(panel, n_hits, genes):
    from cdkbe.source_data import resistance_variants as rv
    h = rv.hits(panel)
    assert len(h) == n_hits
    assert set(h["Gene"]) <= genes
    assert h["key"].is_unique


def test_9b_is_cbe_only_and_colours_by_gene():
    """The compound-mutant screen is CBE-only, so labels mark CDK7 instead."""
    from cdkbe.source_data import resistance_variants as rv
    assert set(rv.hits("S9b")["Editor"]) == {"CBE"}
    assert rv.SPEC["S9b"]["color_by"] == "Gene"
    assert rv.SPEC["S9b"]["label_colors"] == {"CDK7": "tab:red"}
    assert rv.SPEC["S8c"]["color_by"] == "Editor"


# --------------------------------------------------- the two SAV conventions

@pytest.mark.parametrize("sample", ["CBE_subpool1", "ABE_subpool1"])
def test_unfiltered_savs_reproduce_the_shipped_tables(sample):
    from cdkbe import editing, io, paths
    got = editing.savs_from_full_canonical(sample)
    ship = io.read(paths.EDITING / "SAVs_canonical" / f"{sample}_canonical.csv")
    k = ["gRNA_id", "HGVSp"]
    a = got.sort_values(k).reset_index(drop=True)
    b = ship.sort_values(k).reset_index(drop=True)
    assert a[k].equals(b[k])
    assert (a["%Reads"] - b["%Reads"]).abs().max() < 1e-9


def test_the_two_sav_conventions_differ():
    """Filtering to the editor's own edits is not a no-op -- 9b needs it."""
    from cdkbe import editing
    plain = editing.savs_from_full_canonical("CBE_subpool1")
    filtered = editing.savs_from_full_canonical("CBE_subpool1", own_edits_only=True)
    assert len(plain) == 16026
    assert len(filtered) < len(plain)
