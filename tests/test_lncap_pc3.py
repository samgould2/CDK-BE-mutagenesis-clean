"""Supplementary 11c, 11D-11H and 11K.

These seven panels arrived as packages: the values, the script that drew them,
and the published image. The module redraws them from the values, so what is
worth testing is that the redraw still says what the delivered package says --
the DE counts its README prints, the statistics its source-data file carries,
and the gene order its heatmap was built in.

Nothing here re-runs DESeq2 or GSEA. Those are the delivered authors' own fits;
re-running them would need the count matrices, which this repository does not
ship, and would be a different number anyway if a package version moved.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cdkbe.source_data import lncap_pc3 as lp  # noqa: E402

PANELS = ["S11c", "S11D", "S11E", "S11F", "S11G", "S11H", "S11K"]

#: the DE calls the delivered READMEs print, per contrast: (up, down)
DELIVERED_DE = {
    "NT_HQ461_vs_NT_Control": (187, 593),
    "I733V_HQ461_vs_I733V_Control": (12, 27),
    "NT_HQ461_vs_I733V_HQ461": (122, 164),
    "NT_Control_vs_I733V_Control": (0, 1),
    "NT_vs_I733V_Control": (141, 52),
    "NT_vs_G822E_Control": (764, 429),
    "NT_vs_G731K_Control": (696, 363),
}


@pytest.mark.parametrize("panel", PANELS)
def test_each_panel_draws(panel):
    fig, source = lp.FIGURES[panel]()
    assert fig.axes
    assert len(source) > 0
    plt.close(fig)


@pytest.mark.parametrize("panel", PANELS)
def test_each_panel_reads_only_from_this_repository(panel):
    folder = lp.FOLDER[panel]
    assert folder.is_dir(), folder
    assert Path(*folder.parts[:len(folder.parts) - 2]).name == "collaborator_panels" \
        or "collaborator_panels" in folder.parts


@pytest.mark.parametrize("panel", ["S11E", "S11K"])
def test_the_volcanoes_reproduce_the_de_counts_the_readme_prints(panel):
    for _, row in lp.volcano_counts(panel).iterrows():
        assert (row["up"], row["down"]) == DELIVERED_DE[row["comparison"]], \
            row["comparison"]


def test_the_lncap_volcano_panel_draws_all_four_contrasts_in_order():
    stems = [s for s, _ in lp.VOLCANOES["S11E"]]
    assert stems == ["NT_HQ461_vs_NT_Control",
                     "I733V_HQ461_vs_I733V_Control",
                     "NT_HQ461_vs_I733V_HQ461",
                     "NT_Control_vs_I733V_Control"]


def test_the_qpcr_summary_matches_the_replicates_it_came_from():
    import numpy as np

    rep = lp.qpcr_replicates()
    summ = lp.qpcr_summary().set_index(["Gene", "Timepoint_hr"])
    for (gene, tp), row in summ.iterrows():
        x = rep[(rep["Gene"] == gene)
                & (rep["Timepoint_hr"] == tp)]["FoldChange"].to_numpy(float)
        assert len(x) == row["n"]
        assert x.mean() == pytest.approx(row["Mean_FoldChange"], abs=5e-4)
        assert x.std(ddof=1) / np.sqrt(len(x)) == pytest.approx(row["SEM"],
                                                                abs=5e-4)


def test_ercc4_at_zero_hours_has_only_two_replicates():
    """The delivered file says so in its header; the panel's SEM depends on it."""
    rep = lp.qpcr_replicates()
    n = len(rep[(rep["Gene"] == "ERCC4") & (rep["Timepoint_hr"] == 0)])
    assert n == 2
    summ = lp.qpcr_summary().set_index(["Gene", "Timepoint_hr"])
    assert summ.loc[("ERCC4", 0), "n"] == 2
    assert all(summ.loc[(g, t), "n"] == 3 for g in lp.QPCR_GENES
               for t in lp.QPCR_TIMEPOINTS if (g, t) != ("ERCC4", 0))


def test_the_qpcr_stars_are_bonferroni_corrected_across_three_timepoints():
    summ = lp.qpcr_summary()
    treated = summ[summ["Timepoint_hr"] != 0]
    assert treated["p_Bonferroni_vs_0hr"].tolist() == pytest.approx(
        (treated["p_raw_vs_0hr"] * 3).tolist(), rel=2e-3)


def test_the_pca_is_twelve_samples_in_four_conditions():
    d = lp.pca_values()
    assert len(d) == 12
    assert d["condition"].value_counts().tolist() == [3, 3, 3, 3]
    assert set(d["condition"]) == set(lp.PCA_COLOURS)


def test_the_pca_axis_labels_carry_the_variance_explained():
    d = lp.pca_values()
    pcs = [c for c in d.columns if c.startswith("PC")]
    assert pcs[0].startswith("PC1: 31%")
    assert pcs[1].startswith("PC2: 15%")


def test_the_androgen_heatmap_is_ordered_by_log2fc_most_decreased_first():
    d = lp.androgen_table()
    lfc = d["log2FoldChange"].astype(float).tolist()
    assert lfc == sorted(lfc)
    assert d["gene"].iloc[0] == "MAF"
    assert len(d) == 23


def test_the_androgen_heatmap_has_six_nt_samples_in_control_then_hq461_order():
    cols = [c for c in lp.androgen_table().columns if c.endswith("_zscore")]
    assert len(cols) == 6
    assert [c.startswith("NT_Control") for c in cols] == [1, 1, 1, 0, 0, 0]


def test_the_gsea_running_score_starts_and_ends_at_zero():
    d = lp.gsea_table()
    es = d["running_ES"].astype(float)
    assert abs(es.iloc[-1]) < 1e-3
    assert lp.gsea_statistics()["NES"] == pytest.approx(-1.303, abs=1e-3)
    assert lp.gsea_statistics()["nominal p"] == pytest.approx(0.0446, abs=1e-4)


def test_the_gsea_hits_are_the_genes_in_the_set():
    d = lp.gsea_table()
    assert d["in_gene_set"].sum() == 151   # the header says 151 of 164


def test_the_logo_plots_cover_both_amplicons_and_their_own_controls():
    d = lp.logo_table()
    assert set(d["sample"]) == set(lp.LOGO_ORDER)
    amp1 = sorted(d[d["sample"] == "G731K gRNA"]["residue"])
    amp2 = sorted(d[d["sample"] == "G822E/K gRNA"]["residue"])
    assert amp1 == [731, 732, 733]
    assert amp2 == [816, 817, 818, 819, 820, 821, 822]
    assert set(lp.LOGO_WT) == set(amp1) | set(amp2)


def test_every_logo_residue_sums_to_a_hundred_percent():
    d = lp.logo_table().set_index(["sample", "residue"])
    assert d.sum(axis=1).round(6).eq(100.0).all()


def test_the_non_targeting_logo_controls_are_essentially_all_wild_type():
    d = lp.logo_table().set_index(["sample", "residue"])
    for sample in ("NT control (amp1)", "NT control (amp2)"):
        for residue, row in d.loc[sample].iterrows():
            assert row[lp.LOGO_WT[residue]] > 95, (sample, residue)


def test_the_g822e_stop_codon_artifact_is_left_in_its_own_bucket():
    """The authors chose not to fold it back into WT; their README says why."""
    d = lp.logo_table().set_index(["sample", "residue"])
    stop = d.loc[("G822E/K gRNA", 822), "*"]
    assert 4 < stop < 7
    assert d.loc[("NT control (amp2)", 822), "*"] < 0.5


def test_every_panel_is_registered_once():
    assert set(lp.PANELS) == set(lp.FIGURES) == set(lp._STEM) == set(lp.FOLDER)
    assert all(k.startswith("S11") for k in lp.FIGURES)
