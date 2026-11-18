"""Figure 4m and supplementary 13e-13g: the T47D CDK variant ± palbociclib RNA-seq.

One experiment, four panels, and until the expression matrix arrived this
repository had three of them only as published artwork -- `docs/DISCREPANCIES.md`
item 32 is the note that said so. These tests are what closes it, so they are
written against the published SVGs themselves rather than against a stored
redraw: the panels are vendored at `data/raw/collaborator_panels/T47D_MCF7/`,
and nothing is copied into `reference/` because the published file is already
in the repository and one copy is enough.

Three claims, in order of how much they are worth:

1. **Every character of all three panels matches**, including the DE counts the
   volcano legends print and the percent-variance in the PCA axis labels --
   numbers that are computed from the shipped tables, not styling.
2. **The published heatmap's own fills decode back to the shipped z-scores.**
   Item 32 said 4m could not be read back out of its figure because the values
   are in the cell fills. That is still true as a way of *recovering* numbers:
   what the inversion gives is quantised to the colormap's 256 steps. As a
   *check* it is conclusive, and it is the direct evidence that the matrix in
   `data/processed/rnaseq/T47D_CDK_palbo_RNAseq/` is the one the published
   panel was drawn from.
3. **The published PCA's marker positions equal the computed scores exactly**,
   which also pins the sign convention in `PC_ORIENTATION`: with PC2 the other
   way up every point would still land on the axes, and the panel would be a
   mirror image of the published one.
"""
from __future__ import annotations

import base64
import io
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colormaps  # noqa: E402

from cdkbe import paths, svg_values as sv  # noqa: E402
from cdkbe.source_data import rnaseq_t47d_palbo as rp  # noqa: E402

PUBLISHED = paths.RAW / "collaborator_panels" / "T47D_MCF7"
#: panel -> the published file, which is named for the figure it drew
FILES = {"4m": "marker_gene_panel_DMSO_palbo_stacked.svg",
         "S13e": "PCA.svg",
         "S13g": "volcano_variant_vs_control.svg"}
DRAWERS = {"4m": rp.panel_4m, "S13e": rp.panel_S13e, "S13g": rp.panel_S13g}


def published(panel: str) -> str:
    return (PUBLISHED / FILES[panel]).read_text(errors="ignore")


def redrawn(panel: str) -> str:
    """This repository's panel, as SVG.

    Start from matplotlib's own defaults rather than whatever style an earlier
    test left behind, and set the two rcParams the delivered notebooks set that
    change the file: live text, and the 110 dpi they drew at.
    """
    rc = {k: v for k, v in matplotlib.rcParamsDefault.items()
          if k not in {"backend", "interactive"}}
    rc["svg.fonttype"] = "none"
    rc["figure.dpi"] = 110
    with matplotlib.rc_context(rc):
        fig, _ = DRAWERS[panel]()
        buf = io.StringIO()
        fig.savefig(buf, format="svg", bbox_inches="tight")
        plt.close(fig)
    return buf.getvalue()


def texts(svg: str) -> list:
    """Every string the file spells out, in document order.

    `itertext` rather than a regex on `<text>`: matplotlib puts mathtext --
    `log$_2$ Fold Change` -- in nested `<tspan>`s, and a regex that stops at
    the first `<` returns an empty string for every one of them, which would
    make eight panel axis labels compare equal by being equally blank.

    Where a label *is* nested, the indentation between the tspans is dropped:
    matplotlib 3.10 writes `log Fold Change` as one tspan and 3.9 as one tspan
    per character, with a newline and seven spaces between each. Same label,
    and the spaces mathtext itself sets are U+00A0, which survives this.
    """
    root = ET.fromstring(svg)
    out = []
    for t in root.iter(sv.SVG + "text"):
        if len(t):
            out.append(re.sub(r"[ \t\r\n]+", "", "".join(t.itertext())))
        else:
            out.append(t.text or "")
    return out


def heatmap_images(svg: str) -> list:
    """The PNG payloads matplotlib embedded for each `imshow`, decoded."""
    out = []
    for b64 in re.findall(r'xlink:href="data:image/png;base64,([^"]+)"', svg):
        out.append(plt.imread(io.BytesIO(base64.b64decode(b64)), format="png"))
    return out


def marker_positions(svg: str) -> np.ndarray:
    """Every scatter marker in the first axes, in data units.

    `svg_values.scale` inverts the axis from its own tick labels; the markers
    are `<use>` references into a glyph definition, each carrying its centre.
    """
    root = sv.load(io.StringIO(svg))
    fx, fy = sv.scale(root, "x"), sv.scale(root, "y")
    axes = root.find(f".//{sv.SVG}g[@id='axes_1']")
    pts = [(fx(float(u.get("x", 0))), fy(float(u.get("y", 0))))
           for g in axes.iter(sv.SVG + "g")
           if (g.get("id") or "").startswith("PathCollection")
           for u in g.iter(sv.SVG + "use")]
    return np.array(sorted(pts))


# ------------------------------------------------- the experiment's shape --

def test_the_matrix_is_six_lines_in_two_arms_in_triplicate():
    log_expr, metadata, _ = rp.expression()
    assert log_expr.shape[1] == 36
    assert list(dict.fromkeys(metadata["genotype"])) == rp.GENOTYPES
    assert set(metadata["treatment"]) == {"veh", "palbo"}
    assert metadata.groupby(["genotype", "treatment"]).size().unique().tolist() == [3]


def test_it_is_not_the_figure_5_t47d_experiment():
    """The two T47D RNA-seq runs are different experiments (item 32).

    They share three sample *names* and not the libraries behind them, which
    is the trap this pins: `T47D_CBE_NT_veh_1..3` exists in both matrices.
    Within either experiment a pair of CBE NT replicates correlates at 0.994
    or better on log2 TPM; the same name across the two correlates at 0.93,
    which is two runs of one cell line rather than one run read twice.
    """
    from cdkbe.source_data import rnaseq_t47d as fig5

    log_expr, metadata, _ = rp.expression()
    fig5_tpm, fig5_log, fig5_meta, _ = fig5._expression()
    assert len(fig5_meta) == 12 and all("_veh_" in s for s in fig5_meta.index)
    assert set(metadata["genotype"]) & set(fig5_meta) == {"CBE_NT"}

    shared = sorted(set(log_expr.columns) & set(fig5_log.columns))
    assert shared == [f"T47D_CBE_NT_veh_{i}" for i in (1, 2, 3)]
    genes = log_expr.index.intersection(fig5_log.index)
    across = min(np.corrcoef(log_expr.loc[genes, s],
                             fig5_log.loc[genes, s])[0, 1] for s in shared)
    within = np.corrcoef(log_expr.loc[genes, shared[0]],
                         log_expr.loc[genes, shared[1]])[0, 1]
    assert across < 0.95 < within


def test_every_contrast_is_a_variant_against_its_own_arm():
    for variant, control in rp.CONTROL_FOR.items():
        assert control in rp.CONTROLS
        assert variant.startswith("CDK")
        for arm in rp.ARMS:
            # of the 19,664 rows DESeq2 returns, the ones it could call
            assert len(rp.deseq2(variant, arm)) > 10_000


# --------------------------------------------------- against the artwork ---

@pytest.mark.parametrize("panel", list(FILES))
def test_every_character_matches_the_published_panel(panel):
    assert texts(redrawn(panel)) == texts(published(panel))


def test_the_published_volcano_prints_the_counts_this_repository_computes():
    """The legends are the only place the DE counts appear, so pin them here."""
    printed = set(texts(published("S13g")))
    counts = rp.volcano_counts()
    assert len(counts) == 8
    for _, row in counts.iterrows():
        assert f"{row['up']} DE UP" in printed
        assert f"{row['down']} DE DOWN" in printed


def test_the_published_pca_axis_labels_are_the_computed_variance():
    _, var_pct = rp.pca()
    printed = texts(published("S13e"))
    assert f"PC1: {var_pct[0]:.0f}% variance" in printed
    assert f"PC2: {var_pct[1]:.0f}% variance" in printed


def test_the_published_pca_markers_sit_where_the_scores_say():
    """An exact inversion: both sides agree to the SVG's own precision."""
    pca_df, _ = rp.pca()
    computed = np.array(sorted(zip(pca_df["PC1"], pca_df["PC2"])))
    found = marker_positions(published("S13e"))
    assert len(found) == 36
    assert np.abs(found - computed).max() < 1e-3


def test_flipping_pc2_would_not_be_the_published_panel():
    """Why `PC_ORIENTATION` is pinned rather than left to the SVD's own sign."""
    pca_df, _ = rp.pca()
    mirrored = np.array(sorted(zip(pca_df["PC1"], -pca_df["PC2"])))
    found = marker_positions(published("S13e"))
    assert np.abs(found - mirrored).max() > 1.0


def test_the_published_heatmap_fills_decode_to_the_shipped_zscores():
    """Item 32's heatmap, read back out of its own fills.

    Quantised to the colormap's 256 steps and sampled away from the cell
    edges, so the tolerance is a couple of steps -- not a recovery of the
    numbers, but proof that these are the numbers that were drawn.
    """
    arm_z = rp.marker_zscores()
    vmax = np.ceil(max(np.nanmax(np.abs(z.values)) for z in arm_z.values()) * 2) / 2
    step = 2 * vmax / 255
    lut = colormaps["RdBu_r"](np.linspace(0, 1, 256))[:, :3]
    images = heatmap_images(published("4m"))
    assert len(images) == 3                      # two heatmaps and the colour bar

    for image, arm in zip(images, rp.ARMS):
        # matplotlib writes an `imshow` raster bottom-up and flips it with the
        # image's own transform, so the first gene row is the last image row
        image = np.flipud(image[:, :, :3])
        z = arm_z[arm]
        rows, cols = z.shape
        h, w = image.shape[0] / rows, image.shape[1] / cols
        recovered = np.empty((rows, cols))
        for i in range(rows):
            for j in range(cols):
                block = image[int(i * h + h / 4):int(i * h + 3 * h / 4),
                              int(j * w + w / 4):int(j * w + 3 * w / 4)]
                px = np.median(block.reshape(-1, 3), axis=0)
                nearest = int(np.argmin(((lut - px) ** 2).sum(axis=1)))
                recovered[i, j] = nearest / 255 * 2 * vmax - vmax
        assert np.abs(recovered - z.values).max() < 3 * step


# ------------------------------------------- against the delivered tables --

@pytest.mark.parametrize("arm,name", [("veh", "DMSO"), ("palbo", "Palbociclib")])
def test_the_zscores_match_the_tables_the_notebook_wrote(arm, name):
    import pandas as pd

    delivered = pd.read_csv(rp.BASE / "tables" / f"marker_gene_zscore_{name}.csv",
                            index_col=0)
    computed = rp.marker_zscores()[arm]
    assert list(delivered.index) == list(computed.index)
    assert list(delivered.columns) == list(computed.columns)
    assert np.abs(delivered.values - computed.values).max() < 1e-9


def test_the_marker_panel_is_seventeen_genes_in_the_published_order():
    marker, _ = rp.marker_expression()
    assert list(marker.index) == rp.MARKER_GENES
    assert marker.shape == (17, 36)


# ------------------------------------------------------ S13f, source only --

def test_the_gsea_inputs_ship_even_though_the_dot_plot_does_not():
    """S13f is **source only**: its ranks and gene sets, but no redraw."""
    from cdkbe.source_data import PANELS, SUPP_DRAWN, SUPP_PANELS

    assert "S13f" in SUPP_PANELS and "S13f" not in SUPP_DRAWN
    assert "4m" in PANELS
    ranks = rp.gsea_ranks()
    assert set(ranks["variant"]) == set(rp.CONTROL_FOR)
    assert ranks["stat"].notna().all()
    sets = rp.gsea_sets()
    assert list(dict.fromkeys(sets["gene_set"])) == rp.GSEA_SETS


def test_the_five_published_gene_sets_are_the_ones_the_dot_plot_names():
    """The published dot plot spells its five rows out; they are in the cache."""
    cache = rp.hallmark()
    for name in rp.GSEA_SETS:
        assert len(cache[name]) > 100
    assert rp.PROLIF_SETS == rp.GSEA_SETS[:2]
