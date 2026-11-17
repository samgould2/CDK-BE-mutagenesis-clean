"""Figure 5: that each panel is the one that was published, and still reproduces.

The mapping in `cdkbe.source_data.figure5.SPEC` is the fragile part. The two
bodies of code that drew figure 5 letter their panels differently from the
manuscript -- `fig_panels.panel_f` is not figure 5f -- and several stems differ
from their published sibling only in a tick locator or an axis label. Each test
below re-runs the check that told the published panel apart from its near-twins,
so a wrong stem fails here rather than quietly drawing the wrong picture.

WHAT THE REFERENCES ARE. `tests/reference/figure5/` holds the evidence, extracted
once from files that are not part of this repository and vendored so the checks
keep working without them:

    published_figure_text.txt   every text token in `figure5_v2.ai`, the
                                assembled manuscript figure
    panels/<stem>.txt           the same for each source panel PDF, as the
                                BEAT-seq repository and the RNA-seq analysis
                                wrote them
    panel_geometry.csv          each source panel's page size, which is what
                                separates the two stacked severity variants

Comparing a redrawn panel against `panels/<stem>.txt` is deliberately NOT
circular: those extracts came from the original files, not from this repository's
output, so they fail if a port drifts.
"""
from __future__ import annotations

import collections
import csv
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from cdkbe import paths  # noqa: E402
from cdkbe.source_data import figure5, rnaseq_t47d_figures as rna  # noqa: E402

REF = Path(__file__).resolve().parent / "reference" / "figure5"
PUBLISHED = REF / "published_figure_text.txt"

needs_pdftotext = pytest.mark.skipif(
    shutil.which("pdftotext") is None, reason="poppler's pdftotext not installed")


def published_text() -> str:
    return " ".join(PUBLISHED.read_text().split())


def published_numbers() -> set:
    return {t for t in PUBLISHED.read_text().split() if t.isdigit()}


def tokens_of(pdf: Path) -> collections.Counter:
    out = subprocess.run(["pdftotext", "-q", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    return collections.Counter(out.split())


def reference_tokens(stem: str) -> collections.Counter:
    return collections.Counter((REF / "panels" / f"{stem}.txt").read_text().split())


def geometry() -> dict:
    with open(REF / "panel_geometry.csv", newline="") as fh:
        return {r["stem"]: r for r in csv.DictReader(fh)}


# --------------------------------------------------------------- the mapping

def test_every_published_panel_is_specified():
    assert figure5.PANELS == [f"5{c}" for c in "abcdefghijk"]


def test_spec_parts_are_unique_within_a_panel():
    for panel in figure5.PANELS:
        p = figure5.parts(panel)
        assert len(p) == len(set(p)), panel


def test_every_part_carries_its_evidence():
    for panel, spec in figure5.SPEC.items():
        for slug, src, why in spec:
            assert src[0] in ("beatseq", "rnaseq"), (panel, slug)
            assert len(why) > 30, f"{panel}/{slug} has no evidence recorded"


def test_every_beatseq_stem_has_a_vendored_reference():
    for panel, spec in figure5.SPEC.items():
        for slug, src, _ in spec:
            if src[0] == "beatseq":
                assert (REF / "panels" / f"{src[2]}.txt").exists(), f"{panel}/{slug}"


def test_5k_axis_label_rules_out_the_unsuffixed_variants():
    """5k is e1n/e2n, not e1/e2: those differ only in their axis labels."""
    t = published_text()
    assert "max drug log" in t
    assert "maximum drug log" not in t
    assert "in the parent screen" not in t


def test_5j_is_stacked_not_a_single_cell_line_variant():
    """Every non-stacked severity panel carries a y-label the figure does not."""
    assert "SD of the size-matched null" not in published_text()


def test_5j_is_the_wide_stacked_variant_by_aspect_ratio():
    """d12 and d12w plot the same points; only their geometry differs.

    The published plot area is about 1.73 wide to 1 tall, which is d12w. The
    assertion is on the gap between the two candidates, not on a tolerance
    plucked for one of them.
    """
    g = geometry()
    wide = float(g["PANEL_d12w_severity_by_CDK_stacked"]["ratio"])
    tall = float(g["PANEL_d12_severity_by_CDK_stacked"]["ratio"])
    assert 1.70 < wide < 1.78, wide
    assert tall < 1.3, tall
    assert figure5.SPEC["5j"][0][1][2] == "PANEL_d12w_severity_by_CDK_stacked"


def test_5f_to_5i_are_the_mid_layout():
    """The three 4-panel layouts differ only in tick locators; one set matches."""
    nums = published_numbers()
    assert {"80", "160", "240"} <= nums, "the mid layout's volcano ticks are missing"
    assert not {"60", "120", "180"} & nums, "square-layout ticks are present"
    assert rna.PUBLISHED_LAYOUT == "pub_fig_square_4panel_mid"


def test_5a_donut_counts_are_in_the_published_figure():
    nums = published_numbers()
    for n in ("72", "92", "17", "16", "29", "11", "61", "10"):
        assert n in nums, f"{n} missing from the published figure"


def test_5d_is_the_combined_single_axes_variant():
    """Both cell lines' bracket p-values sit on one axes in the published panel."""
    t = published_text()
    assert "0.018" in t and "0.035" in t
    ref = reference_tokens("PANEL_c1s_editing_by_variant_class_combined")
    assert "0.018" in ref and "0.035" in ref


# --------------------------------------------------------------- redrawing

@pytest.mark.parametrize("panel", [f"5{c}" for c in "abcdefghijk"])
def test_panel_draws(panel):
    result = figure5.draw(panel)
    figs = result if isinstance(result, dict) else {figure5.parts(panel)[0]: result}
    assert sorted(figs) == sorted(figure5.parts(panel))
    for f in figs.values():
        assert f.get_axes(), f"{panel} drew an empty figure"
        plt.close(f)


@needs_pdftotext
@pytest.mark.parametrize("panel", ["5a", "5b", "5c", "5d", "5e", "5j", "5k"])
def test_beatseq_parts_match_the_published_panel_files(panel, tmp_path):
    """A redrawn part carries exactly the text its published panel file did."""
    result = figure5.draw(panel)
    figs = result if isinstance(result, dict) else {figure5.parts(panel)[0]: result}
    for slug, fig in figs.items():
        stem = [s for s2, s, _ in figure5.SPEC[panel] if s2 == slug][0][2]
        out = tmp_path / f"{panel}_{slug}.pdf"
        figure5.save(fig, out, panel, slug)
        plt.close(fig)
        assert tokens_of(out) == reference_tokens(stem), \
            f"{panel}/{slug} differs from the published {stem}"


# ------------------------------------------------------- 5f-5i, the RNA-seq

def test_rnaseq_counts_match_the_published_panels():
    fig, (u1, d1, n1, u2, d2, n2) = rna.published()
    plt.close(fig)
    assert (u1, d1) == (343, 477)        # W198* vs CBE NT
    assert (u2, d2) == (1, 9)            # W198* vs D196N
    assert n1 == n2 == 11489             # one gene universe, both contrasts


def test_e2f_gene_count_is_the_published_198():
    assert len(rna._data().e2f_ids) == 198


def test_the_three_layouts_differ_only_in_geometry():
    """If a layout changed the statistics, the mid-layout claim would be unsafe."""
    counts = set()
    for name, *_ in rna.LAYOUTS:
        fig, c = rna.published(name)
        plt.close(fig)
        counts.add(c)
    assert len(counts) == 1, "the layouts disagree on the DEG counts"


@needs_pdftotext
def test_published_layout_reproduces_its_shipped_pdf(tmp_path):
    fig, _ = rna.published()
    out = tmp_path / "mid.pdf"
    rna.save(fig, out)
    plt.close(fig)
    assert tokens_of(out) == reference_tokens("pub_fig_square_4panel_mid")


@needs_pdftotext
def test_a_single_rnaseq_panel_is_the_one_from_the_grid(tmp_path):
    """Each 5f-5i file must carry only its own text, and all of it."""
    grid = reference_tokens("pub_fig_square_4panel_mid")
    seen = collections.Counter()
    for key in ("5f", "5g", "5h", "5i"):
        fig = rna.panel(key)
        out = tmp_path / f"{key}.pdf"
        rna.save(fig, out)
        plt.close(fig)
        t = tokens_of(out)
        assert t, key
        seen += t
    # the four panels together carry the grid's text, less the a/b/c/d letters
    # the grid adds and the per-panel figures do not
    extra = sorted((grid - seen).elements())
    assert extra == ["a", "b", "c", "d"], extra


# ------------------------------------- the manuscript's own panel-source record

#: `docs/source_data_mapping_main.csv` is the mapping kept alongside the
#: manuscript, written independently of anything in this repository. It agrees
#: with `figure5.SPEC` on ten of the eleven panels. The eleventh is a slip in
#: that file, not here -- see docs/DISCREPANCIES.md item 18.
MAPPING_EXCEPTION = {"5k": "PANEL_e1_severity_vs_drug_lfc_A549"}


def mapping_stems() -> dict:
    import re

    import pandas as pd

    x = pd.read_csv(paths.ROOT / "docs" / "source_data_mapping_main.csv")
    x = x.dropna(subset=["figure"])
    x["n"] = x["figure"].astype(str).str.extract(r"(\d+)")
    x = x.dropna(subset=["n"])
    out = {}
    for _, r in x[x["n"] == "5"].iterrows():
        key = "5" + str(r["panel"])
        out[key] = sorted({re.sub(r"\.(pdf|svg|png)$", "", s.split("/")[-1])
                           for s in str(r["source_file"]).split()
                           if "PANEL_" in s or "pub_fig" in s})
    return out


def spec_stems(panel: str) -> list:
    return sorted({src[2] if src[0] == "beatseq" else "pub_fig_square_4panel_mid"
                   for _, src, _ in figure5.SPEC[panel]})


@pytest.mark.parametrize("panel", [f"5{c}" for c in "abcdefghijk"])
def test_spec_agrees_with_the_manuscript_mapping(panel):
    """An independent record of which source file became which panel."""
    want = mapping_stems().get(panel)
    assert want, f"{panel} is not in the mapping file"
    got = spec_stems(panel)
    if panel in MAPPING_EXCEPTION:
        extra = set(want) - set(got)
        assert extra == {MAPPING_EXCEPTION[panel]}, (
            f"{panel}: the disagreement with the mapping file changed; it was "
            f"{MAPPING_EXCEPTION[panel]}, now {sorted(extra)}")
        return
    assert got == want


def test_the_5k_disagreement_is_the_mapping_files_slip():
    """`PANEL_e1`'s axis labels are absent from the published figure; e1n's are.

    If this ever fails, the mapping file was right and `figure5.SPEC` is wrong.
    """
    t = published_text()
    assert "max drug log" in t                      # e1n / e2n
    assert "in the parent screen" not in t          # e1 / e2
    assert "SD of the size-matched null" not in t   # e1 / e2
    assert spec_stems("5k") == ["PANEL_e1n_severity_vs_drug_lfc_A549",
                                "PANEL_e2n_severity_vs_drug_lfc_T47D"]


# ------------------------------------------------- the statistic behind 5j/5k

def test_severity_recomputes_from_the_shipped_tables():
    """5j and 5k plot `severity`; it is derived here, not only read from a CSV."""
    import pandas as pd
    from cdkbe import io as _io
    from cdkbe.beatseq import paths as bp, severity as sev

    m = pd.read_csv(_io.resolve(bp.RESULTS / "master_guide_table.csv"))
    null = pd.read_csv(_io.resolve(bp.RESULTS / "pseudo_guide_null.csv"))
    ship = pd.read_csv(_io.resolve(bp.RESULTS / "severity.csv"))
    out = sev.add_severity(m.copy(), null)

    key = ["arm", "guide_identity"]
    j = ship[key + ["severity", "q"]].merge(out[key + ["severity", "q"]], on=key,
                                            suffixes=("_ship", "_new"))
    assert len(j) == len(ship)
    assert (j.severity_ship - j.severity_new).abs().max() < 1e-9
    assert (j.q_ship - j.q_new).abs().max() < 1e-12
