"""The ChimeraX leader-line panels: figure 6 and supplementary 15e.

These are the two panels where labels are parked in the frame's gutters and
joined to their residue by a leader line. `scripts/chimerax/label_consensus.py`
makes them, and it needs the ChimeraX application, so these tests check the
files a render produced rather than producing them. They skip when it has not
been run.

The published SVGs are vendored beside this file, so the comparison is against
the manuscript's own artwork and not against a previous run of this code.

Two defects are pinned here, because both were silent:

* the SVG took the window's shape from nowhere and always wrote a square
  viewBox, which stretched every vertical gap by 1.76 and distorted the
  structure underneath it;
* `full` mode -- the multi-line variant blocks supplementary 15e is made of --
  was missing from `main()` although `variant_block`, MAX_W_FULL and
  BLOCK_GAP_PER_PT were all still in the file waiting for it.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = Path(__file__).resolve().parent / "reference" / "leader_panels"
RENDERS = ROOT / "figures" / "chimerax"

PANELS = {
    "S15e": (RENDERS / "cxc_core7_peaks" / "peaks_leader_variants.svg",
             REFERENCE / "published_peaks_leader_variants.svg"),
    "6": (RENDERS / "cxc_core7_consensus_strong" / "consensus_leader_variants.svg",
          REFERENCE / "published_consensus_leader_variants.svg"),
}

TEXT = re.compile(r'<text x="([\d.]+)" y="([\d.]+)"[^>]*>([^<]*)</text>')
PATH = re.compile(r'<path d="M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)" '
                  r'stroke="(#[0-9a-f]{6})"')


def labels(path: Path) -> dict:
    """{label text: (x, y)} for one leader panel."""
    return {m.group(3): (float(m.group(1)), float(m.group(2)))
            for m in TEXT.finditer(path.read_text())}


def rendered(key):
    made, published = PANELS[key]
    if not made.exists():
        pytest.skip(f"{made.name} not rendered; run scripts/chimerax/render.py")
    return made, published


@pytest.mark.parametrize("key", sorted(PANELS))
def test_the_viewbox_keeps_the_windows_shape(key):
    """A square viewBox is the bug: it stretches every gap by 1400/796."""
    made, published = rendered(key)
    want = re.search(r'viewBox="([^"]+)"', published.read_text()).group(1)
    got = re.search(r'viewBox="([^"]+)"', made.read_text()).group(1)
    assert got == want == "0 0 1400 796"


def test_15e_is_the_published_panel_label_for_label():
    """Every block, every variant list, every coordinate within half a pixel."""
    made, published = rendered("S15e")
    mine, theirs = labels(made), labels(published)
    assert set(mine) == set(theirs)
    worst = max(max(abs(a - b) for a, b in zip(mine[k], theirs[k]))
                for k in theirs)
    assert worst <= 0.5, f"label positions drift up to {worst:.1f}px"


def test_15e_leader_lines_match_in_position_and_colour():
    made, published = rendered("S15e")
    mine = PATH.findall(made.read_text())
    theirs = PATH.findall(published.read_text())
    assert len(mine) == len(theirs) == 7            # seven detected peaks
    for a, b in zip(mine, theirs):
        assert a[4] == b[4]                         # the ramp colour
        for x, y in zip(a[:4], b[:4]):
            assert abs(float(x) - float(y)) <= 0.5


def test_15e_blocks_are_multi_line_and_headed_by_the_peak():
    """`full` mode: a header line, then one line per contributing gene."""
    made, _ = rendered("S15e")
    text = set(labels(made))
    heads = sorted(t for t in text if t.startswith("P"))
    assert len(heads) == 7
    assert "P4  VII:78  CDK2-135  7 variants / 4 CDKs" in heads
    assert "CDK9 L156F, I157V, T158A, R159C" in text


def test_figure_6_carries_every_published_label_and_one_more():
    """The extra site is a data difference, not a drawing one.

    CDK7 F91L is a strong singleton at the gatekeeper in the current
    `consensus_sites.csv`; the published panel does not show it. Nothing else
    differs. See docs/DISCREPANCIES.md item 22.
    """
    made, published = rendered("6")
    mine, theirs = set(labels(made)), set(labels(published))
    assert theirs <= mine
    assert mine - theirs == {"CDK7 F91L"}
