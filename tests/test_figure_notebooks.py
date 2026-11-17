"""Every panel a figure notebook claims to draw must actually be wired up.

A panel can be registered as a drawer (`screen_figures.FIGURES`) without being
registered as source data (`source_data.PANELS`), or the reverse. The notebook
then fails at run time with a bare KeyError, which is exactly what happened to
2i. These checks are cheap and catch it before anyone opens the notebook.
"""
import json
import re
from pathlib import Path

import pytest

from cdkbe.source_data import PANELS
from cdkbe.source_data import invivo_figures as ivf
from cdkbe.source_data import library_figures as lbf
from cdkbe.source_data import screen_figures as sf

NB_DIR = Path(__file__).resolve().parents[1] / "notebooks" / "08_figures"
DRAWERS = set(ivf.FIGURES) | set(sf.FIGURES) | set(lbf.FIGURES)


def _panels_referenced(path: Path) -> set:
    nb = json.loads(path.read_text())
    src = "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
    return set(re.findall(r'"(\d[a-z])"', src))


@pytest.mark.parametrize("nb", sorted(NB_DIR.glob("figure_*.ipynb")), ids=lambda p: p.stem)
def test_referenced_panels_have_source_data(nb):
    missing = _panels_referenced(nb) - set(PANELS)
    assert not missing, f"{nb.name} asks for {sorted(missing)}, not in source_data.PANELS"


def test_every_drawer_has_source_data():
    missing = DRAWERS - set(PANELS)
    assert not missing, f"drawn but no source data: {sorted(missing)}"


@pytest.mark.parametrize("nb", sorted(NB_DIR.glob("figure_*.ipynb")), ids=lambda p: p.stem)
def test_notebook_is_valid(nb):
    d = json.loads(nb.read_text())
    assert d.get("cells"), f"{nb.name} has no cells"
    assert d.get("nbformat") == 4


SUPP_DIR = Path(__file__).resolve().parents[1] / "notebooks" / "09_supplementary"


@pytest.mark.parametrize("nb", sorted(SUPP_DIR.glob("supp_figure_*.ipynb")),
                         ids=lambda p: p.stem)
def test_supplementary_notebook_is_valid(nb):
    d = json.loads(nb.read_text())
    assert d.get("cells"), f"{nb.name} has no cells"
    # the panel table is the point of these; it must not be empty
    md = "\n".join("".join(c["source"]) for c in d["cells"] if c["cell_type"] == "markdown")
    assert "| panel | source | status | notes |" in md, f"{nb.name} has no panel table"


def test_all_sixteen_supplementary_figures_scaffolded():
    assert len(list(SUPP_DIR.glob("supp_figure_*.ipynb"))) == 16
