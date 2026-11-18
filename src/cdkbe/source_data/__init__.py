"""Per-panel source-data tables for the manuscript.

One function per figure panel. Each returns the dataframe that the panel plots
-- not a re-derivation, but the same object the figure's own code passed to
matplotlib, reached by replaying that code's filters against the processed
tables in `data/`.

Every function carries the original notebook and cell it came from, so a number
in a workbook can be traced back to the code that drew it.

    from cdkbe.source_data import PANELS
    df = PANELS["2m"]()
"""
from __future__ import annotations
from .allele_frequency import source as figure_4f
from .structures import PANELS as _struct_panels
from .allele_specific_mtt import MAIN_PANELS as _mtt_main_panels
from .competition_assays import MAIN_PANELS as _comp_main_panels
from .t47d_cell_cycle import MAIN_PANELS as _cc_main_panels
from .rnaseq_t47d_palbo import MAIN_PANELS as _palbo_main_panels
from .klifs_msa import MAIN_PANELS as _kmsa_main_panels
from .selection_allele_frequency import MAIN_PANELS as _sel_main_panels
from .invivo import figure_2m, figure_2n, figure_4n, figure_4o
from .library import figure_1e, figure_1f
from .screens import (figure_2a, figure_2b, figure_2h, figure_2i, figure_3a,
                      figure_3b, figure_4a, figure_4b, figure_4p)
from .rnaseq import figure_3f, figure_3g, figure_3h
from .validation_mtt import stats as figure_3d
from .rnaseq_t47d import figure_5f, figure_5g, figure_5h, figure_5i
from .singlecell import (figure_5a, figure_5b, figure_5c, figure_5d, figure_5e,
                         figure_5j, figure_5k)

PANELS = {
    "1e": figure_1e, "1f": figure_1f,
    "2a": figure_2a, "2b": figure_2b, "2h": figure_2h, "2i": figure_2i,
    "2c": _struct_panels["2c"],
    "2d": _mtt_main_panels["2d"], "2f": _comp_main_panels["2f"],
    "2j": _struct_panels["2j"],
    "2k": _mtt_main_panels["2k"],
    "2m": figure_2m, "2n": figure_2n,
    "3a": figure_3a, "3b": figure_3b, "3c": _struct_panels["3c"],
    "3d": figure_3d,
    "3f": figure_3f, "3g": figure_3g, "3h": figure_3h,
    "4a": figure_4a, "4b": figure_4b, "4c": _struct_panels["4c"],
    "4d": _struct_panels["4d"], "4e": _struct_panels["4e"],
    "4f": figure_4f, "4g": _sel_main_panels["4g"],
    "4k": _cc_main_panels["4k"],
    "4m": _palbo_main_panels["4m"],
    "4l": _cc_main_panels["4l"], "4p": figure_4p,
    "4n": figure_4n, "4o": figure_4o,
    "5a": figure_5a, "5b": figure_5b, "5c": figure_5c, "5d": figure_5d,
    "5e": figure_5e, "5f": figure_5f, "5g": figure_5g, "5h": figure_5h,
    "5i": figure_5i, "5j": figure_5j, "5k": figure_5k,
    "6": _kmsa_main_panels["6"],
}

#: panel -> the function that *draws* it, merged from the plotting modules. A panel in PANELS but not here has its source data but is not
#: redrawn yet; docs/COVERAGE.md says why for each.
def _drawn():
    from .allele_frequency import FIGURES as _af
    from .allele_specific_mtt import MAIN_FIGURES as _amtt_main
    from .competition_assays import MAIN_FIGURES as _comp_main
    from .t47d_cell_cycle import MAIN_FIGURES as _cc_main
    from .rnaseq_t47d_palbo import MAIN_FIGURES as _palbo_main
    from .klifs_msa import MAIN_FIGURES as _kmsa_main
    from .selection_allele_frequency import MAIN_FIGURES as _sel_main
    from .figure5 import FIGURES as _f5
    from .structures import FIGURES as _struct
    from .invivo_figures import FIGURES as _inv
    from .library_figures import FIGURES as _lib
    from .rnaseq_figures import FIGURES as _rna
    from .screen_figures import FIGURES as _scr
    from .validation_mtt import FIGURES as _vmtt
    merged = {}
    for reg in (_lib, _scr, _rna, _vmtt, _inv, _f5, _amtt_main, _af, _struct,
                _comp_main, _cc_main, _kmsa_main, _sel_main, _palbo_main):
        merged.update(reg)
    unknown = sorted(set(merged) - set(PANELS))
    if unknown:  # a drawer registered under a panel that has no source data
        raise RuntimeError(f"drawn panels missing from PANELS: {unknown}")
    return merged


DRAWN = _drawn()


def _supp():
    """Supplementary panels, merged from the modules that provide them.

    Keys are "S" + the mapping's own figure/panel spelling, e.g. "S12a". A
    panel appears in SUPP_PANELS when its source data reproduces and in
    SUPP_DRAWN when the panel itself is redrawn here.
    """
    from . import cdk9_boltz as _boltz
    from . import competition_assays as _comp
    from . import lncap_pc3 as _lncap
    from . import t47d_cell_cycle as _cc
    from . import selection_allele_frequency as _sel
    from . import cpd as _cpd
    from . import editing_position as _ep
    from . import allele_specific_mtt as _amtt
    from . import beatseq_qc as _bq
    from . import essentiality_panels as _ess
    from . import klifs_msa as _kmsa
    from . import klifs_resistance as _kres
    from . import leukemia_burden as _lb
    from . import resistance_variants as _rvar
    from . import rnaseq_t47d_palbo as _palbo
    from . import rnaseq_volcanoes as _rvol
    from . import screen_dose_response as _sdr
    from . import supp_figures as _sup

    panels, drawn = {}, {}
    for mod in (_cpd, _ep, _lb, _amtt, _sup, _kres, _kmsa, _ess, _bq,
                _rvar, _sdr, _rvol, _boltz, _comp, _lncap, _cc, _sel,
                _palbo):
        panels.update(mod.PANELS)
        drawn.update(mod.FIGURES)
    # Several of these modules also hold main-figure panels (2d and 2k share a
    # plate with supplementary 9d-9f). Those live in their own dicts; anything
    # that leaks into `PANELS`/`FIGURES` would be counted as supplementary.
    stray = sorted(k for k in set(panels) | set(drawn) if not k.startswith("S"))
    if stray:
        raise RuntimeError(f"main-figure keys in the supplementary registry: {stray}")
    return panels, drawn


SUPP_PANELS, SUPP_DRAWN = _supp()

__all__ = ["PANELS", "DRAWN", "SUPP_PANELS", "SUPP_DRAWN"] + sorted(
    f"figure_{k}" for k in PANELS)
