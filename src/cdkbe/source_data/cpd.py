"""Cumulative population doublings during each screen (supplementary 1a-1d).

From `02_library_qc_and_counts/CPDS.ipynb` cells 2, 4, 6, 7.
Every panel is the same plot -- CPD against days, one line per condition, CBE
left and ABE right -- over a different sheet of `CDK_BE_cell_counts.xlsx`, so
the four are specified here as data.

Two things the original did that are easy to lose:

* subpool 1 labels conditions by *drug and concentration* (`KB-0742 LOW` /
  `HIGH`), because several drugs were run at two doses in one screen, and then
  **drops four arms** -- `SEL120 HIGH`, both SY-5609 doses and `KB-0742 HIGH`.
  The SY-5609 arms were rerun as their own screen, which is panel 1b.
* every sheet carries a stale pivoted block to the right of the real columns,
  under blank headers; `_sheet` drops it, so a panel's source table is one
  table rather than two side by side.
* the CDK12/13 and CDK2/4/6 sheets carry a `NONE` drug row that is filtered
  out, and the SY-5609 sheets carry blank rows that are filled with 0 and then
  dropped on `Drug != 0`. Both are sheet artefacts, not conditions.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .. import paths

COUNTS = paths.REFERENCE / "cell_counts" / "CDK_BE_cell_counts.xlsx"

#: subpool 1: (drug, concentration nM) -> condition label
_S1_CONDITIONS = {
    ("KI-CDK9d-32", 100): "KI-CDK9d-32 LOW", ("KI-CDK9d-32", 1000): "KI-CDK9d-32 HIGH",
    ("KI-CDK9d-32N", 1250): "KI-CDK9d-32N LOW", ("KI-CDK9d-32N", 5000): "KI-CDK9d-32N HIGH",
    ("KB-0742", 1500): "KB-0742 LOW", ("KB-0742", 5000): "KB-0742 HIGH",
    ("SY-5609", 1000): "SY-5609 LOW", ("SY-5609", 2500): "SY-5609 HIGH",
    ("Senexin B", 2000): "Senexin B LOW", ("Senexin B", 15000): "Senexin B HIGH",
    ("SEL120", 4000): "SEL120 LOW", ("SEL120", 8000): "SEL120 HIGH",
}
#: arms present in the sheet but not plotted (SY-5609 got its own screen)
_S1_EXCLUDE = ["SEL120 HIGH", "SY-5609 HIGH", "SY-5609 LOW", "KB-0742 HIGH"]
_S1_COLOURS = {
    "DMSO": "black", "KB-0742 LOW": "lightcoral",
    "KI-CDK9d-32 HIGH": "tab:purple", "KI-CDK9d-32 LOW": "thistle",
    "KI-CDK9d-32N HIGH": "crimson", "KI-CDK9d-32N LOW": "palevioletred",
    "SEL120 LOW": "tab:green", "Senexin B HIGH": "darkblue",
    "Senexin B LOW": "cadetblue",
}

_SY_CONDITIONS = {("SY-5609", 10): "SY-5609 LOW", ("SY-5609", 100): "SY-5609 HIGH"}
_SY_COLOURS = {"DMSO": "black", "SY-5609 LOW": "cyan", "SY-5609 HIGH": "darkcyan"}

_S2_COLOURS = {"HQ461": "darkolivegreen", "BSJ-4-116": "tab:green",
               "CDK12-IN-2": "lightgreen", "DMSO": "black"}
_S3_ORDER = ["DMSO", "Ribociclib", "Abemaciclib", "Atirmociclib",
             "Palbociclib", "Tagtociclib", "INX-315"]
_S3_COLOURS = dict(zip(_S3_ORDER, ["black", "lightcoral", "indianred", "firebrick",
                                   "maroon", "darkorange", "goldenrod"]))

#: panel -> (CBE sheet, ABE sheet, hue column, colours, hue order, legend cols, xticks)
SPEC = {
    "S1a": ("CBE_screen", "ABE_screen", "Condition", _S1_COLOURS, None, 2, range(0, 22, 3)),
    "S1b": ("CBE SY-5609", "ABE SY-5609", "Condition", _SY_COLOURS, None, 1, range(0, 31, 3)),
    "S1c": ("CBE_CDK12_13", "ABE_CDK12_13", "Drug", _S2_COLOURS, None, 1, range(0, 22, 3)),
    "S1d": ("CBE_CDK2_4_6", "ABE_CDK2_4_6", "Drug", _S3_COLOURS, _S3_ORDER, 1, range(0, 22, 3)),
}


def _sheet(name: str, panel: str) -> pd.DataFrame:
    d = pd.read_excel(COUNTS, name)
    # Each sheet carries a second, pivoted block off to the right -- days
    # across, one column per replicate -- under blank headers, which pandas
    # reads as `Unnamed: 16` onwards. It is an older working view of the same
    # counts, nothing here reads it, and left in it turns the panel's source
    # table into two tables side by side. The workbook itself is untouched.
    d = d[[c for c in d.columns if not str(c).startswith("Unnamed:")]]
    if panel in ("S1a", "S1b"):
        # the sheet's own Condition column comes back all-NaN (float); assigning
        # labels into it otherwise warns about an incompatible dtype
        d["Condition"] = pd.Series([pd.NA] * len(d), dtype="object")
    if panel == "S1a":
        for (drug, conc), label in _S1_CONDITIONS.items():
            d.loc[(d["Drug"] == drug) & (d["Concentration (nM)"] == conc), "Condition"] = label
        d.loc[d["Drug"] == "DMSO", "Condition"] = "DMSO"
        d = d[~d["Condition"].isin(_S1_EXCLUDE)].dropna(subset="Condition")
        return d.reset_index(drop=True)
    if panel == "S1b":
        for (drug, conc), label in _SY_CONDITIONS.items():
            d.loc[(d["Drug"] == drug) & (d["Concentration (nM)"] == conc), "Condition"] = label
        d.loc[d["Drug"] == "DMSO", "Condition"] = "DMSO"
        d = d.fillna(0)
        return d[d["Drug"] != 0].reset_index(drop=True)
    return d[d["Drug"] != "NONE"].reset_index(drop=True)


def frames(panel: str) -> dict:
    """-> {'CBE': frame, 'ABE': frame} for one panel."""
    cbe, abe, *_ = SPEC[panel]
    return {"CBE": _sheet(cbe, panel), "ABE": _sheet(abe, panel)}


def source(panel: str) -> pd.DataFrame:
    """Both editors in one frame, for the source-data workbook."""
    return pd.concat([f.assign(Editor=e) for e, f in frames(panel).items()],
                     ignore_index=True)


def draw(panel: str, outdir=None, **kw):
    """-> (figure, combined frame)."""
    _, _, hue, colours, order, ncols, xticks = SPEC[panel]
    d = frames(panel)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5), sharey=True, sharex=True)
    for a, editor in zip(ax, ("CBE", "ABE")):
        sns.lineplot(data=d[editor], x="Days Post-Treatment", y="CPDs", hue=hue,
                     marker="o", ax=a, palette=colours, hue_order=order)
        a.spines[["right", "top"]].set_visible(False)
        a.tick_params(axis="both", which="major", labelsize=12)
        a.set_xlabel("Days", fontsize=14)
        a.set_title(editor, fontsize=14)
    ax[0].set_xticks(list(xticks))
    ax[0].set_ylabel("CPD", fontsize=14)
    ax[0].legend(fontsize=9, ncols=ncols, loc="upper left", frameon=False)
    ax[1].get_legend().remove()
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / f"Supp_{panel}_CPD.pdf"), transparent=True)
    return fig, source(panel)


def _wrap(panel):
    cbe, abe, group, *_ = SPEC[panel]
    f = lambda **kw: draw(panel, **kw)          # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    f.__doc__ = (f"Cumulative population doublings against days, one line per "
                 f"{group.lower()}: CBE ({cbe}) left, ABE ({abe}) right.")
    return f


FIGURES = {k: _wrap(k) for k in SPEC}
PANELS = {k: (lambda p=k: source(p)) for k in SPEC}
