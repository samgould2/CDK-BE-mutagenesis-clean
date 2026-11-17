"""Sensor editing by protospacer position (supplementary 3a, 3b).

From `03_editing_analysis/editing_by_CDK_protospacer_location.ipynb`.
One box-and-strip panel per CDK, nine to a figure, one figure per editor.

The per-guide table the panels plot is a reduction of the full canonical
editing calls and is shipped as
`data/processed/editing/editing_by_proto_loc_per_guide_canonical.zip`;
`rebuild()` regenerates it from the `*_full_canonical.zip` files, and
`tests/test_editing_position.py` checks the two agree.

The shaded band on each panel is the library's canonical editing window, taken
from where `ideal_window` sits inside `protospacer` rather than hard-coded.
Note that a guide contributes a row for **every** editable base in its
protospacer, edited or not -- that is what makes the zero values meaningful and
the medians interpretable.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .. import counts as _counts
from .. import io, paths

MIN_READS = 100
EDITORS = {
    "ABE": {"ref": "A", "alt": "G", "color": "tab:purple"},
    "CBE": {"ref": "C", "alt": "T", "color": "tab:blue"},
}
#: editing-file stem per subpool -> the CDKs it covers
SUBPOOLS = {
    "subpool1": ("subpool1", ["CDK7", "CDK8", "CDK9", "CDK19"]),
    "subpool2": ("CDK12_13", ["CDK12", "CDK13"]),
    "subpool3": ("CDK2_4_6", ["CDK2", "CDK4", "CDK6"]),
}
CDK_ORDER = ["CDK2", "CDK4", "CDK6", "CDK7", "CDK8", "CDK9", "CDK12", "CDK13", "CDK19"]
POSITIONS = list(range(1, 21))

TABLE = paths.EDITING / "editing_by_proto_loc_per_guide_canonical.zip"


def window() -> tuple:
    """The library's canonical editing window, in protospacer coordinates."""
    lib = _counts.library()
    w = lib[lib["ideal_window"].notna()]
    start = w.apply(lambda r: r["protospacer"].find(r["ideal_window"]) + 1, axis=1)
    s = int(start.mode()[0])
    return s, s + int(w["ideal_window"].str.len().mode()[0]) - 1


def table() -> pd.DataFrame:
    """The shipped per-guide, per-position editing table."""
    return io.read(TABLE)


def _per_guide(full, library, ref, alt, genes) -> pd.DataFrame:
    full = full.copy()
    edits = full["DNA Change"].str.findall(r"\+(\d+)([ACGT])>([ACGT])")
    # keep only reads whose every edit is the editor's own; a read carrying any
    # other substitution is dropped whole, not partially counted
    keep = edits.apply(lambda e: all(r == ref and a == alt for _, r, a in e))
    full, edits = full[keep], edits[keep]

    totals = full.groupby("gRNA_id")["#Reads"].sum()
    totals = totals[totals >= MIN_READS]

    pos = pd.DataFrame({"gRNA_id": full["gRNA_id"], "#Reads": full["#Reads"],
                        "Location": edits.apply(lambda e: [int(p) for p, _, _ in e])})
    pos = pos.explode("Location").dropna(subset=["Location"])
    pos["Location"] = pos["Location"].astype(int)
    edited = pos.groupby(["gRNA_id", "Location"])["#Reads"].sum().rename("edited_reads")

    lib = library[library["gRNA_id"].isin(totals.index) & library["Gene"].isin(genes)]
    rows = [(g, gene, cls, i + 1)
            for g, gene, cls, proto in lib[["gRNA_id", "Gene", "classification",
                                            "protospacer"]].itertuples(index=False)
            for i, b in enumerate(proto) if b == ref]
    out = pd.DataFrame(rows, columns=["gRNA_id", "Gene", "classification", "Location"])
    out = out.merge(edited.reset_index(), on=["gRNA_id", "Location"], how="left")
    out["edited_reads"] = out["edited_reads"].fillna(0)
    out["total_reads"] = out["gRNA_id"].map(totals)
    out["%Reads"] = 100 * out["edited_reads"] / out["total_reads"]
    return out


def rebuild() -> pd.DataFrame:
    """Recompute the per-guide table from the full canonical editing calls."""
    library = _counts.library()
    out = []
    for editor, e in EDITORS.items():
        for subpool, (stem, genes) in SUBPOOLS.items():
            full = io.read(paths.EDITING / f"{editor}_{stem}_full_canonical.zip")
            d = _per_guide(full, library, e["ref"], e["alt"], genes)
            d.insert(0, "Editor", editor)
            d.insert(1, "Subpool", subpool)
            out.append(d)
    return pd.concat(out, ignore_index=True)


def _box_strip(ax, data, color, win, title=None):
    ax.axvspan(win[0] - 1.5, win[1] - 0.5, color="0.92", zorder=0, lw=0)
    sns.stripplot(data=data, x="Location", y="%Reads", order=POSITIONS, color=color,
                  size=2.2, alpha=0.45, jitter=0.3, ax=ax, zorder=1)
    for c in ax.collections:
        c.set_rasterized(True)      # keep the pdf light; boxes and text stay vector
    line = dict(color="0.15", linewidth=1)
    sns.boxplot(data=data, x="Location", y="%Reads", order=POSITIONS, showfliers=False,
                width=0.6, boxprops=dict(facecolor="none", edgecolor="0.15", linewidth=1),
                medianprops=dict(color="0.15", linewidth=1.8),
                whiskerprops=line, capprops=line, ax=ax, zorder=2)
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_ylim(-2, 102)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="both", which="major", labelsize=10)
    if title:
        ax.set_title(title, fontsize=12)


def draw(editor: str, outdir=None, **kw):
    """-> (figure, the editor's slice of the table): nine CDK panels."""
    d_all = table()
    e = EDITORS[editor]
    win = window()
    fig, axes = plt.subplots(3, 3, figsize=(15, 9), sharex=True, sharey=True)
    for ax, gene in zip(axes.flat, CDK_ORDER):
        d = d_all[(d_all["Editor"] == editor) & (d_all["Gene"] == gene)]
        _box_strip(ax, d, e["color"], win,
                   title=f'{gene} (n = {d["gRNA_id"].nunique():,} guides)')
    for ax in axes[:, 0]:
        ax.set_ylabel(f'{e["ref"]}>{e["alt"]} editing (%)', fontsize=12)
    for ax in axes[-1, :]:
        ax.set_xlabel("Protospacer position", fontsize=12)
    fig.suptitle(f"{editor} sensor editing by protospacer position", fontsize=14)
    fig.tight_layout()
    if outdir is not None:
        fig.savefig(str(outdir / f"editing_by_proto_loc_{editor}_all_CDKs.pdf"),
                    transparent=True, dpi=300)
    return fig, d_all[d_all["Editor"] == editor]


def _wrap(panel, editor):
    f = lambda outdir=None, **kw: draw(editor, outdir)      # noqa: E731
    f.__name__ = f.__qualname__ = f"draw_{panel}"
    f.__doc__ = (f"Sensor editing by protospacer position for {editor}: "
                 "one box-and-strip panel per CDK, nine to the figure.")
    return f


FIGURES = {"S3a": _wrap("S3a", "CBE"), "S3b": _wrap("S3b", "ABE")}
PANELS = {"S3a": lambda: table().query("Editor == 'CBE'"),
          "S3b": lambda: table().query("Editor == 'ABE'")}
