"""Figure 5f-5i: the T47D CDK4 allelic-series RNA-seq panels.

Ported from `data/processed/rnaseq/T47D_RNAseq/cdk4_program_analysis/pub_fig_square.py`,
which the published figure was assembled from. The drawing code is carried over
unchanged -- same constants, same layout table, same statistics -- with three
differences, all mechanical:

    * data loading is deferred into `_data()` and cached, so importing this
      module does not read six tables;
    * the loop that wrote the three layouts at import time is behind
      `if __name__ == "__main__"`;
    * `panels()` cuts the four axes out of the finished figure, so a per-panel
      file is a crop of the published figure rather than a redraw of it.

The published panels are the **mid** layout, `pub_fig_square_4panel_mid.pdf`:

    a -> 5f   CDK4 transcript level per genotype (NMD of the W198* allele)
    b -> 5g   Hallmark E2F target score per sample
    c -> 5h   volcano, W198* vs CBE NT
    d -> 5i   volcano, W198* vs D196N (matched background)

That layout is identified mechanically, not by eye: the three layouts differ in
their y-tick locators, and only the mid one puts {80, 160, 240} on the volcano
axes and 1.0 at the top of the E2F axis. Those are exactly the numbers the
published `figure5_v2.ai` carries; the square layout would need 60/120/180 and
the wide one 100/200/300, and neither appears anywhere in it.
See `tests/test_figure5.py`.

Both volcanoes use contrasts refit with DESeq2 independent filtering disabled,
so the two panels share an identical gene universe and their DEG counts are
comparable. Order throughout: CBE NT, E56K, D196N, W198*.
"""
from __future__ import annotations

import json
from functools import lru_cache

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

from .. import io as _io, paths

BASE = paths.PROCESSED / "rnaseq" / "T47D_RNAseq"
ANALYSIS = BASE / "cdk4_program_analysis"
OUT = ANALYSIS / "results"
FIG = ANALYSIS / "figures"

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
}

BASEMEAN_MIN = 50
LFC_T, PADJ_T = 1.0, 0.05
X_LIM = 6

CTRL = "CBE_NT"
VARIANTS = ["CDK4_E56K", "CDK4_D196N", "CDK4_W198"]
GENOTYPES = [CTRL] + VARIANTS
LABEL = {"CBE_NT": "CBE NT", "CDK4_E56K": "E56K",
         "CDK4_D196N": "D196N", "CDK4_W198": "W198*"}
COLOR = {"CBE_NT": "#9E9E9E", "CDK4_E56K": "#92C5DE",
         "CDK4_D196N": "#2166AC", "CDK4_W198": "#B2182B"}
UP, DOWN, NS = "#B2182B", "#2166AC", "#D9D9D9"

#: (name, figsize, bracket offset, bracket step, headroom, volcano y-ticks, letter x)
LAYOUTS = [
    ("pub_fig_square_4panel",      (6.6, 6.4), .13, .17, .78, 6, -.20),
    ("pub_fig_square_4panel_mid",  (7.6, 5.6), .11, .15, .68, 5, -.17),
    ("pub_fig_square_4panel_wide", (8.8, 4.9), .10, .13, .58, 4, -.13),
]

#: the layout the manuscript used, and the panel letter each of its axes became
PUBLISHED_LAYOUT = "pub_fig_square_4panel_mid"
PANEL_LETTER = {"a": "5f", "b": "5g", "c": "5h", "d": "5i"}


def stars(p: float) -> str:
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < .05 else "n.s."


def holm(p):
    p = np.asarray(p, float)
    o = np.empty(len(p))
    r = 0.
    for k, i in enumerate(np.argsort(p)):
        r = max(r, min(1., (len(p) - k) * p[i]))
        o[i] = r
    return o


class _Data:
    """Everything `build()` reads, loaded once."""

    def __init__(self):
        tpm = _io.read(BASE / "salmon.merged.gene_tpm.tsv", sep="\t", index_col=0)
        gn = tpm.pop("gene_name")
        meta = pd.Series(
            tpm.columns.str.extract(r"^T47D_(.+)_veh_\d+$")[0].values,
            index=tpm.columns)
        samples = [s for g in GENOTYPES for s in meta.index[meta == g]]
        self.tpm, self.meta = tpm[samples], meta[samples]
        self.log_tpm = np.log2(self.tpm + 1)
        order = self.tpm.mean(axis=1).sort_values(ascending=False).index
        s2i = pd.Series(order, index=gn.reindex(order).values)
        self.s2i = s2i[~s2i.index.duplicated()]

        self.res_nt = {v: _io.read(OUT / f"{v}_vs_{CTRL}_veh.csv")
                       .set_index("gene_id") for v in VARIANTS}
        self.volc_nt = _io.read(
            OUT / "CDK4_W198_vs_CBE_NT_veh_nofilt.csv").set_index("gene_id")
        self.volc_wd = _io.read(
            OUT / "CDK4_W198_vs_CDK4_D196N_veh_nofilt.csv").set_index("gene_id")
        hall = json.loads(_io.resolve(BASE / "T47D_CDK4_20260927_results" /
                                      "hallmark_gene_sets_cache.json").read_text())

        self.e2f_ids = [self.s2i[s] for s in hall["E2F Targets"]
                        if s in self.s2i.index
                        and self.log_tpm.loc[self.s2i[s]].mean() > 1]
        m = self.log_tpm.loc[self.e2f_ids]
        self.e2f_score = (m.sub(m.mean(axis=1), axis=0)
                           .div(m.std(axis=1).replace(0, np.nan), axis=0)
                           .mean(axis=0))
        ps = []
        for v in VARIANTS:
            t, p2 = sps.ttest_ind(self.e2f_score[self.grp(v)],
                                  self.e2f_score[self.grp(CTRL)], equal_var=False)
            ps.append(p2 / 2 if t < 0 else 1 - p2 / 2)
        self.e2f_adj = holm(ps)

    def grp(self, g):
        return list(self.meta.index[self.meta == g])


@lru_cache(maxsize=1)
def _data() -> _Data:
    return _Data()


def volcano(ax, r, title, label_genes, yn):
    d = r.dropna(subset=["padj", "log2FoldChange"])
    d = d[d.baseMean >= BASEMEAN_MIN].copy()
    d["y"] = -np.log10(d["padj"].clip(lower=1e-300))
    d["x"] = d["log2FoldChange"].clip(-X_LIM * .98, X_LIM * .98)
    up = (d.padj < PADJ_T) & (d.log2FoldChange > LFC_T)
    dn = (d.padj < PADJ_T) & (d.log2FoldChange < -LFC_T)
    ax.scatter(d.loc[~(up | dn), "x"], d.loc[~(up | dn), "y"], s=2.5,
               color=NS, lw=0, rasterized=True)
    ax.scatter(d.loc[dn, "x"], d.loc[dn, "y"], s=5, color=DOWN, lw=0, rasterized=True)
    ax.scatter(d.loc[up, "x"], d.loc[up, "y"], s=5, color=UP, lw=0, rasterized=True)
    for xv in (LFC_T, -LFC_T):
        ax.axvline(xv, color="#999999", lw=.7, ls="--")
    ax.axhline(-np.log10(PADJ_T), color="#999999", lw=.7, ls="--")
    for sym in label_genes:
        hit = d[d.gene_name == sym]
        if not len(hit):
            continue
        rr = hit.iloc[0]
        ax.scatter(rr["x"], rr["y"], s=26, facecolor="none", edgecolor="black",
                   lw=.9, zorder=5)
        ax.annotate(sym, (rr["x"], rr["y"]), fontsize=6.5, style="italic",
                    zorder=6, xytext=(4, 3), textcoords="offset points")
    ax.set_xlim(-X_LIM, X_LIM)
    ax.set_xlabel("log$_2$ fold change", fontsize=8.5)
    ax.set_ylabel("$-$log$_{10}$ $P_{adj}$", fontsize=8.5)
    ax.set_title(title, fontsize=8.5, pad=3)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(yn))
    ax.tick_params(labelsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(.03, .055, f"{int(dn.sum())} down", transform=ax.transAxes,
            fontsize=7, color=DOWN, va="bottom", ha="left")
    ax.text(.97, .055, f"{int(up.sum())} up", transform=ax.transAxes,
            fontsize=7, color=UP, va="bottom", ha="right")
    return int(up.sum()), int(dn.sum()), len(d)


def _cdk4_tpm(ax):
    """Panel a / figure 5f -- CDK4 transcript level per genotype."""
    D = _data()
    tpm, grp = D.tpm, D.grp
    g = D.s2i["CDK4"]
    rng = np.random.default_rng(1)
    for xi, geno in enumerate(GENOTYPES):
        v = tpm.loc[g, grp(geno)].to_numpy()
        ax.bar(xi, v.mean(), width=.64, color=COLOR[geno], edgecolor="black",
               lw=.6, zorder=2)
        ax.errorbar(xi, v.mean(), yerr=v.std(), color="black", lw=.7,
                    capsize=2.2, zorder=3)
        ax.scatter(xi + rng.uniform(-.1, .1, 3), v, s=9, color="black",
                   zorder=4, lw=0)
    for xi, geno in enumerate(VARIANTS, start=1):
        p = D.res_nt[geno].loc[g, "padj"]
        ax.text(xi, tpm.loc[g, grp(geno)].mean() + tpm.loc[g, grp(geno)].std() + 16,
                stars(p), ha="center", fontsize=8 if p < .05 else 6.2)
    ax.set_ylim(0, 480)
    ax.set_ylabel("CDK4 (TPM)", fontsize=9)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LABEL[x] for x in GENOTYPES], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)


def _e2f_score(ax, br0, brs, head):
    """Panel b / figure 5g -- Hallmark E2F target score per sample.

    `br0`, `brs` and `head` are the significance-bracket geometry, which is the
    only thing that differs between the three layouts in this panel.
    """
    D = _data()
    grp = D.grp
    rng = np.random.default_rng(0)
    vals = {geno: D.e2f_score[grp(geno)].to_numpy() for geno in GENOTYPES}
    for xi, geno in enumerate(GENOTYPES):
        ax.scatter(xi + rng.uniform(-.1, .1, 3), vals[geno], s=28,
                   facecolor=COLOR[geno], edgecolor="black", lw=.5,
                   zorder=3, clip_on=False)
        ax.plot([xi - .27, xi + .27], [vals[geno].mean()] * 2, color="black",
                lw=1.3, zorder=4)
    top = max(v.max() for v in vals.values())
    span = top - min(v.min() for v in vals.values())
    for k, v in enumerate(VARIANTS):
        xv = GENOTYPES.index(v)
        yb = top + span * (br0 + brs * k)
        ax.plot([0, 0, xv, xv], [yb - span * .04, yb, yb, yb - span * .04],
                color="black", lw=.7)
        st = stars(D.e2f_adj[k])
        ax.text(xv / 2, yb + span * .01, st, ha="center", va="bottom",
                fontsize=6.4 if st == "n.s." else 8)
    ax.set_ylim(top=top + span * head)
    ax.set_ylabel(f"Hallmark E2F targets\n(mean $z$-score, $n$ = {len(D.e2f_ids)})",
                  fontsize=8.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LABEL[x] for x in GENOTYPES], rotation=45, ha="right",
                       rotation_mode="anchor", fontsize=8)
    ax.tick_params(labelsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)


#: panel letter -> (drawing function, the arguments it takes from the layout row)
_DRAW = {
    "a": lambda ax, br0, brs, head, yn: _cdk4_tpm(ax),
    "b": lambda ax, br0, brs, head, yn: _e2f_score(ax, br0, brs, head),
    "c": lambda ax, br0, brs, head, yn: volcano(
        ax, _data().volc_nt, "W198* vs CBE NT", ["CDK4"], yn),
    "d": lambda ax, br0, brs, head, yn: volcano(
        ax, _data().volc_wd, "W198* vs D196N", ["CDK4", "PIP", "MGP"], yn),
}


def build(figsize, br0, brs, head, yn, lx):
    """The published 2x2 figure. -> (fig, (up1, down1, n1, up2, down2, n2))."""
    with mpl.rc_context(RC):
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        for ax, L in zip(axes.flat, "abcd"):
            ax.text(lx, 1.07, L, transform=ax.transAxes, fontsize=11,
                    fontweight="bold", va="top")
        counts = {}
        for ax, L in zip(axes.flat, "abcd"):
            counts[L] = _DRAW[L](ax, br0, brs, head, yn)
        fig.tight_layout()
    u1, d1, n1 = counts["c"]
    u2, d2, n2 = counts["d"]
    return fig, (u1, d1, n1, u2, d2, n2)


def panel(key: str, layout: str = PUBLISHED_LAYOUT, figsize=(3.9, 2.9)):
    """One published panel ("5f".."5i") on its own figure.

    The same drawing function the 2x2 uses, on a single-axes figure, with the
    layout's own bracket geometry and tick count so the panel is the one that
    was published rather than a differently-scaled sibling. `tests/test_figure5.py`
    pins it against the 2x2.
    """
    letter = {v: k for k, v in PANEL_LETTER.items()}[key]
    spec = {n: rest for n, *rest in LAYOUTS}[layout]
    br0, brs, head, yn = spec[1], spec[2], spec[3], spec[4]
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=figsize)
        _DRAW[letter](ax, br0, brs, head, yn)
        fig.tight_layout()
    return fig


def save(fig, path, dpi: int = 400, **kw):
    """Write a figure from this module.

    **Save inside the rc context, not just draw inside it.** `pdf.fonttype` and
    `svg.fonttype` are read by the *backend at save time*, so a figure drawn
    under `rc_context(RC)` and saved outside it comes out with Type 3 text
    instead of the TrueType the manuscript needs -- same picture, different
    glyph encoding, and it no longer matches the shipped PDF. Everything that
    writes one of these figures goes through here.
    """
    with mpl.rc_context(RC):
        fig.savefig(path, dpi=dpi, bbox_inches="tight", **kw)


def published(layout: str = PUBLISHED_LAYOUT):
    """The 2x2 figure in the layout the manuscript used."""
    spec = {n: rest for n, *rest in LAYOUTS}
    if layout not in spec:
        raise KeyError(f"{layout!r} is not one of {sorted(spec)}")
    return build(*spec[layout])


def panels(layout: str = PUBLISHED_LAYOUT) -> dict:
    """-> {"5f": fig, "5g": fig, "5h": fig, "5i": fig}, one figure each."""
    return {k: panel(k, layout) for k in PANEL_LETTER.values()}


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    for name, figsize, br0, brs, head, yn, lx in LAYOUTS:
        fig, counts = build(figsize, br0, brs, head, yn, lx)
        for e in ("png", "svg", "pdf"):
            save(fig, FIG / f"{name}.{e}")
        plt.close(fig)
        print(f"{name:30s} {figsize[0]} x {figsize[1]} in")
    u1, d1, n1, u2, d2, n2 = counts
    print(f"\nbaseMean >= {BASEMEAN_MIN}, identical universe ({n1} genes):")
    print(f"  W198* vs CBE NT   {u1} up / {d1} down")
    print(f"  W198* vs D196N    {u2} up / {d2} down")
