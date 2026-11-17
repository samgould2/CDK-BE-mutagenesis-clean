"""Shared-y-axis versions of the combined elongation figure.

Fig10a_elongation_shared_y        all three panels on one y-scale — the baseline effect is
                                  shown at its true size relative to acute drug
Fig10b_elongation_shared_y_inset  same, plus an inset in the baseline panel zoomed onto the
                                  variant-vs-NT curves

Error bars are mean +/- SD across the three biological replicates, from VST-based
replicate-level log2 ratios (see 06_elongation_combined.py for why).

Run:  cd reanalysis/scripts && python 07_elongation_shared_y.py
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from _style import (RES, FIG, GENOS, VARS, DRUGS, GCOL, TLAB, INK, MUTED,
                    save, expressed_pc, meta as load_meta, replicate_profiles)

meta = load_meta()
vst = pd.read_csv(os.path.join(RES, 'vst.csv'), index_col='gene_id')
vst.pop('gene_name')
pc = expressed_pc()
span = pc.span_bp
v = vst.loc[pc.index, meta.index]

bins = pd.qcut(span.rank(method='first'), 10, labels=False)
centre = span.groupby(bins).median()

SHARED_YLIM = (-0.45, 0.85)
ZOOM_YLIM = (-0.10, 0.14)


def samples(g, t):
    return list(meta.index[(meta.genotype == g) & (meta.treatment == t)])


def draw(ax, g, test, ref, ms=3.0, lw=1.3, clip=False):
    R = replicate_profiles(v, bins, test, ref)
    ax.errorbar(centre, R.mean(), yerr=R.std(ddof=1), color=GCOL[g], lw=lw, marker='o',
                ms=ms, mec=INK, mew=0.5, capsize=1.6, capthick=0.6, elinewidth=0.7,
                label=g, zorder=3, clip_on=clip)


def length_axis(ax):
    ax.axhline(0, color=MUTED, lw=0.6, ls='--', zorder=1)
    ax.set_xscale('log')
    ax.set_xlabel('Gene length (bp)')
    ax.set_axisbelow(True)


def build(with_inset):
    fig, axes = plt.subplots(1, 3, figsize=(6.7, 2.05), sharey=True,
                             gridspec_kw={'wspace': 0.14})
    for c_i, d in enumerate(DRUGS):
        ax = axes[c_i]
        for g in GENOS:
            draw(ax, g, samples(g, d), samples(g, 'Control'))
        ax.set_title(f'{TLAB[d]} vs DMSO', fontsize=6.5)
    ax = axes[2]
    for g in VARS:
        draw(ax, g, samples(g, 'Control'), samples('NT', 'Control'))
    ax.set_title('No drug: variant vs NT', fontsize=6.5)

    for a in axes:
        length_axis(a)
    axes[0].set_ylim(*SHARED_YLIM)
    axes[0].set_ylabel('median log$_2$ ratio per decile\n(mean ± SD, n = 3)')
    axes[0].legend(loc='upper right', handlelength=1.2, fontsize=5.5, ncol=2, columnspacing=0.8)

    if with_inset:
        # Inset zooms the y-axis only; x range is the full length axis, so the
        # indicator rectangle spans the panel width.
        ins = axes[2].inset_axes([0.15, 0.52, 0.82, 0.40])
        for g in VARS:
            draw(ins, g, samples(g, 'Control'), samples('NT', 'Control'),
                 ms=2.2, lw=1.0, clip=True)
        ins.axhline(0, color=MUTED, lw=0.5, ls='--', zorder=1)
        ins.set_xscale('log')
        ins.set_ylim(*ZOOM_YLIM)
        ins.set_xlim(axes[2].get_xlim())
        ins.tick_params(labelsize=4.5, width=0.4, length=1.8, pad=1.2)
        ins.set_xticklabels([])
        ins.set_yticks([-0.05, 0, 0.05, 0.10])
        ins.spines[['top', 'right']].set_visible(False)
        for sp in ins.spines.values():
            sp.set_linewidth(0.5)
        ins.set_facecolor('white')
        # Keep the hairline rectangle marking the zoomed y-window, but drop the
        # connector lines: the zoom is y-only across the full x range, so they run
        # from a thin strip to a tall box and read as stray diagonals.
        _rect, conns = axes[2].indicate_inset_zoom(ins, edgecolor=MUTED, lw=0.5, alpha=1.0)
        for c in conns:
            c.set_visible(False)
        axes[2].legend(loc='lower left', handlelength=1.2, fontsize=5.5)
    else:
        axes[2].legend(loc='lower left', handlelength=1.2, fontsize=5.5)
    return fig


save(build(with_inset=False), 'Fig10a_elongation_shared_y')
save(build(with_inset=True), 'Fig10b_elongation_shared_y_inset')
