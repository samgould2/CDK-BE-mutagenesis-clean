#!/usr/bin/env python3
"""
One figure explaining the correction, for a reader who has not read the methods.

Five panels, left to right: what the raw ratio contains, the evidence that the artifact
is real, what each of the two steps removes, and what it buys at the residue level.

Usage:  python3 correction_explained.py
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from . import t0_plasmid_essentiality as E
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
from .essentiality_figures import INK, INK2, MUTED, DEPLETED, ENRICHED, GENES, FIG, TAB

SHOW = ('CBE', 'Subpool2')      # the screen with the largest artifact, so it is visible
say = E.say


def main():
    lib = pd.read_csv(f'{LIBRARY}/CDK_library_final.csv')[
        ['gRNA_id', 'protospacer', 'PAM']]
    lib['GC'] = E.seq_features(lib).GC.values
    gc = dict(zip(lib.gRNA_id, lib.GC))
    g = pd.read_csv(f'{TAB}/guide_essentiality.csv', low_memory=False)
    g['GC'] = g.gRNA_id.map(gc)
    d = g[(g.Editor == SHOW[0]) & (g.Screen_ID == SHOW[1])].dropna(subset=['GC'])
    ctrl = d[d.classification.isin(E.NULL_CLASSES)]
    targ = d[d.classification == 'targeting']

    ladder = pd.read_csv(f'{TAB}/correction_ladder.csv')
    levels, res_gc, frac_depl = [], [], []
    # the unsuffixed directory is the default level, which is now L2
    for lev, suf in [('L1', '_L1'), ('L2', ''), ('L3', '_L3')]:
        t = f'{LFC_DIR}/T0_plasmid_essentiality{suf}'
        if not os.path.exists(f'{t}/residue_guide_rows.csv'):
            continue
        rw = pd.read_csv(f'{t}/residue_guide_rows.csv').rename(columns={'%Reads': 'pct'})
        rw['GC'] = rw.gRNA_id.map(gc)
        r = pd.read_csv(f'{t}/residue_essentiality.csv')
        rg = (rw.groupby(['Gene', 'Codon'])
              .apply(lambda x: np.average(x.GC, weights=x.pct)).rename('res_GC').reset_index())
        m = r.merge(rg, on=['Gene', 'Codon'])
        sig = m[m.FDR_early < 0.05]
        levels.append(lev)
        res_gc.append(scipy.stats.pearsonr(m.res_GC, m.ess_early).statistic)
        frac_depl.append((sig.ess_early < 0).mean())

    fig = plt.figure(figsize=(13.6, 3.5))
    gs = fig.add_gridspec(1, 5, wspace=0.42, left=0.045, right=0.985, top=0.72, bottom=0.17)

    # ---- (a) what the ratio spans
    ax = fig.add_subplot(gs[0]); ax.axis('off')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    for y, lab, sub, c in [(0.84, 'plasmid pool', 'plasmid DNA', '#DADADA'),
                           (0.38, 'T0 cells', 'genomic DNA', '#DADADA')]:
        ax.add_patch(FancyBboxPatch((0.08, y - 0.1), 0.84, 0.2, boxstyle='round,pad=0.012',
                                    fc=c, ec='none'))
        ax.text(0.5, y + 0.03, lab, ha='center', fontsize=8, color=INK)
        ax.text(0.5, y - 0.055, sub, ha='center', fontsize=6.5, color=INK2, style='italic')
    ax.add_patch(FancyArrowPatch((0.33, 0.72), (0.33, 0.50), arrowstyle='-|>',
                                 mutation_scale=11, lw=1.2, color=INK2))
    ax.text(0.40, 0.61, 'packaging · transduction\nintegration · selection\n'
                        'different amplicon', fontsize=6.2, color=ENRICHED, va='center')
    ax.text(0.0, 0.17, 'The ratio spans all of this, and the\n'
                       'plasmid sample cannot report on any\n'
                       'of it — so no amount of dividing by\n'
                       'the plasmid removes it.',
            ha='left', fontsize=6.6, color=INK2, va='top')
    ax.set_title('a   what T0 / plasmid contains', loc='left', fontsize=8.5, color=INK)

    # ---- (b) the evidence
    ax = fig.add_subplot(gs[1])
    ax.scatter(targ.GC, targ.lfc_early_c, s=2.5, color=MUTED, alpha=.25, linewidths=0,
               rasterized=True, label='targeting')
    ax.scatter(ctrl.GC, ctrl.lfc_early_c, s=7, color=ENRICHED, alpha=.8, linewidths=0,
               label='cannot affect fitness')
    xs = np.linspace(ctrl.GC.min(), ctrl.GC.max(), 20)
    ax.plot(xs, np.polyval(np.polyfit(ctrl.GC, ctrl.lfc_early_c, 1), xs),
            color=INK, lw=1.6)
    ax.axhline(0, color=INK2, lw=.6)
    r = scipy.stats.pearsonr(ctrl.GC, ctrl.lfc_early_c).statistic
    ax.set_xlabel('protospacer GC', fontsize=7.5)
    ax.set_ylabel('raw score (L1)', fontsize=7.5)
    ax.set_title(f'b   the artifact is visible in guides\n      with no biology  (r = {r:.2f})',
                 loc='left', fontsize=8.5)
    ax.legend(frameon=False, fontsize=5.8, loc='lower right')
    ax.spines[['top', 'right']].set_visible(False)

    # ---- (c) step 1
    ax = fig.add_subplot(gs[2])
    ax.scatter(targ.GC, targ.step1, s=2.5, color=MUTED, alpha=.25, linewidths=0,
               rasterized=True)
    ax.scatter(ctrl.GC, ctrl.step1, s=7, color=DEPLETED, alpha=.8, linewidths=0)
    ax.plot(xs, np.polyval(np.polyfit(ctrl.GC, ctrl.step1, 1), xs), color=INK, lw=1.6)
    ax.axhline(0, color=INK2, lw=.6)
    r2 = scipy.stats.pearsonr(ctrl.GC, ctrl.step1).statistic
    ax.set_xlabel('protospacer GC', fontsize=7.5)
    ax.set_ylabel('after step 1 (L2)', fontsize=7.5)
    ax.set_title(f'c   step 1: subtract a composition model\n'
                 f'      fit only on those guides  (r = {r2:.2f})', loc='left', fontsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- (d) step 2
    ax = fig.add_subplot(gs[3])
    te = targ.dropna(subset=['edit', 'step1'])
    q = np.quantile(te.edit, np.linspace(0, 1, 13)); q[0] -= 1e-9
    lab = np.digitize(te.edit, q[1:-1])
    med = te.groupby(lab).step1.median(); ctr = te.groupby(lab).edit.median()
    ax.scatter(te.edit, te.step1, s=2, color=MUTED, alpha=.18, linewidths=0, rasterized=True)
    ax.plot(ctr, med, '-o', color=INK, lw=1.8, ms=4.5, label='median per bin — subtracted')
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_ylim(-3, 3)
    ax.set_xlabel('% reads edited', fontsize=7.5)
    ax.set_ylabel('after step 1 (L2)', fontsize=7.5)
    ax.set_title('d   step 2: what is left still rises with\n'
                 '      editing; subtract that curve too', loc='left', fontsize=8.5)
    ax.legend(frameon=False, fontsize=5.8, loc='lower right')
    ax.spines[['top', 'right']].set_visible(False)

    # ---- (e) what it buys
    ax = fig.add_subplot(gs[4])
    x = np.arange(len(levels))
    ax.bar(x - .19, res_gc, .36, color=ENRICHED, label='r (residue GC, score)')
    ax.bar(x + .19, frac_depl, .36, color=DEPLETED, label='fraction of calls depleted')
    for i, (a, b) in enumerate(zip(res_gc, frac_depl)):
        ax.text(i - .19, a + .015, f'{a:.2f}', ha='center', fontsize=6, color=INK2)
        ax.text(i + .19, b + .015, f'{100*b:.0f}%', ha='center', fontsize=6, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(levels, fontsize=7.5)
    ax.set_ylim(0, 1.0)
    ax.set_title('e   and what that buys: the artifact leaves\n'
                 '      the residue scores, and calls become depletions',
                 loc='left', fontsize=8.5)
    ax.legend(frameon=False, fontsize=5.8, loc='upper left')
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('How the essentiality score is corrected — every step is learned from '
                 'guides that cannot change fitness',
                 x=0.045, ha='left', y=0.955, fontsize=11, color=INK)
    fig.text(0.045, 0.875,
             f'Panels b–d show {SHOW[0]} {SHOW[1]}, the screen with the largest artifact. '
             'Step 1 is fit on non-targeting and intron guides only. Step 2 is fit on all '
             'targeting guides and is the one correction that cannot be validated on '
             'controls, because controls have no editing read-out.',
             fontsize=6.4, color=INK2, va='top')
    fig.savefig(f'{FIG}/correction_explained.pdf')
    plt.close(fig)
    say(f'residue-level GC correlation by level: ' +
        ', '.join(f'{l} {v:.3f}' for l, v in zip(levels, res_gc)))
    say(f'fraction of significant calls that are depletions: ' +
        ', '.join(f'{l} {100*v:.0f}%' for l, v in zip(levels, frac_depl)))
    say(f'correction_explained.pdf -> {FIG}')


if __name__ == '__main__':
    main()
