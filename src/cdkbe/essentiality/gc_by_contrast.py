#!/usr/bin/env python3
"""
Does the GC bias cancel in the resistance contrast?

`drug / DMSO` takes both samples from the same infection, the same timepoint and the same
amplicon, so the plasmid-to-genome step cancels exactly and the contrast should be much
cleaner than `T0 / plasmid`. This measures how much of the GC dependence actually goes away.

The test that separates artifact from biology is the same one used throughout: non-targeting
and intron guides cannot respond to a drug, so any GC dependence *they* show is technical.

Usage:  python3 gc_by_contrast.py
Writes: screening_data/.../gc_by_contrast.csv and figures/.../gc_by_contrast.pdf
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from . import t0_plasmid_essentiality as E
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
from .essentiality_figures import INK, INK2, MUTED, DEPLETED, ENRICHED, FIG, TAB

DRUGS = {'SY-5609_10nM', 'SY-5609_100nM', 'SEL120_4000nM', 'Senexin B_2000nM',
         'Senexin B_15000nM', 'KB-0742_1500nM', 'KI-CDK9d-32_100nM', 'KI-CDK9d-32_1000nM',
         'KI-CDK9d-32N_1250nM', 'KI-CDK9d-32N_5000nM', 'HQ461', 'BSJ-4-116', 'CDK12-IN-2',
         'PALBO', 'ABEMA', 'RIBO', 'ATIRMO', 'INX', 'TAGTO'}
say = E.say


def main():
    lib = pd.read_csv(f'{LIBRARY}/CDK_library_final.csv')[
        ['gRNA_id', 'Gene', 'protospacer', 'PAM', 'classification', 'Pool']]
    lib['GC'] = E.seq_features(lib).GC.values
    gc = dict(zip(lib.gRNA_id, lib.GC))
    cls = dict(zip(lib.gRNA_id, lib.classification))
    SAV, BURD = E.canonical_editing()
    ed_map = BURD.groupby(['gRNA_id', 'Editor'])['edit'].mean()

    rows = []
    for ed in ['ABE', 'CBE']:
        m = pd.read_csv(f'{LFC_DIR}/{ed}_master_LFC_table.csv',
                        low_memory=False)
        m['GC'] = m.gRNA_id.map(gc)
        m['cls'] = m.gRNA_id.map(cls)
        m['edit'] = [ed_map.get((g, ed), np.nan) for g in m.gRNA_id]
        for base, nr, kind in [('DMSO', 'base_RAW_DMSO', 'drug / DMSO'),
                               ('T0', 'base_RAW_T0', 'drug / T0'),
                               ('Plasmid', 'base_RAW_Plasmid', 'drug / plasmid')]:
            sub = m[m.Base_Value == base]
            if not len(sub):
                continue
            for c in [c for c in sub.columns
                      if c.startswith('LFC_median_') and c.endswith('_' + base)]:
                cond = c[len('LFC_median_'):-len('_' + base)]
                if cond not in DRUGS:
                    continue
                for scr, d in sub.groupby('Screen_ID'):
                    d = d[d[nr] >= E.MIN_BASE].dropna(subset=[c, 'GC'])
                    t = d[d.cls == 'targeting']
                    n = d[d.cls.isin(E.NULL_CLASSES)]
                    if len(t) < 200 or len(n) < 50:
                        continue
                    te = t.dropna(subset=['edit'])
                    pr = np.nan
                    if len(te) > 200:      # GC vs LFC, holding editing level fixed
                        rx = te.GC - np.polyval(np.polyfit(te.edit, te.GC, 1), te.edit)
                        ry = te[c] - np.polyval(np.polyfit(te.edit, te[c], 1), te.edit)
                        pr = scipy.stats.pearsonr(rx, ry).statistic
                    rows.append(dict(kind=kind, editor=ed, screen=scr, drug=cond,
                                     n_targeting=len(t), n_control=len(n),
                                     r_targeting=scipy.stats.pearsonr(t.GC, t[c]).statistic,
                                     r_control=scipy.stats.pearsonr(n.GC, n[c]).statistic,
                                     r_targeting_given_editing=pr))
    r = pd.DataFrame(rows)
    r.to_csv(f'{TAB}/gc_by_contrast.csv', index=False)

    summ = (r.assign(at=r.r_targeting.abs(), ac=r.r_control.abs())
            .groupby('kind').agg(arms=('drug', 'size'),
                                 median_targeting=('at', 'median'),
                                 max_targeting=('at', 'max'),
                                 median_control=('ac', 'median'),
                                 max_control=('ac', 'max')))
    say('\n|r| between protospacer GC and the LFC, true drug arms only:')
    say(summ.round(3).to_string())
    say('\nper screen, drug / DMSO:')
    per = r[r.kind == 'drug / DMSO'].groupby(['editor', 'screen'])[
        ['r_targeting', 'r_control', 'r_targeting_given_editing']].median()
    say(per.round(3).to_string())

    # -------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(10.6, 3.2))
    order = ['drug / plasmid', 'drug / T0', 'drug / DMSO']
    order = [o for o in order if (r.kind == o).any()]

    a = ax[0]
    for i, k in enumerate(order):
        for j, (col, c, lab) in enumerate([('r_targeting', DEPLETED, 'targeting'),
                                           ('r_control', MUTED, 'controls')]):
            v = r[r.kind == k][col].abs()
            a.scatter(np.full(len(v), i + (j - .5) * .34) +
                      np.random.RandomState(0).uniform(-.06, .06, len(v)),
                      v, s=12, color=c, alpha=.6, linewidths=0)
            a.plot([i + (j - .5) * .34 - .12, i + (j - .5) * .34 + .12],
                   [v.median()] * 2, color=INK, lw=2, zorder=4)
    a.set_xticks(range(len(order)))
    a.set_xticklabels([o.replace(' / ', '\n/ ') for o in order], fontsize=7)
    a.set_ylabel('|r| (protospacer GC, LFC)', fontsize=8)
    a.set_title('a   the same-timepoint ratio halves it,\n      '
                'but does not remove it', loc='left', fontsize=8)
    a.legend(handles=[plt.Line2D([], [], marker='o', ls='none', mfc=c, mec='none', ms=5,
                                 label=l) for l, c in [('targeting', DEPLETED),
                                                       ('controls', MUTED)]],
             frameon=False, fontsize=6.5, loc='upper right')
    a.spines[['top', 'right']].set_visible(False)

    a = ax[1]
    d = r[r.kind == 'drug / DMSO']
    a.scatter(d.r_control, d.r_targeting, s=22, color=DEPLETED, alpha=.75, linewidths=0)
    lim = (-.6, .6)
    a.plot(lim, lim, color=MUTED, lw=.8, ls=(0, (2, 2)))
    a.axhline(0, color=INK2, lw=.6); a.axvline(0, color=INK2, lw=.6)
    a.set_xlim(*lim); a.set_ylim(*lim)
    a.set_xlabel('r among guides that cannot\nrespond to drug', fontsize=7.5)
    a.set_ylabel('r among targeting guides', fontsize=7.5)
    a.set_title('b   whatever is left is still technical\n      (controls show it too)',
                loc='left', fontsize=8)
    a.spines[['top', 'right']].set_visible(False)

    a = ax[2]
    per = per.reset_index()
    ys = np.arange(len(per))[::-1]
    a.barh(ys + .18, per.r_targeting, .34, color=DEPLETED, label='targeting')
    a.barh(ys - .18, per.r_control, .34, color=MUTED, label='controls')
    a.axvline(0, color=INK2, lw=.7)
    a.set_yticks(ys)
    a.set_yticklabels([f'{e} {s}' for e, s in zip(per.editor, per.screen)], fontsize=6.2)
    a.set_xlabel('median r, drug / DMSO', fontsize=7.5)
    a.set_title('c   and it varies a lot by screen', loc='left', fontsize=8)
    a.legend(frameon=False, fontsize=6.5, loc='lower right')
    a.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Does the GC bias cancel in the resistance contrast? '
                 'Mostly, but not entirely',
                 x=0.01, ha='left', y=1.02, fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.savefig(f'{FIG}/gc_by_contrast.pdf')
    plt.close(fig)
    say(f'\ngc_by_contrast.pdf -> {FIG}')


if __name__ == '__main__':
    main()
