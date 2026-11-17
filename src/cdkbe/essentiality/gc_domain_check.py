#!/usr/bin/env python3
"""
Is the GC term the composition model removes just "the kinase domain is GC-rich"?

If it were, the L2 correction would be deleting biology: the catalytic domain is where
the essential residues are, so regressing out anything that tracks domain membership
would flatten the real signal.

Four checks say it is not, and that the bias runs the other way:

  a  the catalytic domain is GC-**poor**, not GC-rich (pooled 0.476 vs 0.548)
  b  guides that cannot affect fitness — non-targeting and intron — show the same
     GC slope as targeting guides (median ratio 1.09), so the slope is not biology
  c  the GC-score correlation is just as strong *inside* the domain as outside, which
     it would not be if GC were a proxy for being in the domain
  d  because the domain is GC-poor and low GC looks depleted, the raw contrast
     **over-states** essentiality in the domain; the correction removes that inflation

Usage:  python3 gc_domain_check.py
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
sys.path.insert(0, f'{E.ROOT}/notebooks/11_ChimeraX')
from .._chimerax_config import STRUCTURES, CDK_COLORS
from .essentiality_figures import INK, INK2, MUTED, DEPLETED, ENRICHED, GENES, FIG, TAB

say = E.say


def main():
    lib = pd.read_csv(f'{LIBRARY}/CDK_library_final.csv')[
        ['gRNA_id', 'Gene', 'protospacer', 'PAM', 'classification', 'Pool']]
    lib['GC'] = E.seq_features(lib).GC.values
    gc = dict(zip(lib.gRNA_id, lib.GC))
    K = pd.read_csv(f'{KLIFS_DIR}/KLIFS_compiled.csv')
    pocket = {g: set(K[g].dropna().astype(int)) for g in GENES}
    kd = {g: STRUCTURES[g]['kd'] for g in GENES}

    g = pd.read_csv(f'{TAB}/guide_essentiality.csv', low_memory=False)
    g['GC'] = g.gRNA_id.map(gc)
    rows = pd.read_csv(f'{TAB}/residue_guide_rows.csv').rename(columns={'%Reads': 'pct'})
    rows['GC'] = rows.gRNA_id.map(gc)
    rows['in_domain'] = [kd[x][0] <= c <= kd[x][1] for x, c in zip(rows.Gene, rows.Codon)]
    rows['in_pocket'] = [c in pocket[x] for x, c in zip(rows.Gene, rows.Codon)]
    uniq = rows.drop_duplicates('gRNA_id')

    # ---- a
    per = []
    for x in GENES:
        d = uniq[uniq.Gene == x]
        a, b = d[d.in_domain].GC, d[~d.in_domain].GC
        per.append(dict(Gene=x, GC_in=a.mean(), GC_out=b.mean(), diff=a.mean() - b.mean(),
                        p=scipy.stats.mannwhitneyu(a, b).pvalue, n_in=len(a), n_out=len(b)))
    per = pd.DataFrame(per)

    # ---- b
    sl = []
    for (e, s), d in g.groupby(['Editor', 'Screen_ID']):
        n = d[d.classification.isin(E.NULL_CLASSES)].dropna(subset=['GC', 'lfc_early_c'])
        t = d[d.classification == 'targeting'].dropna(subset=['GC', 'lfc_early_c'])
        sl.append(dict(Editor=e, Screen_ID=s,
                       slope_controls=np.polyfit(n.GC, n.lfc_early_c, 1)[0],
                       slope_targeting=np.polyfit(t.GC, t.lfc_early_c, 1)[0]))
    sl = pd.DataFrame(sl)
    sl['ratio'] = sl.slope_targeting / sl.slope_controls

    # ---- c
    cc = []
    for lab, col in [('raw (L1)', 'lfc_early_c'), ('corrected (L3)', 'ess_early')]:
        d = rows.dropna(subset=['GC', col])
        for where, m in [('inside domain', d.in_domain), ('outside', ~d.in_domain)]:
            cc.append(dict(level=lab, region=where, n=int(m.sum()),
                           r=scipy.stats.pearsonr(d[m].GC, d[m][col]).statistic))
    cc = pd.DataFrame(cc)

    # ---- d
    rows['shift'] = rows.ess_early - rows.lfc_early_c
    sh = [dict(group=lab, mean_shift=rows.loc[m, 'shift'].mean(), n=int(m.sum()))
          for lab, m in [('catalytic\ndomain', rows.in_domain),
                         ('outside', ~rows.in_domain),
                         ('ATP\npocket', rows.in_pocket),
                         ('non-pocket', ~rows.in_pocket)]]
    sh = pd.DataFrame(sh)

    pd.concat([per.assign(check='GC by domain'), sl.assign(check='GC slope'),
               cc.assign(check='GC correlation'), sh.assign(check='correction shift')],
              ignore_index=True).to_csv(f'{TAB}/qc_gc_domain_check.csv', index=False)
    say(per.round(3).to_string(index=False))
    say(f"\npooled GC: inside {uniq[uniq.in_domain].GC.mean():.3f} vs "
        f"outside {uniq[~uniq.in_domain].GC.mean():.3f}")
    say(f"median targeting/control GC slope ratio: {sl.ratio.median():.2f}")

    # --------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 4, figsize=(12.2, 3.1))

    a = ax[0]
    ys = np.arange(len(per))[::-1]
    a.barh(ys, per['diff'], color=[CDK_COLORS[x] for x in per.Gene], height=.62)
    for y, r in zip(ys, per.itertuples()):
        if r.p < 0.05:
            a.text(r.diff - .004 * np.sign(r.diff) - .004, y, '*', ha='right',
                   va='center', fontsize=9, color=INK)
    a.axvline(0, color=INK2, lw=.7)
    a.set_yticks(ys); a.set_yticklabels(per.Gene, fontsize=6.5)
    for lbl, x in zip(a.get_yticklabels(), per.Gene):
        lbl.set_color(CDK_COLORS[x])
    a.set_xlabel('protospacer GC,\ncatalytic domain − rest', fontsize=7.5)
    a.set_title('a   the domain is GC-POOR,\n      not GC-rich', loc='left', fontsize=8)
    a.spines[['top', 'right']].set_visible(False)

    a = ax[1]
    lim = (0, max(sl.slope_controls.max(), sl.slope_targeting.max()) * 1.12)
    a.plot(lim, lim, color=MUTED, lw=.8, ls=(0, (2, 2)))
    for r in sl.itertuples():
        a.scatter(r.slope_controls, r.slope_targeting, s=38,
                  facecolor=DEPLETED if r.Editor == 'ABE' else ENRICHED,
                  edgecolor='white', linewidths=.6, zorder=3)
    a.set_xlim(*lim); a.set_ylim(*lim)
    a.set_xlabel('slope among guides that\ncannot affect fitness', fontsize=7.5)
    a.set_ylabel('slope among\ntargeting guides', fontsize=7.5)
    a.set_title(f'b   same slope either way\n      (median ratio {sl.ratio.median():.2f})',
                loc='left', fontsize=8)
    a.legend(handles=[plt.Line2D([], [], marker='o', ls='none', mfc=c, mec='none', ms=5,
                                 label=e) for e, c in [('ABE', DEPLETED), ('CBE', ENRICHED)]],
             frameon=False, fontsize=6, loc='lower right')
    a.spines[['top', 'right']].set_visible(False)

    a = ax[2]
    w = .36
    for i, lev in enumerate(['raw (L1)', 'corrected (L3)']):
        for j, where in enumerate(['inside domain', 'outside']):
            v = cc[(cc.level == lev) & (cc.region == where)].r.iloc[0]
            a.bar(i + (j - .5) * w, v, w * .9,
                  color=DEPLETED if j == 0 else MUTED)
            a.text(i + (j - .5) * w, v + .012, f'{v:.2f}', ha='center', fontsize=6,
                   color=INK2)
    a.set_xticks([0, 1]); a.set_xticklabels(['raw (L1)', 'corrected (L3)'], fontsize=7)
    a.set_ylabel('r (GC, score)', fontsize=7.5)
    a.set_title('c   GC acts the same inside the\n      domain as outside it',
                loc='left', fontsize=8)
    a.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=c, ec='none', label=l)
                      for l, c in [('inside domain', DEPLETED), ('outside', MUTED)]],
             frameon=False, fontsize=6, loc='upper right')
    a.spines[['top', 'right']].set_visible(False)

    a = ax[3]
    a.bar(range(len(sh)), sh.mean_shift,
          color=[DEPLETED, MUTED, DEPLETED, MUTED], width=.62)
    for i, r in enumerate(sh.itertuples()):
        a.text(i, r.mean_shift + .012, f'{r.mean_shift:+.2f}', ha='center', va='bottom',
               fontsize=6, color='white' if r.mean_shift < -.2 else INK2)
    a.axhline(0, color=INK2, lw=.7)
    a.set_xticks(range(len(sh))); a.set_xticklabels(sh.group, fontsize=6.2)
    a.set_ylabel('mean change applied\nby the correction', fontsize=7.5)
    a.set_title('d   the raw score over-stated the\n      domain; correction removes it',
                loc='left', fontsize=8)
    a.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Is the GC term just "the kinase domain is GC-rich"?  No — '
                 'the domain is GC-poor, and controls show the same slope',
                 x=0.012, ha='left', y=1.02, fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f'{FIG}/gc_domain_check.pdf')
    plt.close(fig)
    say(f'\ngc_domain_check.pdf -> {FIG}')


if __name__ == '__main__':
    main()
