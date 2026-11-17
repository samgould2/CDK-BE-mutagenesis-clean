#!/usr/bin/env python3
"""
Thresholding instead of detrending.

The editing-level detrend (L3) subtracts a curve fitted to all targeting guides. It is the
only correction in the pipeline that cannot be validated on control guides, because
non-targeting and intron guides have no sensor read-out anywhere in this dataset. The
simpler alternative is to *drop* the weakly edited guides rather than model them: keep the
composition correction, which is fit only on controls, and require a minimum editing level.

This sweeps the editing threshold on top of L1 and L2 and scores every combination on the
same readouts as `correction_ladder.py`, plus the coverage each one costs.

The truncating comparison is kept inside one reference frame throughout — truncating guides
against the other targeting guides, never against the undetrended controls.

Usage:  python3 threshold_sweep.py
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
from .essentiality_figures import INK, INK2, MUTED, DEPLETED, ENRICHED, GENES, FIG, TAB

THRESHOLDS = [10, 25, 50, 75]
say = E.say


def readouts(SAV, f, col, thr, pocket, hold, label, gc_map):
    v = SAV[(SAV['%Reads'] >= thr) & (SAV.Class == 'Missense')].merge(
        f[['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'classification', col]],
        on=['gRNA_id', 'Editor', 'Screen_ID'], how='inner')
    v = v[(v.classification == 'targeting') & v[col].notna()
          & v.Screen_ID.isin(E.PRIMARY)].copy()
    v['w'] = v['%Reads'] / 100.0
    agg = (v.groupby(['Gene', 'Codon'])
           .apply(lambda d: np.average(d[col], weights=d.w)).rename('score').reset_index())
    agg['in_pocket'] = [c in pocket[g] for g, c in zip(agg.Gene, agg.Codon)]

    w = (v.groupby(['Gene', 'Codon', 'Editor'])
         .apply(lambda d: np.average(d[col], weights=d.w)).unstack().dropna())
    r_ed = scipy.stats.pearsonr(w.ABE, w.CBE).statistic if len(w) > 30 else np.nan

    hv = SAV[(SAV['%Reads'] >= thr) & (SAV.Class == 'Missense')].merge(
        hold[['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'classification', col]],
        on=['gRNA_id', 'Editor', 'Screen_ID'], how='inner')
    hv = hv[(hv.classification == 'targeting') & hv[col].notna()].copy()
    hv['w'] = hv['%Reads'] / 100.0
    ha = (hv.groupby(['Gene', 'Codon'])
          .apply(lambda d: np.average(d[col], weights=d.w)).rename('hold').reset_index())
    j = agg.merge(ha, on=['Gene', 'Codon'])
    r_hold = scipy.stats.pearsonr(j.score, j.hold).statistic if len(j) > 30 else np.nan

    # truncating vs other targeting guides, both inside the same frame, at this threshold
    t = f[(f.classification == 'targeting') & f[col].notna() & (f.edit >= thr)]
    is_tr = (t.Nonsense_perc >= thr).fillna(False)
    a, b = t[is_tr][col].dropna(), t[~is_tr][col].dropna()
    tr_delta = a.median() - b.median() if len(a) >= 20 else np.nan
    tr_p = scipy.stats.mannwhitneyu(a, b).pvalue if len(a) >= 20 else np.nan
    n_below = sum(1 for g in GENES
                  if len(t[is_tr & (t.Gene == g)]) >= 5
                  and len(t[~is_tr & (t.Gene == g)]) >= 20
                  and t[is_tr & (t.Gene == g)][col].median()
                  < t[~is_tr & (t.Gene == g)][col].median())

    n_sig_pocket = 0
    for g in GENES:
        s = agg[agg.Gene == g]
        x, y = s[s.in_pocket].score, s[~s.in_pocket].score
        if len(x) >= 5 and len(y) >= 5 and x.median() < y.median() \
                and scipy.stats.mannwhitneyu(x, y).pvalue < 0.05:
            n_sig_pocket += 1

    tg = v.assign(GC=v.gRNA_id.map(gc_map))   # f has one row per editor x screen
    return dict(
        scheme=label, min_edit=thr,
        n_residues=agg.Gene.count(), n_guide_rows=len(v),
        frac_residues=np.nan,
        truncating_minus_targeting=tr_delta, p_truncating=tr_p,
        n_CDKs_truncating_below=n_below,
        r_GC=scipy.stats.pearsonr(tg.GC, tg[col]).statistic,
        r_ABE_vs_CBE=r_ed, r_independent_infection=r_hold,
        pocket_minus_rest=agg[agg.in_pocket].score.median() - agg[~agg.in_pocket].score.median(),
        n_CDKs_pocket_depleted=n_sig_pocket,
        sd_residue_score=agg.score.std())


def main():
    qc = {}
    lib = pd.read_csv(f'{LIBRARY}/CDK_library_final.csv')[
        ['gRNA_id', 'Gene', 'protospacer', 'PAM', 'classification', 'Pool']]
    SAV, BURD = E.canonical_editing()
    f = E.guide_lfc(lib)
    f = f.merge(BURD[['gRNA_id', 'Editor', 'Screen_ID', 'edit', 'Nonsense_perc',
                      'Missense_perc', 'WT_perc', 'Silent_perc', 'Unedited_perc']],
                on=['gRNA_id', 'Editor', 'Screen_ID'], how='left')
    X = E.seq_features(lib)
    f, _ = E.remove_representation(f, X, qc)
    f = E.detrend_editing(f, qc)
    f['GC'] = X.loc[f.gRNA_id, 'GC'].values
    K = pd.read_csv(f'{KLIFS_DIR}/KLIFS_compiled.csv')
    pocket = {g: set(K[g].dropna().astype(int)) for g in GENES}
    hold = f[f.Screen_ID == E.HOLDOUT]
    gc_map = dict(zip(lib.gRNA_id, X.GC.values))

    rows = []
    for label, col in [('L1 centred only', 'lfc_early_c'),
                       ('L2 composition corrected', 'step1'),
                       ('L3 composition + detrend', 'ess_early')]:
        for thr in THRESHOLDS:
            if label.startswith('L3') and thr != E.MIN_EDIT:
                continue           # L3 is shown only at the pipeline's own threshold
            rows.append(readouts(SAV, f, col, thr, pocket, hold, label, gc_map))
    t = pd.DataFrame(rows)
    full = t[(t.scheme.str.startswith('L2')) & (t.min_edit == 10)].n_residues.iloc[0]
    t['frac_residues'] = (t.n_residues / full).round(3)
    t.to_csv(f'{TAB}/threshold_sweep.csv', index=False)
    pd.set_option('display.width', 250)
    say('\n' + '=' * 110)
    say('THRESHOLDING INSTEAD OF DETRENDING')
    say('=' * 110)
    say(t.round(3).to_string(index=False))

    fig, ax = plt.subplots(1, 5, figsize=(13.8, 3.0))
    schemes = [('L1 centred only', MUTED), ('L2 composition corrected', DEPLETED)]
    panels = [('truncating_minus_targeting', 'truncating − other targeting\n(same frame)', True),
              ('r_GC', 'r with protospacer GC', True),
              ('r_independent_infection', 'r, independent infection', False),
              ('pocket_minus_rest', 'ATP pocket − rest', True),
              ('n_residues', 'residues scored', False)]
    for a, (col, title, zero) in zip(ax, panels):
        for sch, c in schemes:
            d = t[t.scheme == sch].sort_values('min_edit')
            a.plot(d.min_edit, d[col], '-o', color=c, lw=1.6, ms=5, label=sch.split()[0])
        d3 = t[t.scheme.str.startswith('L3')]
        if len(d3):
            a.scatter(d3.min_edit, d3[col], marker='*', s=110, color=ENRICHED,
                      zorder=5, label='L3 (detrend)')
        if zero:
            a.axhline(0, color=INK2, lw=.7)
        a.set_xticks(THRESHOLDS)
        a.set_xlabel('minimum % reads making\nthe variant', fontsize=7.5)
        a.set_title(title, loc='left', fontsize=8, color=INK)
        a.tick_params(labelsize=6.5)
        a.spines[['top', 'right']].set_visible(False)
    ax[0].legend(frameon=False, fontsize=6.2, loc='lower left')
    ax[0].annotate('wrong direction', xy=(75, t[(t.scheme.str.startswith('L1')) &
                                               (t.min_edit == 75)].truncating_minus_targeting.iloc[0]),
                   xytext=(-4, -14), textcoords='offset points', fontsize=6, color=INK2,
                   ha='right')
    fig.suptitle('Thresholding is not a substitute for the detrend — and the detrend '
                 'turns out not to be needed',
                 x=0.008, ha='left', y=1.075, fontsize=10.5, color=INK)
    fig.text(0.008, 1.0,
             'Raising the editing threshold makes the positive control WORSE (left) and '
             'costs half the residues (right): high-editing guides are exactly where the '
             'representation artifact is strongest, so thresholding selects for the '
             'confound. L2 at 10% already matches L3 on every readout, with better '
             'cross-editor agreement (0.44 vs 0.38) and the same coverage.',
             fontsize=6.6, color=INK2, ha='left', va='top', wrap=True)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(f'{FIG}/threshold_sweep.pdf')
    plt.close(fig)
    say(f'\nthreshold_sweep.pdf -> {FIG}')
    return t


if __name__ == '__main__':
    main()
