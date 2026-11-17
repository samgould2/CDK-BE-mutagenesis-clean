#!/usr/bin/env python3
"""
What every correction actually does — one readout table, four levels of processing.

The point of this script is auditability. The essentiality score is not the raw
`LFC_median_T0_Plasmid`; three things happen to it first. Each is listed here with
what it buys and what it costs, measured on the same fixed set of readouts, so the
stopping point is a choice made with numbers rather than an assumption.

    L0  raw            `LFC_median_T0_Plasmid` from the master table, untouched
    L1  + centred      minus the median of the non-targeting guides of the same
                       (Editor, Screen). One number per screen. Removes the pool's
                       global dropout offset; nothing else.
    L2  + composition  minus a ridge model of protospacer/PAM composition, fit ONLY on
                       non-targeting and intron guides (out-of-fold), then rescaled per
                       (Editor, Screen) with two parameters fit on those same guides.
                       No targeting guide is ever in the training set.
    L3  + editing      minus the median score of targeting guides in the same editing-level
                       bin (12 quantile bins, per Editor x Screen). This is the aggressive
                       one: it is fit on ALL targeting guides, so any monotone relationship
                       between editing level and score is removed, real signal included,
                       and the reference point becomes "the median targeting guide at the
                       same editing level" rather than "a non-targeting control".

Residue scores are then an editing-weighted mean over guide x residue rows; the
unweighted mean is reported alongside, since weighting is itself a choice.

None of the readouts below was used to pick any of the corrections.

Usage:  python3 correction_ladder.py
Writes: screening_data/.../correction_ladder.csv and
        figures/T0_plasmid_essentiality/correction_ladder.pdf
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from . import t0_plasmid_essentiality as E
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
from .essentiality_figures import (INK, INK2, MUTED, DEPLETED, ENRICHED, GENES, FIG)

TAB = E.TABLE_DIR
LEVELS = [('L0  raw', 'lfc_early'),
          ('L1  + non-targeting centred', 'lfc_early_c'),
          ('L2  + composition corrected', 'step1'),
          ('L3  + editing detrended', 'ess_early')]
say = E.say


def residue_scores(SAV, f, col, weighted=True):
    v = SAV[(SAV['%Reads'] >= E.MIN_EDIT) & (SAV.Class == 'Missense')].merge(
        f[['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'classification', col]],
        on=['gRNA_id', 'Editor', 'Screen_ID'], how='inner')
    v = v[(v.classification == 'targeting') & v[col].notna()]
    v = v[v.Screen_ID.isin(E.PRIMARY)]
    v['w'] = v['%Reads'] / 100.0 if weighted else 1.0
    agg = (v.groupby(['Gene', 'Codon'])
           .apply(lambda d: np.average(d[col], weights=d.w)).rename('score').reset_index())
    return agg, v


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

    # the independent infection, scored the same way at every level
    hold = f[f.Screen_ID == E.HOLDOUT]
    gc_map = dict(zip(lib.gRNA_id, X.GC.values))

    rows = []
    res_by_level = {}
    for label, col in LEVELS:
        t = f[(f.classification == 'targeting') & f[col].notna()]
        nt = f[(f.classification == 'non-targeting control') & f[col].notna()]

        r_gc = [scipy.stats.pearsonr(d.GC, d[col]).statistic
                for _, d in t.groupby(['Editor', 'Screen_ID']) if len(d) > 50]
        r_ed = [scipy.stats.spearmanr(d.edit, d[col]).statistic
                for _, d in t.groupby(['Editor', 'Screen_ID'])
                if d.edit.notna().sum() > 50]

        base = nt[col].median()
        tr = t[t.Nonsense_perc >= E.EDIT_HIGH]
        tr9 = tr[tr.Gene == 'CDK9']
        n_below = sum(1 for g in GENES
                      if len(tr[tr.Gene == g]) >= 5
                      and tr[tr.Gene == g][col].median() < base)

        # L3 detrends targeting guides only — the controls have no editing measurement
        # anywhere in the dataset, so a truncating-vs-non-targeting comparison straddles
        # two reference frames at that level. This one stays inside a single frame:
        # truncating guides against the other targeting guides that edit at least as well.
        hi = t[t.edit >= E.EDIT_HIGH]
        is_tr = (hi.Nonsense_perc >= E.EDIT_HIGH).fillna(False)
        a, b = hi[is_tr][col].dropna(), hi[~is_tr][col].dropna()
        tr_vs_targ = a.median() - b.median() if len(a) >= 20 else np.nan
        tr_vs_targ_p = scipy.stats.mannwhitneyu(a, b).pvalue if len(a) >= 20 else np.nan
        n_below_targ = 0
        for g in GENES:
            sa = hi[is_tr & (hi.Gene == g)][col].dropna()
            sb = hi[~is_tr & (hi.Gene == g)][col].dropna()
            if len(sa) >= 5 and len(sb) >= 20 and sa.median() < sb.median():
                n_below_targ += 1

        agg, v = residue_scores(SAV, f, col)
        v = v.assign(gc=v.gRNA_id.map(gc_map))
        agg_u, _ = residue_scores(SAV, f, col, weighted=False)
        res_by_level[label] = agg
        agg['in_pocket'] = [c in pocket[g] for g, c in zip(agg.Gene, agg.Codon)]

        # cross-editor at residue level
        w = (v.groupby(['Gene', 'Codon', 'Editor'])
             .apply(lambda d: np.average(d[col], weights=d.w)).unstack().dropna())
        r_ed_res = scipy.stats.pearsonr(w.ABE, w.CBE).statistic if len(w) > 30 else np.nan

        # independent infection
        hv = SAV[(SAV['%Reads'] >= E.MIN_EDIT) & (SAV.Class == 'Missense')].merge(
            hold[['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'classification', col]],
            on=['gRNA_id', 'Editor', 'Screen_ID'], how='inner')
        hv = hv[(hv.classification == 'targeting') & hv[col].notna()]
        hv['w'] = hv['%Reads'] / 100.0
        ha = (hv.groupby(['Gene', 'Codon'])
              .apply(lambda d: np.average(d[col], weights=d.w)).rename('hold').reset_index())
        j = agg.merge(ha, on=['Gene', 'Codon'])
        r_hold = scipy.stats.pearsonr(j.score, j.hold).statistic if len(j) > 30 else np.nan

        # the readout that decides this: does the composition artifact survive into the
        # residue scores? Guide-level GC correlations can look fine while the weighted
        # mean per residue still carries the term.
        rg = (v.groupby(['Gene', 'Codon'])
              .apply(lambda d: np.average(d.gc, weights=d.w)).rename('res_GC').reset_index())
        mg = agg.merge(rg, on=['Gene', 'Codon'])
        r_res_gc = scipy.stats.pearsonr(mg.res_GC, mg.score).statistic
        max_gene_gc = max(abs(scipy.stats.pearsonr(d.res_GC, d.score).statistic)
                          for _, d in mg.groupby('Gene') if len(d) > 30)

        pk = agg[agg.in_pocket].score.median() - agg[~agg.in_pocket].score.median()
        n_pocket_sig = 0
        for g in GENES:
            s = agg[agg.Gene == g]
            a, b = s[s.in_pocket].score, s[~s.in_pocket].score
            if len(a) >= 5 and len(b) >= 5 and scipy.stats.mannwhitneyu(a, b).pvalue < 0.05 \
                    and a.median() < b.median():
                n_pocket_sig += 1

        rows.append(dict(
            level=label,
            r_GC_median=np.median(r_gc), r_GC_max=np.max(np.abs(r_gc)),
            rho_editing_median=np.median(r_ed),
            truncating_minus_NT=tr[col].median() - base,
            truncating_minus_other_targeting=tr_vs_targ,
            p_truncating_vs_targeting=tr_vs_targ_p,
            n_CDKs_truncating_below_targeting=n_below_targ,
            CDK9_truncating_minus_NT=tr9[col].median() - base if len(tr9) >= 5 else np.nan,
            n_CDKs_truncating_below_NT=n_below,
            r_ABE_vs_CBE=r_ed_res,
            r_independent_infection=r_hold,
            r_residue_GC_score=r_res_gc, max_per_gene_residue_GC=max_gene_gc,
            pocket_minus_rest=pk,
            n_CDKs_pocket_depleted=n_pocket_sig,
            sd_residue_score=agg.score.std(),
            r_weighted_vs_unweighted=scipy.stats.pearsonr(
                agg.score, agg_u.score).statistic))

    tab = pd.DataFrame(rows)

    # how much does each step move the residue scores?
    labs = [l for l, _ in LEVELS]
    moved = [np.nan] + [scipy.stats.spearmanr(
        res_by_level[labs[i - 1]].merge(res_by_level[labs[i]], on=['Gene', 'Codon'],
                                        suffixes=('_a', '_b')).score_a,
        res_by_level[labs[i - 1]].merge(res_by_level[labs[i]], on=['Gene', 'Codon'],
                                        suffixes=('_a', '_b')).score_b).statistic
        for i in range(1, len(labs))]
    tab['rho_vs_previous_level'] = moved

    tab.to_csv(f'{TAB}/correction_ladder.csv', index=False)
    pd.set_option('display.width', 250)
    say('\n' + '=' * 100)
    say('CORRECTION LADDER — every readout, at every level of processing')
    say('=' * 100)
    say(tab.round(3).to_string(index=False))
    fig_ladder(tab, res_by_level)
    return tab, res_by_level


def fig_ladder(tab, res_by_level):
    labs = list(tab.level)
    short = [l.split()[0] for l in labs]
    panels = [
        ('r_GC_max', 'worst |r| with protospacer GC', 'lower is better', False),
        ('rho_editing_median', 'median ρ with editing level', 'lower is better', False),
        ('truncating_minus_other_targeting',
         'truncating − other targeting guides\n(same reference frame)', 'lower is better', True),
        ('CDK9_truncating_minus_NT', 'CDK9 truncating − non-targeting', 'lower is better', True),
        ('r_independent_infection', 'r, independent infection', 'higher is better', False),
        ('r_ABE_vs_CBE', 'r, ABE vs CBE', 'higher is better', False),
        ('r_residue_GC_score', 'r (residue GC, score)\n— the artifact at residue level',
         'lower is better', True),
        ('n_CDKs_pocket_depleted', 'CDKs with pocket depletion (p<0.05)',
         'higher is better', False),
    ]
    fig, ax = plt.subplots(2, 4, figsize=(11.6, 5.2))
    for a, (col, title, hint, zero) in zip(ax.ravel(), panels):
        v = tab[col].values.astype(float)
        good = ('lower' in hint)
        cols = [DEPLETED if (vv == (np.nanmin(v) if good else np.nanmax(v))) else MUTED
                for vv in v]
        a.bar(range(len(v)), v, color=cols, width=.66)
        if zero:
            a.axhline(0, color=INK2, lw=.7)
        a.set_xticks(range(len(v)))
        a.set_xticklabels(short, fontsize=6.5)
        a.set_title(f'{title}\n({hint})', loc='left', fontsize=7.2, color=INK)
        a.tick_params(labelsize=6.5)
        a.spines[['top', 'right']].set_visible(False)
        for i, vv in enumerate(v):
            if np.isfinite(vv):
                a.text(i, vv, f'{vv:.2f}' if abs(vv) < 10 else f'{vv:.0f}',
                       ha='center', va='bottom' if vv >= 0 else 'top',
                       fontsize=5.8, color=INK2)
    fig.suptitle('What each correction buys, and what it costs  '
                 '(blue = best level for that readout)',
                 x=0.03, ha='left', y=0.985, fontsize=10.5, color=INK)
    fig.text(0.03, 0.005,
             'L0 raw LFC · L1 centred on non-targeting guides of the same screen · '
             'L2 minus a composition model fit only on non-targeting + intron guides · '
             'L3 minus the editing-level trend of all targeting guides',
             fontsize=6.2, color=INK2)
    fig.tight_layout(rect=[0, 0.025, 1, 0.955])
    fig.savefig(f'{FIG}/correction_ladder.pdf')
    plt.close(fig)
    say(f'\ncorrection_ladder.pdf written')


if __name__ == '__main__':
    main()
