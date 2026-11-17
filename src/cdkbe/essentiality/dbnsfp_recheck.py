#!/usr/bin/env python3
"""
dbNSFP variant-effect predictors against the corrected early score.

`dbNSFP_canonical.ipynb` correlated 46 predictors against the raw `LFC_median_T0_Plasmid`
and against `DMSO vs T0`, and found the raw early readout partly artifactual: the sign came
out wrong in CDK2 and CDK6 and the CDK13 correlation was roughly doubled. This re-runs the
same comparison against the corrected score at every processing level, so the question
"does correcting the early contrast fix its agreement with predictors?" gets an answer
rather than an assumption.

A predictor rankscore is 0-1 with higher = more damaging, so against an essentiality score
where more negative = more essential, the expected correlation is **negative**.

The correction makes the agreement *weaker*, which looks bad until the obvious confound is
checked: predictor scores themselves track protospacer GC, because damaging substitutions
sit at particular residues whose codons have characteristic composition. Panel d partials
that out — the raw agreement collapses to exactly where the corrected score already sits,
so the raw agreement was two quantities sharing one artifact.

Usage:  python3 dbnsfp_recheck.py
Writes: screening_data/.../dbnsfp_recheck.csv, figures/.../dbnsfp_recheck.pdf
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
from .essentiality_figures import INK, INK2, MUTED, DEPLETED, ENRICHED, GENES, FIG, TAB
from .._chimerax_config import CDK_COLORS

UNIPROT = {'CDK7': 'P50613', 'CDK8': 'P49336', 'CDK19': 'Q9BWU1', 'CDK9': 'P50750',
           'CDK12': 'Q9NYV4', 'CDK13': 'Q14004', 'CDK2': 'P24941', 'CDK4': 'P11802',
           'CDK6': 'Q00534'}
LEVELS = [('L1 centred only', '_L1'), ('L2 composition', '_L2'), ('L3 default', '')]
say = E.say


def load_dbnsfp(gene, score_cols):
    """Rows for the canonical accession, with the parallel ';' lists exploded together.

    Taking `.split(';')[0]` instead picks whichever transcript happens to be first and
    silently mis-numbers CDK7, CDK9 and CDK13.
    """
    path = f'{DBNSFP}/{gene}_commandline.out'
    base = ['Uniprot_acc', 'aapos', 'aaref', 'aaalt']
    want = set(base + score_cols)
    df = pd.read_csv(path, sep='\t', low_memory=False, usecols=lambda c: c in want)
    n = df['Uniprot_acc'].astype(str).str.count(';')
    par = [c for c in df.columns if df[c].astype(str).str.count(';').equals(n)]
    for c in par:
        df[c] = df[c].astype(str).str.split(';')
    df = df.explode(par)
    df = df[df['Uniprot_acc'] == UNIPROT[gene]].copy()
    df['aapos'] = pd.to_numeric(df['aapos'], errors='coerce')
    for c in df.columns:
        if c not in par and df[c].dtype == object:
            df[c] = df[c].astype(str).str.split(';').str[0]
    for c in score_cols:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df.dropna(subset=['aapos'])


def main():
    head = pd.read_csv(f'{DBNSFP}/CDK7_commandline.out',
                       sep='\t', nrows=0)
    SCORES = [c for c in head.columns if c.endswith('_rankscore')]
    say(f'{len(SCORES)} predictors')

    # per-variant measured score at each level, plus the uncorrected early LFC
    var = {}
    for lab, suf in LEVELS:
        p = f'{LFC_DIR}/T0_plasmid_essentiality{suf}/variant_essentiality.csv'
        if os.path.exists(p):
            var[lab] = pd.read_csv(p)[['Gene', 'SAV', 'Codon', 'ess_early', 'ess_late',
                                       'n_guides', 'max_edit']]
    rows_p = f'{LFC_DIR}/T0_plasmid_essentiality/residue_guide_rows.csv'
    raw = pd.read_csv(rows_p).rename(columns={'%Reads': 'pct'})
    lib = pd.read_csv(f'{LIBRARY}/CDK_library_final.csv')[
        ['gRNA_id', 'protospacer', 'PAM']]
    lib['GC'] = E.seq_features(lib).GC.values
    raw['GC'] = raw.gRNA_id.map(dict(zip(lib.gRNA_id, lib.GC)))
    raw_var = (raw.groupby(['Gene', 'SAV', 'Codon'])
               .apply(lambda d: pd.Series({
                   'ess_raw': np.average(d.lfc_early_c, weights=d.pct),
                   'var_GC': np.average(d.GC, weights=d.pct)})).reset_index())

    out = []
    for gene in GENES:
        db = load_dbnsfp(gene, SCORES)
        db = db.rename(columns={'aapos': 'Codon'})
        db['SAV'] = db.aaref.astype(str) + db.Codon.astype(int).astype(str) + db.aaalt.astype(str)
        db = db.drop_duplicates('SAV').set_index('SAV')

        series = {}
        for lab in var:
            d = var[lab]
            d = d[d.Gene == gene]
            series[lab] = d.set_index('SAV').ess_early
            if lab == 'L3 default':
                series['late (DMSO/T0)'] = d.set_index('SAV').ess_late
        r = raw_var[raw_var.Gene == gene]
        series['L0 raw early LFC'] = r.set_index('SAV').ess_raw
        gcv = r.set_index('SAV').var_GC

        for readout, s in series.items():
            common = db.index.intersection(s.dropna().index)
            if len(common) < 40:
                continue
            for sc in SCORES:
                x = db.loc[common, sc]
                ok = x.notna()
                if ok.sum() < 40:
                    continue
                y = s.loc[common][ok.values]
                rho = scipy.stats.spearmanr(x[ok], y)
                pr = np.nan
                if readout == 'L0 raw early LFC':
                    z = gcv.loc[common][ok.values]
                    rx, ry, rz = [scipy.stats.rankdata(v) for v in (x[ok], y, z)]
                    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
                    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
                    pr = scipy.stats.pearsonr(ex, ey).statistic
                out.append(dict(Gene=gene, readout=readout, predictor=sc,
                                n=int(ok.sum()), rho=rho.statistic, p=rho.pvalue,
                                rho_given_GC=pr,
                                rho_predictor_vs_GC=(
                                    scipy.stats.spearmanr(x[ok], gcv.loc[common][ok.values]).statistic
                                    if readout == 'L0 raw early LFC' else np.nan)))
        say(f'  {gene}: {len(db):,} dbNSFP variants, '
            f'{len(db.index.intersection(series["L3 default"].dropna().index)):,} measured')

    t = pd.DataFrame(out)
    t.to_csv(f'{TAB}/dbnsfp_recheck.csv', index=False)

    ORDER = ['L0 raw early LFC', 'L1 centred only', 'L2 composition', 'L3 default',
             'late (DMSO/T0)']
    ORDER = [o for o in ORDER if (t.readout == o).any()]
    summ = (t.groupby('readout')
            .agg(median_rho=('rho', 'median'),
                 frac_negative=('rho', lambda s: (s < 0).mean()),
                 frac_sig_neg=('rho', lambda s: ((s < 0)).mean()),
                 n=('rho', 'size')).reindex(ORDER))
    say('\nSpearman rho (predictor rankscore vs measured score); expected NEGATIVE:')
    say(summ.round(3).to_string())

    per = t.pivot_table(index='Gene', columns='readout', values='rho',
                        aggfunc='median').reindex(GENES)[ORDER]
    say('\nmedian rho over the 46 predictors, per CDK:')
    say(per.round(3).to_string())
    say('\nCDKs with the WRONG (positive) sign:')
    say((per > 0).sum().to_string())

    # ------------------------------------------------------------- figure
    fig = plt.figure(figsize=(14.2, 4.0))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.1, 1.5, .8, 1.1], wspace=0.34,
                          left=0.05, right=0.99, top=0.80, bottom=0.17)

    ax = fig.add_subplot(gs[0])
    data = [t[t.readout == o].rho.dropna().values for o in ORDER]
    bp = ax.boxplot(data, positions=list(range(len(ORDER))), widths=.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color='white', lw=1.2),
                    whiskerprops=dict(color=INK2, lw=.7), capprops=dict(color=INK2, lw=.7))
    cols = [MUTED, MUTED, '#8FAADC', DEPLETED, ENRICHED][:len(ORDER)]
    for b, c in zip(bp['boxes'], cols):
        b.set_facecolor(c); b.set_edgecolor('none')
    ax.axhline(0, color=INK2, lw=.8)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([o.replace(' ', '\n', 1) for o in ORDER], fontsize=6)
    lo = min(np.percentile(v, 2) for v in data)
    hi = max(np.percentile(v, 98) for v in data)
    pad = 0.12 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_ylabel('Spearman ρ, predictor vs score', fontsize=8)
    ax.set_title('a   46 predictors × 9 CDKs\n      (negative = agreement)',
                 loc='left', fontsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)

    ax = fig.add_subplot(gs[1])
    x = np.arange(len(GENES))
    w = 0.8 / len(ORDER)
    for i, o in enumerate(ORDER):
        ax.bar(x + (i - (len(ORDER) - 1) / 2) * w, per[o].values, w,
               color=cols[i], label=o)
    ax.axhline(0, color=INK2, lw=.8)
    ax.set_xticks(x); ax.set_xticklabels(GENES, fontsize=6.5, rotation=45)
    ax.set_ylabel('median ρ over predictors', fontsize=8)
    ax.set_title('b   per CDK — the raw early readout gets the sign wrong',
                 loc='left', fontsize=8.5)
    ax.legend(frameon=False, fontsize=5.8, ncol=2)
    ax.spines[['top', 'right']].set_visible(False)

    ax = fig.add_subplot(gs[2])
    frac = [(t[t.readout == o].rho < 0).mean() for o in ORDER]
    ax.bar(range(len(ORDER)), frac, color=cols, width=.6)
    for i, v in enumerate(frac):
        ax.text(i, v + .015, f'{100*v:.0f}%', ha='center', fontsize=6.5, color=INK2)
    ax.axhline(0.5, color=MUTED, lw=.7, ls=(0, (2, 2)))
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([o.split()[0] for o in ORDER], fontsize=6.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('fraction of gene × predictor\npairs with the expected sign', fontsize=7.5)
    ax.set_title('c   sign agreement', loc='left', fontsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- (d) the confound
    ax = fig.add_subplot(gs[3])
    g0 = t[t.readout == 'L0 raw early LFC']
    pg = g0.groupby('Gene')[['rho_predictor_vs_GC', 'rho']].median().reindex(GENES).dropna()
    for gene, row in pg.iterrows():
        ax.scatter(row.rho_predictor_vs_GC, row.rho, s=46, color=CDK_COLORS[gene],
                   edgecolor='white', linewidths=.6, zorder=3)
        ax.annotate(gene, (row.rho_predictor_vs_GC, row.rho), textcoords='offset points',
                    xytext=(5, 2), fontsize=5.6, color=CDK_COLORS[gene])
    lim = (-.45, .15)
    ax.plot(lim, lim, color=MUTED, lw=.8, ls=(0, (2, 2)))
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    rr = scipy.stats.pearsonr(pg.rho_predictor_vs_GC, pg.rho).statistic
    ax.set_xlabel('ρ (predictor, protospacer GC)', fontsize=7.5)
    ax.set_ylabel('ρ (predictor, raw score)', fontsize=7.5)
    med_raw = g0.rho.median()
    med_par = g0.rho_given_GC.median()
    ax.set_title(f'd   the raw agreement IS the GC confound\n'
                 f'      r = {rr:.2f};  partialling GC out: '
                 f'{med_raw:.3f} -> {med_par:.3f}', loc='left', fontsize=8.5)
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('dbNSFP predictors against the measured score, before and after correction',
                 x=0.05, ha='left', y=0.955, fontsize=10.5, color=INK)
    fig.savefig(f'{FIG}/dbnsfp_recheck.pdf')
    plt.close(fig)
    say(f'\ndbnsfp_recheck.pdf -> {FIG}')


if __name__ == '__main__':
    main()
