#!/usr/bin/env python3
"""
Early (T0 vs plasmid) residue essentiality across the CDK family.
=================================================================

Remake of `T0_plasmid_analysis.ipynb` / `T0_plasmid_analysis_v2.ipynb`
(the analysis behind supplementary figure 16).

Why the early timepoint
-----------------------
Base editing is essentially complete within days of transduction, so a guide
that destroys an essential residue should already be depleted at T0, before the
long drug-free arm accumulates secondary and clonal effects.  `T0 / plasmid` is
the only early window the screen contains: `plasmid` and `T0` are the only two
samples that bracket it.

What is rebuilt
---------------
1. **Canonical editing.**  Only reads whose every edit is the editor's own
   transition (A>G for ABE, C>T for CBE) are kept; compound `HGVSp` strings are
   exploded so a guide making `S496P_Y502H` counts at *both* residues; a guide is
   attributed to every residue it edits at >= MIN_EDIT% of canonical reads.
   This replaces the single `Top_HGVSp` / `Codon` attribution used before, and it
   also fixes two mask bugs in `T0_plasmid_analysis_v2.ipynb`, where the ABE frame
   was subset with a CBE mask (cells 7, 14) and ABE guides were flagged as
   pocket-targeting using the CBE hit list (cell 8).

2. **Two nuisance terms are removed, both estimated from guides that cannot
   carry a fitness effect.**

   * *Representation.*  `T0 / plasmid` divides a genomic-DNA measurement by a
     plasmid-DNA measurement, so it spans the whole plasmid -> genome step
     (packaging, titre, integration, a different amplicon).  A ridge model on
     protospacer / PAM composition is fit to the non-targeting and
     intron-targeting guides - the latter are edited, cannot change fitness, and
     span the full compositional range of the library - with out-of-fold
     predictions, then rescaled per infection.  This takes the correlation
     between protospacer GC and the early LFC from R ~ 0.53 to R ~ 0.14.
   * *Editing level.*  What is left still rises with how well a guide edits, in
     every infection (Spearman +0.17 to +0.45), which is backwards for a
     dropout readout.  A non-parametric trend in editing level, fit per
     infection across all targeting guides, is subtracted.  This is
     conservative: it also removes the part of the true dropout signal that
     tracks editing efficiency.

3. **Both contrasts are carried through identically** - the early
   `T0 / plasmid` and the late, same-template `DMSO / T0` - so the two can be
   benchmarked against the same controls rather than argued about.

4. **The internal positive control is reported, not assumed.**  Truncating
   (nonsense) variants destroy the protein, so in a working dropout readout they
   must deplete.  `qc_positive_control.csv` reports, per CDK, where they land in
   each contrast.

Outputs
-------
screening_data/03_LFC_FDR_tables/T0_plasmid_essentiality/
    guide_essentiality.csv            guide x editor x screen, both contrasts
    residue_essentiality.csv          gene x codon  (the score; `ess_early`, `ess_late`)
    variant_essentiality.csv          gene x variant
    residue_essentiality_SY5609.csv   same, from the independent 4th infection
    qc_*.csv / parameters.json

Usage:  python3 t0_plasmid_essentiality.py
"""
import os, re, json, warnings
import numpy as np
import pandas as pd
import scipy.stats
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
# ---------------------------------------------------------- processing level
# How much correction is applied, selectable so the whole figure set can be
# regenerated at any level and compared (see correction_ladder.py):
#   L1  centred on the non-targeting guides of the same (Editor, Screen) only
#   L2  L1 + the composition model, fit only on non-targeting and intron guides
#   L3  L2 + the editing-level detrend fit on all targeting guides
#
# L2 is the default: the non-targeting centring plus a composition model fit only on
# guides that cannot change fitness. The editing detrend (L3) is left available but off,
# because it is the one step that cannot be validated on controls and is inconsistent with
# how editing is handled elsewhere in this repo, where a guide is either above the
# attribution threshold or dropped.
#
# What that costs, stated plainly: at L2 the weighted-mean protospacer GC per residue still
# correlates with the score at r = 0.26, and only ~32% of residues clearing FDR < 0.05 are
# depletions (82% at L3). The per-residue significance list is therefore NOT usable at L2 as
# a two-sided call -- everything downstream marks significance as `FDR < 0.05 AND score < 0`,
# i.e. significantly *depleted*, which is the direction this assay is powered for.
# Family-level and per-KLIFS-column results are unaffected: the DFG+2 alanine is rank 1 of 74
# pocket columns at FDR = 0.018 at L2 (0.0074 at L3) in the within-pocket test.
# See correction_ladder.py, threshold_sweep.py and dfg_robustness.py.
LEVEL = os.environ.get('ESS_LEVEL', 'L2').upper()
assert LEVEL in ('L1', 'L2', 'L3'), f'ESS_LEVEL must be L1, L2 or L3, got {LEVEL}'
SUFFIX = '' if LEVEL == 'L2' else f'_{LEVEL}'

TABLE_DIR = os.path.join(LFC_DIR, f'T0_plasmid_essentiality{SUFFIX}')
os.makedirs(TABLE_DIR, exist_ok=True)

# ---------------------------------------------------------------- parameters
MIN_EDIT   = 10.0     # % canonical reads to attribute a guide to a residue
EDIT_HIGH  = 50.0     # high-confidence threshold for the control comparisons
MIN_BASE   = 100      # minimum raw baseline reads
N_EDIT_BINS = 12      # bins for the editing-level detrend
RIDGE_ALPHAS = np.logspace(-1, 5, 30)
N_FOLDS, SEED = 5, 0

GENES  = ['CDK7', 'CDK8', 'CDK19', 'CDK9', 'CDK12', 'CDK13', 'CDK2', 'CDK4', 'CDK6']
EDITS  = {'ABE': ('A', 'G'), 'CBE': ('C', 'T')}
POOLS  = ['subpool1', 'CDK12_13', 'CDK2_4_6']
POOL2SCREEN = {'subpool1': 'Subpool1', 'CDK12_13': 'Subpool2', 'CDK2_4_6': 'Subpool3'}
PRIMARY = ['Subpool1', 'Subpool2', 'Subpool3']
HOLDOUT = 'SY-5609'                      # independent 4th infection of the subpool-1 library
NULL_CLASSES = ['non-targeting control', 'intron']

CONTRASTS = {
    'early': dict(base='Plasmid', lfc='LFC_median_T0_Plasmid', counts='base_RAW_Plasmid',
                  label='T0 vs plasmid (early)'),
    'late':  dict(base='T0',      lfc='LFC_median_DMSO_T0',    counts='base_RAW_T0',
                  label='DMSO vs T0 (late)'),
}

say = lambda *a: print(*a, flush=True)


# ============================================================== 1. editing
def hgvsp_canonical(editor, pool):
    ref, alt = EDITS[editor]
    full = pd.read_csv(f'{EDITING}/{editor}_{pool}_full_canonical.zip')
    e = full['DNA Change'].str.findall(r'\+(\d+)([ACGT])>([ACGT])')
    keep = e.apply(lambda x: all((r == ref) and (a == alt) for _, r, a in x))
    say(f'  {editor}_{pool:<9} dropped {full.loc[~keep, "#Reads"].sum():>11,} of '
        f'{full["#Reads"].sum():>13,} reads carrying a non-{ref}>{alt} edit')
    k = full[keep].copy()
    k['edited_dna'] = k['DNA Change'].fillna('').str.strip().ne('No edit') & \
                      k['DNA Change'].fillna('').str.contains('>')
    tot = k.groupby('gRNA_id')['#Reads'].transform('sum')
    k['%Reads'] = k['#Reads'] / tot * 100
    # protein-silent: the editor did change the DNA, and the protein is unchanged
    silent = (k[(k.HGVSp == 'WT') & k.edited_dna].groupby('gRNA_id')['%Reads'].sum()
              .rename('Silent_perc'))
    unedited = (k[(k.HGVSp == 'WT') & ~k.edited_dna].groupby('gRNA_id')['%Reads'].sum()
                .rename('Unedited_perc'))
    h = k.groupby(['gRNA_id', 'HGVSp'], as_index=False)['#Reads'].sum()
    h['%Reads'] = h['#Reads'] / h.groupby('gRNA_id')['#Reads'].transform('sum') * 100
    return h, pd.concat([silent, unedited], axis=1).fillna(0.0).reset_index()


def canonical_editing():
    say('\n[1] canonical editing calls')
    savs, burd = [], []
    for editor in ['ABE', 'CBE']:
        for pool in POOLS:
            h, quiet = hgvsp_canonical(editor, pool)
            wt = h[h.HGVSp == 'WT'].groupby('gRNA_id')['%Reads'].sum().rename('WT_perc')

            nw = h[h.HGVSp != 'WT'].copy()
            nw['SAV'] = nw['HGVSp'].str.split('_')
            s = nw.explode('SAV')
            s = s[s.SAV.str.contains(r'\d', na=False)]
            s['Codon'] = s['SAV'].str.extract(r'(\d+)').astype(int)
            sav = s.groupby(['gRNA_id', 'SAV', 'Codon'], as_index=False)[['#Reads', '%Reads']].sum()
            sav['Class'] = np.where(sav.SAV.str.endswith('*'), 'Nonsense', 'Missense')

            b = pd.concat([
                sav[sav.Class == 'Nonsense'].groupby('gRNA_id')['%Reads'].sum().rename('Nonsense_perc'),
                sav[sav.Class == 'Missense'].groupby('gRNA_id')['%Reads'].sum().rename('Missense_perc'),
                wt], axis=1).fillna(0.0).reset_index()
            b = b.merge(quiet, on='gRNA_id', how='outer').fillna(0.0)
            b['edit'] = 100 - b.Unedited_perc          # any DNA edit, silent or not
            b['protein_edit'] = 100 - b.WT_perc        # a protein change

            for df in (sav, b):
                df['Editor'] = editor
                df['Screen_ID'] = POOL2SCREEN[pool]
            savs.append(sav); burd.append(b)

    SAV, BURD = pd.concat(savs, ignore_index=True), pd.concat(burd, ignore_index=True)
    # the SY-5609 arm is a second infection of the subpool-1 library: same guides,
    # same constructs, so the subpool-1 editing calls carry over to it
    for frame in ('SAV', 'BURD'):
        d = locals()[frame]
        dup = d[d.Screen_ID == 'Subpool1'].copy(); dup['Screen_ID'] = HOLDOUT
        if frame == 'SAV': SAV = pd.concat([SAV, dup], ignore_index=True)
        else:              BURD = pd.concat([BURD, dup], ignore_index=True)
    say(f'  {len(SAV):,} guide x variant rows over {SAV.gRNA_id.nunique():,} guides')
    return SAV, BURD


# ============================================================== 2. LFC
def guide_lfc(lib):
    say('\n[2] LFC tables')
    out = []
    for e in ['ABE', 'CBE']:
        m = pd.read_csv(f'{LFC_DIR}/{e}_master_LFC_table.csv',
                        low_memory=False)
        parts = []
        for key, c in CONTRASTS.items():
            k = m[m.Base_Value == c['base']][['gRNA_id', 'Screen_ID', c['lfc'], c['counts']]]
            k.columns = ['gRNA_id', 'Screen_ID', f'lfc_{key}', f'n_{key}']
            parts.append(k)
        k = parts[0].merge(parts[1], on=['gRNA_id', 'Screen_ID'], how='outer')
        k['Editor'] = e
        out.append(k)
    f = pd.concat(out, ignore_index=True)
    f = f[f.Screen_ID.isin(PRIMARY + [HOLDOUT])].merge(lib, on='gRNA_id', how='inner')
    for key in CONTRASTS:
        f.loc[f[f'n_{key}'] < MIN_BASE, f'lfc_{key}'] = np.nan
    f = f.dropna(subset=['lfc_early']).reset_index(drop=True)
    say(f'  {len(f):,} guide x editor x screen rows with a usable early LFC')
    return f


# ====================================================== 3. nuisance removal
def seq_features(lib):
    P = lib.protospacer.str.upper().values
    pam = lib.PAM.fillna('').str.upper().values
    X = {}
    for b in 'ACGT':
        X[f'f_{b}'] = np.array([s.count(b) / len(s) for s in P])
    X['GC'] = X['f_G'] + X['f_C']
    for b1 in 'ACGT':
        for b2 in 'ACGT':
            X[f'd_{b1}{b2}'] = np.array([s.count(b1 + b2) / (len(s) - 1) for s in P])
    X['homopolymer'] = np.array([max(len(m.group(0)) for m in re.finditer(r'(.)\1*', s)) for s in P])
    X['pam_GC'] = np.array([(s.count('G') + s.count('C')) / len(s) if len(s) else 0.0 for s in pam])
    return pd.DataFrame(X, index=lib.gRNA_id)


def centre_on_nt(f, col):
    nt = f[f.classification == 'non-targeting control'].groupby(['Editor', 'Screen_ID'])[col].median()
    return f[col].values - pd.MultiIndex.from_arrays([f.Editor, f.Screen_ID]).map(nt).values


def remove_representation(f, X, qc):
    """Ridge on protospacer composition, learned from guides with no possible fitness effect."""
    say(f'\n[3] representation (plasmid -> genome) correction   (ESS_LEVEL = {LEVEL})')
    if LEVEL == 'L1':
        f['lfc_early_c'] = centre_on_nt(f, 'lfc_early')
        f['lfc_late_c'] = centre_on_nt(f, 'lfc_late')
        f['pred'] = 0.0
        f['step1'] = f.lfc_early_c
        qc['representation_calibration'] = pd.DataFrame()
        qc['representation_coefficients'] = pd.DataFrame()
        say('  L1: skipped')
        return f, np.nan
    Z = X.loc[f.gRNA_id].reset_index(drop=True)
    f['lfc_early_c'] = centre_on_nt(f, 'lfc_early')
    f['lfc_late_c'] = centre_on_nt(f, 'lfc_late')

    isnull = f.classification.isin(NULL_CLASSES).values
    say(f'  null guides (non-targeting + intron): {isnull.sum():,} rows, '
        f'{f.loc[isnull, "gRNA_id"].nunique():,} guides')

    g = f.gRNA_id.values
    uniq = np.array(sorted(set(g[isnull])))
    rng = np.random.RandomState(SEED); rng.shuffle(uniq)
    fold = {u: i % N_FOLDS for i, u in enumerate(uniq)}
    folds = np.array([fold.get(x, -1) for x in g])
    oof = np.full(len(f), np.nan)
    for k in range(N_FOLDS):
        tr, te = isnull & (folds != k), isnull & (folds == k)
        mod = make_pipeline(StandardScaler(), RidgeCV(alphas=RIDGE_ALPHAS)).fit(Z[tr], f.lfc_early_c[tr])
        oof[te] = mod.predict(Z[te])
    yn = f.lfc_early_c[isnull].values
    r2 = 1 - np.sum((yn - oof[isnull]) ** 2) / np.sum((yn - yn.mean()) ** 2)
    say(f'  out-of-fold R^2 on null guides: {r2:.3f}')

    full = make_pipeline(StandardScaler(), RidgeCV(alphas=RIDGE_ALPHAS)).fit(Z[isnull], f.lfc_early_c[isnull])
    f['pred'] = full.predict(Z)
    f.loc[isnull, 'pred'] = oof[isnull]

    rows = []
    f['step1'] = np.nan
    for (e, s), idx in f.groupby(['Editor', 'Screen_ID']).groups.items():
        sub = f.loc[idx]; n = sub[sub.classification.isin(NULL_CLASSES)]
        lr = LinearRegression().fit(n[['pred']], n.lfc_early_c)
        f.loc[idx, 'step1'] = sub.lfc_early_c - lr.predict(sub[['pred']])
        res = n.lfc_early_c - lr.predict(n[['pred']])
        rows.append(dict(Editor=e, Screen_ID=s, n_null=len(n), gain=lr.coef_[0],
                         sd_before=n.lfc_early_c.std(), sd_after=res.std(),
                         var_removed=1 - res.var() / n.lfc_early_c.var()))
    f['step1'] = f.step1.values - pd.MultiIndex.from_arrays([f.Editor, f.Screen_ID]).map(
        f[f.classification == 'non-targeting control'].groupby(['Editor', 'Screen_ID']).step1.median()).values
    cal = pd.DataFrame(rows).sort_values(['Screen_ID', 'Editor'])
    qc['representation_calibration'] = cal
    say('\n  per-infection calibration on null guides:')
    say(cal.round(3).to_string(index=False))

    coef = pd.Series(full[-1].coef_, index=Z.columns).sort_values(key=np.abs, ascending=False)
    qc['representation_coefficients'] = coef.rename('coef').reset_index().rename(columns={'index': 'feature'})
    return f, r2


def detrend_editing(f, qc, contributing=None):
    """
    Subtract the dependence of each contrast on how well the guide edits.

    The trend is fit per infection, non-parametrically (medians in 12 editing-level
    bins, linearly interpolated), across *all* targeting guides.  Most residues are
    not essential, so the trend is dominated by guides with nothing to report; what
    survives is a guide's departure from what guides of the same editing level do in
    the same infection, which is the reference the score is read against.

    This is deliberately conservative: the part of a real dropout signal that tracks
    editing efficiency is removed along with the artifact.

    Guides with no editing measurement (non-targeting and intron controls have no
    sensor read-out) keep their un-detrended value and are marked `detrended = False`;
    they are used only to set the null scale, never as a residue measurement.

    The curve is fit only on guides that actually contribute to a residue score — i.e.
    that make some missense variant at >= MIN_EDIT% of canonical reads. Fitting it on all
    targeting guides instead pulls in ~43% that are then discarded, and those sit at
    systematically lower editing (median 73% vs 85%), so they shape bins the analysis
    never uses. In practice it makes almost no difference (residue-level GC correlation
    0.092 vs 0.089, truncating control -0.124 vs -0.121), but the fit set and the
    analysis set should be the same set.
    """
    say(f'\n[4] editing-level detrend   (ESS_LEVEL = {LEVEL})')
    if LEVEL == 'L1':
        say('  L1: no composition correction and no detrend; the score is the '
            'non-targeting-centred LFC')
        f['ess_early'], f['ess_late'] = f.lfc_early_c, f.lfc_late_c
        f['detrended_early'] = f['detrended_late'] = False
        qc['editing_detrend'] = pd.DataFrame()
        return f
    if LEVEL == 'L2':
        say('  L2: composition correction only; the editing-level trend is left in')
        f['ess_early'], f['ess_late'] = f.step1, f.lfc_late_c
        f['detrended_early'] = f['detrended_late'] = False
        qc['editing_detrend'] = pd.DataFrame()
        return f
    rows = []
    for key, src in [('early', 'step1'), ('late', 'lfc_late_c')]:
        col = f'ess_{key}'
        f[col] = f[src]
        f[f'detrended_{key}'] = False
        for (e, s), idx in f.groupby(['Editor', 'Screen_ID']).groups.items():
            d = f.loc[idx]
            d = d[(d.classification == 'targeting') & d.edit.notna() & d[src].notna()]
            fit = d[d.contributes] if contributing is not None else d
            if len(fit) < 200:
                continue
            q = np.quantile(fit.edit, np.linspace(0, 1, N_EDIT_BINS + 1)); q[0] -= 1e-9
            lab = np.digitize(fit.edit, q[1:-1])
            med = fit.groupby(lab)[src].median(); ctr = fit.groupby(lab)['edit'].median()
            f.loc[d.index, col] = d[src] - np.interp(d.edit, ctr.values, med.values)
            f.loc[d.index, f'detrended_{key}'] = True
            rows.append(dict(contrast=key, Editor=e, Screen_ID=s,
                             n_fit=len(fit), n_applied=len(d),
                             rho_before=scipy.stats.spearmanr(d.edit, d[src]).statistic,
                             rho_after=scipy.stats.spearmanr(d.edit, f.loc[d.index, col]).statistic))
    tab = pd.DataFrame(rows)
    qc['editing_detrend'] = tab
    say(tab.round(3).to_string(index=False))
    return f


def null_scale(f, qc):
    """
    How wide is the score when nothing is happening?

    Three reference sets, all measured on the final score:
      * unedited targeting guides (canonical WT >= 90%) - the matched null: the same
        amplicon, the same detrend, no protein change;
      * intron-targeting guides - edited, but in an intron;
      * non-targeting controls.
    The matched null sets the scale used for the per-residue z and FDR.
    """
    rows = {}
    sets = {
        'unedited targeting (WT>=90%)': (f.classification == 'targeting') & (f.WT_perc >= 90),
        'intron':                        f.classification == 'intron',
        'non-targeting':                 f.classification == 'non-targeting control',
    }
    for lab, m in sets.items():
        rows[lab] = {k: dict(n=int((m & f[f'ess_{k}'].notna()).sum()),
                             median=f.loc[m, f'ess_{k}'].median(),
                             sd=f.loc[m, f'ess_{k}'].std()) for k in ('early', 'late')}
    tab = pd.DataFrame([dict(reference=lab, contrast=k, **v)
                        for lab, d in rows.items() for k, v in d.items()])
    qc['null_scale'] = tab
    say('\n  null scale of the final score:')
    say(tab.round(3).to_string(index=False))
    ref = 'unedited targeting (WT>=90%)'
    return {k: float(rows[ref][k]['sd']) for k in ('early', 'late')}


# ============================================================ 4. validation
def validate(f, X, qc):
    say('\n[5] validation')
    # (a) composition
    rows = []
    for (e, s), d in f.groupby(['Editor', 'Screen_ID']):
        gc = X.loc[d.gRNA_id, 'GC'].values
        for cls, lab in [(NULL_CLASSES, 'null (NT+intron)'), (['targeting'], 'targeting')]:
            m = d.classification.isin(cls).values
            if m.sum() < 30: continue
            ok = m & d.ess_early.notna().values
            rows.append(dict(Editor=e, Screen_ID=s, guide_class=lab, n=int(m.sum()),
                             R_GC_raw=scipy.stats.pearsonr(gc[m], d.lfc_early_c[m]).statistic,
                             R_GC_step1=scipy.stats.pearsonr(gc[m], d.step1[m]).statistic,
                             R_GC_final=(scipy.stats.pearsonr(gc[ok], d.ess_early[ok]).statistic
                                         if ok.sum() > 30 else np.nan)))
    qc['gc_dependence'] = pd.DataFrame(rows)
    say('\n  (a) protospacer GC vs the early score, through the two corrections:')
    say(qc['gc_dependence'].round(3).to_string(index=False))

    # (b) cross-editor reproducibility
    rows = []
    for stage, col in [('raw', 'lfc_early_c'), ('corrected', 'ess_early'), ('late', 'ess_late')]:
        w = f.pivot_table(index=['gRNA_id', 'Screen_ID', 'classification'],
                          columns='Editor', values=col).dropna().reset_index()
        for s in PRIMARY + [HOLDOUT]:
            for cls in NULL_CLASSES + ['targeting']:
                d = w[(w.Screen_ID == s) & (w.classification == cls)]
                if len(d) < 30: continue
                rows.append(dict(stage=stage, Screen_ID=s, guide_class=cls, n=len(d),
                                 R=scipy.stats.pearsonr(d.ABE, d.CBE).statistic))
    qc['cross_editor'] = pd.DataFrame(rows)
    say('\n  (b) ABE vs CBE, two independent infections of one plasmid pool:')
    say(qc['cross_editor'].pivot_table(index=['Screen_ID', 'guide_class'], columns='stage',
                                       values='R').round(3).to_string())
    return f


GUIDE_CLASSES = ['Non-targeting', 'Intron', 'Unedited', 'Silent', 'Missense', 'Truncating']


def classify_guides(f, qc):
    """
    One mutational class per guide, for the distribution panels.

    Mutually exclusive, in this order:
      Truncating    >= 25% of canonical reads carry a stop
      Missense      < 25% truncating, >= 25% missense
      Silent        protein unchanged in >= 90% of reads, yet >= 50% of reads carry a
                    DNA edit - the editor fired and the protein did not change
      Unedited      >= 90% of reads carry no edit at all
      Intron / Non-targeting   from the library annotation
    Guides that fall in none of these (low or mixed editing) are left unlabelled and
    are not plotted.
    """
    cls = pd.Series(np.nan, index=f.index, dtype=object)
    t = f.classification == 'targeting'
    cls[t & (f.Unedited_perc >= 90)] = 'Unedited'
    cls[t & (f.WT_perc >= 90) & (f.Silent_perc >= 50)] = 'Silent'
    cls[t & (f.Missense_perc >= 25) & (f.Nonsense_perc < 25)] = 'Missense'
    cls[t & (f.Nonsense_perc >= 25)] = 'Truncating'
    cls[f.classification == 'intron'] = 'Intron'
    cls[f.classification == 'non-targeting control'] = 'Non-targeting'
    f['GuideClass'] = cls
    tab = (f.dropna(subset=['GuideClass'])
           .groupby('GuideClass')[['ess_early', 'ess_late']]
           .agg(['count', 'median', 'std']).reindex(GUIDE_CLASSES))
    qc['guide_class_distribution'] = tab.round(3).reset_index()
    say('\n  guide classes:')
    say(tab.round(3).to_string())
    return f


def positive_control(f, qc):
    """Truncating variants destroy the protein: in a dropout readout they must deplete."""
    say('\n[6] positive control - truncating variants')
    g = f[(f.classification == 'targeting') & f.Nonsense_perc.notna()]
    rows = []
    for gene in GENES:
        s = g[(g.Gene == gene) & (g.Nonsense_perc >= EDIT_HIGH)]
        if len(s) < 5: continue
        o = g[(g.Gene != gene) & (g.Nonsense_perc >= EDIT_HIGH)]
        rows.append(dict(
            Gene=gene, n=len(s),
            early_raw=s.lfc_early_c.median(), early=s.ess_early.median(),
            late=s.ess_late.median(),
            p_early=scipy.stats.mannwhitneyu(s.ess_early.dropna(), o.ess_early.dropna()).pvalue,
            p_late=scipy.stats.mannwhitneyu(s.ess_late.dropna(), o.ess_late.dropna()).pvalue))
    pc = pd.DataFrame(rows)
    qc['positive_control'] = pc
    say(pc.round(3).to_string(index=False))

    rows = []
    for lo, hi, lab in [(25, 50, '25-50%'), (50, 75, '50-75%'), (75, 101, '>=75%')]:
        s = g[(g.Nonsense_perc >= lo) & (g.Nonsense_perc < hi)]
        rows.append(dict(bin=lab, n=len(s), early_raw=s.lfc_early_c.median(),
                         early=s.ess_early.median(), late=s.ess_late.median()))
    qc['truncating_dose'] = pd.DataFrame(rows)
    say('\n  pooled dose response (should fall as editing rises):')
    say(qc['truncating_dose'].round(3).to_string(index=False))

    # per-replicate consistency of the CDK9 inversion
    rows = []
    for (e, s), d in g[g.Gene == 'CDK9'].groupby(['Editor', 'Screen_ID']):
        d = d[d.Nonsense_perc >= EDIT_HIGH]
        if len(d) < 5: continue
        rows.append(dict(Editor=e, Screen_ID=s, n=len(d),
                         early=d.ess_early.median(), late=d.ess_late.median()))
    qc['cdk9_by_infection'] = pd.DataFrame(rows)
    if len(qc['cdk9_by_infection']):
        say('\n  CDK9 truncating guides, per infection:')
        say(qc['cdk9_by_infection'].round(3).to_string(index=False))
    return pc


# ======================================================= 5. residue scores
def residue_scores(f, SAV, qc):
    say('\n[7] residue scores')
    v = SAV[(SAV['%Reads'] >= MIN_EDIT) & (SAV.Class == 'Missense')].merge(
        f[['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'classification',
           'ess_early', 'ess_late', 'lfc_early_c', 'edit']],
        on=['gRNA_id', 'Editor', 'Screen_ID'], how='inner')
    v = v[v.classification == 'targeting']
    main = v[v.Screen_ID.isin(PRIMARY)].copy()
    hold = v[v.Screen_ID == HOLDOUT].copy()
    say(f'  {len(main):,} guide x variant rows, {main.gRNA_id.nunique():,} guides, '
        f'{main.SAV.nunique():,} variants')

    null_sd = null_scale(f, qc)
    say(f'  null sd used for z / FDR: early {null_sd["early"]:.3f}, late {null_sd["late"]:.3f}')

    def score(df):
        df = df.copy(); df['w'] = df['%Reads'] / 100.0
        def one(d):
            out = {}
            for k in ('early', 'late'):
                y = d[f'ess_{k}']; ok = y.notna()
                if ok.sum():
                    w = d.w[ok]
                    out[f'ess_{k}'] = float(np.average(y[ok], weights=w))
                    out[f'n_eff_{k}'] = float(w.sum() ** 2 / np.sum(w ** 2))
                else:
                    out[f'ess_{k}'] = np.nan; out[f'n_eff_{k}'] = 0.0
            out['ess_early_raw'] = float(np.average(d.lfc_early_c, weights=d.w))
            out['n_guides'] = d.gRNA_id.nunique()
            out['n_obs'] = len(d)
            out['n_editors'] = d.Editor.nunique()
            out['max_edit'] = d['%Reads'].max()
            out['mean_edit'] = d['%Reads'].mean()
            out['variants'] = ';'.join(sorted(set(d.SAV)))
            return pd.Series(out)
        agg = df.groupby(['Gene', 'Codon']).apply(one).reset_index()
        for k in ('early', 'late'):
            agg[f'z_{k}'] = agg[f'ess_{k}'] / (null_sd[k] / np.sqrt(agg[f'n_eff_{k}'].clip(lower=1e-9)))
            p = 2 * scipy.stats.norm.sf(np.abs(agg[f'z_{k}']))
            agg[f'p_{k}'] = p
            agg[f'FDR_{k}'] = np.nan
            for gene, idx in agg.groupby('Gene').groups.items():
                pv = agg.loc[idx, f'p_{k}'].values
                m = np.isfinite(pv)
                q = np.full(len(pv), np.nan)
                if m.sum():
                    pp = pv[m]; o = np.argsort(pp); qq = np.empty_like(pp)
                    qq[o] = np.minimum.accumulate(
                        (pp[o] * len(pp) / (np.arange(len(pp)) + 1))[::-1])[::-1]
                    q[m] = np.clip(qq, 0, 1)
                agg.loc[idx, f'FDR_{k}'] = q
        return agg

    res = score(main)
    hol = score(hold) if len(hold) else pd.DataFrame()

    def wmean(y, w):
        ok = y.notna()
        return float(np.average(y[ok], weights=w[ok])) if ok.any() else np.nan

    var = (main.assign(w=lambda d: d['%Reads'] / 100.0)
           .groupby(['Gene', 'SAV', 'Codon'])
           .apply(lambda d: pd.Series({
               'ess_early': wmean(d.ess_early, d.w),
               'ess_late': wmean(d.ess_late, d.w),
               'n_guides': d.gRNA_id.nunique(), 'max_edit': d['%Reads'].max(),
               'editors': ','.join(sorted(set(d.Editor)))})).reset_index())

    rep = {}
    for k in ('early', 'late'):
        w = (main.assign(w=lambda d: d['%Reads'] / 100.0)
             .groupby(['Gene', 'Codon', 'Editor'])
             .apply(lambda d: wmean(d[f'ess_{k}'], d.w))
             .unstack().dropna())
        if len(w) > 30:
            rep[f'{k}: ABE vs CBE'] = dict(n=len(w), R=scipy.stats.pearsonr(w.ABE, w.CBE).statistic,
                                           rho=scipy.stats.spearmanr(w.ABE, w.CBE).statistic)
    if len(hol):
        j = res.merge(hol, on=['Gene', 'Codon'], suffixes=('', '_h')).dropna(subset=['ess_early', 'ess_early_h'])
        if len(j) > 30:
            rep['early: primary vs SY-5609'] = dict(
                n=len(j), R=scipy.stats.pearsonr(j.ess_early, j.ess_early_h).statistic,
                rho=scipy.stats.spearmanr(j.ess_early, j.ess_early_h).statistic)
    jj = res.dropna(subset=['ess_early', 'ess_late'])
    rep['early vs late'] = dict(n=len(jj), R=scipy.stats.pearsonr(jj.ess_early, jj.ess_late).statistic,
                                rho=scipy.stats.spearmanr(jj.ess_early, jj.ess_late).statistic)
    # ---- matched empirical null for the residue score.
    # The same guides, the same editing weights and the same number of guides per
    # residue; only the guide -> residue assignment is permuted, within each gene.
    # Repeating it gives every residue its own null distribution, so a residue
    # measured by one weak guide is judged against a null built the same way.
    N_PERM = 200
    rng = np.random.RandomState(SEED)
    keys = main.groupby(['Gene', 'Codon']).ngroup().values
    w = (main['%Reads'] / 100.0).values
    y = main.ess_early.values
    ok = np.isfinite(y)
    nk = keys.max() + 1
    perm_vals = np.full((N_PERM, nk), np.nan)
    gene_idx = {g: np.where(main.Gene.values == g)[0] for g in main.Gene.unique()}
    for b in range(N_PERM):
        yp = y.copy()
        for g, idx in gene_idx.items():
            sub = idx[ok[idx]]
            yp[sub] = y[sub][rng.permutation(len(sub))]
        num = np.bincount(keys[ok], weights=(yp * w)[ok], minlength=nk)
        den = np.bincount(keys[ok], weights=w[ok], minlength=nk)
        perm_vals[b] = np.where(den > 0, num / den, np.nan)
    key_tab = (main.assign(_k=keys).groupby(['Gene', 'Codon'])['_k'].first().reset_index())
    null_sd_res = pd.Series(np.nanstd(perm_vals, axis=0), name='null_sd_res')
    key_tab['null_sd_res'] = null_sd_res.values[key_tab._k.values]
    res = res.merge(key_tab[['Gene', 'Codon', 'null_sd_res']], on=['Gene', 'Codon'], how='left')
    qc['residue_null'] = pd.DataFrame({'ess_early_permuted': perm_vals[np.isfinite(perm_vals)]})
    say(f'  empirical null over {N_PERM} permutations: sd {np.nanstd(perm_vals):.3f} '
        f'vs observed {res.ess_early.std():.3f}')

    # re-derive the early z / FDR from the matched empirical null
    res['z_early'] = res.ess_early / res.null_sd_res.replace(0, np.nan)
    res['p_early'] = 2 * scipy.stats.norm.sf(np.abs(res.z_early))
    res['FDR_early'] = np.nan
    for gene, idx in res.groupby('Gene').groups.items():
        pv = res.loc[idx, 'p_early'].values
        m = np.isfinite(pv); q = np.full(len(pv), np.nan)
        if m.sum():
            pp = pv[m]; o = np.argsort(pp); qq = np.empty_like(pp)
            qq[o] = np.minimum.accumulate((pp[o] * len(pp) / (np.arange(len(pp)) + 1))[::-1])[::-1]
            q[m] = np.clip(qq, 0, 1)
        res.loc[idx, 'FDR_early'] = q
    say(f'  residues at FDR < 0.05 (early): {(res.FDR_early < 0.05).sum()}')

    qc['reproducibility'] = pd.DataFrame(rep).T.reset_index().rename(columns={'index': 'comparison'})
    say('\n  reproducibility of the residue score:')
    say(qc['reproducibility'].round(3).to_string(index=False))

    # the guide x residue rows behind every score, so downstream work (split-half
    # reliability, structural tests) does not have to re-derive the attribution
    keep = ['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'SAV', 'Codon', '%Reads',
            'edit', 'ess_early', 'ess_late', 'lfc_early_c']
    main[keep].to_csv(f'{TABLE_DIR}/residue_guide_rows.csv', index=False)
    say(f'  wrote {len(main):,} guide x residue rows')
    return res, hol, var, null_sd


def annotate(res):
    K = pd.read_csv(f'{KLIFS_DIR}/KLIFS_compiled.csv')
    pocket, idx = {}, {}
    for g in GENES:
        s = K[['KLIFS_idx', g]].dropna()
        pocket[g] = set(s[g].astype(int)); idx[g] = dict(zip(s[g].astype(int), s.KLIFS_idx))
    res['in_pocket'] = [c in pocket[g] for g, c in zip(res.Gene, res.Codon)]
    res['KLIFS_idx'] = [idx[g].get(c) for g, c in zip(res.Gene, res.Codon)]
    return res


def main():
    qc = {}
    lib = pd.read_csv(f'{LIBRARY}/CDK_library_final.csv')[
        ['gRNA_id', 'Gene', 'protospacer', 'PAM', 'classification', 'Pool']]
    SAV, BURD = canonical_editing()
    f = guide_lfc(lib)
    f = f.merge(BURD[['gRNA_id', 'Editor', 'Screen_ID', 'edit', 'protein_edit',
                      'Nonsense_perc', 'Missense_perc', 'WT_perc', 'Silent_perc',
                      'Unedited_perc']],
                on=['gRNA_id', 'Editor', 'Screen_ID'], how='left')

    X = seq_features(lib)
    f, r2 = remove_representation(f, X, qc)
    # guides that make a missense variant above the attribution threshold — the only
    # ones that can ever enter a residue score, and so the only ones the detrend is fit on
    contrib = (SAV[(SAV['%Reads'] >= MIN_EDIT) & (SAV.Class == 'Missense')]
               [['gRNA_id', 'Editor', 'Screen_ID']].drop_duplicates().assign(contributes=True))
    f = f.merge(contrib, on=['gRNA_id', 'Editor', 'Screen_ID'], how='left')
    f['contributes'] = f.contributes.fillna(False)
    say(f'\n  {int(f.contributes.sum()):,} of {int((f.classification == "targeting").sum()):,} '
        'targeting guide rows contribute to a residue score; the detrend is fit on those')
    f = detrend_editing(f, qc, contributing=True)
    f = validate(f, X, qc)
    f = classify_guides(f, qc)
    pc = positive_control(f, qc)
    res, hol, var, null_sd = residue_scores(f, SAV, qc)
    res = annotate(res)

    rows = []
    for gene in GENES:
        s = res[res.Gene == gene]
        for k in ('early', 'late'):
            a, b = s[s.in_pocket][f'ess_{k}'].dropna(), s[~s.in_pocket][f'ess_{k}'].dropna()
            if len(a) >= 5 and len(b) >= 5:
                rows.append(dict(Gene=gene, contrast=k, n_pocket=len(a), n_other=len(b),
                                 med_pocket=a.median(), med_other=b.median(),
                                 delta=a.median() - b.median(),
                                 p=scipy.stats.mannwhitneyu(a, b).pvalue))
    qc['pocket_contrast'] = pd.DataFrame(rows)
    say('\n[8] ATP pocket vs the rest of the protein:')
    say(qc['pocket_contrast'].pivot_table(index='Gene', columns='contrast',
                                          values=['delta', 'p']).round(4).to_string())

    f.to_csv(f'{TABLE_DIR}/guide_essentiality.csv', index=False)
    res.to_csv(f'{TABLE_DIR}/residue_essentiality.csv', index=False)
    var.to_csv(f'{TABLE_DIR}/variant_essentiality.csv', index=False)
    if len(hol): hol.to_csv(f'{TABLE_DIR}/residue_essentiality_SY5609.csv', index=False)
    for k, t in qc.items():
        t.to_csv(f'{TABLE_DIR}/qc_{k}.csv', index=False)
    json.dump(dict(level=LEVEL, MIN_EDIT=MIN_EDIT, EDIT_HIGH=EDIT_HIGH, MIN_BASE=MIN_BASE,
                   null_classes=NULL_CLASSES, null_sd={k: float(v) for k, v in null_sd.items()},
                   representation_oof_r2=float(r2), genes=GENES,
                   primary_screens=PRIMARY, holdout=HOLDOUT),
              open(f'{TABLE_DIR}/parameters.json', 'w'), indent=2)
    say(f'\nwrote tables to {TABLE_DIR}')
    say(res.groupby('Gene').agg(residues=('Codon', 'nunique'),
                                med_early=('ess_early', 'median'),
                                med_late=('ess_late', 'median'),
                                sig_late=('FDR_late', lambda s: int((s < 0.05).sum())),
                                sig_early=('FDR_early', lambda s: int((s < 0.05).sum()))
                                ).reindex(GENES).round(3).to_string())


if __name__ == '__main__':
    main()
