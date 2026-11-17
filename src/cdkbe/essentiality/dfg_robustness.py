#!/usr/bin/env python3
"""
Does the DFG finding survive without the editing detrend?

The detrend (L3) is the one correction that cannot be validated on control guides and is
inconsistent with how the rest of this repo treats editing — everywhere else a guide is
either above the attribution threshold or it is dropped. The alternative the data allows is
L2 (composition correction only, fit on non-targeting + intron guides) plus a higher
attribution threshold, which is the convention used elsewhere.

This sweeps (processing level) x (minimum % reads making the variant) and reports, for each
combination, the things the figure actually claims:

  * the DFG / activation-loop peak: height, permutation p, FDR
  * the two columns that carry it (a.l:85 = DFG+2 Ala, xDFG:83 = DFG Gly): rank and p
  * residue-level protospacer GC correlation — the artifact check that chose L3 before
  * how many residues clear FDR < 0.05 and what fraction are depletions
  * the truncating control, compared inside one reference frame
  * coverage

Usage:  python3 dfg_robustness.py
Writes: screening_data/.../dfg_robustness.csv
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import scipy.stats
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from . import t0_plasmid_essentiality as E
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
from .essentiality_peaks import (alignment_maps, column_annotations, kinase_domain_window,
                                SIGMA, DISTANCE, HEIGHT, PROMINENCE, DFG_BLOCK, GENES, TAB)

N_PERM_RES, N_PERM_PEAK, SEED = 200, 1000, 0
SCHEMES = [('L1 centred only', 'lfc_early_c'),
           ('L2 composition', 'step1'),
           ('L3 composition + detrend', 'ess_early')]
THRESHOLDS = [10, 25, 50, 75]
say = E.say


def residue_table(SAV, f, col, thr, gc_map, rng):
    v = SAV[(SAV['%Reads'] >= thr) & (SAV.Class == 'Missense')].merge(
        f[['gRNA_id', 'Editor', 'Screen_ID', 'Gene', 'classification', col]],
        on=['gRNA_id', 'Editor', 'Screen_ID'], how='inner')
    v = v[(v.classification == 'targeting') & v[col].notna()
          & v.Screen_ID.isin(E.PRIMARY)].copy()
    v['w'] = v['%Reads'] / 100.0
    v['gc'] = v.gRNA_id.map(gc_map)
    keys = v.groupby(['Gene', 'Codon']).ngroup().values
    w, y = v.w.values, v[col].values
    nk = keys.max() + 1
    num = np.bincount(keys, weights=y * w, minlength=nk)
    den = np.bincount(keys, weights=w, minlength=nk)
    score = num / den
    gcs = np.bincount(keys, weights=v.gc.values * w, minlength=nk) / den

    # matched permutation null: shuffle guide -> residue within gene
    gene_idx = {g: np.where(v.Gene.values == g)[0] for g in v.Gene.unique()}
    perm = np.empty((N_PERM_RES, nk))
    for b in range(N_PERM_RES):
        yp = y.copy()
        for g, idx in gene_idx.items():
            yp[idx] = y[idx][rng.permutation(len(idx))]
        perm[b] = np.bincount(keys, weights=yp * w, minlength=nk) / den
    sd = np.nanstd(perm, axis=0)

    key_tab = v.assign(_k=keys).groupby(['Gene', 'Codon'])['_k'].first().reset_index()
    t = key_tab.assign(score=score[key_tab._k], res_GC=gcs[key_tab._k],
                       null_sd=sd[key_tab._k])
    t['z'] = t.score / t.null_sd.replace(0, np.nan)
    t['p'] = 2 * scipy.stats.norm.sf(np.abs(t.z))
    t['FDR'] = np.nan
    for gene, idx in t.groupby('Gene').groups.items():
        pv = t.loc[idx, 'p'].values
        m = np.isfinite(pv); q = np.full(len(pv), np.nan)
        if m.sum():
            pp = pv[m]; o = np.argsort(pp); qq = np.empty_like(pp)
            qq[o] = np.minimum.accumulate((pp[o] * len(pp) / (np.arange(len(pp)) + 1))[::-1])[::-1]
            q[m] = np.clip(qq, 0, 1)
        t.loc[idx, 'FDR'] = q
    return t


def peak_test(t, res2col, col2res, col_klifs, cons, aln_len, kd, rng):
    """The DFG peak, and the two columns that carry it, under this scheme."""
    calls = np.zeros((len(GENES), aln_len))
    cont = np.full((len(GENES), aln_len), np.nan)
    for i, g in enumerate(GENES):
        d = t[t.Gene == g]
        for codon, s, q in zip(d.Codon, d.score, d.FDR):
            c = res2col[g].get(int(codon))
            if c is None:
                continue
            cont[i, c] = -s
            if q < 0.05 and s < 0:
                calls[i, c] = 1

    y = calls.sum(0)
    ys = gaussian_filter1d(y, SIGMA)
    idx, _ = find_peaks(ys, height=HEIGHT, distance=DISTANCE, prominence=PROMINENCE)
    idx = idx[(idx >= kd[0]) & (idx <= kd[1])]
    _, _, left, right = peak_widths(ys, idx, rel_height=0.5)

    h_null = []
    for _ in range(N_PERM_PEAK):
        sh = np.vstack([np.roll(calls[i], rng.randint(aln_len)) for i in range(len(GENES))])
        yp = gaussian_filter1d(sh.sum(0), SIGMA)
        ip, _ = find_peaks(yp, height=HEIGHT, distance=DISTANCE, prominence=PROMINENCE)
        ip = ip[(ip >= kd[0]) & (ip <= kd[1])]
        h_null.extend(yp[ip])
    h_null = np.array(h_null)

    out = dict(n_peaks=len(idx), dfg_height=np.nan, dfg_p=np.nan, dfg_FDR=np.nan,
               dfg_found=False)
    ps = []
    for apex, lo, hi in zip(idx, left, right):
        lo, hi = int(np.floor(lo)), int(np.ceil(hi))
        labels = [col_klifs[c] for c in range(lo, hi + 1) if c in col_klifs]
        p = (h_null >= ys[apex]).mean() if len(h_null) else np.nan
        ps.append(p)
        if any(str(l).startswith(DFG_BLOCK) for l in labels):
            out.update(dfg_found=True, dfg_height=float(ys[apex]), dfg_p=float(p))
    if len(ps):
        order = scipy.stats.rankdata(ps)
        for (apex, lo, hi), p, r in zip(zip(idx, left, right), ps, order):
            lo, hi = int(np.floor(lo)), int(np.ceil(hi))
            labels = [col_klifs[c] for c in range(lo, hi + 1) if c in col_klifs]
            if any(str(l).startswith(DFG_BLOCK) for l in labels):
                out['dfg_FDR'] = float(min(p * len(ps) / r, 1))

    # column-level: a.l:85 and xDFG:83
    obs = np.nanmean(cont, axis=0)
    null = np.empty((N_PERM_PEAK, aln_len))
    for b in range(N_PERM_PEAK):
        sh = np.vstack([np.roll(cont[i], rng.randint(aln_len)) for i in range(len(GENES))])
        with np.errstate(invalid='ignore'):
            null[b] = np.nan_to_num(np.nanmean(sh, axis=0), nan=0.0)
    pcol = ((null >= obs).sum(0) + 1) / (N_PERM_PEAK + 1)
    ncov = (~np.isnan(cont)).sum(0)
    rows = [dict(klifs=col_klifs[c], dep=obs[c], p=pcol[c])
            for c in range(aln_len) if ncov[c] >= 7 and col_klifs.get(c)]
    ct = pd.DataFrame(rows).sort_values('dep', ascending=False).reset_index(drop=True)
    ct['rank'] = np.arange(1, len(ct) + 1)
    for k, tag in [('a.l:85', 'A85'), ('xDFG:83', 'G83')]:
        r = ct[ct.klifs == k]
        out[f'rank_{tag}'] = int(r['rank'].iloc[0]) if len(r) else np.nan
        out[f'p_{tag}'] = float(r.p.iloc[0]) if len(r) else np.nan
    out['n_cols_tested'] = len(ct)
    out['dfg_in_top7'] = int(ct.head(7).klifs.str.startswith(DFG_BLOCK).sum())
    return out


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
    contrib = (SAV[(SAV['%Reads'] >= E.MIN_EDIT) & (SAV.Class == 'Missense')]
               [['gRNA_id', 'Editor', 'Screen_ID']].drop_duplicates().assign(contributes=True))
    f = f.merge(contrib, on=['gRNA_id', 'Editor', 'Screen_ID'], how='left')
    f['contributes'] = f.contributes.fillna(False)
    f = E.detrend_editing(f, qc, contributing=True)
    gc_map = dict(zip(lib.gRNA_id, X.GC.values))

    aln_len, res2col, col2res, chars = alignment_maps()
    col_klifs, cons, modal, ngap = column_annotations(aln_len, res2col, chars)
    kd = kinase_domain_window(res2col)

    rows = []
    for label, col in SCHEMES:
        for thr in THRESHOLDS:
            rng = np.random.RandomState(SEED)
            t = residue_table(SAV, f, col, thr, gc_map, rng)
            sig = t[t.FDR < 0.05]
            pk = peak_test(t, res2col, col2res, col_klifs, cons, aln_len, kd, rng)

            # truncating control inside one frame, at this threshold
            tg = f[(f.classification == 'targeting') & f[col].notna() & (f.edit >= thr)]
            istr = (tg.Nonsense_perc >= thr).fillna(False)
            a, b = tg[istr][col].dropna(), tg[~istr][col].dropna()

            rows.append(dict(
                scheme=label, min_edit=thr, n_residues=len(t), n_sig=len(sig),
                frac_sig_depleted=(sig.score < 0).mean() if len(sig) else np.nan,
                r_residue_GC=scipy.stats.pearsonr(t.res_GC, t.score).statistic,
                trunc_minus_targeting=a.median() - b.median() if len(a) >= 20 else np.nan,
                p_trunc=scipy.stats.mannwhitneyu(a, b).pvalue if len(a) >= 20 else np.nan,
                **pk))
            r = rows[-1]
            say(f"{label:26s} >= {thr:2d}%  residues {r['n_residues']:5d}  sig {r['n_sig']:4d} "
                f"({100*(r['frac_sig_depleted'] or 0):3.0f}% depleted)  GC {r['r_residue_GC']:+.3f}  "
                f"DFG p {r['dfg_p'] if np.isfinite(r['dfg_p']) else float('nan'):.3f}  "
                f"A85 rank {r['rank_A85']}")
    t = pd.DataFrame(rows)
    t.to_csv(f'{TAB}/dfg_robustness.csv', index=False)
    say(f'\nwrote {TAB}/dfg_robustness.csv')
    return t


if __name__ == '__main__':
    main()
