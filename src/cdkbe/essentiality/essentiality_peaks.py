#!/usr/bin/env python3
"""
Peak finding on the early essentiality score, and orthogonal checks that the
score is measuring protein function.

The peak analysis is the mirror image of the resistance peak analysis in
`CDK_resistance_patterns.ipynb` / `canonical_resistance_figures.py`: every CDK's
residues are placed on the family-wide Clustal alignment, the per-column signal is
summed across CDKs, smoothed with the same Gaussian, and `scipy.signal.find_peaks`
is run with the same settings — except the signal counted here is a **depleted**
residue (FDR < 0.05, score < 0) rather than an enriched one.

Everything is tested against a null that keeps each gene's own distribution and
its own coverage, and only breaks the link between a residue and its position:
the per-gene column profile is **circularly shifted** by a random offset, 2,000
times. That preserves the number of calls, the autocorrelation of the track and
the alignment's gap structure, so the only thing destroyed is whether CDKs agree
on *where* the depleted residues are.

Orthogonal checks (none of them uses the screen):
  * alignment conservation  - how many of the nine CDKs share the modal residue
  * catalytic motifs        - the invariant VAIK / HRD / DFG columns
  * predicted disorder      - metapredict, source_data/20_metapredict_disorder
  * burial                  - relative solvent accessibility, computed with
                              Shrake-Rupley on each CDK's own structure
  * AlphaMissense           - source_data/12_alphamissense_scores (CDK7/8/9/19)

Usage:  python3 essentiality_peaks.py
Writes: figures/T0_plasmid_essentiality/{essentiality_peaks_MSA,peak_significance,
        orthogonal_validation}.pdf and the matching tables.
"""
import os, re, json, warnings, collections
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
from Bio import AlignIO
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
ALN = f'{KLIFS_DIR}/clustalo-I20250722-154001-0730-44590675-p1m.aln-clustal_num'

import sys
from .._chimerax_config import STRUCTURES

from .essentiality_figures import (CMAP, NORM, VLIM, DEPLETED, ENRICHED, NO_DATA,
                                  INK, INK2, MUTED, POCKET_FILL, GENES, res, KLIFS,
                                  TAB, FIG, LEVEL)   # FIG/TAB follow ESS_LEVEL
os.makedirs(FIG, exist_ok=True)

# same smoothing / detection settings as the resistance peaks
SIGMA, DISTANCE = 5, 10
HEIGHT, PROMINENCE = 0.2, 0.3
N_PERM = 2000
SEED = 0
REFERENCE = 'CDK2'
FDR_CUT = 0.05

say = lambda *a: print(*a, flush=True)


# ------------------------------------------------------------- alignment
def alignment_maps():
    aln = AlignIO.read(ALN, 'clustal')
    res2col, col2res, chars = {}, {}, {}
    for rec in aln:
        d, inv, n = {}, {}, 0
        for c, ch in enumerate(str(rec.seq)):
            if ch != '-':
                n += 1
                d[n] = c; inv[c] = n
        res2col[rec.id], col2res[rec.id], chars[rec.id] = d, inv, str(rec.seq)
    return aln.get_alignment_length(), res2col, col2res, chars


def column_annotations(aln_len, res2col, chars):
    """KLIFS label, conservation and the family's modal residue, per alignment column."""
    col_klifs = {}
    for _, row in KLIFS.iterrows():
        for g in GENES:
            if pd.notna(row[g]):
                c = res2col[g].get(int(row[g]))
                if c is not None:
                    col_klifs.setdefault(c, row['KLIFS_idx'])
    cons, modal, ngap = np.zeros(aln_len), [''] * aln_len, np.zeros(aln_len, int)
    for c in range(aln_len):
        col = [chars[g][c] for g in GENES]
        aa = [x for x in col if x != '-']
        ngap[c] = len(aa)
        if aa:
            m, k = collections.Counter(aa).most_common(1)[0]
            modal[c] = m
            cons[c] = k / len(GENES)          # fraction of the nine sharing it
    return col_klifs, cons, modal, ngap


# ---------------------------------------------------------- peak finding
def build_profiles(res2col, aln_len):
    """Per-gene column profiles: the depleted-call count, and the mean depletion."""
    calls = np.zeros((len(GENES), aln_len))
    cont = np.full((len(GENES), aln_len), np.nan)
    for i, g in enumerate(GENES):
        d = res[res.Gene == g].dropna(subset=['ess_early'])
        for codon, e, q in zip(d.Codon, d.ess_early, d.FDR_early):
            c = res2col[g].get(int(codon))
            if c is None:
                continue
            cont[i, c] = -e                   # positive = depleted
            if (q < FDR_CUT) and (e < 0):
                calls[i, c] = 1
    return calls, cont


def detect(y, height=HEIGHT, prominence=PROMINENCE):
    ys = gaussian_filter1d(y, sigma=SIGMA)
    idx, _ = find_peaks(ys, height=height, distance=DISTANCE, prominence=prominence)
    return ys, idx


def continuous_profile(cont):
    """
    Mean depletion per alignment column, over the CDKs that measured it.

    The call-count profile (the resistance analysis's signal) throws away everything
    below FDR 0.05; 211 calls over 1,749 columns is thin. Averaging the score itself
    uses every residue, and a column is only counted where a CDK actually has a
    measurement, so coverage differences do not create peaks.
    """
    with np.errstate(invalid='ignore'):
        m = np.nanmean(cont, axis=0)
    return np.nan_to_num(m, nan=0.0)


def permutation_null(calls, rng):
    """Circularly shift each gene's own profile; keeps n, coverage and autocorrelation."""
    n_peaks, max_h, heights = [], [], []
    aln_len = calls.shape[1]
    for _ in range(N_PERM):
        shifted = np.vstack([np.roll(calls[i], rng.randint(aln_len))
                             for i in range(len(GENES))])
        ys, idx = detect(shifted.sum(0))
        n_peaks.append(len(idx))
        max_h.append(ys[idx].max() if len(idx) else 0.0)
        heights.extend(ys[idx])
    return np.array(n_peaks), np.array(max_h), np.array(heights)


DOMAIN_PAD = 30      # alignment columns of margin on each side of the catalytic domain


def kinase_domain_window(res2col, aln_len=None, pad=DOMAIN_PAD):
    """
    Alignment columns spanned by the catalytic domain of every CDK, plus a margin.

    The margin matters. `cdk_config`'s catalytic ranges stop at the end of the kinase fold
    proper — CDK2's is 8-149, and 149 *is* the DFG+2 alanine — so an unpadded window puts
    the activation-loop signal exactly on the boundary, where Gaussian smoothing is distorted
    by edge padding and `find_peaks` cannot see a maximum at all because there is no
    right-hand descent. 30 columns of margin clears the activation segment while still
    excluding the long non-catalytic extensions of CDK12 and CDK13, which is the point of
    restricting at all.
    """
    lo, hi = [], []
    for g in GENES:
        a, b = STRUCTURES[g]['kd']
        cols = [res2col[g][r] for r in range(a, b + 1) if r in res2col[g]]
        if cols:
            lo.append(min(cols)); hi.append(max(cols))
    a, b = min(lo) - pad, max(hi) + pad
    return max(a, 0), (min(b, aln_len - 1) if aln_len else b)


def run_peaks(y, label, rng, aln_len, calls, cont, col_klifs, col2res, cons,
              window=None, height=HEIGHT, prominence=PROMINENCE):
    """
    Detect peaks in one profile and test every peak against circular shifts.

    When `window` is given, everything happens **inside** it: the profile is cropped to the
    catalytic domain first, smoothed there, and the null shifts circularly within the same
    window. Smoothing the full-length profile and then discarding the outside is not the
    same thing — CDK12 and CDK13 carry most of their calls in their long non-catalytic
    extensions, and with sigma = 5 that signal bleeds across the domain boundary, while a
    whole-alignment null is dominated by regions the test never looks at.
    """
    lo_w, hi_w = (window if window else (0, aln_len - 1))
    sl = slice(lo_w, hi_w + 1)
    n_w = hi_w - lo_w + 1

    ys_w = gaussian_filter1d(y[sl], sigma=SIGMA)
    idx_w, _ = find_peaks(ys_w, height=height, distance=DISTANCE, prominence=prominence)
    _, _, left_w, right_w = peak_widths(ys_w, idx_w, rel_height=0.5)
    ys = np.full(aln_len, np.nan); ys[sl] = ys_w
    idx = idx_w + lo_w

    n_null, h_null = [], []
    for _ in range(N_PERM):
        if label.startswith('calls'):
            sh = np.vstack([np.roll(calls[i, sl], rng.randint(n_w)) for i in range(len(GENES))])
            yp = sh.sum(0)
        else:
            sh = np.vstack([np.roll(cont[i, sl], rng.randint(n_w)) for i in range(len(GENES))])
            with np.errstate(invalid='ignore'):
                yp = np.nan_to_num(np.nanmean(sh, axis=0), nan=0.0)
        yps = gaussian_filter1d(yp, sigma=SIGMA)
        ip, _ = find_peaks(yps, height=height, distance=DISTANCE, prominence=prominence)
        n_null.append(len(ip))
        h_null.extend(yps[ip])
    n_null, h_null = np.array(n_null), np.array(h_null)

    rows = []
    for rank, (apex, lo, hi) in enumerate(zip(idx, left_w + lo_w, right_w + lo_w), 1):
        lo, hi = int(np.floor(lo)), int(np.ceil(hi))
        cols = list(range(lo, hi + 1))
        contributing = sorted({GENES[i] for i in range(len(GENES))
                               for c in cols if calls[i, c] > 0})
        labels = [col_klifs[c] for c in cols if c in col_klifs]
        near = sorted((abs(c - apex), c) for c in cols if c in col2res[REFERENCE])
        ref_resi = col2res[REFERENCE][near[0][1]] if near else None
        rows.append(dict(
            profile=label, rank=rank, apex_col=int(apex), start_col=lo, end_col=hi,
            height=float(ys[apex]), p_height=float((h_null >= ys[apex]).mean()),
            n_cdks_called=len(contributing), cdks=';'.join(contributing),
            n_calls=int(calls[:, lo:hi + 1].sum()),
            n_cdks_measured=int((~np.isnan(cont[:, lo:hi + 1])).any(1).sum()),
            mean_depletion=float(np.nanmean(cont[:, lo:hi + 1])),
            klifs=';'.join(sorted(set(labels))) or '', in_pocket=bool(labels),
            conservation=float(np.nanmean(cons[lo:hi + 1])),
            ref_gene=REFERENCE, ref_residue=ref_resi,
            ref_residues=';'.join(str(col2res[REFERENCE][c]) for c in cols
                                  if c in col2res[REFERENCE])))
    t = pd.DataFrame(rows)
    if len(t):
        t['FDR'] = np.minimum(t.p_height * len(t) / scipy.stats.rankdata(t.p_height), 1)
    return ys, idx, t, n_null, h_null


SUBDOMAIN_NAME = {
    'I': 'β1', 'g.l': 'G-loop', 'II': 'β2', 'III': 'β3 (VAIK)', 'αC': 'αC helix',
    'b.l': 'β4–β5', 'IV': 'β4', 'V': 'β5', 'GK': 'gatekeeper', 'hinge': 'hinge',
    'linker': 'linker', 'αD': 'αD helix', 'αE': 'αE helix', 'VI': 'β6',
    'c.l': 'catalytic loop (HRD)', 'VII': 'β7', 'VIII': 'β8',
    'xDFG': 'DFG motif', 'a.l': 'activation loop',
}



DFG_BLOCK = ('xDFG', 'a.l:84', 'a.l:85')
DFG_SUBDOMAINS = ('xDFG', 'a.l')   # the same block at subdomain granularity


def within_pocket_test(cont, col_klifs, col2res, modal, rng, min_cdks=7, n_perm=20000):
    """
    Family-level depletion per pocket position, against a **within-pocket** null.

    `klifs_column_test` shifts each gene's profile around a window of the alignment. This
    null instead permutes each CDK's values **among the pocket columns only**, so a column
    is judged against the other pocket positions of the same kinases. Any background common
    to the whole window — the catalytic domain is GC-poor, so at lower correction levels it
    reads as uniformly depleted — cancels exactly. This is the test to quote.
    """
    ks = sorted(set(col_klifs.values()),
                key=lambda k: (k.split(':')[0], int(k.split(':')[1])))
    kcol = {}
    for c, k in col_klifs.items():
        kcol.setdefault(k, c)
    M = np.full((len(GENES), len(ks)), np.nan)
    for j, k in enumerate(ks):
        c = kcol[k]
        for i in range(len(GENES)):
            M[i, j] = cont[i, c]
    keep = (~np.isnan(M)).sum(0) >= min_cdks
    sub, names = M[:, keep], [k for k, m in zip(ks, keep) if m]
    obs = np.nanmean(sub, axis=0)
    null = np.empty((n_perm, int(keep.sum())))
    null_full = np.empty((n_perm, sub.shape[0], sub.shape[1]))
    for b in range(n_perm):
        pm = np.vstack([np.where(np.isnan(r), np.nan, r[rng.permutation(len(r))])
                        for r in sub])
        null_full[b] = pm
        with np.errstate(invalid='ignore'):
            null[b] = np.nanmean(pm, axis=0)
    p = ((null >= obs).sum(0) + 1) / (n_perm + 1)

    # the same permutations, aggregated by KLIFS subdomain: αC as a unit, the catalytic
    # loop as a unit, and so on. A subdomain is the level most readers think in, and
    # pooling its columns is better powered than any one of them alone.
    sd = np.array([k.split(':')[0] for k in names])
    uniq_sd = sorted(set(sd), key=lambda x: min(int(names[i].split(':')[1])
                                                for i in range(len(names)) if sd[i] == x))
    obs_sd = np.array([np.nanmean(sub[:, sd == u]) for u in uniq_sd])
    null_sd_ = np.empty((n_perm, len(uniq_sd)))
    for b in range(n_perm):
        null_sd_[b] = [np.nanmean(null_full[b][:, sd == u]) for u in uniq_sd]
    p_sd = ((null_sd_ >= obs_sd).sum(0) + 1) / (n_perm + 1)
    sdt = pd.DataFrame({'subdomain': uniq_sd, 'name': [SUBDOMAIN_NAME.get(u, u) for u in uniq_sd],
                        'n_columns': [int((sd == u).sum()) for u in uniq_sd],
                        'n_obs': [int((~np.isnan(sub[:, sd == u])).sum()) for u in uniq_sd],
                        'mean_depletion': obs_sd, 'p': p_sd})
    sdt['FDR'] = np.minimum(sdt.p * len(sdt) / scipy.stats.rankdata(sdt.p), 1)
    sdt = sdt.sort_values('mean_depletion', ascending=False).reset_index(drop=True)
    sdt['rank'] = np.arange(1, len(sdt) + 1)
    say(f'\nby KLIFS subdomain ({len(sdt)} subdomains):')
    say(sdt.head(6)[['rank', 'subdomain', 'name', 'n_columns', 'mean_depletion', 'p', 'FDR']]
        .round(4).to_string(index=False))

    t = pd.DataFrame({'klifs': names, 'modal_aa': [modal[kcol[k]] for k in names],
                      'n_cdks': (~np.isnan(sub)).sum(0),
                      'mean_depletion': obs, 'p': p})
    t = t.sort_values('mean_depletion', ascending=False).reset_index(drop=True)
    t['rank'] = np.arange(1, len(t) + 1)
    t['FDR'] = np.minimum(t.p * len(t) / scipy.stats.rankdata(t.p), 1)
    t['DFG_block'] = t.klifs.str.startswith(DFG_BLOCK)
    say(f'\nwithin-pocket test ({len(t)} columns, {n_perm:,} permutations among pocket positions):')
    say(t.head(6)[['rank', 'klifs', 'modal_aa', 'n_cdks', 'mean_depletion', 'p', 'FDR']]
        .round(4).to_string(index=False))
    sig = t[t.FDR < 0.05]
    say(f'  columns at FDR < 0.05: {len(sig)}' +
        (f"  ({', '.join(sig.klifs)})" if len(sig) else ''))
    return t, sdt


def klifs_column_test(cont, col_klifs, modal, rng, aln_len, min_cdks=7):
    """
    Family-level depletion at each pocket position, against the circular-shift null.

    Kept for continuity with the peak analysis, which uses the same null. Where the two
    disagree, `within_pocket_test` is the one to trust: its null holds the window's own
    background fixed.
    """
    obs = np.nanmean(cont, axis=0)
    null = np.empty((N_PERM, aln_len))
    for b in range(N_PERM):
        sh = np.vstack([np.roll(cont[i], rng.randint(aln_len)) for i in range(len(GENES))])
        with np.errstate(invalid='ignore'):
            null[b] = np.nan_to_num(np.nanmean(sh, axis=0), nan=0.0)
    p = ((null >= obs).sum(0) + 1) / (N_PERM + 1)
    ncov = (~np.isnan(cont)).sum(0)
    rows = [dict(aln_col=c, klifs=col_klifs[c], modal_aa=modal[c], n_cdks=int(ncov[c]),
                 mean_depletion=float(obs[c]), p=float(p[c]))
            for c in range(aln_len) if ncov[c] >= min_cdks and col_klifs.get(c)]
    t = pd.DataFrame(rows)
    t['FDR'] = np.minimum(t.p * len(t) / scipy.stats.rankdata(t.p), 1)
    t = t.sort_values('mean_depletion', ascending=False).reset_index(drop=True)
    t['rank'] = np.arange(1, len(t) + 1)
    t['DFG_block'] = t.klifs.str.startswith(DFG_BLOCK)
    K = int(t.DFG_block.sum())
    top = max(int(t.DFG_block.head(7).sum()), 1)
    hg = scipy.stats.hypergeom.sf(top - 1, len(t), K, 7)
    say(f'\nKLIFS column test (circular-shift null): {len(t)} columns with >= {min_cdks} CDKs')
    say(f'  the DFG / activation-loop block is {K} of them, and holds {top} of the top 7 '
        f'(hypergeometric p = {hg:.2g})')
    return t, hg


def main():
    rng = np.random.RandomState(SEED)
    aln_len, res2col, col2res, chars = alignment_maps()
    col_klifs, cons, modal, ngap = column_annotations(aln_len, res2col, chars)
    calls, cont = build_profiles(res2col, aln_len)
    kd = kinase_domain_window(res2col, aln_len)
    say(f'catalytic domain spans alignment columns {kd[0]}-{kd[1]}')
    say(f'alignment {aln_len} columns; '
        f'{int(calls.sum())} depleted calls (FDR < {FDR_CUT}) over {len(GENES)} CDKs')

    y_calls = calls.sum(0)
    y_cont = continuous_profile(cont)

    results = {}
    for label, prof, h, pr in [
            ('calls (FDR < 0.05)', y_calls, HEIGHT, PROMINENCE),
            ('continuous (mean depletion)', y_cont, 0.05, 0.05)]:
        ys, idx, t, n_null, h_null = run_peaks(prof, label, rng, aln_len, calls, cont,
                                               col_klifs, col2res, cons,
                                               window=kd, height=h, prominence=pr)
        p_count = (n_null >= len(idx)).mean()
        say(f'\n[{label}]  profile, smoothing and null all restricted to the '
            f'catalytic domain (columns {kd[0]}-{kd[1]})')
        say(f'  {len(idx)} peaks; null {n_null.mean():.1f} ± {n_null.std():.1f}  '
            f'p(count) = {p_count:.4f}')
        if len(t):
            say(t.sort_values('height', ascending=False)
                [['rank', 'height', 'p_height', 'FDR', 'n_cdks_measured', 'klifs',
                  'ref_residue']].head(8).round(3).to_string(index=False))
        results[label] = dict(ys=ys, idx=idx, peaks=t, n_null=n_null, h_null=h_null,
                              profile=prof, p_count=p_count)

    peaks = pd.concat([r['peaks'] for r in results.values()], ignore_index=True)
    peaks.to_csv(f'{TAB}/essentiality_peaks.csv', index=False)
    # the discrete profile is the one that beats its null; the continuous variant is
    # kept as a sensitivity check and written to the same table
    PRIMARY = 'calls (FDR < 0.05)'
    R = results[PRIMARY]
    idx, ys = R['idx'], R['ys']
    n_null, h_null = R['n_null'], R['h_null']

    cols_t, hg_p = klifs_column_test(cont, col_klifs, modal, rng, aln_len)
    cols_t.to_csv(f'{TAB}/klifs_column_tests.csv', index=False)
    wp, wsd = within_pocket_test(cont, col_klifs, col2res, modal, rng)
    wp.to_csv(f'{TAB}/klifs_within_pocket_tests.csv', index=False)
    wsd.to_csv(f'{TAB}/klifs_subdomain_tests.csv', index=False)

    # ------------------------------------------------------- orthogonal data
    ortho = orthogonal_table(res2col, cons, modal, col_klifs)
    ortho.to_csv(f'{TAB}/orthogonal_annotations.csv', index=False)
    stats = orthogonal_stats(ortho)
    stats.to_csv(f'{TAB}/qc_orthogonal_validation.csv', index=False)
    say('\northogonal checks:')
    say(stats.round(4).to_string(index=False))

    json.dump(dict(sigma=SIGMA, height=HEIGHT, distance=DISTANCE,
                   prominence=PROMINENCE, n_perm=N_PERM, fdr_cut=FDR_CUT,
                   kinase_domain_columns=list(kd), primary_profile=PRIMARY,
                   runs={k: dict(n_peaks=int(len(v['idx'])),
                                 null_mean=float(v['n_null'].mean()),
                                 null_sd=float(v['n_null'].std()),
                                 p_count=float(v['p_count']))
                         for k, v in results.items()}),
              open(f'{TAB}/peak_parameters.json', 'w'), indent=2)

    fig_peaks(aln_len, calls, cont, R['profile'], ys, idx, R['peaks'], col_klifs, kd)
    sig_peaks = R['peaks'][R['peaks'].FDR < 0.05]
    if len(sig_peaks):
        write_peak_residues(sig_peaks, col2res, col_klifs)
    write_depleted_residues(res)
    fig_klifs_significance(wp, wsd)
    fig_klifs_significance(wp, wsd, fname='peak_significance_wide.pdf',
                           figsize=(17.0, 7.0), orientation='h', n_cols=30)
    fig_klifs_significance(wp, wsd, fname='peak_significance_wide_flipped.pdf',
                           figsize=(17.0, 6.6), orientation='v', n_cols=30)
    fig_orthogonal(ortho, stats)
    say(f'\nfigures -> {FIG}')


# ------------------------------------------------------ orthogonal evidence
def peak_label(klifs_field, ref_residue):
    """A name a reader recognises, not a KLIFS index."""
    if not klifs_field:
        return (f'non-pocket ({REFERENCE} {int(ref_residue)})'
                if pd.notna(ref_residue) else 'non-pocket')
    subs, seen = [], set()
    for k in klifs_field.split(';'):
        sd = k.split(':')[0]
        if sd not in seen:
            seen.add(sd); subs.append(SUBDOMAIN_NAME.get(sd, sd))
    name = ' / '.join(subs[:3]) + ('…' if len(subs) > 3 else '')
    return f'{name} ({REFERENCE} {int(ref_residue)})' if pd.notna(ref_residue) else name


def rsa_table():
    """Relative solvent accessibility per residue, from each CDK's own structure."""
    MAX_ASA = dict(A=129, R=274, N=195, D=193, C=167, E=223, Q=225, G=104, H=224,
                   I=197, L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285,
                   Y=263, V=174)
    three = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLU='E', GLN='Q',
                 GLY='G', HIS='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F',
                 PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
    sr = ShrakeRupley()
    out = []
    for g in GENES:
        cfg = STRUCTURES[g]
        path = f'{ROOT}/notebooks/11_ChimeraX/structures/{cfg["file"]}'
        if not os.path.exists(path):
            continue
        parser = MMCIFParser(QUIET=True) if path.endswith('.cif') else PDBParser(QUIET=True)
        model = parser.get_structure(g, path)[0]
        chain = cfg['kinase_chain']
        if chain not in model:
            continue
        # SASA of the isolated kinase chain, so the cyclin does not bury its surface
        class _S:  # minimal container holding just the kinase chain
            pass
        for other in [c.id for c in model if c.id != chain]:
            model.detach_child(other)
        sr.compute(model, level='R')
        for r in model[chain]:
            aa = three.get(r.get_resname())
            if aa and r.id[0] == ' ':
                out.append(dict(Gene=g, Codon=r.id[1],
                                RSA=min(r.sasa / MAX_ASA[aa], 1.5)))
    return pd.DataFrame(out)


def orthogonal_table(res2col, cons, modal, col_klifs):
    d = res.copy()
    d['aln_col'] = [res2col[g].get(int(c)) for g, c in zip(d.Gene, d.Codon)]
    d['conservation'] = [cons[c] if c is not None else np.nan for c in d.aln_col]
    d['modal_aa'] = [modal[c] if c is not None else '' for c in d.aln_col]
    d['invariant'] = d.conservation >= 1.0
    kd = {g: STRUCTURES[g]['kd'] for g in GENES}
    d['in_domain'] = [kd[g][0] <= c <= kd[g][1] for g, c in zip(d.Gene, d.Codon)]

    dis = []
    for g in GENES:
        p = f'{DISORDER}/{g}_disorder_scores.csv'
        if os.path.exists(p):
            t = pd.read_csv(p)
            t.columns = [c.strip() for c in t.columns]
            dis.append(pd.DataFrame(dict(Gene=g, Codon=t['Index'],
                                         disorder=t['Disorder'],
                                         ppLDDT=t.get('ppLDDT'))))
    if dis:
        d = d.merge(pd.concat(dis, ignore_index=True), on=['Gene', 'Codon'], how='left')

    am = []
    for g in GENES:
        p = f'{ALPHAMISSENSE}/AlphaMissense-{g}.csv'
        if os.path.exists(p):
            t = pd.read_csv(p, comment=None)
            t.columns = [c.strip() for c in t.columns]
            col = [c for c in t.columns if 'pathogenicity score' in c]
            if col:
                a = (t.groupby('position')[col[0]].mean().reset_index()
                     .rename(columns={'position': 'Codon', col[0]: 'alphamissense'}))
                a['Gene'] = g
                am.append(a)
    if am:
        d = d.merge(pd.concat(am, ignore_index=True), on=['Gene', 'Codon'], how='left')

    r = rsa_table()
    if len(r):
        d = d.merge(r, on=['Gene', 'Codon'], how='left')
        d['buried'] = d.RSA < 0.20
    return d


def cliffs_delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) == 0 or len(b) == 0:
        return np.nan
    u = scipy.stats.mannwhitneyu(a, b).statistic
    return 2 * u / (len(a) * len(b)) - 1


def orthogonal_stats(d):
    rows = []

    def corr(col, label, direction):
        s = d.dropna(subset=['ess_early', col])
        if len(s) < 50: return
        rho = scipy.stats.spearmanr(s[col], s.ess_early)
        rows.append(dict(check=label, kind='correlation', n=len(s),
                         statistic=rho.statistic, cliffs_delta=np.nan,
                         p=rho.pvalue, expected=direction))

    def split(mask, label, direction):
        s = d.dropna(subset=['ess_early'])
        a, b = s[mask.reindex(s.index).fillna(False)], s[~mask.reindex(s.index).fillna(False)]
        if len(a) < 10 or len(b) < 10: return
        u = scipy.stats.mannwhitneyu(a.ess_early, b.ess_early)
        delta = a.ess_early.median() - b.ess_early.median()
        rows.append(dict(check=label, kind='group difference', n=len(a),
                         statistic=delta, cliffs_delta=cliffs_delta(a.ess_early, b.ess_early),
                         p=u.pvalue, expected=direction))

    corr('conservation', 'family conservation (alignment)', 'negative')
    if 'in_domain' in d:
        split(d.in_domain, 'inside the catalytic domain', 'negative')
    split(d.invariant, 'invariant alignment column', 'negative')
    split(d.in_pocket, 'ATP pocket (KLIFS)', 'negative')
    if 'disorder' in d:
        corr('disorder', 'predicted disorder (metapredict)', 'positive')
        split(d.disorder < 0.3, 'ordered (disorder < 0.3)', 'negative')
    if 'alphamissense' in d:
        corr('alphamissense', 'AlphaMissense pathogenicity', 'negative')
    if 'RSA' in d:
        corr('RSA', 'solvent accessibility', 'positive')
        split(d.get('buried'), 'buried (RSA < 0.20)', 'negative')
    t = pd.DataFrame(rows)
    # for a 'positive' expectation the sign flips
    t['passes'] = [(r.statistic > 0 if r.expected == 'positive' else r.statistic < 0)
                   and r.p < 0.05 for r in t.itertuples()]
    return t


def write_depleted_residues(res_tab):
    """Per-CDK residues that are significantly depleted, for the structural panels."""
    d = res_tab[(res_tab.FDR_early < 0.05) & (res_tab.ess_early < 0)].copy()
    out = f'{ROOT}/notebooks/11_ChimeraX/essentiality_depleted_residues.csv'
    d[['Gene', 'Codon', 'ess_early', 'z_early', 'FDR_early', 'n_guides', 'variants']].rename(
        columns={'Codon': 'position'}).sort_values(['Gene', 'position']).to_csv(out, index=False)
    say(f'\n{len(d)} significantly depleted residues -> {out}')
    say(d.groupby('Gene').size().to_string())
    return d


def write_peak_residues(sig, col2res, col_klifs):
    """
    The reference kinase's residues under each significant peak, in the
    `hit_residues.csv` schema the 11_ChimeraX scripts already read.
    """
    rows = []
    for r in sig.itertuples():
        for c in range(r.start_col, r.end_col + 1):
            resi = col2res[REFERENCE].get(c)
            if resi is None:
                continue
            rows.append(dict(gene=REFERENCE, position=resi,
                             variant=f'peak{r.rank}', variants_at_position=f'peak{r.rank}',
                             klifs=col_klifs.get(c, ''), peak_rank=r.rank,
                             peak_height=round(r.height, 3), peak_FDR=round(r.FDR, 4)))
    t = pd.DataFrame(rows).drop_duplicates(subset=['gene', 'position'])
    out = f'{ROOT}/notebooks/11_ChimeraX/essentiality_peak_residues.csv'
    t.to_csv(out, index=False)
    say(f'\n{len(t)} {REFERENCE} residues under {len(sig)} significant peak(s) -> {out}')
    return t


# -------------------------------------------------------------- figures
def fig_peaks(aln_len, calls, cont, y, ys, idx, peaks, col_klifs, kd):
    n = len(GENES)
    fig, ax = plt.subplots(n + 1, 1, figsize=(11, 7.6),
                           height_ratios=[2.6] + [1] * n, sharex=True)
    x = np.arange(aln_len)
    kl_cols = sorted(col_klifs)
    spans = []
    if kl_cols:
        start = prev = kl_cols[0]
        for c in kl_cols[1:]:
            if c > prev + 3:
                spans.append((start, prev)); start = c
            prev = c
        spans.append((start, prev))

    a0 = ax[0]
    for lo, hi in spans:
        a0.axvspan(lo, hi, color=POCKET_FILL, lw=0, zorder=0)
    a0.plot(x, y, color=MUTED, lw=.7, label='depleted residues per column')
    a0.plot(x, ys, color=INK, lw=1.4,
            label=f'smoothed (σ = {SIGMA}), catalytic domain only')
    sig = peaks[peaks.FDR < 0.05]
    a0.plot(x[idx], ys[idx], 'x', ms=5, mew=1.1, color=MUTED, label='peak')
    if len(sig):
        a0.plot(sig.apex_col, ys[sig.apex_col.values], 'x', ms=6.5, mew=1.5,
                color=DEPLETED, label='peak, FDR < 0.05')
    for r in sig.itertuples():
        a0.annotate(peak_label(r.klifs, r.ref_residue), (r.apex_col, ys[r.apex_col]),
                    textcoords='offset points', xytext=(6, 4), ha='left',
                    fontsize=6, color=DEPLETED)
    a0.set_ylabel('CDKs with a depleted\nresidue (FDR < 0.05)', fontsize=8)
    a0.legend(frameon=False, fontsize=6.5, ncol=4, loc='upper left')
    a0.spines[['top', 'right']].set_visible(False)
    a0.set_ylim(min(np.nanmin(ys), y.min()) * 1.1,
                max(np.nanmax(ys), y.max()) * 1.55)
    for a in ax:
        a.axvspan(0, kd[0], color='white', alpha=.55, lw=0, zorder=6)
        a.axvspan(kd[1], aln_len, color='white', alpha=.55, lw=0, zorder=6)
    a0.text(np.mean(kd), a0.get_ylim()[1] * .99, 'catalytic domain', ha='center',
            va='top', fontsize=6, color=INK2)

    for i, g in enumerate(GENES):
        a = ax[i + 1]
        for lo, hi in spans:
            a.axvspan(lo, hi, color=POCKET_FILL, lw=0, zorder=0)
        g_prof = gaussian_filter1d(calls[i], SIGMA)
        a.plot(x, g_prof, color=INK, lw=.9)
        a.fill_between(x, 0, g_prof, color=DEPLETED, alpha=.35, lw=0)
        a.set_yticks([]); a.set_ylim(0, None)
        a.set_ylabel(g, rotation=0, fontsize=7.5, ha='right', va='center', color=INK)
        a.spines[['top', 'right', 'left']].set_visible(False)
    ax[-1].set_xlabel('alignment column', fontsize=9)
    ax[-1].set_xlim(0, aln_len)
    fig.suptitle('Peaks in the early essentiality signal across the CDK alignment',
                 y=0.985, fontsize=10.5, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f'{FIG}/essentiality_peaks_MSA_full.pdf'); plt.close(fig)

    # the same profiles, zoomed on the window the test was actually run in
    pad = 20
    fig, ax = plt.subplots(n + 1, 1, figsize=(8.2, 7.2),
                           height_ratios=[2.6] + [1] * n, sharex=True)
    a0 = ax[0]
    for lo, hi in spans:
        a0.axvspan(lo, hi, color=POCKET_FILL, lw=0, zorder=0)
    a0.plot(x, y, color=MUTED, lw=.8)
    a0.plot(x, ys, color=INK, lw=1.5)
    a0.plot(x[idx], ys[idx], 'x', ms=5, mew=1.1, color=MUTED)
    if len(sig):
        a0.plot(sig.apex_col, ys[sig.apex_col.values], 'x', ms=7, mew=1.6, color=DEPLETED)
    for r in peaks.itertuples():
        c = DEPLETED if r.FDR < 0.05 else INK2
        a0.annotate(f'{peak_label(r.klifs, r.ref_residue)}\np = {r.p_height:.3f}',
                    (r.apex_col, ys[r.apex_col]), textcoords='offset points',
                    xytext=(0, 9), ha='center', fontsize=5.8, color=c)
    a0.set_ylabel('CDKs with a depleted\nresidue (FDR < 0.05)', fontsize=8)
    a0.set_ylim(-0.1, max(y.max(), np.nanmax(ys)) * 1.7)
    a0.spines[['top', 'right']].set_visible(False)
    for i, g in enumerate(GENES):
        a = ax[i + 1]
        for lo, hi in spans:
            a.axvspan(lo, hi, color=POCKET_FILL, lw=0, zorder=0)
        gp = gaussian_filter1d(calls[i], SIGMA)
        a.plot(x, gp, color=INK, lw=.9)
        a.fill_between(x, 0, gp, color=DEPLETED, alpha=.35, lw=0)
        a.set_yticks([]); a.set_ylim(0, None)
        a.set_ylabel(g, rotation=0, fontsize=7.5, ha='right', va='center', color=INK)
        a.spines[['top', 'right', 'left']].set_visible(False)
    # KLIFS subdomain ticks instead of raw alignment columns
    ticks, labels = [], []
    seen = set()
    for c in sorted(col_klifs):
        sd = col_klifs[c].split(':')[0]
        if sd not in seen and kd[0] - pad <= c <= kd[1] + pad:
            seen.add(sd); ticks.append(c); labels.append(sd)
    ax[-1].set_xticks(ticks); ax[-1].set_xticklabels(labels, fontsize=5.8, rotation=90)
    ax[-1].set_xlabel('KLIFS subdomain', fontsize=9)
    ax[-1].set_xlim(kd[0] - pad, kd[1] + pad)
    fig.suptitle('Essentiality peaks inside the catalytic domain', y=0.985,
                 fontsize=10.5, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f'{FIG}/essentiality_peaks_kinase_domain.pdf'); plt.close(fig)
    say('peak MSA figures written (full + kinase domain)')


def fig_klifs_significance(wp, wsd, fname='peak_significance.pdf',
                           figsize=(11.8, 7.4), orientation='h', n_cols=24):
    """
    Significance of the KLIFS column calls — the test this analysis rests on.

    The MSA peak test is not used: once the null is restricted to the same window the test
    looks at, no peak is significant (see the README). Subdomain first, because pooling a
    site's columns is better powered than any one of them and is the level readers think in;
    the per-column ranking shows what sits inside those sites.

    `orientation='h'` draws horizontal bars (labels down the side, good for a tall figure);
    `'v'` draws vertical bars with the labels on the x axis, which suits a wide one.
    """
    horiz = orientation == 'h'
    fig = plt.figure(figsize=figsize)
    if horiz:
        gs = fig.add_gridspec(1, 2, width_ratios=[1.28, 1.0], wspace=0.34,
                              left=0.135, right=0.985, top=0.855, bottom=0.085)
    else:
        gs = fig.add_gridspec(1, 2, width_ratios=[len(wsd), n_cols], wspace=0.13,
                              left=0.052, right=0.995, top=0.80, bottom=0.235)

    nsd = int((wsd.FDR < 0.05).sum())
    panels = [
        dict(ax=fig.add_subplot(gs[0]), data=wsd[::-1] if horiz else wsd,
             label=lambda r: f'{r.name}  ({r.n_columns} col)',
             colour=lambda r: DEPLETED if r.subdomain in DFG_SUBDOMAINS else MUTED,
             bold=lambda r: r.subdomain in DFG_SUBDOMAINS,
             title=f'a   by KLIFS subdomain — {nsd} clear FDR < 0.05',
             pad=2.05 if horiz else 1.72, show_p=True),
        dict(ax=fig.add_subplot(gs[1]),
             data=wp.head(n_cols)[::-1] if horiz else wp.head(n_cols),
             label=lambda r: f'{r.klifs} ({r.modal_aa})',
             colour=lambda r: DEPLETED if r.DFG_block else MUTED,
             bold=lambda r: r.DFG_block,
             title=f'b   the {len(wp)} pocket columns, top {n_cols}   * = FDR < 0.05',
             pad=1.18, show_p=False),
    ]

    for p_ in panels:
        ax, d = p_['ax'], p_['data']
        xs = np.arange(len(d))
        cols = [p_['colour'](r) for r in d.itertuples()]
        if horiz:
            ax.barh(xs, d.mean_depletion, color=cols, height=.74)
            ax.axvline(0, color=INK2, lw=.8)
            ax.set_yticks(xs)
            ax.set_yticklabels([p_['label'](r) for r in d.itertuples()],
                               fontsize=9.5 if p_['show_p'] else 8.5)
            lbls = ax.get_yticklabels()
            ax.set_xlabel('mean depletion across CDKs', fontsize=11)
            ax.tick_params(axis='x', labelsize=9.5)
            ax.set_xlim(min(0, d.mean_depletion.min() * 1.15),
                        d.mean_depletion.max() * p_['pad'])
            for y, r in zip(xs, d.itertuples()):
                if p_['show_p']:
                    mark = '**' if r.FDR < 0.01 else ('*' if r.FDR < 0.05 else '')
                    pstr = 'p < 0.001' if r.p < 0.001 else f'p = {r.p:.3f}'
                    txt = f'{mark} {pstr}'.strip() + (f'   FDR {r.FDR:.3f}'
                                                      if r.FDR < 0.05 else '')
                    ax.text(max(r.mean_depletion, 0) + .035, y, txt, va='center',
                            fontsize=8, color=INK if r.FDR < 0.05 else MUTED)
                elif r.FDR < 0.05:
                    ax.text(r.mean_depletion + .04, y, '*', va='center',
                            fontsize=14, color=INK)
        else:
            ax.bar(xs, d.mean_depletion, color=cols, width=.74)
            ax.axhline(0, color=INK2, lw=.8)
            ax.set_xticks(xs)
            ax.set_xticklabels([p_['label'](r) for r in d.itertuples()],
                               fontsize=9.5 if p_['show_p'] else 8.5,
                               rotation=90, ha='center')
            lbls = ax.get_xticklabels()
            ax.set_ylabel('mean depletion across CDKs', fontsize=11)
            ax.tick_params(axis='y', labelsize=9.5)
            top = d.mean_depletion.max()
            ax.set_ylim(min(0, d.mean_depletion.min() * 1.2), top * p_['pad'])
            for x, r in zip(xs, d.itertuples()):
                base = max(r.mean_depletion, 0) + top * 0.03
                mark = '**' if r.FDR < 0.01 else ('*' if r.FDR < 0.05 else '')
                if p_['show_p']:
                    # one rotated string per bar: a separate marker glyph collides with it
                    pstr = 'p < 0.001' if r.p < 0.001 else f'p = {r.p:.3f}'
                    ax.text(x, base, f'{mark} {pstr}'.strip(), ha='center', va='bottom',
                            rotation=90, fontsize=8,
                            color=INK if r.FDR < 0.05 else MUTED,
                            fontweight='bold' if r.FDR < 0.05 else 'normal')
                elif mark:
                    ax.text(x, base, mark, ha='center', va='bottom', fontsize=13,
                            color=INK)
        for lbl, r in zip(lbls, d.itertuples()):
            if p_['bold'](r):
                lbl.set_color(DEPLETED); lbl.set_fontweight('bold')
        ax.set_title(p_['title'], loc='left', fontsize=11.5, color=INK, pad=8)
        ax.spines[['top', 'right']].set_visible(False)

    x0 = 0.135 if horiz else 0.052
    fig.suptitle('Significance of the KLIFS column calls',
                 x=x0, ha='left', y=0.965, fontsize=13.5, color=INK)
    fig.text(x0, 0.925,
             'Null: each CDK permuted among the pocket positions, so the catalytic '
             'domain\'s own background cancels. Blue = DFG motif / activation loop.',
             fontsize=9, color=INK2, va='top')
    fig.savefig(f'{FIG}/{fname}')
    say(f'{fname} written ({orientation}, {figsize[0]}x{figsize[1]}in)')
    # returned rather than closed, so a notebook can render it inline;
    # the pdf above is still written exactly as before
    return fig


def fig_orthogonal(d, stats):
    fig = plt.figure(figsize=(9.6, 5.6))
    gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.34,
                          left=0.08, right=0.97, top=0.88, bottom=0.1)

    def by_decile(ax, col, ylab, title):
        '''Bin residues by their score and read the orthogonal property off each bin.

        Averaging ~500 residues per bin is far steadier than binning the other way
        round, and it answers the question directly: do the residues the screen calls
        most essential look more important by a measure the screen never saw?'''
        s = d.dropna(subset=['ess_early', col])
        if len(s) < 200:
            ax.axis('off'); return
        q = pd.qcut(s.ess_early, 10, labels=False, duplicates='drop')
        g = s.groupby(q)[col].agg(['mean', 'count', 'sem'])
        xs = np.arange(len(g))
        ax.errorbar(xs, g['mean'], yerr=g['sem'], fmt='o-', color=DEPLETED,
                    lw=1.5, ms=4.5, capsize=2, elinewidth=.8)
        ax.set_xticks([0, len(g) - 1])
        ax.set_xticklabels(['most\ndepleted', 'most\nenriched'], fontsize=6)
        ax.set_xlabel('residue score decile', fontsize=7)
        ax.set_ylabel(ylab, fontsize=7)
        rho = scipy.stats.spearmanr(s[col], s.ess_early)
        ax.set_title(f'{title}\nρ = {rho.statistic:+.3f}, p = {rho.pvalue:.1e}, '
                     f'n = {len(s):,}', loc='left', fontsize=7.5)
        ax.spines[['top', 'right']].set_visible(False)

    by_decile(fig.add_subplot(gs[0, 0]), 'conservation',
              'fraction of CDKs sharing\nthe residue', 'a   family conservation')
    if 'disorder' in d:
        by_decile(fig.add_subplot(gs[0, 1]), 'disorder', 'metapredict disorder',
                  'b   predicted disorder')
    if 'RSA' in d:
        by_decile(fig.add_subplot(gs[0, 2]), 'RSA', 'relative solvent accessibility',
                  'c   burial')

    ax = fig.add_subplot(gs[1, 0])
    groups = [('invariant\nacross CDKs', d.invariant.fillna(False)),
              ('ATP pocket', d.in_pocket.fillna(False))]
    if 'buried' in d:
        groups.append(('buried\n(RSA < 0.2)', d.buried.fillna(False)))
    if 'disorder' in d:
        groups.append(('disordered\n(> 0.5)', (d.disorder > 0.5).fillna(False)))
    xs, labs = [], []
    for i, (lab, m) in enumerate(groups):
        for j, (sub, col) in enumerate([(d[m], DEPLETED), (d[~m], MUTED)]):
            v = sub.ess_early.dropna()
            if len(v) < 10: continue
            pos = i + (j - .5) * .36
            ax.add_patch(plt.Rectangle((pos - .16, v.quantile(.25)), .32,
                                       v.quantile(.75) - v.quantile(.25), fc=col,
                                       ec='none', alpha=.85))
            ax.plot([pos - .16, pos + .16], [v.median()] * 2, color='white', lw=1.1, zorder=3)
        labs.append(lab)
    ax.axhline(0, color=INK2, lw=.6)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(labs, fontsize=6)
    ax.set_xlim(-.6, len(groups) - .4); ax.set_ylim(-1.3, 0.9)
    ax.set_ylabel('residue score', fontsize=7)
    ax.set_title('d   in the class vs outside it', loc='left', fontsize=7.5)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=DEPLETED, ec='none', label='in class'),
                       plt.Rectangle((0, 0), 1, 1, fc=MUTED, ec='none', label='rest')],
              frameon=False, fontsize=6, loc='lower right', ncol=2)
    ax.spines[['top', 'right']].set_visible(False)

    if 'alphamissense' in d:
        ax = fig.add_subplot(gs[1, 1])
        s = d.dropna(subset=['ess_early', 'alphamissense'])
        q = pd.qcut(s.ess_early, 10, labels=False, duplicates='drop')
        g = s.groupby(q).alphamissense.agg(['mean', 'sem'])
        ax.errorbar(np.arange(len(g)), g['mean'], yerr=g['sem'], fmt='o-',
                    color=DEPLETED, lw=1.5, ms=4.5, capsize=2, elinewidth=.8)
        ax.set_xticks([0, len(g) - 1])
        ax.set_xticklabels(['most\ndepleted', 'most\nenriched'], fontsize=6)
        ax.set_xlabel('residue score decile', fontsize=7)
        ax.set_ylabel('AlphaMissense (mean)', fontsize=7)
        rho = scipy.stats.spearmanr(s.alphamissense, s.ess_early)
        ax.set_title(f'e   AlphaMissense, {s.Gene.nunique()} CDKs\n'
                     f'ρ = {rho.statistic:+.3f}, p = {rho.pvalue:.1e}, n = {len(s):,}',
                     loc='left', fontsize=7.5)
        ax.spines[['top', 'right']].set_visible(False)

    ax = fig.add_subplot(gs[1, 2])
    t = stats.copy()[::-1]
    ax.barh(range(len(t)), -np.log10(t.p.clip(lower=1e-300)),
            color=[DEPLETED if p else MUTED for p in t.passes])
    ax.set_yticks(range(len(t))); ax.set_yticklabels(t.check, fontsize=5.6)
    ax.set_xlabel('−log10 p', fontsize=7)
    ax.set_title('f   every check, and whether it points\nthe predicted way (blue)',
                 loc='left', fontsize=7.5)
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Does the early essentiality score track protein function? '
                 '(no panel uses the screen twice)',
                 y=0.965, fontsize=10, color=INK, x=0.08, ha='left')
    fig.savefig(f'{FIG}/orthogonal_validation.pdf'); plt.close(fig)
    say('orthogonal validation figure written')


if __name__ == '__main__':
    main()
