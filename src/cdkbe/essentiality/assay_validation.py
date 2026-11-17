#!/usr/bin/env python3
"""
One figure asking whether this assay measures residue essentiality at all.

Seven panels, each a different kind of evidence, none of which can be satisfied by the
representation artifact the corrections remove:

  a  specificity      where every class of guide lands, controls included. Truncating
                      variants are the only class that leaves the non-targeting band.
  b  reliability      split-half: the guides hitting each residue are split in two at
                      random and scored separately. This is the assay's own ceiling —
                      no external data, no annotation, just: do independent guides at the
                      same residue agree?
  c  reproducibility  that split-half (Spearman-Brown corrected), the two editors — which
                      make *different* substitutions at the same residue — and an
                      independent infection of the same library.
  d  artifact removed what the corrections do to the two readouts that matter most.
  e  orthogonality    AlphaMissense and predicted disorder across score deciles.
  f  3D coherence     depleted residues are closer together in the folded protein than
                      random residues of the same count, per CDK, against 1,000 draws.
  g  landmark         the one family-wide peak, against its circular-shift null.

Usage:  python3 assay_validation.py
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio.PDB import MMCIFParser, PDBParser

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)

from .essentiality_figures import (CMAP, NORM, DEPLETED, ENRICHED, NO_DATA, INK, INK2,
                                  MUTED, POCKET_FILL, GENES, res, KLIFS, TAB, FIG,
                                  CLASSES, CLASS_COLOR, guide)
from .._chimerax_config import STRUCTURES, CDK_COLORS

N_SPLIT, N_PERM3D, SEED = 200, 1000, 0
say = lambda *a: print(*a, flush=True)

rows = pd.read_csv(f'{TAB}/residue_guide_rows.csv').rename(columns={'%Reads': 'pct'})
ladder = pd.read_csv(f'{TAB}/correction_ladder.csv')
peaks = pd.read_csv(f'{TAB}/essentiality_peaks.csv')
ortho = pd.read_csv(f'{TAB}/orthogonal_annotations.csv')
hold = pd.read_csv(f'{TAB}/residue_essentiality_SY5609.csv')
qc_rep = pd.read_csv(f'{TAB}/qc_reproducibility.csv')


# ---------------------------------------------------------- split-half
def split_half():
    """Score each residue twice, from disjoint halves of its guides."""
    d = rows.dropna(subset=['ess_early']).copy()
    d['key'] = d.Gene + '_' + d.Codon.astype(str)
    n = d.groupby('key').size()
    keep = set(n[n >= 2].index)
    d = d[d.key.isin(keep)]
    rng = np.random.RandomState(SEED)
    rs, pairs = [], None
    for it in range(N_SPLIT):
        half = rng.rand(len(d)) < .5
        a = d[half].groupby('key').apply(
            lambda x: np.average(x.ess_early, weights=x.pct))
        b = d[~half].groupby('key').apply(
            lambda x: np.average(x.ess_early, weights=x.pct))
        j = pd.concat([a.rename('a'), b.rename('b')], axis=1).dropna()
        if len(j) > 50:
            rs.append(scipy.stats.pearsonr(j.a, j.b).statistic)
            if pairs is None:
                pairs = j
    r = float(np.mean(rs))
    sb = 2 * r / (1 + r)                       # Spearman-Brown, half -> full length
    say(f'split-half over {N_SPLIT} draws: r = {r:.3f} (sd {np.std(rs):.3f}), '
        f'Spearman-Brown = {sb:.3f}, residues = {len(keep):,}')
    return r, sb, pairs, len(keep)


# ------------------------------------------------------ 3D coherence
def ca_coords():
    """CA coordinates per CDK, restricted to the confidently modelled range."""
    out = {}
    for g in GENES:
        cfg = STRUCTURES[g]
        path = f'{ROOT}/notebooks/11_ChimeraX/structures/{cfg["file"]}'
        if not os.path.exists(path):
            continue
        par = MMCIFParser(QUIET=True) if path.endswith('.cif') else PDBParser(QUIET=True)
        ch = par.get_structure(g, path)[0][cfg['kinase_chain']]
        dr = cfg.get('display_range')      # CDK19 is an AlphaFold model whose tail is
        out[g] = {r.id[1]: r['CA'].coord for r in ch   # disordered spaghetti; exclude it
                  if r.id[0] == ' ' and 'CA' in r
                  and (dr is None or dr[0] <= r.id[1] <= dr[1])}
    return out


def morans_i(coords, cutoff=10.0):
    """
    Spatial autocorrelation of the score itself, per CDK.

    `clustering_3d` only looks at the significant tail. Moran's I uses every measured
    residue, so it asks the better-powered question: do residues that sit near each other
    in the folded protein get similar scores? Weights are 1 for CA pairs within `cutoff`
    angstroms, 0 otherwise; the null shuffles the scores over the same residues, which
    keeps the score distribution and the structure and breaks only their pairing.
    """
    rng = np.random.RandomState(SEED)
    out = []
    for g in GENES:
        if g not in coords:
            continue
        d = res[(res.Gene == g) & res.ess_early.notna()]
        d = d[d.Codon.isin(coords[g])]
        if len(d) < 40:
            continue
        P = np.array([coords[g][c] for c in d.Codon])
        D = np.linalg.norm(P[:, None] - P[None], axis=-1)
        W = ((D > 0) & (D <= cutoff)).astype(float)
        if W.sum() == 0:
            continue
        x = d.ess_early.values

        def I(v):
            z = v - v.mean()
            return len(v) / W.sum() * (z @ W @ z) / (z @ z)

        obs = I(x)
        null = np.array([I(rng.permutation(x)) for _ in range(N_PERM3D)])
        out.append(dict(Gene=g, n=len(d), I=obs, null_mean=null.mean(),
                        null_sd=null.std(), z=(obs - null.mean()) / null.std(),
                        p=((null >= obs).sum() + 1) / (N_PERM3D + 1)))
    t = pd.DataFrame(out)
    say(f'\nMoran\'s I of the residue score (CA within {cutoff:.0f} A):')
    say(t.round(3).to_string(index=False))
    return t


def clustering_3d(coords):
    """Are the depleted residues closer together than random residues of the same count?"""
    rng = np.random.RandomState(SEED)
    out = []
    for g in GENES:
        if g not in coords:
            continue
        d = res[(res.Gene == g) & res.ess_early.notna()]
        d = d[d.Codon.isin(coords[g])]
        hit = d[(d.FDR_early < 0.05) & (d.ess_early < 0)].Codon.values
        pool = d.Codon.values
        if len(hit) < 4 or len(pool) < 4 * len(hit):
            continue
        def meand(cs):
            P = np.array([coords[g][c] for c in cs])
            D = np.linalg.norm(P[:, None] - P[None], axis=-1)
            return D[np.triu_indices(len(cs), 1)].mean()
        obs = meand(hit)
        null = np.array([meand(rng.choice(pool, len(hit), replace=False))
                         for _ in range(N_PERM3D)])
        out.append(dict(Gene=g, n_hits=len(hit), n_pool=len(pool), observed=obs,
                        null_mean=null.mean(), null_sd=null.std(),
                        z=(obs - null.mean()) / null.std(),
                        p=((null <= obs).sum() + 1) / (N_PERM3D + 1)))
    t = pd.DataFrame(out)
    say('\n3D clustering of depleted residues:')
    say(t.round(3).to_string(index=False))
    return t


# ------------------------------------------------------------- figure
def main():
    r_half, sb, pairs, n_res = split_half()
    coords = ca_coords()
    clus = clustering_3d(coords)
    clus.to_csv(f'{TAB}/qc_3d_clustering.csv', index=False)
    moran = morans_i(coords)
    moran.to_csv(f'{TAB}/qc_spatial_autocorrelation.csv', index=False)

    fig = plt.figure(figsize=(12.2, 7.2))
    gs = fig.add_gridspec(2, 4, hspace=0.58, wspace=0.42,
                          left=0.065, right=0.985, top=0.875, bottom=0.085)

    # ---- a  specificity
    ax = fig.add_subplot(gs[0, 0])
    gd = guide.dropna(subset=['GuideClass'])
    band = gd[gd.GuideClass == 'Non-targeting'].ess_early.dropna()
    ax.axhspan(band.quantile(.25), band.quantile(.75), color=MUTED, alpha=.22, lw=0)
    ax.axhline(0, color=INK2, lw=.7)
    data = [gd[gd.GuideClass == k].ess_early.dropna().values for k in CLASSES]
    bp = ax.boxplot(data, positions=range(len(CLASSES)), widths=.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color='white', lw=1.1),
                    whiskerprops=dict(color=INK2, lw=.7), capprops=dict(color=INK2, lw=.7))
    for b, k in zip(bp['boxes'], CLASSES):
        b.set_facecolor(CLASS_COLOR[k]); b.set_edgecolor('none')
    for i, k in enumerate(CLASSES):
        ax.text(i, 2.35, f'{len(data[i]):,}', ha='center', fontsize=5.2, color=INK2)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels([k.replace('Non-targeting', 'Non-targ.') for k in CLASSES],
                       fontsize=5.8, rotation=40, ha='right')
    ax.set_ylim(-3.2, 2.6)
    ax.set_ylabel('guide score', fontsize=8)
    tr = gd[gd.GuideClass == 'Truncating'].ess_early.dropna()
    pv = scipy.stats.mannwhitneyu(tr, band).pvalue
    ax.text(len(CLASSES) - 1, -3.0, f'median {tr.median() - band.median():+.2f}\n'
            f'p = {pv:.1e}', ha='right', fontsize=5.6, color=INK)
    ax.set_title('a   truncating guides are the only\n      class shifted off the controls',
                 loc='left', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- b  split-half
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(pairs.a, pairs.b, s=3, color=DEPLETED, alpha=.25, linewidths=0,
               rasterized=True)
    lim = (-4.5, 3)
    ax.plot(lim, lim, color=MUTED, lw=.7, ls=(0, (2, 2)))
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel('score from half the guides', fontsize=8)
    ax.set_ylabel('score from the other half', fontsize=8)
    ax.set_title(f'b   independent guides agree\n      r = {r_half:.2f}  '
                 f'(n = {n_res:,} residues)', loc='left', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- c  reproducibility bars
    ax = fig.add_subplot(gs[0, 2])
    def rep(name):
        m = qc_rep[qc_rep.comparison == name]
        return float(m.R.iloc[0]) if len(m) else np.nan
    bars = [('split-half\n(Spearman–Brown)', sb),
            ('ABE vs CBE\n(different variants)', rep('early: ABE vs CBE')),
            ('independent\ninfection', rep('early: primary vs SY-5609'))]
    ax.bar(range(len(bars)), [b for _, b in bars], color=DEPLETED, width=.6)
    for i, (_, b) in enumerate(bars):
        ax.text(i, b + .02, f'{b:.2f}', ha='center', fontsize=7, color=INK)
    ax.set_xticks(range(len(bars)))
    ax.set_xticklabels([n for n, _ in bars], fontsize=5.8)
    ax.set_ylim(0, 1); ax.set_ylabel('Pearson r', fontsize=8)
    ax.set_title('c   the same residue, measured\n      independently', loc='left',
                 fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- d  what the corrections do
    ax = fig.add_subplot(gs[0, 3])
    x = np.arange(len(ladder))
    ax.plot(x, ladder.r_GC_max, '-o', color=ENRICHED, lw=1.5, ms=5,
            label='worst |r| with GC')
    ax.plot(x, ladder.truncating_minus_NT, '-o', color=DEPLETED, lw=1.5, ms=5,
            label='truncating − control')
    ax.axhline(0, color=INK2, lw=.7)
    ax.set_xticks(x); ax.set_xticklabels([l.split()[0] for l in ladder.level], fontsize=6.5)
    ax.set_ylabel('value', fontsize=8)
    ax.set_title('d   artifact out, positive\n      control in', loc='left', fontsize=8)
    ax.legend(frameon=False, fontsize=5.8, loc='upper right')
    ax.spines[['top', 'right']].set_visible(False)

    # ---- e  orthogonal measures by score decile
    ax = fig.add_subplot(gs[1, 0])
    for col, c, lab in [('alphamissense', DEPLETED, 'AlphaMissense'),
                        ('disorder', ENRICHED, 'disorder (metapredict)')]:
        if col not in ortho:
            continue
        s = ortho.dropna(subset=['ess_early', col])
        q = pd.qcut(s.ess_early, 10, labels=False, duplicates='drop')
        g = s.groupby(q)[col].agg(['mean', 'sem'])
        z = (g['mean'] - g['mean'].mean()) / g['mean'].std()
        ax.errorbar(np.arange(len(g)), z, yerr=g['sem'] / g['mean'].std(), fmt='-o',
                    color=c, lw=1.4, ms=4, capsize=2, elinewidth=.7, label=lab)
    ax.axhline(0, color=INK2, lw=.7)
    ax.set_xticks([0, 9]); ax.set_xticklabels(['most\ndepleted', 'most\nenriched'], fontsize=6)
    ax.set_xlabel('residue score decile', fontsize=8)
    ax.set_ylabel('z-scored mean', fontsize=8)
    ax.set_title('e   two measures that never\n      saw the screen', loc='left', fontsize=8)
    ax.legend(frameon=False, fontsize=5.8)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- f  spatial autocorrelation
    ax = fig.add_subplot(gs[1, 1])
    if len(moran):
        c = moran.sort_values('z', ascending=False)
        ys = np.arange(len(c))[::-1]
        ax.barh(ys, c.z, color=[CDK_COLORS[g] for g in c.Gene], height=.6)
        ax.axvline(0, color=INK2, lw=.7)
        for y, r in zip(ys, c.itertuples()):
            if r.p < 0.05:
                ax.text(r.z + .12, y, '*', va='center', fontsize=9, color=INK)
        ax.set_yticks(ys); ax.set_yticklabels(c.Gene, fontsize=6.5)
        for lbl, g in zip(ax.get_yticklabels(), c.Gene):
            lbl.set_color(CDK_COLORS[g])
        ax.set_xlabel('z, Moran\'s I vs 1,000 shuffles', fontsize=7.5)
        n_sig = int((c.p < 0.05).sum())
        ax.set_title(f'f   neighbouring residues get similar\n      scores '
                     f'({n_sig}/{len(c)} CDKs, p < 0.05)', loc='left', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- g  the landmark peak
    ax = fig.add_subplot(gs[1, 2])
    pk = peaks[peaks.profile == 'calls (FDR < 0.05)'].sort_values('height', ascending=False)
    ys = np.arange(len(pk))[::-1]
    ax.barh(ys, pk.height, color=[DEPLETED if q < .05 else MUTED for q in pk.FDR], height=.6)
    for y, r in zip(ys, pk.itertuples()):
        ax.text(r.height + .01, y, f'p={r.p_height:.3f}', va='center', fontsize=5.6,
                color=INK2)
    ax.set_yticks(ys)
    from essentiality_peaks import peak_label
    ax.set_yticklabels([peak_label(r.klifs if isinstance(r.klifs, str) else '',
                                   r.ref_residue).replace(' / ', ' /\n')
                        for r in pk.itertuples()], fontsize=5.6)
    ax.set_xlabel('peak height', fontsize=8)
    ax.set_xlim(0, pk.height.max() * 1.45)
    ax.set_title('g   the family-wide peak lands on\n      the catalytic motif',
                 loc='left', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    # ---- h  per-CDK pocket contrast
    ax = fig.add_subplot(gs[1, 3])
    pc = pd.read_csv(f'{TAB}/qc_pocket_contrast.csv')
    pc = pc[pc.contrast == 'early'].set_index('Gene').reindex(GENES).dropna(how='all')
    ys = np.arange(len(pc))[::-1]
    ax.barh(ys, pc.delta, color=[CDK_COLORS[g] for g in pc.index], height=.6)
    for y, g in zip(ys, pc.index):
        if pc.loc[g, 'p'] < 0.05:
            ax.text(pc.loc[g, 'delta'] - .02, y, '*', ha='right', va='center',
                    fontsize=9, color=INK)
    ax.axvline(0, color=INK2, lw=.7)
    ax.set_yticks(ys); ax.set_yticklabels(pc.index, fontsize=6.5)
    for lbl, g in zip(ax.get_yticklabels(), pc.index):
        lbl.set_color(CDK_COLORS[g])
    ax.set_xlabel('ATP pocket − rest of protein', fontsize=7.5)
    ax.set_title('h   the pocket is more essential\n      than the rest', loc='left',
                 fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    fig.suptitle('Does this assay measure residue essentiality?',
                 x=0.065, ha='left', y=0.965, fontsize=12, color=INK)
    fig.text(0.065, 0.925,
             'Eight independent lines of evidence. None of them can be produced by the '
             'plasmid-to-genome representation artifact, which is sequence-driven and '
             'carries no residue-level structure.',
             fontsize=7, color=INK2)
    fig.savefig(f'{FIG}/assay_validation.pdf')
    plt.close(fig)
    say(f'\nassay_validation.pdf -> {FIG}')


if __name__ == '__main__':
    main()
