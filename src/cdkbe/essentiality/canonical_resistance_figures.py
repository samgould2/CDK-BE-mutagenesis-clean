"""
Resistance figures built on the canonical editing calls.

Two figures, both driven by `resistance_calls_canonical.csv` (written by
CDK_resistance_patterns.ipynb, one enriched guide x residue x drug test per row):

1. MSA_resistance_peaks_* -- the smoothed enrichment profile along the CDK
   alignment with peak detection, one strip per CDK underneath.  Same figure the
   notebook draws, but the panel count follows the analysis set instead of being
   fixed at nine, so `core7` no longer renders two empty axes at the bottom.

2. clinical_alignment_* -- the same CDK profile against the clinically observed
   resistance mutations (OncoKB R1/R2) of ten other kinases, on the wider
   kinase-domain alignment from build_kinase_domain_msa.py.  On the CDK-only
   alignment a clinical mutation can only be placed if it sits on one of the 85
   numbered KLIFS positions; here every kinase has its own row, so mutations
   anywhere in the domain land on a shared column.

    python3 canonical_resistance_figures.py              # core7 and all9
    python3 canonical_resistance_figures.py --set core7
"""
import argparse, json, os, re, textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths
from Bio import AlignIO
from adjustText import adjust_text
from .. import paths

plt.rc('font', family='Helvetica')
plt.rcParams['pdf.fonttype'] = 42

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
# paths resolved from the repository root -- see _paths.py
from ._paths import (LIBRARY, KLIFS_DIR, DISORDER, ALPHAMISSENSE, DBNSFP, LFC_DIR, EDITING, FIG_DIR)  # noqa: E402
from ._paths import KLIFS_DIR, LFC_DIR, FIG_DIR  # noqa: E402  repo-root paths

GENES_ALL = ['CDK7', 'CDK8', 'CDK19', 'CDK9', 'CDK12', 'CDK13', 'CDK2', 'CDK4', 'CDK6']
EXCLUDE = {'core7': ['CDK8', 'CDK19'], 'all9': []}

CLINICAL = ['ABL1', 'ALK', 'BTK', 'EGFR', 'KIT', 'MET', 'NTRK1', 'NTRK3', 'PDGFRA', 'ROS1']

# the structure the ChimeraX peak panel is built on (notebooks/11_ChimeraX)
REFERENCE = 'CDK2'

color_dict = {'CDK7': 'cyan', 'CDK8': 'tab:blue', 'CDK19': 'cornflowerblue',
              'CDK9': 'darkviolet', 'CDK12': 'tab:green', 'CDK13': 'darkolivegreen',
              'CDK2': 'goldenrod', 'CDK4': 'tab:red', 'CDK6': 'tab:orange'}

KLIFS_ELEMENTS = [((1, 3), 'I', 'goldenrod'), ((4, 9), 'g.l', 'green'), ((10, 13), 'II', 'goldenrod'),
                  ((14, 19), 'III', 'goldenrod'), ((20, 30), 'aC', 'red'), ((31, 37), 'b.l', 'green'),
                  ((38, 41), 'IV', 'goldenrod'), ((42, 44), 'V', 'goldenrod'), ((45, 45), 'GK', 'orange'),
                  ((46, 48), 'Hinge', 'purple'), ((49, 52), 'Linker', 'blue'), ((53, 59), 'aD', 'red'),
                  ((60, 64), 'aE', 'red'), ((65, 67), 'VI', 'goldenrod'), ((68, 75), 'c.l', 'orange'),
                  ((76, 78), 'VII', 'goldenrod'), ((79, 79), 'VIII', 'goldenrod'),
                  ((80, 83), 'xDFG', 'blue'), ((84, 85), 'a.l', 'blue')]

# peak detection, as in the notebook
SIGMA, HEIGHT, DISTANCE, PROMINENCE = 5, 0.2, 10, 0.3


def find_ranges(numbers):
    """Collapse a list of ints into contiguous (start, end) ranges."""
    s = sorted(set(numbers))
    if not s:
        return []
    out, start = [], s[0]
    for i in range(1, len(s)):
        if s[i] != s[i - 1] + 1:
            out.append((start, s[i - 1]))
            start = s[i]
    out.append((start, s[-1]))
    return out


SAV_RE = re.compile(r'^([A-Z])(\d+)([A-Z*])$')


def variant_lines(sub_calls):
    """['CDK7 L317P/S, K318R, E319G, ...', ...] -- one line per CDK, genes in
    GENES_ALL order and residues in sequence order.

    Substitutions at a shared residue collapse onto one stem, the same compression
    notebooks/11_ChimeraX/label_consensus.py applies to the structural panel, so the
    two figures name a peak's variants identically.
    """
    out = []
    for gene in [g for g in GENES_ALL if g in set(sub_calls.Gene)]:
        stems = {}
        for sav in sub_calls[sub_calls.Gene == gene].SAV:
            m = SAV_RE.match(sav)
            if m is None:                      # anything unparsable is kept verbatim
                stems.setdefault((sav, 0), set())
                continue
            stems.setdefault((m.group(1), int(m.group(2))), set()).add(m.group(3))
        body = ', '.join(f'{wt}{pos}{"/".join(sorted(alts))}' if alts else f'{wt}'
                         for (wt, pos), alts in sorted(stems.items(),
                                                       key=lambda kv: kv[0][1]))
        out.append(f'{gene} {body}')
    return out


def klifs_json(gene):
    for sub in ('', 'oncokb_klifs'):
        p = os.path.join(KLIFS_DIR, sub, f'{gene}.json')
        if os.path.exists(p):
            return json.load(open(p))
    raise FileNotFoundError(f'no KLIFS json for {gene}')


# --------------------------------------------------------------- figure 1


def peak_figure(calls, genes, aln_len, klifs_ranges, xlim, fname, title,
                ref_cols=None, reference='CDK2'):
    """Smoothed enrichment along the CDK alignment, with one strip per CDK.

    Peaks are numbered P1..Pn in detection order, which is the same order
    notebooks/11_ChimeraX/generate_peaks_cxc.py uses: both run find_peaks on the
    same smoothed profile over the whole alignment, so a peak carries the same
    rank here and on the structure.  `ref_cols` is the reference kinase's
    column -> residue map; peaks it cannot place are flagged, because those are
    exactly the ones missing from the structural panel.
    """
    n = len(genes)
    fig, ax = plt.subplots(n + 1, 1, figsize=(14, 9),
                           height_ratios=[n] + [1] * n, sharex=True)
    x = np.arange(aln_len)

    up = np.zeros((n, aln_len))
    for i, gene in enumerate(genes):
        sub = calls[calls.Gene == gene]
        np.add.at(up[i], sub.aln_col.values, 1)

    for i, gene in enumerate(genes):
        y = up[i]
        ys = gaussian_filter1d(y, sigma=SIGMA)
        a = ax[i + 1]
        if y.max() > 0:
            a.plot(x, y / y.max(), color='black', alpha=0.3)
        a.plot(x, ys, color=color_dict[gene])
        a.set_yticks([])
        a.spines[['right', 'top', 'left']].set_visible(False)
        a.set_ylabel(gene, rotation=0, fontsize=13, ha='right', va='center',
                     color=color_dict[gene])
        a.set_facecolor('none')
        for lo, hi in klifs_ranges:
            a.axvspan(lo, hi, alpha=0.12, color='tab:blue', zorder=0)

    y = up.sum(axis=0)
    ys = gaussian_filter1d(y, sigma=SIGMA)
    peaks, _ = find_peaks(ys, height=HEIGHT, distance=DISTANCE, prominence=PROMINENCE)

    ax[0].plot(x, y, color='black', alpha=0.3, label='Observed')
    ax[0].plot(x, ys, color='black', alpha=0.7, label='Smoothed (combined)')
    ax[0].plot(x[peaks], ys[peaks], 'x', color='red', label='Detected peaks')

    # P1..Pn against each marker, so a peak can be followed between this figure and
    # the CDK2 panel.  A peak whose half-prominence window holds no reference
    # residue cannot be drawn on the structure and is marked with an asterisk.
    _, _, lo_w, hi_w = peak_widths(ys, peaks, rel_height=0.5)
    called = set(calls.aln_col.unique())
    off = []
    for rank, (apex, lo, hi) in enumerate(zip(peaks, lo_w, hi_w), start=1):
        cols = [c for c in range(int(np.floor(lo)), int(np.ceil(hi)) + 1) if c in called]
        placed = ref_cols is None or any(c in ref_cols for c in cols)
        if not placed:
            off.append((rank, apex, lo, hi))

    # headroom for the P-labels, without cropping the observed trace, whose
    # single-column spikes run far above the smoothed profile.  Call-outs get an
    # empty band of their own above the data rather than being laid over it.
    data_top = max(y.max(), ys.max() * 1.18) * 1.02
    top = data_top * (1 + 0.34 * len(off)) if off else data_top
    ax[0].set_ylim(-0.2, top)

    off_ranks = {r for r, *_ in off}
    for rank, apex in zip(range(1, len(peaks) + 1), peaks):
        on = rank not in off_ranks
        ax[0].annotate(f'P{rank}' + ('' if on else '*'),
                       (x[apex], ys[apex]), textcoords='offset points',
                       xytext=(0, 7), ha='center', fontsize=11, fontweight='bold',
                       color='red' if on else 'dimgrey')

    # A peak with no reference residue gets no panel on the structure, so this figure
    # is the only place its variants can be read: spell them out in a call-out.
    lo_x, hi_x = xlim if xlim else (0, aln_len)
    for i, (rank, apex, lo, hi) in enumerate(off):
        win = calls[calls.aln_col.between(int(np.floor(lo)), int(np.ceil(hi)))]
        # distinct variants, not table rows: a variant can be hit by several guides
        # and scored in several drug arms, and the structural panel counts it once
        nv = len({(g, v) for g, v in zip(win.Gene, win.SAV)})
        ng = win.Gene.nunique()
        lines = [f'P{rank}  col {apex}  {nv} variant{"s" if nv > 1 else ""} / '
                 f'{ng} CDK{"s" if ng > 1 else ""}  (no {reference} residue)']
        for ln in variant_lines(win):
            lines += textwrap.wrap(ln, width=58,
                                   subsequent_indent=' ' * (len(ln.split(' ')[0]) + 1))
        ax[0].annotate('\n'.join(lines),
                       (x[apex], ys[apex]),
                       xytext=(lo_x + 0.012 * (hi_x - lo_x),
                               top - (top - data_top) * (0.06 + 0.52 * i)),
                       textcoords='data', fontsize=8, color='black',
                       ha='left', va='top', multialignment='left',
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                                 edgecolor='dimgrey', linewidth=0.6, alpha=0.9),
                       arrowprops=dict(arrowstyle='-', color='dimgrey', linewidth=0.8,
                                       shrinkB=4))

    ax[0].set_ylabel('No. of enriched gRNAs\n(FDR < 0.01 & LFC > 0)', fontsize=13)
    handles, labels = ax[0].get_legend_handles_labels()
    if off:
        handles.append(Line2D([], [], linestyle='none'))
        labels.append(f'* no {reference} residue (gap in the reference)')
    ax[0].legend(handles, labels, fontsize=11, frameon=False)
    ax[0].spines[['right', 'top']].set_visible(False)
    ax[0].tick_params(axis='both', labelsize=12)
    for lo, hi in klifs_ranges:
        ax[0].axvspan(lo, hi, alpha=0.2, color='tab:blue')
    ax[0].set_title(title, fontsize=14)

    ax[n].set_xlabel('Alignment position', fontsize=14)
    ax[n].tick_params(axis='both', labelsize=12)
    if xlim:
        ax[0].set_xlim(*xlim)
    fig.tight_layout()
    fig.savefig(fname, transparent=True)
    # the figure is returned rather than closed so a notebook can render it
    # inline; the pdf above is still written exactly as before
    return x[peaks], ys[peaks], [r for r, *_ in off], fig


# --------------------------------------------------------------- figure 2


def load_wide_alignment():
    """res2col / col2res for the 9 CDKs + 10 clinical kinases, full-protein numbering."""
    path = os.path.join(KLIFS_DIR, 'kinase_domain_with_clinical.aln.fa')
    if not os.path.exists(path):
        raise SystemExit('run build_kinase_domain_msa.py first')
    aln = AlignIO.read(path, 'fasta')
    off = pd.read_csv(os.path.join(
        KLIFS_DIR, 'kinase_domain_with_clinical_offsets.csv')).set_index('gene')
    res2col, col2res = {}, {}
    for rec in aln:
        lo = int(off.loc[rec.id, 'domain_start'])
        d, inv, n = {}, {}, 0
        for c, ch in enumerate(str(rec.seq)):
            if ch != '-':
                d[lo + n], inv[c] = c, lo + n
                n += 1
        res2col[rec.id], col2res[rec.id] = d, inv
    return aln.get_alignment_length(), res2col, col2res


def clinical_figure(calls, genes, fname, title, tablename):
    """CDK resistance profile against OncoKB resistance in ten other kinases."""
    aln_len, res2col, col2res = load_wide_alignment()

    # CDK calls -> wide-alignment columns (a few sites fall outside the padded domain)
    cols = [res2col[g].get(int(c), -1) for g, c in zip(calls.Gene, calls.Codon)]
    calls = calls.assign(wide_col=cols)
    dropped = (calls.wide_col < 0).sum()
    calls = calls[calls.wide_col >= 0]

    # OncoKB mutations -> wide-alignment columns
    onco = pd.read_csv(os.path.join(KLIFS_DIR, 'OncoKB_resistance.csv'))
    onco = onco[onco.Gene.isin(CLINICAL)].drop_duplicates(['Gene', 'Mutation']).copy()
    onco['wide_col'] = [res2col[g].get(int(c), -1) for g, c in zip(onco.Gene, onco.Codon)]
    onco_out = (onco.wide_col < 0).sum()
    onco = onco[onco.wide_col >= 0]

    # KLIFS pocket on this coordinate, taken from every kinase so the band is complete.
    # A label is 'element:n' with n the 1-85 pocket index; entries with no residue
    # (the kinase has a gap at that pocket position) come back as None.
    pocket_cols, klifs_label, cols_of_pos = set(), {}, {}
    for g in genes + CLINICAL:
        for lbl, resi in klifs_json(g)['KLIFS2UniProtIdx'].items():
            if resi is None:
                continue
            c = res2col[g].get(int(resi))
            if c is not None:
                pocket_cols.add(c)
                klifs_label.setdefault(c, lbl)
                cols_of_pos.setdefault(int(lbl.split(':')[-1]), set()).add(c)
    klifs_ranges = find_ranges(pocket_cols)

    x = np.arange(aln_len)
    cdk_y = np.zeros(aln_len)
    np.add.at(cdk_y, calls.wide_col.values, 1)
    cdk_s = gaussian_filter1d(cdk_y, sigma=SIGMA)
    peaks, _ = find_peaks(cdk_s, height=HEIGHT, distance=DISTANCE, prominence=PROMINENCE)

    onco_y = np.zeros(aln_len)
    np.add.at(onco_y, onco.wide_col.values, 1)

    # columns where a CDK resistance variant and a clinical resistance mutation meet
    shared_cols = sorted(set(calls.wide_col) & set(onco.wide_col))

    fig, ax = plt.subplots(5, 1, figsize=(15, 10.4),
                           height_ratios=[0.3, 3, 2.4, 0.22, len(CLINICAL) * 0.42],
                           sharex=True)

    # --- KLIFS subdomain strip ---
    a = ax[0]
    for (lo, hi), label, col in KLIFS_ELEMENTS:
        cs = sorted({c for p in range(lo, hi + 1) for c in cols_of_pos.get(p, ())})
        if not cs:
            continue
        a.axvspan(min(cs) - 0.5, max(cs) + 0.5, color=col, alpha=0.75)
        # the pocket elements are only a few columns wide on this coordinate, so a
        # label is drawn only where it actually fits
        if max(cs) - min(cs) >= 8:
            a.text((min(cs) + max(cs)) / 2, 0.5, label, ha='center', va='center',
                   fontsize=8, color='white', fontweight='bold')
    a.set_yticks([])
    # the axes share an x, so the strip's ticks are hidden rather than cleared --
    # set_xticks([]) here would strip the labels off the bottom panel too
    a.tick_params(axis='x', bottom=False, labelbottom=False)
    a.spines[:].set_visible(False)
    a.set_ylabel('KLIFS\npocket', rotation=0, fontsize=9, ha='right', va='center')

    # --- CDK resistance profile ---
    a = ax[1]
    a.plot(x, cdk_y, color='black', alpha=0.3, label='Observed')
    a.plot(x, cdk_s, color='black', alpha=0.75, label='Smoothed (combined)')
    a.plot(x[peaks], cdk_s[peaks], 'x', color='red', label='Detected peaks')
    for lo, hi in klifs_ranges:
        a.axvspan(lo, hi, alpha=0.18, color='tab:blue', zorder=0)
    a.set_ylabel('CDK screen\nenriched gRNAs', fontsize=12)
    a.legend(fontsize=10, frameon=False, loc='upper right')
    a.spines[['right', 'top']].set_visible(False)
    a.set_ylim(-0.2, )

    # --- clinical resistance density, with the best-known mutations named ---
    a = ax[2]
    a.vlines(x[onco_y > 0], 0, onco_y[onco_y > 0], color='dimgrey', linewidth=1.6)
    a.plot(x[onco_y > 0], onco_y[onco_y > 0], 'o', color='dimgrey', markersize=3.5)
    for lo, hi in klifs_ranges:
        a.axvspan(lo, hi, alpha=0.18, color='tab:blue', zorder=0)
    texts = []
    for c in np.argsort(onco_y)[::-1][:9]:
        if onco_y[c] == 0:
            continue
        muts = onco[onco.wide_col == c]
        lbl = ','.join(sorted({f'{g}-{m}' for g, m in zip(muts.Gene, muts.Mutation)}))
        if len(lbl) > 34:
            lbl = f'{len(muts)} muts, {muts.Gene.nunique()} kinases'
        texts.append(a.text(c, onco_y[c], lbl, fontsize=7.5, ha='center'))
    adjust_text(texts, ax=a, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
    a.set_ylabel('OncoKB R1/R2\n(other kinases)', fontsize=12)
    a.spines[['right', 'top']].set_visible(False)
    a.set_ylim(0, onco_y.max() * 1.55)

    # --- shared columns, on their own rug so nothing is drawn over the data ---
    a = ax[3]
    a.plot(shared_cols, [0] * len(shared_cols), '|', color='firebrick',
           markersize=11, mew=1.4)
    a.set_yticks([])
    a.set_ylim(-0.5, 0.5)
    a.spines[:].set_visible(False)
    a.set_ylabel('CDK +\nclinical', rotation=0, fontsize=8, ha='right', va='center',
                 color='firebrick')

    # --- one row per clinical kinase ---
    a = ax[4]
    for i, g in enumerate(CLINICAL):
        sub = onco[onco.Gene == g]
        a.plot([0, aln_len], [i, i], color='lightgrey', lw=0.6, zorder=0)
        a.plot(sub.wide_col, [i] * len(sub), '|', color='black', markersize=8, mew=1.2)
    for lo, hi in klifs_ranges:
        a.axvspan(lo, hi, alpha=0.18, color='tab:blue', zorder=0)
    a.set_yticks(range(len(CLINICAL)))
    a.set_yticklabels(CLINICAL, fontsize=10)
    a.set_ylim(-0.7, len(CLINICAL) - 0.3)
    a.spines[['right', 'top']].set_visible(False)
    a.set_xlabel('Kinase-domain alignment position (9 CDKs + 10 clinically resistant kinases)',
                 fontsize=13)

    ax[1].legend(handles=ax[1].get_legend().legend_handles +
                 [Line2D([], [], color='firebrick', lw=0, marker='|', mew=1.4,
                         markersize=9, label='CDK + clinical at the same column')],
                 fontsize=10, frameon=False, loc='upper right')

    for a in ax:
        a.set_xlim(0, aln_len)
    ax[0].set_title(title, fontsize=14)
    fig.tight_layout()
    fig.savefig(fname, transparent=True)
    plt.close(fig)

    # --- where the two coincide ---
    rows = []
    for c in sorted(set(calls.wide_col) & set(onco.wide_col)):
        cd = calls[calls.wide_col == c]
        ok = onco[onco.wide_col == c]
        rows.append(dict(
            aln_col=c, KLIFS_idx=klifs_label.get(c, ''), in_KLIFS=c in pocket_cols,
            n_cdks=cd.Gene.nunique(),
            cdks=','.join(sorted(cd.Gene.unique())),
            cdk_variants=','.join(sorted({f'{g}-{s}' for g, s in zip(cd.Gene, cd.SAV)})),
            n_clinical_kinases=ok.Gene.nunique(),
            clinical_variants=','.join(sorted({f'{g}-{m}' for g, m in zip(ok.Gene, ok.Mutation)}))))
    shared = pd.DataFrame(rows).sort_values(['n_cdks', 'n_clinical_kinases'], ascending=False)
    shared.to_csv(tablename, index=False)

    print(f'  wide alignment {aln_len} columns; {dropped} CDK calls and {onco_out} OncoKB '
          f'mutations fell outside the padded kinase domain')
    print(f'  {len(shared)} alignment columns carry both a CDK resistance variant and an '
          f'OncoKB resistance mutation ({shared.in_KLIFS.sum()} inside the KLIFS pocket)')
    return shared


# --------------------------------------------------------------- driver


def run(aset):
    print(f'\n=== {aset} ===')
    genes = [g for g in GENES_ALL if g not in EXCLUDE[aset]]
    tdir = os.path.join(LFC_DIR, 'resistance_patterns', aset)
    fdir = os.path.join(FIG_DIR, 'resistance_patterns', aset)
    os.makedirs(fdir, exist_ok=True)

    calls = pd.read_csv(os.path.join(tdir, 'resistance_calls_canonical.csv'))
    assert set(calls.Gene) <= set(genes), 'calls table holds genes outside the analysis set'
    print(f'{len(calls):,} enriched guide x residue x drug tests, {len(genes)} CDKs')

    # KLIFS pocket on the CDK-only alignment, for the shaded band
    KL = pd.read_csv(os.path.join(KLIFS_DIR, 'KLIFS_compiled.csv'))
    aln = AlignIO.read(os.path.join(
        KLIFS_DIR, 'clustalo-I20250722-154001-0730-44590675-p1m.aln-clustal_num'), 'clustal')
    aln_len = aln.get_alignment_length()
    res2col, col2res = {}, {}
    for rec in aln:
        d, inv, n = {}, {}, 0
        for c, ch in enumerate(str(rec.seq)):
            if ch != '-':
                n += 1
                d[n], inv[c] = c, n
        res2col[rec.id], col2res[rec.id] = d, inv
    klifs_cols = {res2col[g][int(c)] for g in genes for c in KL[g].dropna()
                  if int(c) in res2col[g]}
    klifs_ranges = find_ranges(klifs_cols)

    # REFERENCE is the kinase the structural panel is drawn on, so the two figures
    # agree about which peaks could be projected
    ref_cols = col2res[REFERENCE]
    pk, _, off, _fig = peak_figure(calls, genes, aln_len, klifs_ranges, (700, 1400),
                             f'{fdir}/MSA_resistance_peaks_kinase_domain.pdf',
                             'Resistance enrichment across the CDK alignment '
                             '(kinase domain; KLIFS pocket shaded)', ref_cols)
    peak_figure(calls, genes, aln_len, klifs_ranges, None,
                f'{fdir}/MSA_resistance_peaks_full.pdf',
                'Resistance enrichment across the full CDK alignment', ref_cols)
    print(f'  {len(pk)} peaks: ' +
          ', '.join(f'P{i}=col {c}' + ('*' if i in off else '')
                    for i, c in enumerate(pk, start=1)))
    if off:
        print(f'  * peaks {off} have no {REFERENCE} residue and are absent from the '
              f'structural panel')

    clinical_figure(calls, genes,
                    f'{fdir}/clinical_alignment_comparison.pdf',
                    'CDK base-editing resistance vs clinically observed kinase resistance '
                    '(OncoKB R1/R2)',
                    os.path.join(tdir, 'clinical_alignment_shared_columns.csv'))
    print(f'  wrote figures to {paths.rel(fdir)}')
    # the kinase-domain panel is supplementary 15d; returned so a
    # notebook can render it inline
    return _fig


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--set', dest='sets', action='append', choices=['core7', 'all9'])
    a = ap.parse_args()
    for s in a.sets or ['core7', 'all9']:
        run(s)
