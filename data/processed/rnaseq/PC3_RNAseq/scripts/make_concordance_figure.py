"""
Standalone: regenerate lfc_concordance_de_upset.pdf from existing DESeq2 CSVs.
Run from any env that has pandas + matplotlib + scipy.
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as _pearsonr

plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica Neue', 'Helvetica', 'DejaVu Sans']
plt.rcParams['pdf.fonttype']    = 42
plt.rcParams['ps.fonttype']     = 42

# ── Paths ──────────────────────────────────────────────────────────────────────
_HERE   = os.path.dirname(os.path.abspath(__file__))
_ROOT   = os.path.dirname(_HERE)
RESDIR  = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results')
OUTDIR  = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '03_hypomorphic')
os.makedirs(OUTDIR, exist_ok=True)

_FDR = 0.05

# ── Load DESeq2 results ────────────────────────────────────────────────────────
def _load(fname, lfc_col, padj_col):
    df = pd.read_csv(os.path.join(RESDIR, fname))
    df = df.dropna(subset=['log2FoldChange', 'padj', 'gene_name'])
    df = df.drop_duplicates(subset='gene_name')
    return df.set_index('gene_name')[['log2FoldChange', 'padj']].rename(
        columns={'log2FoldChange': lfc_col, 'padj': padj_col})

_r_g822e = _load('NT_vs_G822E_Control.csv', 'lfc_G822E', 'padj_G822E')
_r_g731k = _load('NT_vs_G731K_Control.csv', 'lfc_G731K', 'padj_G731K')
_r_i733v = _load('NT_vs_I733V_Control.csv', 'lfc_I733V', 'padj_I733V')

_conc = _r_g822e.join(_r_g731k, how='inner').join(_r_i733v, how='inner').dropna()
print(f"Genes in concordance table: {len(_conc):,}")

# ── Panel helper ───────────────────────────────────────────────────────────────
def _conc_panel(ax, df, xcol, ycol, xpadj, ypadj,
                xlabel, ylabel, dot_color, title, de_only=False):
    both = (df[xpadj] < _FDR) & (df[ypadj] < _FDR)
    one  = ((df[xpadj] < _FDR) | (df[ypadj] < _FDR)) & ~both

    if de_only:
        df   = df[both | one].copy()
        both = (df[xpadj] < _FDR) & (df[ypadj] < _FDR)
        one  = ((df[xpadj] < _FDR) | (df[ypadj] < _FDR)) & ~both
        none = pd.Series(np.zeros(len(df), dtype=bool), index=df.index)
    else:
        none = ~(both | one)

    x, y = df[xcol].values, df[ycol].values
    if len(x) == 0:
        ax.text(0.5, 0.5, 'No DE genes', transform=ax.transAxes,
                ha='center', va='center')
        return 0.0, 0

    ax.scatter(x[none.values], y[none.values], s=3,  c='#CCCCCC',  alpha=0.35,
               linewidths=0, rasterized=True)
    ax.scatter(x[one.values],  y[one.values],  s=7,  c=dot_color,  alpha=0.45,
               linewidths=0, rasterized=True)
    ax.scatter(x[both.values], y[both.values], s=14, c=dot_color,  alpha=0.85,
               linewidths=0)

    lim = max(np.abs(x).max(), np.abs(y).max()) * 1.08
    ax.plot([-lim, lim], [-lim, lim], color='#888888', lw=0.8, ls='--', zorder=1)
    ax.axhline(0, color='#CCCCCC', lw=0.5, zorder=0)
    ax.axvline(0, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

    fin  = np.isfinite(x) & np.isfinite(y)
    r_all, _ = _pearsonr(x[fin], y[fin])

    de_m = fin & (both.values | one.values)
    if de_m.sum() > 2:
        r_de, _ = _pearsonr(x[de_m], y[de_m])
        r_de_str = f'R_DE = {r_de:.3f}'
    else:
        r_de_str = 'R_DE = n/a'

    ax.text(0.04, 0.97,
            f'R_all = {r_all:.3f}\n{r_de_str}\nn(both FDR<0.05) = {both.sum():,}',
            transform=ax.transAxes, fontsize=8, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC', alpha=0.9))

    if both.sum() > 0:
        _bs = df[both].copy()
        _bs['_d'] = np.sqrt(_bs[xcol]**2 + _bs[ycol]**2)
        for gene, row in _bs.nlargest(14, '_d').iterrows():
            ax.text(row[xcol] + lim * 0.02, row[ycol], gene,
                    fontsize=6.5, color='#222222', clip_on=True, va='center')

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=9, fontweight='bold', pad=4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=8)
    return r_all, int(both.sum())

# ── UpSet set membership ───────────────────────────────────────────────────────
_sig_G822E = set(_conc.index[_conc['padj_G822E'] < _FDR])
_sig_G731K = set(_conc.index[_conc['padj_G731K'] < _FDR])
_sig_I733V = set(_conc.index[_conc['padj_I733V'] < _FDR])

_upset_cats = [
    ('G822E ∩ G731K\n(not I733V)',
     len((_sig_G822E & _sig_G731K) - _sig_I733V),   '#c0392b'),
    ('All three',
     len(_sig_G822E & _sig_G731K & _sig_I733V),      '#7f8c8d'),
    ('G822E only',
     len(_sig_G822E - _sig_G731K - _sig_I733V),      '#e74c3c'),
    ('G731K only',
     len(_sig_G731K - _sig_G822E - _sig_I733V),      '#c0392b'),
    ('G822E ∩ I733V\n(not G731K)',
     len((_sig_G822E & _sig_I733V) - _sig_G731K),   '#2980b9'),
    ('G731K ∩ I733V\n(not G822E)',
     len((_sig_G731K & _sig_I733V) - _sig_G822E),   '#2980b9'),
    ('I733V only',
     len(_sig_I733V - _sig_G822E - _sig_G731K),      '#3498db'),
]
_upset_cats.sort(key=lambda c: -c[1])

for label, count, _ in _upset_cats:
    print(f"  {label.replace(chr(10), ' ')}: {count:,}")

# ── Build figure ───────────────────────────────────────────────────────────────
_fig = plt.figure(figsize=(20, 5.2))
_gs  = _fig.add_gridspec(1, 4, width_ratios=[3.2, 3.2, 3.2, 3.0], wspace=0.36)
_ax_de = [_fig.add_subplot(_gs[i]) for i in range(3)]
_ax_up  = _fig.add_subplot(_gs[3])

_conc_panel(_ax_de[0], _conc, 'lfc_G822E', 'lfc_G731K', 'padj_G822E', 'padj_G731K',
            'log₂FC  (NT vs G822E)', 'log₂FC  (NT vs G731K)',
            '#c0392b', 'G822E vs G731K\n(both confirmed hypomorphs)\n[DE genes only]',
            de_only=True)
_conc_panel(_ax_de[1], _conc, 'lfc_G822E', 'lfc_I733V', 'padj_G822E', 'padj_I733V',
            'log₂FC  (NT vs G822E)', 'log₂FC  (NT vs I733V)',
            '#2980b9', 'G822E vs I733V\n(confirmed vs unconfirmed)\n[DE genes only]',
            de_only=True)
_conc_panel(_ax_de[2], _conc, 'lfc_G731K', 'lfc_I733V', 'padj_G731K', 'padj_I733V',
            'log₂FC  (NT vs G731K)', 'log₂FC  (NT vs I733V)',
            '#2980b9', 'G731K vs I733V\n(confirmed vs unconfirmed)\n[DE genes only]',
            de_only=True)

# UpSet bar chart
_labels = [c[0] for c in _upset_cats]
_counts = [c[1] for c in _upset_cats]
_colors = [c[2] for c in _upset_cats]
_ybars  = range(len(_upset_cats))

_bars = _ax_up.barh(_ybars, _counts, color=_colors,
                    edgecolor='white', height=0.65, zorder=3)
for _bar, _n in zip(_bars, _counts):
    _ax_up.text(_bar.get_width() + max(_counts) * 0.02,
                _bar.get_y() + _bar.get_height() / 2,
                f'{_n:,}', va='center', ha='left', fontsize=8.5, color='#333333')

_ax_up.set_yticks(list(_ybars))
_ax_up.set_yticklabels(_labels, fontsize=8)
_ax_up.set_xlabel('Number of DE genes (FDR < 0.05)', fontsize=9)
_ax_up.set_title('DE Gene Overlap\nbetween Variants', fontsize=9, fontweight='bold', pad=6)
_ax_up.set_xlim(0, max(_counts) * 1.22)
_ax_up.spines[['top', 'right']].set_visible(False)
_ax_up.tick_params(axis='y', length=0)
_ax_up.grid(axis='x', lw=0.4, alpha=0.4, zorder=0)

_fig.suptitle(
    'LFC Concordance — DE Genes Only (FDR < 0.05 in ≥ 1 comparison)  |  Control condition',
    fontsize=10.5, y=1.02)

_out_base = os.path.join(OUTDIR, 'lfc_concordance_de_upset')
plt.savefig(f'{_out_base}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{_out_base}.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved: {_out_base}.pdf / .png")
