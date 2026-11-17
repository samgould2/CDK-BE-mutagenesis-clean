"""
Shared hypomorphic CDK12 transcriptional program (G822E ∩ G731K, not I733V)
Three analyses:
  1. ORA (Enrichr) — up-regulated vs down-regulated gene subsets
  2. Gene length CDF — shared hypomorphs vs I733V vs background
  3. Pre-ranked GSEA — all genes ranked by mean(lfc_G822E, lfc_G731K)

Inputs:
  pc3_rnaseq_results/deseq2_results/shared_hypomorphic_DEGs.csv
  pc3_rnaseq_results/deseq2_results/NT_vs_G822E_Control.csv
  pc3_rnaseq_results/deseq2_results/NT_vs_G731K_Control.csv
  pc3_rnaseq_results/deseq2_results/NT_vs_I733V_Control.csv

Requires: pandas numpy matplotlib scipy gseapy  (pip install gseapy if missing)
"""

import os
import urllib.parse
import urllib.request
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Liberation Sans', 'Arial', 'Helvetica Neue',
                          'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype':      42,
    'ps.fonttype':       42,
    'axes.linewidth':    0.8,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
})

# ── Paths ──────────────────────────────────────────────────────────────────────
_HERE  = os.path.dirname(os.path.abspath(__file__))
_ROOT  = os.path.dirname(_HERE)
RESDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results')
OUTDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '03_hypomorphic')
os.makedirs(OUTDIR, exist_ok=True)

FDR = 0.05

# ── Load data ─────────────────────────────────────────────────────────────────
shared = pd.read_csv(os.path.join(RESDIR, 'shared_hypomorphic_DEGs.csv'), index_col=0)
print(f"Shared hypomorphic DEGs: {len(shared):,}")

def _load_lfc(fname):
    df = pd.read_csv(os.path.join(RESDIR, fname))
    df = df.dropna(subset=['log2FoldChange', 'padj', 'gene_name'])
    df = df.drop_duplicates('gene_name').set_index('gene_name')
    return df[['log2FoldChange', 'padj']]

_g822e = _load_lfc('NT_vs_G822E_Control.csv')
_g731k = _load_lfc('NT_vs_G731K_Control.csv')
_i733v = _load_lfc('NT_vs_I733V_Control.csv')

# All genes with both hypomorph comparisons (for prerank)
_all = _g822e.join(_g731k, how='inner', lsuffix='_G822E', rsuffix='_G731K').dropna()
_all['lfc_mean']  = (_all['log2FoldChange_G822E'] + _all['log2FoldChange_G731K']) / 2
_all['padj_G822E'] = _all['padj_G822E']
_all['padj_G731K'] = _all['padj_G731K']

# I733V significant genes (for gene length grouping)
_sig_i733v = set(_i733v.index[_i733v['padj'] < FDR])

print(f"All genes (both hypomorph comparisons): {len(_all):,}")

# ── Gene sub-lists ────────────────────────────────────────────────────────────
up_genes = shared.index[shared['lfc_mean_hypomorphs'] > 0].tolist()
dn_genes = shared.index[shared['lfc_mean_hypomorphs'] < 0].tolist()
print(f"  Up: {len(up_genes):,}   Down: {len(dn_genes):,}")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 — ORA via Enrichr
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[1/3] Running ORA (Enrichr)...")

try:
    import gseapy as gp
    _HAS_GSEAPY = True
except ImportError:
    print("  gseapy not found — skipping ORA and prerank. pip install gseapy")
    _HAS_GSEAPY = False

ORA_SETS   = ['MSigDB_Hallmark_2020', 'Reactome_2022', 'GO_Biological_Process_2023']
TOP_N_TERMS = 15

def _run_enrichr(gene_list, label):
    if not gene_list:
        return pd.DataFrame()
    try:
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=ORA_SETS,
                         organism='human',
                         outdir=None,
                         verbose=False)
        df = enr.results.copy()
        df = df[df['Adjusted P-value'] < 0.05].copy()
        # Shorten term names
        df['Term'] = (df['Term']
                      .str.replace(r'_\d{4}.*', '', regex=True)
                      .str.replace('_', ' ', regex=False)
                      .str.replace('HALLMARK ', '', regex=False)
                      .str.replace('GOBP ', '', regex=False))
        df['-log10q'] = -np.log10(df['Adjusted P-value'].clip(lower=1e-30))
        # Overlap count
        df['overlap_n'] = df['Overlap'].apply(lambda s: int(s.split('/')[0]))
        df['gene_set']  = df['Gene_set']
        print(f"  {label}: {len(df)} significant terms")
        return df
    except Exception as e:
        print(f"  Enrichr failed ({e})")
        return pd.DataFrame()

if _HAS_GSEAPY:
    _enr_up = _run_enrichr(up_genes, 'Up')
    _enr_dn = _run_enrichr(dn_genes, 'Down')

    def _top_terms(df, n=TOP_N_TERMS):
        if df.empty:
            return df
        return (df.sort_values('-log10q', ascending=False)
                  .drop_duplicates('Term')
                  .head(n)
                  .sort_values('-log10q'))

    _top_up = _top_terms(_enr_up)
    _top_dn = _top_terms(_enr_dn)

    # ── ORA figure ──────────────────────────────────────────────────────────
    _n_up = len(_top_up)
    _n_dn = len(_top_dn)

    if _n_up > 0 or _n_dn > 0:
        _max_n  = max(_n_up, _n_dn, 5)
        _fig_h  = max(5.5, _max_n * 0.38 + 1.8)
        _fig, (_ax_up, _ax_dn) = plt.subplots(1, 2, figsize=(14, _fig_h))

        _DOT_SCALE = 8  # pt² per gene

        for ax, df, label, color in [
            (_ax_up, _top_up, f'Up-regulated ({len(up_genes):,} genes)', '#c0392b'),
            (_ax_dn, _top_dn, f'Down-regulated ({len(dn_genes):,} genes)', '#2980b9'),
        ]:
            if df.empty:
                ax.text(0.5, 0.5, 'No significant terms', transform=ax.transAxes,
                        ha='center', va='center', fontsize=10)
            else:
                ys  = range(len(df))
                sc  = ax.scatter(df['-log10q'].values, list(ys),
                                 s=df['overlap_n'].values * _DOT_SCALE,
                                 c=color, alpha=0.8, linewidths=0.5,
                                 edgecolors='white', zorder=3)
                ax.set_yticks(list(ys))
                ax.set_yticklabels(df['Term'].tolist(), fontsize=8.5)
                ax.axvline(-np.log10(0.05), color='#AAAAAA', lw=0.8,
                           ls='--', zorder=1, label='FDR = 0.05')

                # Dot size legend
                _sizes = [10, 50, 100]
                for _sz in _sizes:
                    ax.scatter([], [], s=_sz * _DOT_SCALE, c='#AAAAAA',
                               label=f'n = {_sz}', alpha=0.7, linewidths=0)
                ax.legend(frameon=False, fontsize=7.5, loc='lower right',
                          handletextpad=0.3, labelspacing=0.4)

            ax.set_xlabel('-log₁₀(adj. p-value)', fontsize=10)
            ax.set_title(label, fontsize=10, fontweight='bold', pad=8, color=color)
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(labelsize=8.5, axis='x')
            ax.tick_params(length=3, axis='x')
            ax.tick_params(length=0, axis='y')
            ax.xaxis.grid(True, lw=0.35, color='#EEEEEE', zorder=0)
            ax.set_axisbelow(True)

        _fig.suptitle(
            'ORA — Shared CDK12 Hypomorphic Program (G822E ∩ G731K, not I733V)',
            fontsize=11, fontweight='bold', y=1.01)
        plt.tight_layout()
        _out = os.path.join(OUTDIR, 'hypomorph_program_ora')
        plt.savefig(f'{_out}.pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f'{_out}.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {_out}.pdf")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 — Gene length CDF
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[2/3] Fetching gene lengths from Ensembl BioMart...")

_GL_CACHE = os.path.join(RESDIR, '_gene_lengths_cache.csv')

def _fetch_biomart_lengths(gene_symbols, chunk=400):
    """Query Ensembl BioMart for max transcript length per HGNC symbol."""
    all_rows = []
    unique_genes = list(set(gene_symbols))
    for i in range(0, len(unique_genes), chunk):
        batch = unique_genes[i:i + chunk]
        val_str = ','.join(batch)
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<!DOCTYPE Query>'
            '<Query virtualSchemaName="default" formatter="TSV" header="1" '
            'uniqueRows="0" count="">'
            '<Dataset name="hsapiens_gene_ensembl" interface="default">'
            f'<Filter name="hgnc_symbol" value="{val_str}"/>'
            '<Attribute name="hgnc_symbol"/>'
            '<Attribute name="transcript_length"/>'
            '</Dataset></Query>'
        )
        url = ('https://www.ensembl.org/biomart/martservice?query='
               + urllib.parse.quote(xml))
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                raw = resp.read().decode('utf-8')
            rows = [r.split('\t') for r in raw.strip().split('\n')[1:] if r]
            all_rows.extend(rows)
            print(f"  BioMart: fetched {i+len(batch)}/{len(unique_genes)} genes", end='\r')
        except Exception as e:
            print(f"\n  BioMart chunk failed: {e}")
    print()
    if not all_rows:
        return pd.Series(dtype=float)
    df = pd.DataFrame(all_rows, columns=['gene_name', 'tx_length'])
    df['tx_length'] = pd.to_numeric(df['tx_length'], errors='coerce')
    return df.groupby('gene_name')['tx_length'].max()

# Load or fetch gene lengths
all_gene_symbols = _all.index.tolist()
if os.path.exists(_GL_CACHE):
    print("  Using cached gene lengths.")
    _gl = pd.read_csv(_GL_CACHE, index_col=0)['tx_length']
else:
    _gl = _fetch_biomart_lengths(all_gene_symbols)
    if len(_gl) > 0:
        _gl.to_frame('tx_length').to_csv(_GL_CACHE)
        print(f"  Cached {len(_gl):,} gene lengths → {_GL_CACHE}")

if len(_gl) > 100:
    # Assign group membership
    _shared_set  = set(shared.index)
    _i733v_only  = _sig_i733v - _shared_set
    _g822e_only  = set(_all.index[(_all['padj_G822E'] < FDR) &
                                   (_all['padj_G731K'] >= FDR)]) - _sig_i733v
    _g731k_only  = set(_all.index[(_all['padj_G731K'] < FDR) &
                                   (_all['padj_G822E'] >= FDR)]) - _sig_i733v
    _background  = set(_all.index) - _shared_set - _i733v_only - _g822e_only - _g731k_only

    GROUPS = {
        'Shared hypomorphs\n(G822E ∩ G731K)': (_shared_set, '#c0392b'),
        'G822E only':                          (_g822e_only, '#e74c3c'),
        'G731K only':                          (_g731k_only, '#8e0152'),
        'I733V only':                          (_i733v_only, '#2980b9'),
        'Background\n(not DE in any)':         (_background, '#AAAAAA'),
    }

    from scipy.stats import mannwhitneyu

    _fig_gl, _ax_gl = plt.subplots(figsize=(6.5, 4.8))

    _ref_lengths = None
    for label, (gene_set, color) in GROUPS.items():
        lengths = _gl.reindex(list(gene_set)).dropna().values
        if len(lengths) < 5:
            continue
        lengths_log = np.log10(lengths + 1)
        xs = np.sort(lengths_log)
        ys = np.arange(1, len(xs) + 1) / len(xs)
        lw = 2.2 if 'Shared' in label else 1.5
        ls = '-'  if 'Shared' in label or 'Background' in label else '--'
        _ax_gl.plot(xs, ys, color=color, lw=lw, ls=ls,
                    label=f'{label} (n={len(lengths):,})', alpha=0.9)
        if 'Background' in label:
            _ref_lengths = lengths

    # Mann-Whitney vs background for shared hypomorphs
    if _ref_lengths is not None:
        _sh_len = _gl.reindex(list(_shared_set)).dropna().values
        _stat, _pval = mannwhitneyu(_sh_len, _ref_lengths, alternative='greater')
        _ax_gl.text(0.97, 0.06,
                    f'Shared vs background\nMann-Whitney p = {_pval:.2e}',
                    transform=_ax_gl.transAxes, fontsize=8,
                    ha='right', va='bottom', color='#c0392b',
                    bbox=dict(boxstyle='round,pad=0.3', fc='white',
                              ec='#DDDDDD', alpha=0.9))

    _ax_gl.set_xlabel('Max transcript length  (log₁₀ bp)', fontsize=10)
    _ax_gl.set_ylabel('Cumulative fraction', fontsize=10)
    _ax_gl.set_title(
        'Gene Length Distribution — CDK12 Variant DE Gene Groups',
        fontsize=10, fontweight='bold', pad=8)
    _ax_gl.legend(frameon=False, fontsize=8, loc='upper left',
                  handlelength=1.5, labelspacing=0.4)
    _ax_gl.spines[['top', 'right']].set_visible(False)
    _ax_gl.tick_params(labelsize=9, length=3)
    _ax_gl.xaxis.grid(True, lw=0.35, color='#EEEEEE', zorder=0)
    _ax_gl.set_axisbelow(True)
    _ax_gl.set_ylim(0, 1.05)

    _out_gl = os.path.join(OUTDIR, 'hypomorph_gene_length_cdf')
    plt.tight_layout()
    plt.savefig(f'{_out_gl}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{_out_gl}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {_out_gl}.pdf")
else:
    print("  Gene length fetch returned too few entries — skipping CDF.")

# ═══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 — Pre-ranked GSEA (mean hypomorph LFC)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[3/3] Pre-ranked GSEA...")

if _HAS_GSEAPY:
    # Ranking metric: mean LFC across both confirmed hypomorphs
    _rnk = _all['lfc_mean'].sort_values(ascending=False).dropna()
    _rnk.index.name = 'gene'
    print(f"  Ranking {len(_rnk):,} genes by mean(lfc_G822E, lfc_G731K)")

    PRERANK_SETS = ['MSigDB_Hallmark_2020', 'Reactome_2022']

    try:
        _pr = gp.prerank(
            rnk=_rnk,
            gene_sets=PRERANK_SETS,
            min_size=15,
            max_size=500,
            permutation_num=1000,
            seed=42,
            outdir=None,
            verbose=False,
        )
        _pr_df = _pr.res2d.copy()
        _pr_df = _pr_df[_pr_df['FDR q-val'] < 0.25].copy()
        _pr_df['Term'] = (_pr_df.index
                          .str.replace(r'_\d{4}.*', '', regex=True)
                          .str.replace('_', ' ', regex=False)
                          .str.replace('HALLMARK ', '', regex=False))
        _pr_df['NES']   = pd.to_numeric(_pr_df['NES'],   errors='coerce')
        _pr_df['FDR']   = pd.to_numeric(_pr_df['FDR q-val'], errors='coerce')
        _pr_df['-log10fdr'] = -np.log10(_pr_df['FDR'].clip(lower=1e-10))
        print(f"  {len(_pr_df)} significant terms (FDR < 0.25)")

        # Top 15 up + top 15 down
        _up_pr = (_pr_df[_pr_df['NES'] > 0]
                  .sort_values('NES', ascending=False)
                  .drop_duplicates('Term').head(15))
        _dn_pr = (_pr_df[_pr_df['NES'] < 0]
                  .sort_values('NES')
                  .drop_duplicates('Term').head(15))
        _plot_pr = pd.concat([_dn_pr, _up_pr]).sort_values('NES')

        if len(_plot_pr) > 0:
            _fig_h   = max(5, len(_plot_pr) * 0.36 + 1.8)
            _fig_pr, _ax_pr = plt.subplots(figsize=(8, _fig_h))

            _colors = ['#2980b9' if n < 0 else '#c0392b'
                       for n in _plot_pr['NES'].values]
            _bars = _ax_pr.barh(
                range(len(_plot_pr)),
                _plot_pr['NES'].values,
                color=_colors,
                edgecolor='black', linewidth=0.5,
                height=0.65, zorder=3
            )

            # FDR stars
            def _stars(fdr):
                if fdr < 0.01: return '***'
                if fdr < 0.05: return '**'
                if fdr < 0.25: return '*'
                return ''

            for i, (_, row) in enumerate(_plot_pr.iterrows()):
                st = _stars(row['FDR'])
                if st:
                    xpos = row['NES'] + (0.06 if row['NES'] > 0 else -0.06)
                    ha   = 'left' if row['NES'] > 0 else 'right'
                    _ax_pr.text(xpos, i, st, va='center', ha=ha,
                                fontsize=8, color='#333333')

            _ax_pr.set_yticks(range(len(_plot_pr)))
            _ax_pr.set_yticklabels(_plot_pr['Term'].tolist(), fontsize=8.5)
            _ax_pr.axvline(0, color='black', lw=0.8, zorder=4)
            _ax_pr.set_xlabel('Normalized Enrichment Score (NES)', fontsize=10)
            _ax_pr.set_title(
                'Pre-ranked GSEA — mean(lfc_G822E, lfc_G731K)\n'
                'Hallmark + Reactome   |   FDR < 0.25 shown   |   * <0.25  ** <0.05  *** <0.01',
                fontsize=9, fontweight='bold', pad=8)
            _ax_pr.spines[['top', 'right']].set_visible(False)
            _ax_pr.tick_params(labelsize=8.5, length=3, axis='x')
            _ax_pr.tick_params(length=0, axis='y')
            _ax_pr.xaxis.grid(True, lw=0.35, color='#EEEEEE', zorder=0)
            _ax_pr.set_axisbelow(True)

            plt.tight_layout()
            _out_pr = os.path.join(OUTDIR, 'hypomorph_program_gsea_prerank')
            plt.savefig(f'{_out_pr}.pdf', dpi=300, bbox_inches='tight')
            plt.savefig(f'{_out_pr}.png', dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  Saved: {_out_pr}.pdf")

    except Exception as e:
        print(f"  Pre-ranked GSEA failed: {e}")

print("\nDone.")
