"""
nt_gsea_preservation_score.py
──────────────────────────────
Computes a GSEA-based Signature Preservation Score (SPS) for each CDK12
variant under each drug — an orthogonal complement to the Pearson r score
in nt_resistance_score.py.

Method
------
1. Build NT drug-response gene sets from NT DESeq2 results:
       NT_{drug}_UP  (padj < 0.05, LFC > 1)
       NT_{drug}_DN  (padj < 0.05, LFC < -1)

2. For each variant × drug, run prerank GSEA of these NT signatures
   against the variant's drug-response Wald stat ranking.

3. Compute raw SPS:
       raw_SPS = (NES_UP + |NES_DN|) / 2
   A drug-sensitive variant should have NES_UP > 0 and NES_DN < 0
   (genes move in the same direction as NT), so both terms add positively.

4. Normalize by NT self-score (NT signatures scored against NT ranking)
   so NT = 1.0, directly comparable to the Pearson r plot.

Why orthogonal to Pearson r
---------------------------
    Pearson r  — parametric, genome-wide, equal weight to all ~15k genes
    SPS (GSEA) — rank-based, non-parametric, focused on NT-defined DEGs,
                 weighted toward rank extremes

Outputs
-------
    figures/04_gsea_cdk12/GSEA_preservation_score.{pdf,png,svg}
    figures/04_gsea_cdk12/GSEA_preservation_score_raw.csv

Runtime: ~3-5 min (8 prerank GSEA calls × 1000 permutations)

Run from the project root
-------------------------
    cd ~/Desktop/"Comprehensive PC3 RNA-seq Analysis"
    python3 scripts/nt_gsea_preservation_score.py
"""

import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import gseapy as gs

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams['font.family']      = 'sans-serif'
plt.rcParams['font.sans-serif']  = ['Liberation Sans', 'Arial', 'Helvetica Neue',
                                     'Helvetica', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype']     = 42
plt.rcParams['ps.fonttype']      = 42

# ── Paths ─────────────────────────────────────────────────────────────────────
_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'deseq2_results')
FIGDIR = os.path.join(_ROOT, 'pc3_rnaseq_results', 'figures', '04_gsea_cdk12')
os.makedirs(FIGDIR, exist_ok=True)

# ── Load DESeq2 CSVs ──────────────────────────────────────────────────────────
def load_csv(fname):
    path = os.path.join(RESDIR, f'{fname}.csv')
    df = pd.read_csv(path, index_col=0)
    df = df.dropna(subset=['stat', 'log2FoldChange', 'padj', 'gene_name'])
    df = df.drop_duplicates(subset='gene_name').set_index('gene_name')
    return df

print("Loading DESeq2 results...")
nt_hq461 = load_csv('NT_HQ461_vs_NT_Control')
nt_zly   = load_csv('NT_ZLY025_vs_NT_Control')
print(f"  NT HQ461:   {len(nt_hq461):,} genes")
print(f"  NT ZLY-025: {len(nt_zly):,} genes")

# ── Build universal NT CDK inhibition signature (intersection across both drugs)
def de_set(df, direction='down'):
    mask = (df['padj'] < 0.05) & (
        (df['log2FoldChange'] < -1) if direction == 'down' else (df['log2FoldChange'] > 1)
    )
    return set(df[mask].index)

core_up = sorted(de_set(nt_hq461, 'up') & de_set(nt_zly, 'up'))
core_dn = sorted(de_set(nt_hq461, 'down') & de_set(nt_zly, 'down'))
print(f"  Core NT CDK signature: {len(core_up)} UP, {len(core_dn)} DN genes")
print(f"  (intersection of HQ461 and ZLY-025 NT drug responses)")

CORE_SIGS = {'NT_CORE_UP': core_up, 'NT_CORE_DN': core_dn}

# ── GSEA helper ───────────────────────────────────────────────────────────────
def run_prerank(df, gene_sets_dict, permutations=1000, seed=42):
    """Run prerank GSEA; returns res2d indexed by Term."""
    ranking = df[['stat']].sort_values('stat', ascending=False)
    pre = gs.prerank(
        rnk=ranking,
        gene_sets=gene_sets_dict,
        permutation_num=permutations,
        seed=seed,
        min_size=5,
        max_size=5000,
        verbose=False,
        outdir=None,
    )
    return pre.res2d.set_index('Term')

def get_nes(res, term):
    """Safely retrieve NES, FDR, NOM p-val from GSEA result table."""
    if term not in res.index:
        return np.nan, np.nan, np.nan
    row = res.loc[term]
    return float(row['NES']), float(row['FDR q-val']), float(row['NOM p-val'])

# ── NT self-scores (normalization reference) ──────────────────────────────────
# ── NT self-scores: score each NT drug response against the core signature ────
print("\nComputing NT self-scores (normalization reference)...")
nt_self    = {}   # combined SPS (kept for reference)
nt_nes_ref = {}   # per-component NT NES for normalization
for drug, ref_df in [('HQ461', nt_hq461), ('ZLY-025', nt_zly)]:
    res = run_prerank(ref_df, CORE_SIGS)
    nes_up, _, _ = get_nes(res, 'NT_CORE_UP')
    nes_dn, _, _ = get_nes(res, 'NT_CORE_DN')
    self_sps = (nes_up + abs(nes_dn)) / 2
    nt_self[drug]    = self_sps
    nt_nes_ref[drug] = {'UP': nes_up, 'DN': abs(nes_dn)}
    print(f"  NT {drug}: NES_UP={nes_up:.3f}  NES_DN={nes_dn:.3f}  SPS={self_sps:.3f}  (→ normalizes to 1.0)")

# ── Score each variant ────────────────────────────────────────────────────────
COMPARISONS = [
    ('I733V_HQ461_vs_I733V_Control',  'I733V', 'HQ461'),
    ('G822E_HQ461_vs_G822E_Control',  'G822E', 'HQ461'),
    ('G731K_HQ461_vs_G731K_Control',  'G731K', 'HQ461'),
    ('I733V_ZLY025_vs_I733V_Control', 'I733V', 'ZLY-025'),
    ('G822E_ZLY025_vs_G822E_Control', 'G822E', 'ZLY-025'),
    ('G731K_ZLY025_vs_G731K_Control', 'G731K', 'ZLY-025'),
]

print("\nScoring variants...")
records = []
for fname, variant, drug in COMPARISONS:
    var_df    = load_csv(fname)
    res = run_prerank(var_df, CORE_SIGS)

    nes_up, fdr_up, nom_up = get_nes(res, 'NT_CORE_UP')
    nes_dn, fdr_dn, nom_dn = get_nes(res, 'NT_CORE_DN')

    # SPS: average only components with NOM p < 0.05
    ref = nt_nes_ref[drug]
    parts = []
    if nom_up < 0.05: parts.append(nes_up   / ref['UP'])
    if nom_dn < 0.05: parts.append(abs(nes_dn) / ref['DN'])
    norm_sps  = np.mean(parts) if parts else 0.0
    asymmetric = len(parts) == 1   # only one component contributed
    raw_sps   = (nes_up + abs(nes_dn)) / 2

    fdr_min = max(fdr_up, fdr_dn)
    nom_min = max(nom_up, nom_dn)

    records.append({
        'variant':    variant,
        'drug':       drug,
        'NES_UP':     nes_up,
        'NES_DN':     nes_dn,
        'FDR_UP':     fdr_up,
        'FDR_DN':     fdr_dn,
        'NOM_UP':     nom_up,
        'NOM_DN':     nom_dn,
        'raw_SPS':    raw_sps,
        'SPS':        norm_sps,
        'FDR_min':    fdr_min,
        'NOM_min':    nom_min,
        'asymmetric': asymmetric,
    })
    print(f"  {variant} {drug:7s}  NES_UP={nes_up:+.3f}  NES_DN={nes_dn:+.3f}  "
          f"SPS={norm_sps:.3f}  FDR(max)={fdr_min:.3f}")

scores = pd.DataFrame(records)
csv_out = os.path.join(FIGDIR, 'GSEA_preservation_score_raw.csv')
scores.to_csv(csv_out, index=False)
print(f"\nRaw scores saved: {csv_out}")

# ── Plot ──────────────────────────────────────────────────────────────────────
VARIANT_COLORS = {
    'I733V': '#4393C3',   # blue
    'G822E': '#F4A582',   # orange
    'G731K': '#D6604D',   # red
}
DRUGS = ['HQ461', 'ZLY-025']

fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.5), sharey=False)
fig.subplots_adjust(wspace=0.38)

def stars(fdr):
    return '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else ''

for ax, drug in zip(axes, DRUGS):
    # Sort ascending by SPS: most resistant on left
    sub = (scores[scores['drug'] == drug]
           .sort_values('SPS', ascending=True)
           .reset_index(drop=True))
    n     = len(sub)
    x_pos = np.arange(n)
    colors = [VARIANT_COLORS[v] for v in sub['variant']]

    ax.bar(x_pos, sub['SPS'], color=colors, width=0.55,
           edgecolor='black', linewidth=0.7, zorder=3)

    # FDR stars (conservative: weaker of FDR_UP / FDR_DN) + NES text
    for xi, (_, row) in enumerate(sub.iterrows()):
        fdr_star = stars(row['FDR_min'])
        asym_flag = '†' if row['asymmetric'] else ''
        label = fdr_star + asym_flag + "\n" + f"\u2191{row['NES_UP']:+.2f}  \u2193{row['NES_DN']:+.2f}"
        ax.text(xi, row['SPS'] + 0.02, label,
                ha='center', va='bottom', fontsize=6.5, color='#333333',
                linespacing=1.4)

    # NT reference at 1.0
    ax.axhline(1.0, color='#555555', lw=1.0, linestyle='--', zorder=2)

    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1.45)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(list(sub['variant']), fontsize=9)
    ax.set_title(drug, fontsize=10, fontweight='bold', pad=6)
    ax.set_ylabel('Signature Preservation Score\n(GSEA-based, NT = 1.0)', fontsize=8.5)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.tick_params(axis='both', labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, color='#e0e0e0', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

fig.text(
    0.5, 0.01,
    'SPS = (NES↑ + |NES↓|) / 2, normalized to NT self-score (NT = 1.0).  '
    'Gene sets: core NT CDK inhibition signature (intersection of HQ461 ∩ ZLY-025 NT DEGs, padj<0.05, |LFC|>1).  '
    'Stars: FDR of weaker GSEA component (***<0.001, **<0.01, *<0.05).  '
    '†: SPS computed from one significant component only (other NOM p≥0.05).',
    ha='center', va='bottom', fontsize=6.5, color='#555555', style='italic'
)

fig.suptitle('Core NT CDK Inhibition Signature Preservation — PC3 CDK12 variants',
             fontsize=10, y=1.02)

outpath = os.path.join(FIGDIR, 'GSEA_preservation_score')
plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{outpath}.svg', dpi=300, bbox_inches='tight')
plt.close()

print(f"\nSaved: {outpath}.{{pdf,png,svg}}")
print("\nFinal scores:")
print(scores[['variant', 'drug', 'NES_UP', 'NES_DN', 'SPS', 'FDR_min']].to_string(index=False))
