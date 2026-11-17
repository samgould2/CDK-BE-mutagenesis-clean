"""
PC3 CDK12 v2 RNA-seq: DESeq2 + GSEA
4 genotypes (NT, I733V, G822E, G731K) × 3 treatments (Control, HQ461, ZLY025) × 3 replicates = 36 samples

Inputs:
    - salmon.merged.gene_counts.tsv  (nf-core/rnaseq output)
    - deseq2.pca.vals.txt            (nf-core/rnaseq QC output, optional)
Outputs:
    - pc3_rnaseq_results/deseq2_results/  (per-comparison CSVs)
    - pc3_rnaseq_results/figures/         (volcano plots, PCA, GSEA heatmaps)

Install dependencies:
    pip install pydeseq2 gseapy adjustText pandas numpy matplotlib
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings('ignore')

# ── Style ─────────────────────────────────────────────────────────────────────
# Font: try Liberation Sans (Linux), fall back to Arial / Helvetica / DejaVu Sans
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'Helvetica Neue',
                                    'Helvetica', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype'] = 42   # editable text in Illustrator
plt.rcParams['ps.fonttype']      = 42
plt.rcParams['svg.fonttype']     = 'none'

# ── Paths ─────────────────────────────────────────────────────────────────────
_ROOT       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COUNTS_FILE = os.path.join(_ROOT, 'data', 'salmon.merged.gene_counts.tsv')
PCA_FILE    = os.path.join(_ROOT, 'data', 'deseq2.pca.vals.txt')
OUTDIR      = os.path.join(_ROOT, 'pc3_rnaseq_results')

# Create output subfolders
for _d in [
    'deseq2_results',
    'figures/01_qc',
    'figures/02_volcano',
    'figures/03_hypomorphic',
    'figures/04_gsea_cdk12',
    'figures/05_hallmark_gsea',
    'figures/06_reactome',
    'figures/07_dna_damage_stress',
]:
    os.makedirs(os.path.join(OUTDIR, _d), exist_ok=True)

# ── Load counts ───────────────────────────────────────────────────────────────
print("Loading count matrix...")
counts_df = pd.read_csv(COUNTS_FILE, sep='\t', index_col=0)
gene_names = counts_df['gene_name'].copy()
counts_df  = counts_df.drop(columns=['gene_name']).round().astype(int)
print(f"  {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
print(f"  Columns: {list(counts_df.columns)}\n")

# ── Sample metadata ───────────────────────────────────────────────────────────
# Column format: {num}_{genotype}_{treatment}Rep{rep}
# e.g. 1_NT_ControlRep1, 15_NT_HQ461Rep3, 27_NT_ZLYRep3
def parse_sample(col):
    m = re.match(r'X?\d+_(\w+?)_(Control|HQ461|ZLY)Rep\d+', col)
    if not m:
        raise ValueError(f"Cannot parse sample name: {col}")
    genotype  = m.group(1)
    treatment = m.group(2)
    if treatment == 'ZLY':
        treatment = 'ZLY025'
    return genotype, treatment

metadata = pd.DataFrame(
    [parse_sample(c) for c in counts_df.columns],
    columns=['genotype', 'treatment'],
    index=counts_df.columns
)
metadata['condition'] = metadata['genotype'] + '_' + metadata['treatment']
print("Sample metadata:")
print(metadata.to_string())
print()

GENOTYPES  = ['NT', 'I733V', 'G822E', 'G731K']
TREATMENTS = ['Control', 'HQ461', 'ZLY025']

# ── PyDESeq2 ──────────────────────────────────────────────────────────────────
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

keep        = (counts_df >= 10).sum(axis=1) >= 3
counts_filt = counts_df.loc[keep]
print(f"After filtering: {counts_filt.shape[0]} genes retained\n")

dds = DeseqDataSet(
    counts=counts_filt.T,
    metadata=metadata[['condition']],
    design_factors='condition',
    refit_cooks=True,
)
dds.deseq2()

# ── Comparisons ───────────────────────────────────────────────────────────────
# (test, ref, filename, plot_title)
comparisons = [
    # Drug effects per genotype — HQ461
    ('NT_HQ461',    'NT_Control',    'NT_HQ461_vs_NT_Control',          'NT: HQ461 vs Control'),
    ('I733V_HQ461', 'I733V_Control', 'I733V_HQ461_vs_I733V_Control',    'I733V: HQ461 vs Control'),
    ('G822E_HQ461', 'G822E_Control', 'G822E_HQ461_vs_G822E_Control',    'G822E: HQ461 vs Control'),
    ('G731K_HQ461', 'G731K_Control', 'G731K_HQ461_vs_G731K_Control',    'G731K: HQ461 vs Control'),
    # Drug effects per genotype — ZLY025
    ('NT_ZLY025',    'NT_Control',    'NT_ZLY025_vs_NT_Control',         'NT: ZLY-025 vs Control'),
    ('I733V_ZLY025', 'I733V_Control', 'I733V_ZLY025_vs_I733V_Control',  'I733V: ZLY-025 vs Control'),
    ('G822E_ZLY025', 'G822E_Control', 'G822E_ZLY025_vs_G822E_Control',  'G822E: ZLY-025 vs Control'),
    ('G731K_ZLY025', 'G731K_Control', 'G731K_ZLY025_vs_G731K_Control',  'G731K: ZLY-025 vs Control'),
    # Resistance under HQ461 (NT is sensitive/numerator)
    ('NT_HQ461', 'I733V_HQ461', 'NT_vs_I733V_HQ461', 'HQ461: NT vs I733V'),
    ('NT_HQ461', 'G822E_HQ461', 'NT_vs_G822E_HQ461', 'HQ461: NT vs G822E'),
    ('NT_HQ461', 'G731K_HQ461', 'NT_vs_G731K_HQ461', 'HQ461: NT vs G731K'),
    # Resistance under ZLY025
    ('NT_ZLY025', 'I733V_ZLY025', 'NT_vs_I733V_ZLY025', 'ZLY-025: NT vs I733V'),
    ('NT_ZLY025', 'G822E_ZLY025', 'NT_vs_G822E_ZLY025', 'ZLY-025: NT vs G822E'),
    ('NT_ZLY025', 'G731K_ZLY025', 'NT_vs_G731K_ZLY025', 'ZLY-025: NT vs G731K'),
    # Intrinsic differences — no drug (Control)
    ('NT_Control', 'I733V_Control', 'NT_vs_I733V_Control', 'Control: NT vs I733V'),
    ('NT_Control', 'G822E_Control', 'NT_vs_G822E_Control', 'Control: NT vs G822E'),
    ('NT_Control', 'G731K_Control', 'NT_vs_G731K_Control', 'Control: NT vs G731K'),
]

results = {}
for test, ref, name, title in comparisons:
    print(f"Running: {name}")
    stat_res = DeseqStats(dds, contrast=['condition', test, ref], alpha=0.05)
    stat_res.summary()

    coeff     = f'condition[T.{test}]'
    available = list(dds.varm['LFC'].columns)
    if coeff in available:
        stat_res.lfc_shrink(coeff=coeff)
    else:
        print(f"  Note: {coeff} not available, skipping shrinkage")

    res = stat_res.results_df.copy().join(gene_names)
    res = res.reset_index().rename(columns={'index': 'gene_id'}).sort_values('padj')

    sig = res[(res['padj'] < 0.05) & (res['log2FoldChange'].abs() > 1)].shape[0]
    print(f"  Significant (padj<0.05, |LFC|>1): {sig}\n")

    res.to_csv(f'{OUTDIR}/deseq2_results/{name}.csv', index=False)
    results[name] = (res, title)

print("All DESeq2 comparisons done.\n")

# ── Volcano plot helper ────────────────────────────────────────────────────────
LFC_THRESH  = 1.0
PADJ_THRESH = 0.05

PADJ_FLOOR = 1e-300  # floor before log10 to avoid Inf

def compute_axis_limits(res_subset):
    all_lfc, all_y = [], []
    for name, (res, _) in res_subset.items():
        rc = res.dropna(subset=['padj', 'log2FoldChange'])
        all_lfc.extend(rc['log2FoldChange'].values)
        y_vals = -np.log10(rc['padj'].clip(lower=PADJ_FLOOR).values)
        all_y.extend(y_vals[np.isfinite(y_vals)])
    x_lim = np.ceil(max(abs(np.nanmin(all_lfc)), abs(np.nanmax(all_lfc)))) + 0.5
    y_lim = np.ceil(np.nanmax(all_y)) + 2
    return x_lim, y_lim

def draw_volcano(ax, res, title, x_lim, y_lim):
    rc = res.dropna(subset=['padj', 'log2FoldChange']).copy()
    rc['padj'] = rc['padj'].clip(lower=PADJ_FLOOR)
    up   = (rc['padj'] < PADJ_THRESH) & (rc['log2FoldChange'] >  LFC_THRESH)
    down = (rc['padj'] < PADJ_THRESH) & (rc['log2FoldChange'] < -LFC_THRESH)
    ns   = ~(up | down)
    ax.scatter(rc.loc[ns,   'log2FoldChange'], -np.log10(rc.loc[ns,   'padj']),
               s=4, color="#AAAAAA", alpha=0.4, linewidths=0, rasterized=True)
    ax.scatter(rc.loc[up,   'log2FoldChange'], -np.log10(rc.loc[up,   'padj']),
               s=25, color="tab:red",  alpha=0.7, linewidths=0.3,
               edgecolors='black', rasterized=True, label=f'{up.sum()} UP')
    ax.scatter(rc.loc[down, 'log2FoldChange'], -np.log10(rc.loc[down, 'padj']),
               s=25, color="tab:blue", alpha=0.7, linewidths=0.3,
               edgecolors='black', rasterized=True, label=f'{down.sum()} DOWN')
    ax.set_xlim(-x_lim, x_lim)
    ax.set_ylim(0, y_lim)
    ax.legend(fontsize=10, frameon=True, framealpha=0.9, edgecolor='#CCCCCC')
    ax.axvline( LFC_THRESH, color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.axvline(-LFC_THRESH, color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.axhline(-np.log10(PADJ_THRESH), color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.set_xlabel('log₂ Fold Change', fontsize=12)
    ax.set_ylabel('-log₁₀(padj)', fontsize=12)
    ax.set_title(title, fontsize=12)

# ── Figure 1: Drug effects (2 rows × 4 cols) ─────────────────────────────────
drug_keys = [
    ('NT_HQ461_vs_NT_Control',       'I733V_HQ461_vs_I733V_Control',
     'G822E_HQ461_vs_G822E_Control', 'G731K_HQ461_vs_G731K_Control'),
    ('NT_ZLY025_vs_NT_Control',      'I733V_ZLY025_vs_I733V_Control',
     'G822E_ZLY025_vs_G822E_Control','G731K_ZLY025_vs_G731K_Control'),
]
drug_subset = {k: results[k] for row in drug_keys for k in row}
x_lim, y_lim = compute_axis_limits(drug_subset)

fig, axes = plt.subplots(2, 4, figsize=(22, 10))
for row_idx, row_keys in enumerate(drug_keys):
    for col_idx, key in enumerate(row_keys):
        draw_volcano(axes[row_idx, col_idx], results[key][0], results[key][1], x_lim, y_lim)

fig.suptitle('PC3 CDK12 v2 — Drug Effects', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/02_volcano/volcano_drug_effects.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/02_volcano/volcano_drug_effects.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: volcano_drug_effects")

# ── Figure 2: Resistance comparisons (2 rows × 3 cols) ───────────────────────
resist_keys = [
    ('NT_vs_I733V_HQ461', 'NT_vs_G822E_HQ461', 'NT_vs_G731K_HQ461'),
    ('NT_vs_I733V_ZLY025','NT_vs_G822E_ZLY025','NT_vs_G731K_ZLY025'),
]
resist_subset = {k: results[k] for row in resist_keys for k in row}
x_lim_r, y_lim_r = compute_axis_limits(resist_subset)

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
for row_idx, row_keys in enumerate(resist_keys):
    for col_idx, key in enumerate(row_keys):
        draw_volcano(axes[row_idx, col_idx], results[key][0], results[key][1], x_lim_r, y_lim_r)

fig.suptitle('PC3 CDK12 v2 — Resistance Comparisons (NT = sensitive)', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/02_volcano/volcano_resistance.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/02_volcano/volcano_resistance.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: volcano_resistance")

# ── Figure 3: Intrinsic differences — No Drug (1 row × 3 cols) ───────────────
nodrug_keys = [('NT_vs_I733V_Control', 'NT_vs_G822E_Control', 'NT_vs_G731K_Control')]
nodrug_subset = {k: results[k] for row in nodrug_keys for k in row}
x_lim_nd, y_lim_nd = compute_axis_limits(nodrug_subset)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for col_idx, key in enumerate(nodrug_keys[0]):
    draw_volcano(axes[col_idx], results[key][0], results[key][1], x_lim_nd, y_lim_nd)

fig.suptitle('PC3 CDK12 v2 — Intrinsic Differences (No Drug)', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/02_volcano/volcano_nodrug.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/02_volcano/volcano_nodrug.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: volcano_nodrug")

# ── PCA ───────────────────────────────────────────────────────────────────────
if os.path.exists(PCA_FILE):
    pca = pd.read_csv(PCA_FILE, sep='\t')
    pc1_col = [c for c in pca.columns if 'PC1' in c][0]
    pc2_col = [c for c in pca.columns if 'PC2' in c][0]

    def label_pca(s):
        m = re.match(r'X?\d+_(\w+?)_(Control|HQ461|ZLY)Rep\d+', str(s))
        if not m:
            return 'other'
        g, t = m.group(1), m.group(2)
        return f'{g}_{t}' if t != 'ZLY' else f'{g}_ZLY025'

    pca['condition'] = pca.iloc[:, 0].apply(label_pca)

    colors_map = {
        'NT_Control':    '#2166ac', 'NT_HQ461':    '#92c5de', 'NT_ZLY025':    '#4393c3',
        'I733V_Control': '#d6604d', 'I733V_HQ461': '#f4a582', 'I733V_ZLY025': '#e08060',
        'G822E_Control': '#4dac26', 'G822E_HQ461': '#b8e186', 'G822E_ZLY025': '#7fbc41',
        'G731K_Control': '#8e0152', 'G731K_HQ461': '#de77ae', 'G731K_ZLY025': '#c51b7d',
    }
    markers_map = {
        'NT': 'o', 'I733V': 's', 'G822E': '^', 'G731K': 'D'
    }

    fig, ax = plt.subplots(figsize=(11, 6))
    for cond, grp in pca.groupby('condition'):
        parts = cond.rsplit('_', 1)
        geno = parts[0] if len(parts) == 2 else cond
        marker = markers_map.get(geno, 'o')
        ax.scatter(grp[pc1_col], grp[pc2_col],
                   label=cond.replace('_', ' '),
                   color=colors_map.get(cond, 'gray'),
                   marker=marker, s=100, zorder=3, linewidths=0.5,
                   edgecolors='white')

    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel(pc1_col, fontsize=12)
    ax.set_ylabel(pc2_col, fontsize=12)
    ax.set_title('PCA — PC3 CDK12 v2 RNA-seq', fontsize=13)
    # Legend outside plot to the right, sorted for readability
    handles, labels = ax.get_legend_handles_labels()
    sorted_pairs = sorted(zip(labels, handles), key=lambda x: x[0])
    sorted_labels, sorted_handles = zip(*sorted_pairs)
    ax.legend(sorted_handles, sorted_labels,
              fontsize=9, frameon=True, framealpha=0.9,
              edgecolor='#CCCCCC', ncol=1,
              bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/figures/01_qc/PCA.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{OUTDIR}/figures/01_qc/PCA.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: PCA")

# ── Hypomorphic Scoring ───────────────────────────────────────────────────────
# Two approaches:
# 1. ssGSEA: per-sample CDK12/DDR activity from log-normalized counts
# 2. Drug-response attenuation: LFC of CDK12 targets in mutants vs NT under drug
import gseapy as gs

# CDK12 signature sets (defined here, reused by GSEA section below)
CDK12_LOSS_DN      = ["BRCA1","BRCA2","ATR","ATRIP","TOPBP1","ETAA1","FANCA","FANCC","FANCD2","FANCI","FANCL","FANCM","RAD51","RAD51C","RAD51D","PALB2","BARD1","NBN","MRE11","RAD50","CHEK1","RPA1","RPA2","RPA3","RFC1","RFC2","RFC3","RFC4","RFC5","PCNA","POLD1","POLD2","POLD3","POLE","POLE2","LIG1","LIG3","XRCC1","XRCC2","XRCC3","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC45","CDT1","ORC1","EME1","MUS81","SLX4","GEN1","BLM","WRN","EXO1","DNA2"]
CDK12_LOSS_UP      = ["CGAS","STING1","IRF3","IRF7","IFNB1","IFNA1","IFNA2","MX1","MX2","OAS1","OAS2","OAS3","OASL","ISG15","ISG20","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","STAT1","STAT2","IRF9","DDX58","IFIH1","TRIM56","TRIM32","RSAD2","HERC5","XAF1"]
CDK13_LOSS_DN      = ["MCL1","CCND1","CDK4","CDK6","MYC","MYCN","SRSF1","SRSF2","SRSF3","SRSF7","U2AF1","U2AF2","SF3B1","ELL","ELL2","AFF1","AFF4","CDC7","TIPIN","CLASPIN","E2F1","E2F2","E2F3"]
CDK12_CDK13_SHARED = ["CCNK","BRCA1","BRCA2","ATR","FANCD2","FANCI","CHEK1","MCL1","CCND1","MYC","WEE1","CDC25A","POLR2A","SUPT5H","SUPT4H1","PAF1","CTR9","LEO1","CDC73","SPT6"]
HQ461_DN           = sorted(set(CDK12_LOSS_DN + CDK13_LOSS_DN + CDK12_CDK13_SHARED))

HYPO_SIGS = {
    "CDK12_LOSS_DN":         CDK12_LOSS_DN,
    "CDK12_LOSS_UP":         CDK12_LOSS_UP,
    "CDK12_CDK13_SHARED_DN": CDK12_CDK13_SHARED,
}

# -- 1. ssGSEA per-sample activity --------------------------------------------
print("Running ssGSEA for hypomorphic scoring...")
# Build gene-symbol expression matrix from size-factor-normalized counts
normed_counts = pd.DataFrame(
    dds.layers['normed_counts'],
    index=dds.obs_names,   # samples (rows in DDS)
    columns=dds.var_names  # gene IDs
).T                         # → gene_id × samples

# Map gene IDs → symbols (drop unmapped / duplicated symbols)
normed_sym = normed_counts.copy()
normed_sym.index = normed_sym.index.map(lambda x: gene_names.get(x, ''))
normed_sym = normed_sym[normed_sym.index != '']
normed_sym = normed_sym[~normed_sym.index.duplicated(keep='first')]
log_normed  = np.log1p(normed_sym)

ssgsea_out = gs.ssgsea(
    data=log_normed,
    gene_sets=HYPO_SIGS,
    outdir=None,
    no_plot=True,
    processes=1,
)
ss = ssgsea_out.res2d.pivot(index='Term', columns='Name', values='NES').T
ss.index = ss.index.str.replace(r'^X?\d+_', '', regex=True)   # strip leading num_
# re-attach sample names via metadata
ss_meta = metadata[['genotype', 'treatment']].copy()
ss_meta.index = ss_meta.index.str.replace(r'^X?\d+_', '', regex=True, )
ss = ss.join(ss_meta, how='left')

geno_order  = ['NT', 'I733V', 'G822E', 'G731K']
treat_colors = {'Control': '#555555', 'HQ461': '#d6604d', 'ZLY025': '#4393c3'}
treat_markers = {'Control': 'o', 'HQ461': 's', 'ZLY025': '^'}

sig_names = list(HYPO_SIGS.keys())
fig, axes = plt.subplots(1, len(sig_names), figsize=(5 * len(sig_names), 4.5), sharey=False)
if len(sig_names) == 1:
    axes = [axes]

for ax, sig in zip(axes, sig_names):
    x_tick_pos = np.arange(len(geno_order))
    for gi, geno in enumerate(geno_order):
        for ti, treat in enumerate(TREATMENTS):
            grp = ss[(ss['genotype'] == geno) & (ss['treatment'] == treat)]
            if grp.empty:
                continue
            x_off = (ti - 1) * 0.22
            vals  = grp[sig].values
            ax.scatter(
                [gi + x_off] * len(vals), vals,
                color=treat_colors[treat], marker=treat_markers[treat],
                s=55, zorder=4, alpha=0.9,
                label=treat if gi == 0 else ''
            )
            ax.plot(
                [gi + x_off - 0.07, gi + x_off + 0.07],
                [np.mean(vals)] * 2,
                color=treat_colors[treat], lw=2.5, zorder=5
            )
    ax.set_xticks(x_tick_pos)
    ax.set_xticklabels(geno_order, fontsize=10)
    ax.set_title(sig.replace('_', ' '), fontsize=9)
    ax.set_ylabel('ssGSEA NES', fontsize=10)
    ax.axhline(0, color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.spines[['top', 'right']].set_visible(False)
    if sig == sig_names[0]:
        ax.legend(fontsize=8, frameon=False, title='Treatment', title_fontsize=8)

fig.suptitle('CDK12 Pathway Activity by Genotype — ssGSEA', fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_ssgsea.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_ssgsea.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: hypomorphic_ssgsea")

# -- 2. Drug-response attenuation ---------------------------------------------
# Define the NT CDK12-inhibition signature: genes significantly changed
# in NT by HQ461 or ZLY025 (union, |LFC|>1, padj<0.05)
def get_sig_genes(comp_name):
    res = results[comp_name][0].dropna(subset=['padj', 'log2FoldChange', 'gene_name'])
    return set(res[(res['padj'] < 0.05) & (res['log2FoldChange'].abs() > 1)]['gene_name'])

nt_drug_genes = get_sig_genes('NT_HQ461_vs_NT_Control') | get_sig_genes('NT_ZLY025_vs_NT_Control')
print(f"NT CDK12-drug signature: {len(nt_drug_genes)} genes")

# For each mutant × drug, compute mean LFC of NT-drug genes
attn_rows = []
for treat_short, treat_full in [('HQ461', 'HQ461'), ('ZLY025', 'ZLY025')]:
    nt_res = results[f'NT_{treat_full}_vs_NT_Control'][0].dropna(subset=['log2FoldChange','gene_name'])
    nt_lfc = nt_res[nt_res['gene_name'].isin(nt_drug_genes)].set_index('gene_name')['log2FoldChange']
    for geno in ['NT', 'I733V', 'G822E', 'G731K']:
        key = f'{geno}_{treat_full}_vs_{geno}_Control'
        res_g = results[key][0].dropna(subset=['log2FoldChange', 'gene_name'])
        res_g = res_g[res_g['gene_name'].isin(nt_drug_genes)].set_index('gene_name')['log2FoldChange']
        common = nt_lfc.index.intersection(res_g.index)
        if len(common) == 0:
            continue
        attenuation = res_g[common].mean() / (nt_lfc[common].mean() + 1e-9)
        correlation  = float(np.corrcoef(nt_lfc[common].values, res_g[common].values)[0, 1])
        attn_rows.append({
            'genotype': geno, 'treatment': treat_full,
            'mean_lfc': res_g[common].mean(),
            'nt_mean_lfc': nt_lfc[common].mean(),
            'attenuation_ratio': attenuation,
            'pearson_r': correlation,
            'n_genes': len(common),
        })

attn_df = pd.DataFrame(attn_rows)
attn_df.to_csv(f'{OUTDIR}/deseq2_results/hypomorphic_attenuation.csv', index=False)

# Plot: mean LFC per genotype/treatment compared to NT
fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
for ax, treat in zip(axes, ['HQ461', 'ZLY025']):
    sub = attn_df[attn_df['treatment'] == treat]
    nt_val = sub[sub['genotype'] == 'NT']['mean_lfc'].values[0]
    mutants = sub[sub['genotype'] != 'NT']
    genos  = mutants['genotype'].tolist()
    colors = ['#d6604d', '#4dac26', '#8e0152'][:len(genos)]
    bars = ax.bar(genos, mutants['mean_lfc'].values, color=colors, width=0.5, zorder=3, alpha=0.85)
    ax.axhline(nt_val, color='#2166ac', lw=1.8, linestyle='--', label=f'NT ({nt_val:.2f})', zorder=4)
    ax.axhline(0, color='black', lw=0.5, alpha=0.4)
    # annotate Pearson r
    for bar, (_, row) in zip(bars, mutants.iterrows()):
        y = row['mean_lfc']
        ax.text(bar.get_x() + bar.get_width()/2, y + (0.05 if y >= 0 else -0.1),
                f"r={row['pearson_r']:.2f}", ha='center', va='bottom' if y >= 0 else 'top',
                fontsize=8.5)
    ax.set_title(f'{treat} drug response\n(NT CDK12-drug genes, n={sub["n_genes"].iloc[0]})', fontsize=10)
    ax.set_ylabel('Mean log₂FC vs own Control', fontsize=10)
    ax.set_xlabel('Genotype', fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)

fig.suptitle('Drug-Response Attenuation — Hypomorphic Score', fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_attenuation.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_attenuation.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: hypomorphic_attenuation")

# -- 3. Empirical hypomorphic score — do mutant Controls look pre-treated? -----
# Strategy: build a drug-response gene signature from NT data (top-N genes most
# strongly down-regulated in NT by each drug), then ask whether Control samples
# from the mutant genotypes already show suppression of these genes at baseline.
# If yes (negative NES relative to NT), the mutant CDK12 behaves like partial
# pharmacological inhibition even without drug.
print("Building empirical drug-response signatures and running ssGSEA on Controls...")

def top_down_genes(comp_name, n=200):
    """Return top-N most down-regulated gene symbols (padj<0.05) from comparison."""
    res_df = results[comp_name][0].dropna(subset=['padj', 'log2FoldChange', 'gene_name'])
    sig    = res_df[(res_df['padj'] < 0.05) & (res_df['log2FoldChange'] < -0.5)]
    return sig.nsmallest(n, 'log2FoldChange')['gene_name'].tolist()

emp_sigs = {}
for drug in ['HQ461', 'ZLY025']:
    key   = f'NT_{drug}_vs_NT_Control'
    genes = top_down_genes(key, n=200)
    emp_sigs[f'NT_{drug}_Down'] = genes
    print(f"  NT {drug} down-signature: {len(genes)} genes")

# ssGSEA restricted to Control samples (no drug) — genes are already mapped to symbols
ctrl_cols        = [c for c in log_normed.columns if 'Control' in c]
log_normed_ctrl  = log_normed[ctrl_cols]

ssgsea_emp = gs.ssgsea(
    data=log_normed_ctrl,
    gene_sets=emp_sigs,
    outdir=None,
    no_plot=True,
    processes=1,
)
ss_emp = ssgsea_emp.res2d.pivot(index='Term', columns='Name', values='NES').T
ss_emp.index = ss_emp.index.str.replace(r'^X?\d+_', '', regex=True)

ss_emp_meta = metadata[['genotype']].copy()
ss_emp_meta.index = ss_emp_meta.index.str.replace(r'^X?\d+_', '', regex=True)
ss_emp = ss_emp.join(ss_emp_meta, how='left')
ss_emp.to_csv(f'{OUTDIR}/deseq2_results/hypomorphic_empirical_ssgsea.csv')

# Plot: 1 row × 2 panels (HQ461 signature | ZLY025 signature)
emp_sig_names = list(emp_sigs.keys())
geno_colors_emp = {'NT': '#2166ac', 'I733V': '#d6604d', 'G822E': '#4dac26', 'G731K': '#8e0152'}

fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=False)
for ax, sig in zip(axes, emp_sig_names):
    geno_means = {}
    for gi, geno in enumerate(geno_order):
        grp = ss_emp[ss_emp['genotype'] == geno][sig].values
        if len(grp) == 0:
            continue
        geno_means[geno] = np.mean(grp)
        rng = np.random.RandomState(42 + gi)
        jitter = rng.uniform(-0.12, 0.12, len(grp))
        ax.scatter([gi + j for j in jitter], grp,
                   color=geno_colors_emp[geno], s=65, zorder=4, alpha=0.9,
                   edgecolors='white', linewidths=0.5)
        ax.plot([gi - 0.18, gi + 0.18], [geno_means[geno]] * 2,
                color=geno_colors_emp[geno], lw=2.5, zorder=5)

    # dotted reference at NT mean
    nt_mean = geno_means.get('NT', 0)
    ax.axhline(nt_mean, color='#2166ac', lw=1.0, linestyle=':', alpha=0.5,
               label=f'NT mean ({nt_mean:.2f})')
    ax.axhline(0, color='black', lw=0.5, linestyle='--', alpha=0.35)

    drug_label = sig.replace('NT_', '').replace('_Down', '')
    ax.set_xticks(range(len(geno_order)))
    ax.set_xticklabels(geno_order, fontsize=11)
    ax.set_ylabel('ssGSEA NES (Control samples only)', fontsize=10)
    ax.set_title(f'NT {drug_label} response signature\n'
                 f'(top-200 down-regulated genes, n={len(emp_sigs[sig])})', fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)

fig.suptitle(
    'Empirical Hypomorphic Score — Baseline CDK12-inhibition in mutant Controls\n'
    '(NES < NT baseline → genes already suppressed without drug)',
    fontsize=10, y=1.05
)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_empirical.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_empirical.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: hypomorphic_empirical")

# Scatter: NT drug-effect LFC vs mutant drug-effect LFC (one panel per drug × mutant)
# Correlation slope < 1 → attenuation; slope > 1 → amplification
print("Plotting LFC correlation scatter (NT drug response vs mutant)...")
fig, axes = plt.subplots(2, 3, figsize=(13, 9))
for row_idx, drug in enumerate(['HQ461', 'ZLY025']):
    nt_res = results[f'NT_{drug}_vs_NT_Control'][0].dropna(subset=['log2FoldChange', 'gene_name'])
    nt_res = nt_res[nt_res['gene_name'].isin(nt_drug_genes)].set_index('gene_name')
    for col_idx, geno in enumerate(['I733V', 'G822E', 'G731K']):
        ax = axes[row_idx, col_idx]
        mut_res = results[f'{geno}_{drug}_vs_{geno}_Control'][0].dropna(subset=['log2FoldChange', 'gene_name'])
        mut_res = mut_res[mut_res['gene_name'].isin(nt_drug_genes)].set_index('gene_name')
        common = nt_res.index.intersection(mut_res.index)
        if len(common) < 5:
            ax.set_visible(False)
            continue
        x = nt_res.loc[common, 'log2FoldChange'].values
        y = mut_res.loc[common, 'log2FoldChange'].values
        r = float(np.corrcoef(x, y)[0, 1])
        # OLS slope
        slope = np.polyfit(x, y, 1)[0]
        ax.scatter(x, y, s=12, alpha=0.5, color=geno_colors_emp[geno],
                   edgecolors='none', rasterized=True)
        lim = np.ceil(max(abs(x).max(), abs(y).max())) + 0.5
        ax.plot([-lim, lim], [-lim, lim], color='black', lw=0.7, linestyle='--',
                alpha=0.5, label='slope = 1')
        ax.plot([-lim, lim], [slope * -lim, slope * lim], color=geno_colors_emp[geno],
                lw=1.5, label=f'slope = {slope:.2f}')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel('NT log₂FC', fontsize=10)
        ax.set_ylabel(f'{geno} log₂FC', fontsize=10)
        drug_disp = drug.replace('ZLY025', 'ZLY-025')
        ax.set_title(f'{drug_disp}: NT vs {geno}\nr={r:.2f}, slope={slope:.2f}', fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[['top', 'right']].set_visible(False)

fig.suptitle('Drug-Response LFC Correlation — NT vs Mutants\n'
             '(slope < 1 = attenuated response; r = concordance)',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_lfc_scatter.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/hypomorphic_lfc_scatter.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: hypomorphic_lfc_scatter")

# -- 4. Long-gene elongation index — CDK12 (P)-Ser2 proxy ---------------------
# CDK12 phosphorylates RNAPII CTD-Ser2, which is required for elongation through
# long introns. Hypomorphic CDK12 → reduced (P)-Ser2 → long DDR/HR genes are
# disproportionately suppressed relative to short housekeeping genes.
# Metric: mean log1p(normed) of CDK12-sensitive long genes /
#         mean log1p(normed) of short housekeeping genes.
# Lower ratio in Control mutants → elongation defect consistent with hypomorphism.
print("Computing long-gene elongation index (CDK12 Ser2 proxy)...")

CDK12_LONG_GENES = [
    # Published CDK12-sensitive long genes (>50 kb) — DDR, HR, replication
    # Harrod et al. 2020, Krajewska et al. 2019, Blazek et al. 2011
    "ATR", "ATRIP", "TOPBP1", "ETAA1",
    "BRCA1", "BRCA2", "PALB2", "BARD1",
    "FANCD2", "FANCA", "FANCM", "FANCL", "FANCC",
    "NBN", "MRE11", "RAD50",
    "RAD51", "RAD51C",
    "EXO1", "BLM", "WRN", "DNA2",
    "CHEK1", "RPA1", "RFC1", "POLD1",
    "SLX4", "ERCC4",
]

SHORT_HK_GENES = [
    # Short (<10 kb), highly expressed, CDK12-independent housekeeping genes
    "ACTB", "GAPDH", "B2M", "RPL19", "RPS18", "RPLP0",
    "UBB", "UBC", "LDHA", "PGK1", "TPI1", "PPIA", "TUBB",
]

def mean_expr_genes(gene_list, matrix):
    present = [g for g in gene_list if g in matrix.index]
    print(f"    {len(present)}/{len(gene_list)} genes found in matrix")
    if not present:
        return pd.Series(np.nan, index=matrix.columns)
    return matrix.loc[present].mean(axis=0)

long_mean  = mean_expr_genes(CDK12_LONG_GENES, log_normed)
short_mean = mean_expr_genes(SHORT_HK_GENES, log_normed)
elon_index = long_mean / (short_mean + 1e-6)

elon_df = pd.DataFrame({
    'elongation_index': elon_index,
    'long_mean': long_mean,
    'short_mean': short_mean,
})
elon_df.index = elon_df.index.str.replace(r'^X?\d+_', '', regex=True)
elon_meta = metadata[['genotype', 'treatment']].copy()
elon_meta.index = elon_meta.index.str.replace(r'^X?\d+_', '', regex=True)
elon_df = elon_df.join(elon_meta, how='left')
elon_df.to_csv(f'{OUTDIR}/deseq2_results/elongation_index.csv')

ctrl_elon = elon_df[elon_df['treatment'] == 'Control']

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel 1: Control only — cleanest test of hypomorphism at baseline
ax = axes[0]
for gi, geno in enumerate(geno_order):
    vals = ctrl_elon[ctrl_elon['genotype'] == geno]['elongation_index'].values
    rng  = np.random.RandomState(gi + 10)
    jit  = rng.uniform(-0.12, 0.12, len(vals))
    ax.scatter([gi + j for j in jit], vals,
               color=geno_colors_emp[geno], s=70, zorder=4, alpha=0.9,
               edgecolors='white', linewidths=0.5)
    ax.plot([gi - 0.18, gi + 0.18], [np.mean(vals)] * 2,
            color=geno_colors_emp[geno], lw=2.5, zorder=5)

nt_mean_ei = np.mean(ctrl_elon[ctrl_elon['genotype'] == 'NT']['elongation_index'].values)
ax.axhline(nt_mean_ei, color='#2166ac', lw=1.0, linestyle=':', alpha=0.5,
           label=f'NT mean ({nt_mean_ei:.3f})')
ax.set_xticks(range(len(geno_order)))
ax.set_xticklabels(geno_order, fontsize=11)
ax.set_ylabel('Long / Short gene expression ratio', fontsize=10)
ax.set_title('CDK12 Elongation Index — Control only\n'
             f'(CDK12-sensitive long genes, n={len([g for g in CDK12_LONG_GENES if g in log_normed.index])})',
             fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.spines[['top', 'right']].set_visible(False)

# Panel 2: All conditions — does drug further suppress long genes in NT but not mutants?
ax = axes[1]
for gi, geno in enumerate(geno_order):
    for ti, treat in enumerate(TREATMENTS):
        grp = elon_df[(elon_df['genotype'] == geno) & (elon_df['treatment'] == treat)]
        if grp.empty:
            continue
        x_off = (ti - 1) * 0.22
        vals  = grp['elongation_index'].values
        ax.scatter([gi + x_off] * len(vals), vals,
                   color=treat_colors[treat], marker=treat_markers[treat],
                   s=55, zorder=4, alpha=0.9,
                   label=treat if gi == 0 else '')
        ax.plot([gi + x_off - 0.07, gi + x_off + 0.07],
                [np.mean(vals)] * 2,
                color=treat_colors[treat], lw=2.5, zorder=5)

ax.set_xticks(range(len(geno_order)))
ax.set_xticklabels(geno_order, fontsize=11)
ax.set_ylabel('Long / Short gene expression ratio', fontsize=10)
ax.set_title('CDK12 Elongation Index — All conditions\n'
             '(drug effect in NT vs no further suppression in mutants)', fontsize=10)
ax.legend(fontsize=8, frameon=False, title='Treatment', title_fontsize=8)
ax.spines[['top', 'right']].set_visible(False)

fig.suptitle(
    'CDK12 Elongation Index — Proxy for RNAPII (P)-Ser2 Activity\n'
    '(lower ratio = impaired Pol II elongation through long DDR/HR genes)',
    fontsize=10, y=1.04
)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/elongation_index.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/elongation_index.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: elongation_index")

# ── LFC Concordance Scatter ────────────────────────────────────────────────────
# Key question: do G822E and G731K (confirmed hypomorphs) drive concordant
# transcriptional changes that are absent in I733V (unconfirmed)?
# High R between G822E and G731K + low R vs I733V = mechanistic proof.
print("\nPlotting LFC concordance scatter...")

from scipy.stats import pearsonr as _pearsonr

def _get_lfc(key):
    """Extract gene_name + log2FoldChange + padj from a DESeq2 result."""
    df = results[key][0].copy()
    df = df.dropna(subset=['log2FoldChange', 'padj', 'gene_name'])
    df = df.drop_duplicates(subset='gene_name')
    return df.set_index('gene_name')[['log2FoldChange', 'padj']]

_r_g822e = _get_lfc('NT_vs_G822E_Control').rename(
    columns={'log2FoldChange': 'lfc_G822E', 'padj': 'padj_G822E'})
_r_g731k = _get_lfc('NT_vs_G731K_Control').rename(
    columns={'log2FoldChange': 'lfc_G731K', 'padj': 'padj_G731K'})
_r_i733v = _get_lfc('NT_vs_I733V_Control').rename(
    columns={'log2FoldChange': 'lfc_I733V', 'padj': 'padj_I733V'})

_conc = _r_g822e.join(_r_g731k, how='inner').join(_r_i733v, how='inner').dropna()

_FDR = 0.05

def _conc_panel(ax, df, xcol, ycol, xpadj, ypadj, xlabel, ylabel, dot_color, title):
    both  = (df[xpadj] < _FDR) & (df[ypadj] < _FDR)
    one   = ((df[xpadj] < _FDR) | (df[ypadj] < _FDR)) & ~both
    none  = ~(both | one)

    x, y = df[xcol].values, df[ycol].values

    ax.scatter(x[none], y[none], s=3,  c='#CCCCCC', alpha=0.35, linewidths=0, rasterized=True)
    ax.scatter(x[one],  y[one],  s=7,  c=dot_color, alpha=0.35, linewidths=0, rasterized=True)
    ax.scatter(x[both], y[both], s=14, c=dot_color, alpha=0.85, linewidths=0)

    lim = max(np.abs(x).max(), np.abs(y).max()) * 1.08
    ax.plot([-lim, lim], [-lim, lim], color='#888888', lw=0.8, ls='--', zorder=1)
    ax.axhline(0, color='#CCCCCC', lw=0.5, zorder=0)
    ax.axvline(0, color='#CCCCCC', lw=0.5, zorder=0)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

    fin = np.isfinite(x) & np.isfinite(y)

    # R on all genes (inflated by ~15k zeros — less discriminating)
    r_all, _ = _pearsonr(x[fin], y[fin])

    # R on DE genes only: significant in at least one comparison (the informative subset)
    de_mask = fin & (both | one)
    if de_mask.sum() > 2:
        r_de, _ = _pearsonr(x[de_mask], y[de_mask])
        r_de_str = f'R_DE = {r_de:.3f}  (FDR<0.05 in ≥1)'
    else:
        r_de_str = 'R_DE = n/a'

    ax.text(0.04, 0.97,
            f'R_all = {r_all:.3f}\n{r_de_str}\nn(shared FDR<0.05) = {both.sum():,}',
            transform=ax.transAxes, fontsize=8, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC', alpha=0.9))

    # Label top concordantly-changed genes in the shared-sig set
    if both.sum() > 0:
        _bs = df[both].copy()
        _bs['_d'] = np.sqrt(_bs[xcol]**2 + _bs[ycol]**2)
        for gene, row in _bs.nlargest(14, '_d').iterrows():
            ax.text(row[xcol] + lim * 0.02, row[ycol], gene,
                    fontsize=6.5, color='#222222', clip_on=True, va='center')

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=9, fontweight='bold', pad=6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=8)
    return r

_fig, _axes = plt.subplots(1, 3, figsize=(14, 4.8))

r_hm = _conc_panel(
    _axes[0], _conc, 'lfc_G822E', 'lfc_G731K', 'padj_G822E', 'padj_G731K',
    'log₂FC  (NT vs G822E)', 'log₂FC  (NT vs G731K)',
    '#c0392b', 'G822E vs G731K\n(both confirmed hypomorphs)')

_conc_panel(
    _axes[1], _conc, 'lfc_G822E', 'lfc_I733V', 'padj_G822E', 'padj_I733V',
    'log₂FC  (NT vs G822E)', 'log₂FC  (NT vs I733V)',
    '#2980b9', 'G822E vs I733V\n(confirmed vs unconfirmed)')

_conc_panel(
    _axes[2], _conc, 'lfc_G731K', 'lfc_I733V', 'padj_G731K', 'padj_I733V',
    'log₂FC  (NT vs G731K)', 'log₂FC  (NT vs I733V)',
    '#2980b9', 'G731K vs I733V\n(confirmed vs unconfirmed)')

_fig.suptitle('LFC Concordance — CDK12 Variant Comparisons at Baseline (Control)',
              fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/lfc_concordance.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/03_hypomorphic/lfc_concordance.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: lfc_concordance")

# Export the shared hypomorphic DEG list (G822E ∩ G731K, NOT I733V)
_both_hm  = (_conc['padj_G822E'] < _FDR) & (_conc['padj_G731K'] < _FDR)
_not_i733 = _conc['padj_I733V'] >= _FDR
_shared   = _conc[_both_hm & _not_i733].copy()
_shared['lfc_mean_hypomorphs'] = (_shared['lfc_G822E'] + _shared['lfc_G731K']) / 2
_shared.sort_values('lfc_mean_hypomorphs').to_csv(
    f'{OUTDIR}/deseq2_results/shared_hypomorphic_DEGs.csv')
print(f"  Shared hypomorphic DEGs (G822E ∩ G731K, not I733V): {len(_shared)} genes")

# ── GSEA ──────────────────────────────────────────────────────────────────────

DDR_SETS = {
    'GOBP_DSB_REPAIR_VIA_HR': [
        "BRCA1","BRCA2","RAD51","RAD51C","RAD51D","RAD51B","PALB2","BARD1",
        "NBN","MRE11","RAD50","ATM","RPA1","RPA2","RPA3","FANCD2","FANCI",
        "FANCA","FANCC","FANCF","FANCG","FANCL","FANCM","BRIP1","SLX4",
        "MUS81","EME1","GEN1","BLM","WRN","RTEL1","EXO1","DNA2",
        "TOPBP1","ATRIP","RAD54L","RAD54B","XRCC2","XRCC3","DMC1",
    ],
    'KEGG_FANCONI_ANEMIA_PATHWAY': [
        "FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI",
        "FANCJ","FANCL","FANCM","FANCN","FANCO","FANCP","FANCQ","FANCR",
        "FANCS","FANCT","BRCA1","BRCA2","PALB2","RAD51","RAD51C","SLX4",
        "BRIP1","UBE2T","FAAP100","FAAP24","FAAP20","RFWD3","MAD2L2",
        "REV1","REV3L","POLH","RPA1","RPA2","RPA3",
    ],
    'GOBP_DOUBLE_STRAND_BREAK_REPAIR': [
        "BRCA1","BRCA2","RAD51","PALB2","BARD1","NBN","MRE11","RAD50",
        "ATM","RPA1","RPA2","XRCC4","LIG4","NHEJ1","PRKDC","XRCC5","XRCC6",
        "FANCD2","FANCI","FANCA","FANCC","FANCL","SLX4","MUS81","GEN1",
        "BLM","EXO1","DNA2","TOPBP1","XRCC2","XRCC3","POLL","POLM",
    ],
    'GOBP_DNA_REPAIR': [
        "BRCA1","BRCA2","RAD51","RAD51C","RAD51D","PALB2","BARD1","NBN",
        "MRE11","RAD50","ATM","ATR","CHEK1","CHEK2","RPA1","RPA2","RPA3",
        "FANCD2","FANCI","FANCA","FANCC","FANCF","FANCG","FANCL","FANCM",
        "BRIP1","SLX4","MUS81","EME1","GEN1","BLM","WRN","RECQL4","RECQL5",
        "EXO1","DNA2","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5",
        "XPC","XPA","DDB1","DDB2","PCNA","RFC1","RFC2","RFC3","RFC4","RFC5",
        "LIG1","LIG3","LIG4","XRCC1","XRCC2","XRCC3","XRCC4","XRCC5","XRCC6",
        "PRKDC","POLB","POLD1","POLD2","POLD3","POLE","POLE2",
        "MLH1","MSH2","MSH3","MSH6","PMS1","PMS2",
        "UNG","APEX1","OGG1","MUTYH","NEIL1","NEIL2",
    ],
}

CUSTOM_SETS = {
    "CDK12_LOSS_DN":         CDK12_LOSS_DN,
    "CDK12_LOSS_UP":         CDK12_LOSS_UP,
    "CDK13_LOSS_DN":         CDK13_LOSS_DN,
    "CDK12_CDK13_SHARED_DN": CDK12_CDK13_SHARED,
    "HQ461_CDK12_CDK13_DN":  HQ461_DN,
}

def run_gsea(results_dict, gene_sets, label):
    gsea_res = {}
    for name, (res, title) in results_dict.items():
        rc = res.dropna(subset=['stat', 'gene_name']).drop_duplicates(subset='gene_name')
        ranking = rc[['gene_name', 'stat']].set_index('gene_name').sort_values('stat', ascending=False)
        pre = gs.prerank(
            rnk=ranking, gene_sets=gene_sets,
            permutation_num=1000, seed=42,
            min_size=5, max_size=500, verbose=False, outdir=None,
        )
        gsea_res[title] = pre.res2d.set_index('Term')
        print(f"  {label} GSEA done: {title}")
    return gsea_res

def plot_gsea_heatmap(gsea_res, sig_order, outpath, title, group_separators=None):
    comp_labels = list(gsea_res.keys())
    nes_df = pd.DataFrame({l: gsea_res[l]['NES']       for l in comp_labels}).loc[sig_order]
    fdr_df = pd.DataFrame({l: gsea_res[l]['FDR q-val'] for l in comp_labels}).loc[sig_order]
    nom_df = pd.DataFrame({l: gsea_res[l]['NOM p-val'] for l in comp_labels}).loc[sig_order]

    vmax = max(abs(nes_df.values.max()), abs(nes_df.values.min()))
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu_r

    n_rows = len(comp_labels)
    n_cols = len(sig_order)
    cell_w, cell_h = 1.7, 0.7

    fig, ax = plt.subplots(figsize=(n_cols * cell_w + 2, n_rows * cell_h + 1.5))
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_aspect('equal')
    ax.axis('off')

    for i, gene_set in enumerate(sig_order):
        for j, comp in enumerate(comp_labels):
            nes   = nes_df.loc[gene_set, comp]
            fdr   = fdr_df.loc[gene_set, comp]
            nom   = nom_df.loc[gene_set, comp]
            color = cmap(norm(nes))
            rect  = plt.Rectangle((i, n_rows - j - 1), 1, 1, color=color, zorder=1)
            ax.add_patch(rect)
            stars = '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else ('+' if nom < 0.05 else '')
            text_color = 'white' if abs(nes) > vmax * 0.6 else 'black'
            ax.text(i + 0.5, n_rows - j - 0.5, f'{nes:.2f}{stars}',
                    ha='center', va='center', fontsize=8.5, color=text_color, zorder=2)

    for i, gs_name in enumerate(sig_order):
        ax.text(i + 0.5, -0.15, gs_name, ha='right', va='top', fontsize=9, rotation=30)
    for j, comp in enumerate(comp_labels):
        ax.text(-0.05, n_rows - j - 0.5, comp, ha='right', va='center', fontsize=9)

    # Draw group separator lines if requested
    if group_separators:
        for sep in group_separators:
            y = n_rows - sep  # sep = number of rows from top
            ax.axhline(y=y, xmin=0, xmax=1, color='white', linewidth=2, zorder=3)
            ax.axhline(y=y, xmin=0, xmax=1, color='black', linewidth=0.8, linestyle='--', zorder=4)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, aspect=15)
    cbar.set_label('NES', fontsize=9)

    ax.set_title(title, fontsize=11, pad=20)
    plt.tight_layout(pad=2.0)
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

print("\nRunning GSEA — DDR (published gene sets)...")
ddr_gsea = run_gsea(results, DDR_SETS, 'DDR')
_DDR_SETS_ORDER = ['GOBP_DSB_REPAIR_VIA_HR', 'KEGG_FANCONI_ANEMIA_PATHWAY',
                   'GOBP_DOUBLE_STRAND_BREAK_REPAIR', 'GOBP_DNA_REPAIR']
_DDR_SUFFIX = '\n(negative NES = enriched in second/denominator condition)'

# Panel A: Drug effect — HQ461
_ddr_hq461 = {k: v for k, v in ddr_gsea.items() if 'HQ461 vs Control' in k}
plot_gsea_heatmap(_ddr_hq461, _DDR_SETS_ORDER,
    f'{OUTDIR}/figures/04_gsea_cdk12/GSEA_DDR_heatmap_drug_HQ461',
    f'DDR GSEA — Drug Effect (HQ461 vs Control){_DDR_SUFFIX}')

# Panel B: Drug effect — ZLY-025
_ddr_zly = {k: v for k, v in ddr_gsea.items() if 'ZLY-025 vs Control' in k}
plot_gsea_heatmap(_ddr_zly, _DDR_SETS_ORDER,
    f'{OUTDIR}/figures/04_gsea_cdk12/GSEA_DDR_heatmap_drug_ZLY025',
    f'DDR GSEA — Drug Effect (ZLY-025 vs Control){_DDR_SUFFIX}')

# Panel C: Variant effect across conditions (HQ461 / ZLY-025 / Control)
_ddr_var_keys = [
    'HQ461: NT vs I733V', 'HQ461: NT vs G822E', 'HQ461: NT vs G731K',
    'ZLY-025: NT vs I733V', 'ZLY-025: NT vs G822E', 'ZLY-025: NT vs G731K',
    'Control: NT vs I733V', 'Control: NT vs G822E', 'Control: NT vs G731K',
]
_ddr_var = {k: ddr_gsea[k] for k in _ddr_var_keys}
plot_gsea_heatmap(_ddr_var, _DDR_SETS_ORDER,
    f'{OUTDIR}/figures/04_gsea_cdk12/GSEA_DDR_heatmap_variant_effect',
    f'DDR GSEA — Variant Effect (NT vs variant){_DDR_SUFFIX}',
    group_separators=[3, 6])  # dashed lines after HQ461 trio and ZLY-025 trio

print("\nRunning GSEA — CDK12/13 custom signatures...")
custom_gsea = run_gsea(results, CUSTOM_SETS, 'CDK12/13')
plot_gsea_heatmap(
    custom_gsea,
    sig_order=["CDK12_LOSS_DN", "CDK13_LOSS_DN", "CDK12_CDK13_SHARED_DN",
               "HQ461_CDK12_CDK13_DN", "CDK12_LOSS_UP"],
    outpath=f'{OUTDIR}/figures/04_gsea_cdk12/GSEA_CDK12_13_heatmap',
    title='CDK12/13 Signature GSEA — PC3 CDK12 v2\n(negative NES = enriched in second/denominator condition)',
)

# ── Empirical NT drug sensitivity signature ────────────────────────────────────
print("\nBuilding empirical NT drug sensitivity signatures...")
_nt_hq461_res = results['NT_HQ461_vs_NT_Control'][0].dropna(subset=['padj', 'gene_name'])
_nt_zly_res   = results['NT_ZLY025_vs_NT_Control'][0].dropna(subset=['padj', 'gene_name'])

NT_DRUG_SIGS = {
    'NT_HQ461_DN':  _nt_hq461_res.query('padj < 0.05 and log2FoldChange < -1')['gene_name'].dropna().tolist(),
    'NT_HQ461_UP':  _nt_hq461_res.query('padj < 0.05 and log2FoldChange >  1')['gene_name'].dropna().tolist(),
    'NT_ZLY025_DN': _nt_zly_res.query('padj < 0.05 and log2FoldChange < -1')['gene_name'].dropna().tolist(),
    'NT_ZLY025_UP': _nt_zly_res.query('padj < 0.05 and log2FoldChange >  1')['gene_name'].dropna().tolist(),
}
for _sig_name, _genes in NT_DRUG_SIGS.items():
    print(f"  {_sig_name}: {len(_genes)} genes")

# Test against variant drug responses only (I733V, G822E, G731K × HQ461, ZLY-025)
# Excludes NT (circular) and Control-only comparisons
_var_drug_keys = [k for k in results
                  if ('HQ461_vs' in k or 'ZLY025_vs' in k) and not k.startswith('NT_')]
_var_drug_results = {k: results[k] for k in _var_drug_keys}

print("\nRunning GSEA — NT sensitivity signatures vs variant drug responses...")
_nt_sig_gsea = run_gsea(_var_drug_results, NT_DRUG_SIGS, 'NT_sigs')

# Order rows: HQ461 variants first, then ZLY-025 variants
_nt_sig_row_order = [
    'I733V: HQ461 vs Control', 'G822E: HQ461 vs Control', 'G731K: HQ461 vs Control',
    'I733V: ZLY-025 vs Control', 'G822E: ZLY-025 vs Control', 'G731K: ZLY-025 vs Control',
]
_nt_sig_gsea_ordered = {k: _nt_sig_gsea[k] for k in _nt_sig_row_order}

plot_gsea_heatmap(
    _nt_sig_gsea_ordered,
    sig_order=['NT_HQ461_DN', 'NT_HQ461_UP', 'NT_ZLY025_DN', 'NT_ZLY025_UP'],
    outpath=f'{OUTDIR}/figures/04_gsea_cdk12/GSEA_NT_sensitivity_signature',
    title='NT Drug Sensitivity Signature — scored across variants\n'
          '(negative NES = drug response preserved; near 0 = resistant)',
    group_separators=[3],  # dashed line between HQ461 and ZLY-025 variant groups
)

# ── Hallmark GSEA dot plots — NT: HQ461 vs Control & NT: ZLY-025 vs Control ──
print("\nRunning GSEA — MSigDB Hallmarks...")

HALLMARKS = {
    'HALLMARK_E2F_TARGETS': ["CCNE1","CCNE2","CDC20","CDC25A","CDC45","CDKN2C","CHEK1","DHFR","E2F1","E2F2","E2F3","EZH2","FOXM1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MKI67","MYBL2","PCNA","PLK1","POLE","PTTG1","RAD51","RAD51AP1","RRM1","RRM2","TYMS","UBE2C","BRCA1","BRCA2","BUB1","CCNA2","CCNB1","CDC6","CDK1","CDK2","GINS1","GINS2"],
    'HALLMARK_G2M_CHECKPOINT': ["AURKA","AURKB","BUB1","BUB1B","BUB3","CCNA2","CCNB1","CCNB2","CDC20","CDC25C","CDK1","CDKN2C","CENPA","CENPE","CENPF","CHEK1","CHEK2","ESCO2","FOXM1","H2AX","KIF11","KIF23","MAD1L1","MAD2L1","MCM10","MKI67","NDC80","NEK2","NUF2","NUSAP1","PLK1","PLK4","PTTG1","RAD51","RRM2","SPC25","TACC3","TOP2A","TPX2","TTK","UBE2C"],
    'HALLMARK_MYC_TARGETS_V1': ["AP4","CAD","CCT3","CCT4","CCT5","CCT6A","CCT7","CCT8","CDK4","DCTPP1","DUT","EIF4E","EXOSC5","FASN","G3BP1","GNL1","HSPA4","HSPE1","LDHA","MCM5","MCM6","MCM7","MRPS18B","MYC","MYCN","NME1","NOP16","NOP58","NPM1","PRDX3","PTMA","RAN","RPL11","RPL14","RPL18A","RPL23","RPL26","RPL27","RPL3","RPL30","RPL37A","RPS11","RPS14","RPS15","RPS18","RPS19","RPS20","RPS3","RPS5","RPS6","SET","SRM"],
    'HALLMARK_MYC_TARGETS_V2': ["AIMP2","CDC5L","EIF2S1","EIF4B","ERAL1","FARSA","FARSB","FBL","HARS","IARS","KARS","MARS","MCM4","MRPS12","MRTO4","NARS","NOP2","NOP56","NPM3","POLR1B","POLR2B","POLR2D","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR2L","PPAT","PRIM1","RFC2","RFC3","RFC4","RFC5","RPS28","SNRPD1","SNRPD2","SNRPD3","VARS"],
    'HALLMARK_DNA_REPAIR': ["APEX1","ATM","ATR","BRCA1","BRCA2","CHEK1","CHEK2","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","EXO1","FANCD2","FANCI","LIG1","LIG3","MLH1","MLH3","MSH2","MSH3","MSH6","MUTYH","NBN","OGG1","PARP1","PCNA","PLK3","POLD1","POLD2","POLD3","POLE","POLE2","PMS1","PMS2","RAD17","RAD18","RAD21","RAD23A","RAD23B","RAD50","RAD51","RAD51AP1","RAD52","RAD54L","REV1","REV3L","RFC1","RFC2","RFC3","RFC4","RFC5","RPA1","RPA2","RRM1","RRM2","TIPIN","TOPBP1","TP53","UNG","XPA","XPC"],
    'HALLMARK_P53_PATHWAY': ["AEN","APAF1","BAX","BBC3","CCNG1","CCNG2","CDKN1A","DDIT3","DDB2","DINB","EI24","FAS","GADD45A","GADD45B","GDF15","HIPK2","LRDD","MDM2","PIDD1","PLK2","PLK3","PMAIP1","PPM1D","PTEN","RRM2B","SERPINE1","SFN","STEAP3","TP53","TP53INP1","TP53I3","TRIAP1","ZMAT3","CDKN2A","SIAH1","ATM","CHEK1","CHEK2"],
    'HALLMARK_APOPTOSIS': ["ADD1","AKT3","APAF1","BAD","BAX","BCL2","BCL2L1","BCL2L11","BID","BIRC2","BIRC3","BIRC5","BNIP3","CASP1","CASP3","CASP6","CASP7","CASP8","CASP9","CASP10","CFLAR","CIDEA","CYCS","DFFA","DFFB","DIABLO","ENDOG","FADD","FAS","FASLG","GADD45A","GADD45B","HRK","MCL1","PMAIP1","TNFRSF10A","TNFRSF10B","TNFRSF25","TNFSF10","TRADD","TRAF2","XIAP","BBC3","BIK","BOK"],
    'HALLMARK_ANDROGEN_RESPONSE': ["ABCC4","AR","B4GALT1","BMPR1B","CNTN1","DRG1","EAF2","ELL2","FKBP5","GNMT","GUCY1A1","IGF1R","INPP4B","KLK2","KLK3","KRT19","MED28","MYO6","NKX3-1","NRCAM","PDLIM5","PMEPA1","PTPRF","RAB4A","SEC14L2","SNAI2","SRD5A2","STEAP4","TMPRSS2","ZBTB10","HERC3","MAF","OPRK1","PPAP2A","PRKD1","PTGER4"],
    'HALLMARK_MTORC1_SIGNALING': ["ACACA","ACLY","ACSL3","ACTB","AKT1S1","AMD1","CAD","CANX","CBX3","CCNF","CDK7","DDIT4","DEPTOR","EEF1A1","EEF2","EIF2S1","EIF3A","EIF4B","EIF4E","EIF4EBP1","FASN","GAPDH","HMGCR","HMGCS1","HSPA5","IDH1","LDHA","MDM2","MRPS15","MYC","NDRG1","NOP56","NUDT3","NUDT5","PPAT","PRKAA1","RACK1","RRAGC","RRAGD","RPTOR","RPS6KB1","SLC3A2","SLC7A5","SQSTM1","TALDO1","TKT","TRAP1","ULK1","VIM","WEE1","YTHDF1"],
    'HALLMARK_UNFOLDED_PROTEIN_RESPONSE': ["ARFGAP1","ARF4","ATF3","ATF4","ATF6","CALR","CANX","CCND1","CEBPB","CEBPG","CEBPZ","CUL2","DDIT3","DNAJB1","DNAJB9","DNAJC3","EIF2AK3","ERN1","GADD45A","HERPUD1","HMOX1","HSP90B1","HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HYOU1","MAP3K5","NFYA","NFYB","NFYC","NUCB2","PDIA3","PDIA4","PDIA6","PPP1R15A","SEC23A","SEC23B","SEC24A","SLC1A4","SLC7A5","TRIB3","UBE2J1","XBP1"],
    'HALLMARK_INFLAMMATORY_RESPONSE': ["ADORA2B","AHR","AREG","BCL6","CCL2","CCL20","CEBPB","CEBPD","CXCL1","CXCL2","CXCL3","CXCL8","DUSP1","DUSP4","DUSP5","EBI3","EGR1","EGR3","F3","FOSL1","FOSL2","FPR1","HIF1A","HMOX1","ICAM1","IL1A","IL1B","IL1R1","IL1RN","IL6","IL6R","IL8","INHBA","IRF1","KYNU","LIF","NFKB1","NFKB2","NFKBIA","NFKBIZ","NLRP3","PTGS2","RELA","RELB","RIPK2","S100A8","S100A9","SELE","TRAF1","VCAM1"],
    'HALLMARK_INTERFERON_GAMMA_RESPONSE': ["AIM2","APOL6","ARID5B","BST2","BTG1","C1S","CASP1","CASP4","CASP7","CCL2","CD274","CIITA","CXCL10","CXCL9","DDX58","GBP1","GBP2","GBP4","GBP5","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","ICAM1","IDO1","IFIH1","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IRF1","IRF2","IRF7","IRF9","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","OASL","STAT1","STAT2","TAP1","TAP2","TAPBP","TRIM21","TRIM25","UBE2L6","USP18","WARS"],
    'HALLMARK_INTERFERON_ALPHA_RESPONSE': ["ADAR","BST2","DDX58","DDX60","HERC5","HERC6","IFI6","IFI27","IFI35","IFI44","IFI44L","IFIT1","IFIT2","IFIT3","IFIT5","IFITM1","IFITM2","IFITM3","IRF7","IRF9","ISG15","ISG20","LAMP3","MX1","MX2","OAS1","OAS2","OAS3","OASL","PLSCR1","RSAD2","RTP4","SAMD9","SAMD9L","STAT1","STAT2","TRIM22","TRIM5","UBE2L6","USP18","XAF1","ZNFX1"],
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB': ["ABCA1","ACKR3","AREG","ATF3","BCL2A1","BCL3","BIRC2","BIRC3","BTG1","BTG2","CCL2","CCL20","CCL4","CCND1","CD44","CDKN1A","CFLAR","CXCL1","CXCL2","CXCL3","CXCL6","DUSP1","DUSP2","DUSP4","DUSP5","EGR1","EGR3","F3","FOSL1","GADD45A","GADD45B","HBEGF","ICAM1","IER3","IL1A","IL1B","IL6","IL6R","INHBA","IRF1","JUN","JUNB","KLF2","KLF4","KLF6","KYNU","LITAF","MAFF","MCL1","NFKB1","NFKB2","NFKBIA","NFKBIB","NFKBIZ","NR4A1","NR4A2","NR4A3","PTGS2","RELB","RIPK2","SELE","SOD2","SQSTM1","TANK","TIMP1","TNFAIP2","TNFAIP3","TNFAIP6","TRAF1","VEGFA","XIAP"],
    'HALLMARK_CELL_CYCLE': ["ANAPC1","ANAPC2","BUB1","BUB1B","BUB3","CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2","CDC14A","CDC20","CDC25A","CDC25B","CDC6","CDC7","CDK1","CDK2","CDK4","CDK6","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","E2F1","E2F2","E2F3","MAD1L1","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PLK1","PLK4","PRKDC","RB1","SKP2","WEE1"],
    'HALLMARK_HYPOXIA': ["ACKR3","ADM","ALDOA","ALDOB","ALDOC","ANGPTL4","BNIP3","CA9","CDKN1A","COL5A1","DDIT4","EGFR","ENO1","ENO2","ENO3","ERO1A","GAPDH","GPC1","GPI","HMOX1","HRAS","HIF1A","IGFBP1","IGFBP3","JMJD6","LDHA","LOX","MXI1","MYH9","NDRG1","NFIL3","P4HA1","P4HA2","PDK1","PFKFB3","PFKL","PFKP","PGAM1","PGK1","PLOD1","PLOD2","PRKCA","SLC2A1","SLC2A3","SLC6A6","TGFB3","TPI1","VEGFA","VEGFB","VHL"],
}

def hallmark_dotplot(res_key, res_label, outpath):
    res_r = results[res_key][0]
    rc    = res_r.dropna(subset=['stat', 'gene_name']).drop_duplicates(subset='gene_name')
    rank  = rc[['gene_name', 'stat']].set_index('gene_name').sort_values('stat', ascending=False)
    pre   = gs.prerank(rnk=rank, gene_sets=HALLMARKS, permutation_num=1000, seed=42,
                       min_size=5, max_size=500, verbose=False, outdir=None)
    hdf   = pre.res2d.copy()
    hdf['Term_clean']    = hdf['Term'].str.replace('HALLMARK_', '').str.replace('_', ' ').str.title()
    hdf['neg_log10_fdr'] = -np.log10(hdf['FDR q-val'].clip(1e-3))
    hdf_plot = hdf[hdf['FDR q-val'] < 0.25].sort_values('NES').reset_index(drop=True)

    if hdf_plot.empty:
        print(f"  No significant gene sets (FDR < 0.25) for {res_label} — skipping")
        return

    max_y = hdf_plot['neg_log10_fdr'].max() or 1
    sizes = 40 + 310 * (hdf_plot['neg_log10_fdr'] / max_y)

    fig, ax = plt.subplots(figsize=(6.5, 8))
    for i, row in hdf_plot.iterrows():
        color = 'tab:red' if row['NES'] > 0 else 'tab:blue'
        alpha = min(0.9, 0.3 + 0.6 * (row['neg_log10_fdr'] / max_y))
        ax.scatter(row['NES'], i, s=sizes[i], color=color,
                   alpha=alpha, edgecolors='black', linewidths=0.4, zorder=3)

    ax.set_yticks(range(len(hdf_plot)))
    ax.set_yticklabels(hdf_plot['Term_clean'], fontsize=11)
    ax.axvline(0, color='black', lw=0.8, zorder=2)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=13)
    ax.set_title(f'PC3 {res_label}\nMSigDB Hallmark GSEA', fontsize=13)
    ax.tick_params(axis='x', labelsize=12)
    ax.set_xlim(-3, 3)

    for fdr_val, label in [(0.3, 'FDR=0.30'), (0.1, 'FDR=0.10'), (0.01, 'FDR=0.01')]:
        s = 40 + 310 * (-np.log10(fdr_val) / max_y)
        ax.scatter([], [], s=s, c='grey', alpha=0.7, edgecolors='black', linewidths=0.4, label=label)
    ax.legend(title='Dot size\n= -log₁₀(FDR)', fontsize=9, title_fontsize=9,
              frameon=True, framealpha=0.9, edgecolor='#CCCCCC',
              loc='lower right', handletextpad=1.0)

    plt.tight_layout()
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

_HM = f'{OUTDIR}/figures/05_hallmark_gsea'
hallmark_dotplot('NT_HQ461_vs_NT_Control',  'NT: HQ461 vs Control',   f'{_HM}/hallmark_gsea_NT_HQ461')
hallmark_dotplot('NT_ZLY025_vs_NT_Control', 'NT: ZLY-025 vs Control', f'{_HM}/hallmark_gsea_NT_ZLY025')

# No-drug mutant vs NT (positive NES = higher in NT, negative NES = higher in mutant)
hallmark_dotplot('NT_vs_I733V_Control', 'NT vs I733V — no drug', f'{_HM}/hallmark_gsea_NT_vs_I733V_Control')
hallmark_dotplot('NT_vs_G822E_Control', 'NT vs G822E — no drug', f'{_HM}/hallmark_gsea_NT_vs_G822E_Control')
hallmark_dotplot('NT_vs_G731K_Control', 'NT vs G731K — no drug', f'{_HM}/hallmark_gsea_NT_vs_G731K_Control')

# ── Reactome GSEA — NT vs hypomorphic mutants at baseline ─────────────────────
print("\nRunning GSEA — Reactome (downloading from Enrichr)...")

def reactome_dotplot(res_key, res_label, outpath, gene_sets, top_n=30):
    """Reactome preranked GSEA dot-plot. Shows top_n significant terms by |NES|."""
    res_r = results[res_key][0]
    rc    = res_r.dropna(subset=['stat', 'gene_name']).drop_duplicates(subset='gene_name')
    rank  = rc[['gene_name', 'stat']].set_index('gene_name').sort_values('stat', ascending=False)

    pre  = gs.prerank(rnk=rank, gene_sets=gene_sets, permutation_num=1000, seed=42,
                      min_size=5, max_size=500, verbose=False, outdir=None)
    rdf  = pre.res2d.copy()

    # Clean Reactome term names (strip 'REACTOME_' prefix if present)
    rdf['Term_clean'] = (rdf['Term']
                         .str.replace(r'^REACTOME[_ ]', '', regex=True)
                         .str.replace('_', ' ')
                         .str.title())

    rdf['neg_log10_fdr'] = -np.log10(rdf['FDR q-val'].clip(1e-4))

    # Filter to significant and pick top_n by |NES|
    sig = rdf[rdf['FDR q-val'] < 0.25].copy()
    if sig.empty:
        print(f"  No significant Reactome terms (FDR<0.25) for {res_label} — skipping")
        return
    sig = sig.reindex(sig['NES'].abs().sort_values(ascending=False).index)
    sig = sig.head(top_n).sort_values('NES').reset_index(drop=True)

    max_y  = sig['neg_log10_fdr'].max() or 1
    sizes  = 40 + 310 * (sig['neg_log10_fdr'] / max_y)

    fig_h  = max(6, len(sig) * 0.38)
    fig, ax = plt.subplots(figsize=(7.5, fig_h))

    for i, row in sig.iterrows():
        color = 'tab:red' if row['NES'] > 0 else 'tab:blue'
        alpha = min(0.9, 0.3 + 0.6 * (row['neg_log10_fdr'] / max_y))
        ax.scatter(row['NES'], i, s=sizes[i], color=color,
                   alpha=alpha, edgecolors='black', linewidths=0.4, zorder=3)

    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(sig['Term_clean'], fontsize=9)
    ax.axvline(0, color='black', lw=0.8, zorder=2)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    ax.set_title(f'PC3 {res_label}\nReactome GSEA (top {top_n}, FDR<0.25)', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)

    for fdr_val, label in [(0.25, 'FDR=0.25'), (0.10, 'FDR=0.10'), (0.01, 'FDR=0.01')]:
        s = 40 + 310 * (-np.log10(fdr_val) / max_y)
        ax.scatter([], [], s=s, c='grey', alpha=0.7, edgecolors='black', linewidths=0.4, label=label)
    ax.legend(title='Dot size\n= -log10(FDR)', fontsize=9, title_fontsize=9,
              frameon=True, framealpha=0.9, edgecolor='#CCCCCC',
              loc='lower right', handletextpad=1.0)

    plt.tight_layout()
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")

try:
    REACTOME = gs.get_library('Reactome_2022')
    print(f"  Loaded Reactome_2022: {len(REACTOME)} gene sets")
    _RT = f'{OUTDIR}/figures/06_reactome'
    reactome_dotplot('NT_vs_G822E_Control', 'NT vs G822E — no drug', f'{_RT}/reactome_NT_vs_G822E', REACTOME)
    reactome_dotplot('NT_vs_G731K_Control', 'NT vs G731K — no drug', f'{_RT}/reactome_NT_vs_G731K', REACTOME)
    reactome_dotplot('NT_vs_I733V_Control', 'NT vs I733V — no drug', f'{_RT}/reactome_NT_vs_I733V', REACTOME)
except Exception as _e:
    print(f"  Reactome download failed ({_e}) — skipping Reactome GSEA")

# ── DNA Damage & Stress Signatures — ssGSEA heatmap ──────────────────────────
# Curated gene sets covering DDR sub-pathways at higher resolution than Hallmark DNA_REPAIR
# Scored per-sample with ssGSEA on control-condition samples only
# Positive score = pathway active; compare hypomorphic mutants vs NT
print("\nRunning DNA Damage & Stress ssGSEA...")

DDR_STRESS_SIGS = {
    # ── Replication stress & fork protection ──────────────────────────────────
    'Replication_Stress_Response': [
        'ATR','ATRIP','TOPBP1','ETAA1','CLASPIN','TIPIN','TIMELESS',
        'RPA1','RPA2','RPA3','RAD17','RFC2','RFC3','RFC4','RFC5',
        'CHEK1','WEE1','CDC25A','HLTF','SHPRH','ZRANB3','SMARCAL1',
        'RAD51','BRCA2','PALB2','FANCD2','FANCI','BRCA1','BARD1',
    ],
    # ── ATR-CHK1 signaling ────────────────────────────────────────────────────
    'ATR_CHK1_Signaling': [
        'ATR','ATRIP','TOPBP1','ETAA1','RAD9A','RAD1','HUS1',
        'RAD17','CLASPIN','TIPIN','TIMELESS','CHEK1',
        'CDC25A','CDC25B','WEE1','CDK1','CDK2','PKMYT1',
        'MCM10','FANCD2','FANCI','UHRF1','RIF1',
    ],
    # ── ATM-CHK2 / DSB signaling ──────────────────────────────────────────────
    'ATM_CHK2_DSB_Signaling': [
        'ATM','CHEK2','MDM2','TP53','H2AX','MDC1','RNF8','RNF168',
        'BRCA1','53BP1','RBBP8','CtIP','NBN','MRE11','RAD50',
        'PRKDC','XRCC5','XRCC6','KAP1','PALB2','ABRAXAS1','FAM175A',
    ],
    # ── Homologous recombination (HR) ─────────────────────────────────────────
    'Homologous_Recombination': [
        'BRCA1','BRCA2','PALB2','BARD1','RAD51','RAD51C','RAD51D',
        'XRCC2','XRCC3','NBN','MRE11','RAD50','ATM','RPA1','RPA2',
        'RAD54L','RAD54B','RTEL1','BLM','EXO1','DNA2','RECQL5',
        'ABRAXAS1','FAM175A','SLX4','MUS81','EME1','GEN1','SLX1A',
    ],
    # ── Fanconi Anemia / ICL repair ───────────────────────────────────────────
    'Fanconi_Anemia_ICL_Repair': [
        'FANCA','FANCB','FANCC','FANCD2','FANCE','FANCF','FANCG',
        'FANCI','FANCJ','FANCL','FANCM','FANCN','FANCP','FANCQ',
        'FANCS','FANCT','BRCA1','BRCA2','PALB2','RAD51','SLX4',
        'BRIP1','UBE2T','RFWD3','REV1','REV3L','MAD2L2',
    ],
    # ── Transcription-replication conflicts / R-loops ─────────────────────────
    'R_Loop_Resolution': [
        'SETX','AQR','DDX23','DHX9','ADAR1','RNH1','RNASEH1',
        'RNASEH2A','RNASEH2B','RNASEH2C','TOP1','TOP2A','TOP2B',
        'BRCA1','BRCA2','FANCD2','FANCA','FANCM','RECQL5',
        'SRSF1','SRSF2','HNRNPA1','HNRNPA2B1','THO1',
    ],
    # ── CDK12-dependent elongation targets (long DDR genes) ───────────────────
    'CDK12_Elongation_Targets': [
        'BRCA1','BRCA2','ATR','ATRIP','TOPBP1','FANCD2','FANCI',
        'FANCA','FANCC','FANCF','FANCG','FANCL','FANCM',
        'RAD51','RAD51C','RAD51D','PALB2','BARD1','NBN','MRE11',
        'RAD50','CHEK1','DNA2','EXO1','BLM','WRN','RTEL1',
        'RPA1','RPA2','RFC1','LIG1','POLD1','POLE',
    ],
    # ── cGAS-STING / innate immune sensing of DNA damage ─────────────────────
    'cGAS_STING_Innate_Sensing': [
        'CGAS','STING1','TBK1','IRF3','IRF7','IFNB1',
        'IFNA1','IFNA2','MX1','MX2','OAS1','OAS2','OAS3','OASL',
        'ISG15','ISG20','IFIT1','IFIT2','IFIT3','IFITM1','IFITM2',
        'IFITM3','STAT1','STAT2','IRF9','DDX58','IFIH1',
        'TRIM56','TRIM32','HERC5','RSAD2','XAF1',
    ],
    # ── CHK1-mediated G2/M arrest ─────────────────────────────────────────────
    'G2M_Checkpoint_CHK1': [
        'CHEK1','CDC25A','CDC25B','CDC25C','CDK1','CCNB1','CCNB2',
        'WEE1','PKMYT1','PLK1','AURKA','AURKB','BUB1','BUB1B',
        'MAD1L1','MAD2L1','BUBR1','ESPL1','PTTG1','SGOL1','SGOL2',
    ],
    # ── p21 / p53 / senescence response ──────────────────────────────────────
    'p53_p21_Stress_Response': [
        'TP53','CDKN1A','MDM2','MDM4','PPM1D','PLK3','HIPK2',
        'GADD45A','GADD45B','BBC3','PMAIP1','FAS','CASP3',
        'CASP7','CASP9','BAX','BAD','BCL2','BCL2L1','MCL1',
        'CDKN2A','RB1','E2F1','TP53INP1','TIGAR','DRAM1',
    ],
}

# Control-condition samples only
ctrl_meta = metadata[metadata['treatment'] == 'Control']
ctrl_sample_ids = ctrl_meta.index.tolist()

# Build log-normed expression matrix for control samples
ctrl_log = log_normed[ctrl_sample_ids]

print(f"  ssGSEA on {ctrl_log.shape[1]} control samples × {ctrl_log.shape[0]} genes")

ddrstress_out = gs.ssgsea(
    data=ctrl_log,
    gene_sets=DDR_STRESS_SIGS,
    outdir=None,
    no_plot=True,
    processes=1,
)

ss_ddr = ddrstress_out.res2d.pivot(index='Term', columns='Name', values='NES')
# Reorder columns so they're grouped by genotype
ctrl_geno_order = []
for g in ['NT', 'I733V', 'G822E', 'G731K']:
    cols = [c for c in ss_ddr.columns if f'_{g}_Control' in c or c.endswith(f'{g}_ControlRep1'.split('Rep')[0])]
    # simpler: match via metadata
    g_samples = ctrl_meta[ctrl_meta['genotype'] == g].index.tolist()
    ctrl_geno_order.extend([s for s in g_samples if s in ss_ddr.columns])

# Fall back to original order if matching failed
if len(ctrl_geno_order) == 0:
    ctrl_geno_order = list(ss_ddr.columns)
else:
    # add any unmatched columns
    ctrl_geno_order += [c for c in ss_ddr.columns if c not in ctrl_geno_order]

ss_ddr = ss_ddr[ctrl_geno_order]

# ── Shared helpers for pooled ssGSEA heatmaps ────────────────────────────────
_GENO_ORDER  = ['NT', 'I733V', 'G822E', 'G731K']
geno_palette = {'NT': '#2166ac', 'I733V': '#d6604d', 'G822E': '#4dac26', 'G731K': '#8e0152'}

def _pool_replicates(ss_raw, meta_df, geno_col_order):
    """Mean NES across replicates per genotype → DataFrame[term × genotype]."""
    pooled = {}
    for g in _GENO_ORDER:
        g_cols = [c for c in geno_col_order
                  if c in meta_df.index and meta_df.loc[c, 'genotype'] == g]
        if g_cols:
            pooled[g] = ss_raw[g_cols].mean(axis=1)
    present = [g for g in _GENO_ORDER if g in pooled]
    return pd.DataFrame(pooled)[present]

def _render_heatmap(mat_z, title, outpath_base, cmap):
    """Render an elegant pooled ssGSEA heatmap (rows=gene sets, cols=genotypes)."""
    n_rows, n_cols = mat_z.shape
    fig_w = n_cols * 1.25 + 2.8
    fig_h = n_rows * 0.56 + 1.6
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    arr  = mat_z.values.astype(float)
    flat = arr[~np.isnan(arr)]
    vmax = float(np.percentile(np.abs(flat), 95)) if len(flat) else 2.0
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im   = ax.imshow(arr, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

    # Annotate each cell with its z-score value
    for r in range(n_rows):
        for c in range(n_cols):
            val = arr[r, c]
            if not np.isnan(val):
                fg = 'white' if abs(val) > vmax * 0.62 else '#1a1a1a'
                ax.text(c, r, f'{val:+.2f}', ha='center', va='center',
                        fontsize=8.5, color=fg, fontweight='semibold')

    # x-axis: genotype names in their palette colours (no rotation needed — 4 cols)
    ax.set_xticks(range(n_cols))
    tick_lbls = ax.set_xticklabels(mat_z.columns, fontsize=12, fontweight='bold', rotation=0)
    for lbl, g in zip(tick_lbls, mat_z.columns):
        lbl.set_color(geno_palette.get(g, '#333333'))
    ax.xaxis.set_tick_params(length=0)

    # y-axis: gene set names
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([k.replace('_', ' ') for k in mat_z.index], fontsize=9.5)
    ax.yaxis.set_tick_params(length=0)

    ax.spines[:].set_visible(False)

    # Fine white grid lines separating cells
    for xi in np.arange(-0.5, n_cols, 1):
        ax.axvline(xi, color='white', lw=1.5, zorder=5)
    for yi in np.arange(-0.5, n_rows, 1):
        ax.axhline(yi, color='white', lw=0.8, zorder=5)

    ax.set_title(title, fontsize=11, pad=10)

    cbar = plt.colorbar(im, ax=ax, shrink=0.55, pad=0.03, aspect=18)
    cbar.set_label('z-score (mean NES)', fontsize=9)
    cbar.ax.tick_params(labelsize=8, length=0)
    cbar.outline.set_visible(False)

    fig.tight_layout()
    plt.savefig(f'{outpath_base}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath_base}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {outpath_base}')

# Pool DDR replicates → z-score → render
ss_ddr_pooled = _pool_replicates(ss_ddr, ctrl_meta, ctrl_geno_order)
_p_mean  = ss_ddr_pooled.mean(axis=1)
_p_std   = ss_ddr_pooled.std(axis=1).replace(0, 1)
ss_ddr_z = ss_ddr_pooled.subtract(_p_mean, axis=0).divide(_p_std, axis=0).dropna(how='all')

_DS = f'{OUTDIR}/figures/07_dna_damage_stress'
_render_heatmap(
    ss_ddr_z,
    'DNA Damage & Replication Stress — Control, Replicates Pooled\n(ssGSEA, z-scored per gene set)',
    f'{_DS}/ddr_stress_ssgsea_heatmap',
    plt.cm.RdBu_r,
)

# ── Prostate Lineage Identity — ssGSEA heatmap ────────────────────────────────
# Asks: do CDK12 hypomorphic mutations shift lineage state?
# Luminal identity, AR/FOXA1/NKX3-1 programs, basal/stem, NE, and plasticity
# Scored per-sample on control-condition samples only
print("\nRunning Prostate Lineage ssGSEA...")

PROSTATE_LINEAGE_SIGS = {
    # ── Luminal / AR-driven identity ──────────────────────────────────────────
    'AR_Canonical_Targets': [
        'KLK3','KLK2','TMPRSS2','NKX3-1','FKBP5','STEAP1','STEAP2',
        'PMEPA1','ABCC4','SLC45A3','MSMB','ALDH1A3','SORD','BMPR1B',
        'ELL2','HERC3','GNMT','MAF','ZBTB10','SEC14L2','DRG1','CNTN1',
        'PDLIM5','KRT19','GUCY1A1','INPP4B',
    ],
    'FOXA1_Luminal_Program': [
        'FOXA1','AR','KRT8','KRT18','KRT19','GATA3','TFF1','TFF3',
        'SPDEF','NKX3-1','HOXB13','EHF','CLDN3','CLDN4','CLDN7',
        'TACSTD2','MUC1','AGR2','PSCA','SLC45A3','TMPRSS2',
    ],
    'NKX3-1_Targets': [
        'NKX3-1','KLK3','KLK2','MSMB','SPDEF','STEAP1','LTF',
        'TMPRSS2','SLC45A3','PSCA','PMEPA1','ALDH1A3','PCA3',
    ],
    'HOXB13_Targets': [
        'HOXB13','AR','FOXA1','NKX3-1','KLK3','MEIS2','ETS2',
        'TMPRSS2','SPDEF','ALDH1A3','ABCC4','SLC45A3',
    ],
    # ── Basal / stem-like identity ────────────────────────────────────────────
    'Basal_Stem_Markers': [
        'TP63','KRT5','KRT14','KRT17','KRT6A','CD44','ITGA6',
        'ALDH1A1','LAMB3','COL17A1','CDH3','S100A2','S100A14',
        'SOX2','MYC','PROM1','CD49F','NGFR',
    ],
    # ── Neuroendocrine differentiation ───────────────────────────────────────
    'Neuroendocrine_Program': [
        'CHGA','CHGB','SYP','ENO2','NCAM1','AURKA','EZH2',
        'SOX2','MYCN','ASCL1','DLL3','NEUROD1','PEG10',
        'ONECUT2','BRN2','REST','INSM1','SCG3','SCG5',
        'SYT1','SNAP25','VGF','NPY',
    ],
    # ── Epithelial–mesenchymal / plasticity ───────────────────────────────────
    'EMT_Plasticity': [
        'VIM','CDH2','FN1','SNAI1','SNAI2','TWIST1','TWIST2',
        'ZEB1','ZEB2','FOXC2','ITGB6','MMP2','MMP9',
        'TGFB1','TGFB2','WNT5A','CDH1',   # CDH1 is E-cad (loss marker)
    ],
    # ── Androgen receptor splice variants / AR-V context ─────────────────────
    'AR_Coregulators': [
        'SRC','EP300','CREBBP','MED1','NCOA1','NCOA2','NCOA3',
        'NCOR1','NCOR2','SMRT','HDAC1','HDAC2','HDAC3',
        'BRD4','BRD2','CDK9','CCNT1','MYC','MAX','E2F1',
    ],
    # ── WNT / stem-cell signaling ─────────────────────────────────────────────
    'Wnt_Stemness': [
        'CTNNB1','APC','AXIN1','AXIN2','LRP5','LRP6','FZD1','FZD7',
        'TCF4','TCF7','LEF1','MYC','CCND1','CD44','LGR5','PROM1',
        'SOX9','ALDH1A1','WNT5A','WNT3A','RNF43','ZNRF3',
    ],
    # ── EZH2 / Polycomb — epigenetic reprogramming ───────────────────────────
    'EZH2_Polycomb_Targets': [
        'EZH2','EED','SUZ12','BMI1','RING1','RNF2','CBX7','CBX8',
        'PHC1','PHC2','RYBP','JARID2','ARID1A','SMARCB1',
        'CDKN2A','CDKN2B','HOXA','HOXB','HOXC','HOXD',
        'CDH1','PTEN','DAB2IP','IRF8','RUNX3',
    ],
}

lineage_out = gs.ssgsea(
    data=ctrl_log,
    gene_sets=PROSTATE_LINEAGE_SIGS,
    outdir=None,
    no_plot=True,
    processes=1,
)

ss_lin = lineage_out.res2d.pivot(index='Term', columns='Name', values='NES')
ss_lin = ss_lin[ctrl_geno_order]   # reuse column order from DDR heatmap

# Pool lineage replicates → z-score → render
ss_lin_pooled = _pool_replicates(ss_lin, ctrl_meta, ctrl_geno_order)
_lp_mean  = ss_lin_pooled.mean(axis=1)
_lp_std   = ss_lin_pooled.std(axis=1).replace(0, 1)
ss_lin_z  = ss_lin_pooled.subtract(_lp_mean, axis=0).divide(_lp_std, axis=0).dropna(how='all')

if ss_lin_z.empty:
    print("  No valid lineage gene sets after ssGSEA — skipping prostate lineage heatmap")
else:
    _PL = f'{OUTDIR}/figures/08_prostate_lineage'
    os.makedirs(_PL, exist_ok=True)
    _render_heatmap(
        ss_lin_z,
        'Prostate Lineage Identity — Control, Replicates Pooled\n(ssGSEA, z-scored per gene set)',
        f'{_PL}/prostate_lineage_ssgsea_heatmap',
        plt.cm.PuOr_r,
    )

print(f"\n✓ All done. Results in: {OUTDIR}")
