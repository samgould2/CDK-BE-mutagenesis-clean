"""
regen_ddr_heatmaps.py
─────────────────────
Standalone script: re-generates all DDR GSEA heatmaps
from the existing per-comparison DESeq2 CSVs (no need to re-run
DESeq2 or the full pipeline).

Run from the project root:
    cd ~/Desktop/"Comprehensive PC3 RNA-seq Analysis"
    pip3 install gseapy   # if not already installed
    python3 scripts/regen_ddr_heatmaps.py 2>&1 | grep -v findfont

Outputs (overwrite in place):
    pc3_rnaseq_results/figures/04_gsea_cdk12/GSEA_DDR_heatmap_drug_HQ461.{pdf,png,svg}
    pc3_rnaseq_results/figures/04_gsea_cdk12/GSEA_DDR_heatmap_drug_ZLY025.{pdf,png,svg}
    pc3_rnaseq_results/figures/04_gsea_cdk12/GSEA_DDR_heatmap_variant_effect.{pdf,png,svg}
    pc3_rnaseq_results/figures/04_gsea_cdk12/GSEA_DDR_heatmap_variant_G731K.{pdf,png,svg}
    pc3_rnaseq_results/figures/04_gsea_cdk12/GSEA_DDR_heatmap_variant_G822E.{pdf,png,svg}
    pc3_rnaseq_results/figures/04_gsea_cdk12/GSEA_DDR_heatmap_variant_I733V.{pdf,png,svg}
"""

import os
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import gseapy as gs

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams['font.family']      = 'sans-serif'
plt.rcParams['font.sans-serif']  = ['Liberation Sans', 'Arial', 'Helvetica Neue',
                                     'Helvetica', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype']     = 42
plt.rcParams['ps.fonttype']      = 42
plt.rcParams['svg.fonttype']     = 'none'

# ── Paths ─────────────────────────────────────────────────────────────────────
_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(_ROOT, 'pc3_rnaseq_results')
RESDIR = os.path.join(OUTDIR, 'deseq2_results')
FIGDIR = os.path.join(OUTDIR, 'figures', '04_gsea_cdk12')
os.makedirs(FIGDIR, exist_ok=True)

# ── CSV lookup: searches subdirectories ───────────────────────────────────────
def find_csv(fname, base_dir):
    """Return path to {fname}.csv searching base_dir and one level of subdirs."""
    direct = os.path.join(base_dir, f'{fname}.csv')
    if os.path.exists(direct):
        return direct
    for entry in os.scandir(base_dir):
        if entry.is_dir():
            candidate = os.path.join(entry.path, f'{fname}.csv')
            if os.path.exists(candidate):
                return candidate
    raise FileNotFoundError(f"Cannot find {fname}.csv under {base_dir}. "
                            f"Run pc3_deseq2_gsea.py first to generate CSVs.")

# ── Comparison map: CSV filename → display label ──────────────────────────────
COMPARISONS = [
    ('NT_HQ461_vs_NT_Control',         'NT: HQ461 vs Control'),
    ('I733V_HQ461_vs_I733V_Control',   'I733V: HQ461 vs Control'),
    ('G822E_HQ461_vs_G822E_Control',   'G822E: HQ461 vs Control'),
    ('G731K_HQ461_vs_G731K_Control',   'G731K: HQ461 vs Control'),
    ('NT_ZLY025_vs_NT_Control',        'NT: ZLY-025 vs Control'),
    ('I733V_ZLY025_vs_I733V_Control',  'I733V: ZLY-025 vs Control'),
    ('G822E_ZLY025_vs_G822E_Control',  'G822E: ZLY-025 vs Control'),
    ('G731K_ZLY025_vs_G731K_Control',  'G731K: ZLY-025 vs Control'),
    ('NT_vs_I733V_HQ461',              'HQ461: NT vs I733V'),
    ('NT_vs_G822E_HQ461',              'HQ461: NT vs G822E'),
    ('NT_vs_G731K_HQ461',              'HQ461: NT vs G731K'),
    ('NT_vs_I733V_ZLY025',             'ZLY-025: NT vs I733V'),
    ('NT_vs_G822E_ZLY025',             'ZLY-025: NT vs G822E'),
    ('NT_vs_G731K_ZLY025',             'ZLY-025: NT vs G731K'),
    ('NT_vs_I733V_Control',            'Control: NT vs I733V'),
    ('NT_vs_G822E_Control',            'Control: NT vs G822E'),
    ('NT_vs_G731K_Control',            'Control: NT vs G731K'),
]

# ── DDR gene sets ─────────────────────────────────────────────────────────────
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

# ── Load DESeq2 CSVs ──────────────────────────────────────────────────────────
print("Loading DESeq2 results...")
results = {}
for fname, title in COMPARISONS:
    path = find_csv(fname, RESDIR)
    df = pd.read_csv(path, index_col=0)
    results[fname] = (df, title)
print(f"  Loaded {len(results)} comparisons.\n")

# ── GSEA function ─────────────────────────────────────────────────────────────
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

# ── Heatmap function ──────────────────────────────────────────────────────────
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
            stars = ('***' if fdr < 0.001 else '**' if fdr < 0.01
                     else '*' if fdr < 0.05 else ('+' if nom < 0.05 else ''))
            text_color = 'white' if abs(nes) > vmax * 0.6 else 'black'
            ax.text(i + 0.5, n_rows - j - 0.5, f'{nes:.2f}{stars}',
                    ha='center', va='center', fontsize=8.5, color=text_color, zorder=2)

    for i, gs_name in enumerate(sig_order):
        ax.text(i + 0.5, -0.15, gs_name, ha='right', va='top', fontsize=9, rotation=30)
    for j, comp in enumerate(comp_labels):
        ax.text(-0.05, n_rows - j - 0.5, comp, ha='right', va='center', fontsize=9)

    if group_separators:
        for sep in group_separators:
            y = n_rows - sep
            ax.axhline(y=y, xmin=0, xmax=1, color='white', linewidth=2, zorder=3)
            ax.axhline(y=y, xmin=0, xmax=1, color='black', linewidth=0.8,
                       linestyle='--', zorder=4)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, aspect=15)
    cbar.set_label('NES', fontsize=9)

    ax.set_title(title, fontsize=11, pad=20)
    plt.tight_layout(pad=2.0)
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.{{pdf,png,svg}}")

# ── Run DDR GSEA ──────────────────────────────────────────────────────────────
print("Running DDR GSEA (1000 permutations × 17 comparisons)...")
ddr_gsea = run_gsea(results, DDR_SETS, 'DDR')

_DDR_ORDER  = ['GOBP_DSB_REPAIR_VIA_HR', 'KEGG_FANCONI_ANEMIA_PATHWAY',
               'GOBP_DOUBLE_STRAND_BREAK_REPAIR', 'GOBP_DNA_REPAIR']
_DDR_SUFFIX = '\n(negative NES = enriched in second/denominator condition)'

# Panel A — Drug effect: HQ461 + ZLY-025 combined (8 rows, separator at 4)
print("\nPlotting Panel A: Drug Effect — combined...")
_drug_keys = [
    'NT: HQ461 vs Control', 'I733V: HQ461 vs Control',
    'G822E: HQ461 vs Control', 'G731K: HQ461 vs Control',
    'NT: ZLY-025 vs Control', 'I733V: ZLY-025 vs Control',
    'G822E: ZLY-025 vs Control', 'G731K: ZLY-025 vs Control',
]
_drug = {k: ddr_gsea[k] for k in _drug_keys}
plot_gsea_heatmap(_drug, _DDR_ORDER,
    os.path.join(FIGDIR, 'GSEA_DDR_heatmap_drug_effect'),
    f'DDR GSEA — Drug Effect (vs Control){_DDR_SUFFIX}',
    group_separators=[4])

# Panel C — Variant effect: NT vs variant, across all conditions (9 rows)
print("\nPlotting Panel C: Variant Effect (all 9 rows)...")
_var_keys = [
    'HQ461: NT vs I733V', 'HQ461: NT vs G822E', 'HQ461: NT vs G731K',
    'ZLY-025: NT vs I733V', 'ZLY-025: NT vs G822E', 'ZLY-025: NT vs G731K',
    'Control: NT vs I733V', 'Control: NT vs G822E', 'Control: NT vs G731K',
]
_var = {k: ddr_gsea[k] for k in _var_keys}
plot_gsea_heatmap(_var, _DDR_ORDER,
    os.path.join(FIGDIR, 'GSEA_DDR_heatmap_variant_effect'),
    f'DDR GSEA — Variant Effect (NT vs variant){_DDR_SUFFIX}',
    group_separators=[3, 6])

# Panels C1–C3 — Variant effect split by drug condition (3 rows each) ──────
print("\nPlotting split variant heatmaps (one per drug condition)...")

# HQ461 + ZLY-025 combined (6 rows, separator at 3)
_drug_combined_keys = [
    'HQ461: NT vs G731K', 'HQ461: NT vs G822E', 'HQ461: NT vs I733V',
    'ZLY-025: NT vs G731K', 'ZLY-025: NT vs G822E', 'ZLY-025: NT vs I733V',
]
_drug_combined = {k: ddr_gsea[k] for k in _drug_combined_keys}
plot_gsea_heatmap(_drug_combined, _DDR_ORDER,
    os.path.join(FIGDIR, 'GSEA_DDR_heatmap_variant_drug'),
    f'DDR GSEA — Variant Effect (HQ461 + ZLY-025){_DDR_SUFFIX}',
    group_separators=[3])

# Control (3 rows)
_ctrl_keys = ['Control: NT vs G731K', 'Control: NT vs G822E', 'Control: NT vs I733V']
_ctrl = {k: ddr_gsea[k] for k in _ctrl_keys}
plot_gsea_heatmap(_ctrl, _DDR_ORDER,
    os.path.join(FIGDIR, 'GSEA_DDR_heatmap_variant_Control'),
    f'DDR GSEA — Variant Effect (Control){_DDR_SUFFIX}')

print("\nDone. All panels written to:")
print(f"  {FIGDIR}/")
