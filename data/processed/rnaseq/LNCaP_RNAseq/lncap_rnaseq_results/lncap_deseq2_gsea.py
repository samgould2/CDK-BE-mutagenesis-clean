"""
LNCaP CDK12 I733V RNA-seq: DESeq2 + GSEA
Inputs:
    - salmon.merged.gene_counts.tsv  (nf-core/rnaseq output)
    - deseq2.pca.vals.txt            (nf-core/rnaseq QC output, optional)
Outputs:
    - lncap_rnaseq_results/deseq2_results/  (per-comparison CSVs)
    - lncap_rnaseq_results/figures/         (volcano plots, PCA, GSEA heatmaps)

Install dependencies:
    pip install pydeseq2 gseapy adjustText pandas numpy matplotlib
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings('ignore')

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rc('font', family='Liberation Sans')
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ── Paths ─────────────────────────────────────────────────────────────────────
COUNTS_FILE = os.path.expanduser('~/Desktop/LNCAP RNA-Seq Analysis/salmon.merged.gene_counts.tsv')
PCA_FILE    = os.path.expanduser('~/Desktop/LNCAP RNA-Seq Analysis/deseq2.pca.vals.txt')
OUTDIR      = os.path.expanduser('~/Desktop/LNCAP RNA-Seq Analysis/lncap_rnaseq_results')
os.makedirs(f'{OUTDIR}/deseq2_results', exist_ok=True)
os.makedirs(f'{OUTDIR}/figures', exist_ok=True)

# ── Load counts ───────────────────────────────────────────────────────────────
print("Loading count matrix...")
counts_df = pd.read_csv(COUNTS_FILE, sep='\t', index_col=0)
gene_names = counts_df['gene_name'].copy()
counts_df  = counts_df.drop(columns=['gene_name']).round().astype(int)
print(f"  {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")

# ── Sample metadata ───────────────────────────────────────────────────────────
metadata = pd.DataFrame({
    'genotype':  ['I733V']*6 + ['NT']*6,
    'treatment': ['Control']*3 + ['HQ461']*3 + ['Control']*3 + ['HQ461']*3,
}, index=counts_df.columns)
metadata['condition'] = metadata['genotype'] + '_' + metadata['treatment']

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
    ('NT_HQ461',    'NT_Control',    'NT_HQ461_vs_NT_Control',          'NT: HQ461 vs Control'),
    ('I733V_HQ461', 'I733V_Control', 'I733V_HQ461_vs_I733V_Control',    'Cdk12 I733V: HQ461 vs Control'),
    ('NT_HQ461',    'I733V_HQ461',   'NT_HQ461_vs_I733V_HQ461',         'HQ461: NT vs Cdk12 I733V'),
    ('NT_Control',  'I733V_Control', 'NT_Control_vs_I733V_Control',     'Control: NT vs Cdk12 I733V'),
]

highlight_genes = {}  # No gene labels on volcano plots

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

# ── Volcano plots ─────────────────────────────────────────────────────────────
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    print("adjustText not installed — labels may overlap. Install with: pip install adjustText")

LFC_THRESH  = 1.0
PADJ_THRESH = 0.05

# Compute shared axis limits across all comparisons
all_lfc, all_y = [], []
for name, (res, title) in results.items():
    rc = res.dropna(subset=['padj', 'log2FoldChange'])
    all_lfc.extend(rc['log2FoldChange'].values)
    all_y.extend(-np.log10(rc['padj'].values))

x_lim = np.ceil(max(abs(np.nanmin(all_lfc)), abs(np.nanmax(all_lfc)))) + 0.5
y_lim = np.ceil(np.nanmax(all_y)) + 2

fig, axes = plt.subplots(1, 4, figsize=(20, 5))

for ax, (name, (res, title)) in zip(axes, results.items()):
    rc = res.dropna(subset=['padj', 'log2FoldChange'])
    up   = (rc['padj'] < PADJ_THRESH) & (rc['log2FoldChange'] >  LFC_THRESH)
    down = (rc['padj'] < PADJ_THRESH) & (rc['log2FoldChange'] < -LFC_THRESH)
    ns   = ~(up | down)

    ax.scatter(rc.loc[ns,   'log2FoldChange'], -np.log10(rc.loc[ns,   'padj']),
               s=4, color="#AAAAAA", alpha=0.4, linewidths=0, rasterized=True)
    ax.scatter(rc.loc[up,   'log2FoldChange'], -np.log10(rc.loc[up,   'padj']),
               s=25, color="tab:red",  alpha=0.7, linewidths=0.3,
               edgecolors='black', rasterized=True,
               label=f'{up.sum()} DE UP')
    ax.scatter(rc.loc[down, 'log2FoldChange'], -np.log10(rc.loc[down, 'padj']),
               s=25, color="tab:blue", alpha=0.7, linewidths=0.3,
               edgecolors='black', rasterized=True,
               label=f'{down.sum()} DE DOWN')

    ax.set_xlim(-x_lim, x_lim)
    ax.set_ylim(0, y_lim)

    ax.legend(fontsize=12, frameon=True, framealpha=0.9, edgecolor='#CCCCCC')
    ax.axvline( LFC_THRESH, color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.axvline(-LFC_THRESH, color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.axhline(-np.log10(PADJ_THRESH), color='black', lw=0.5, linestyle='--', alpha=0.4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xlabel('log₂ Fold Change', fontsize=13)
    ax.set_ylabel('-log₁₀(padj)', fontsize=13)
    ax.set_title(title, fontsize=13)

plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/volcano_all_comparisons.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/volcano_all_comparisons.png', dpi=300, bbox_inches='tight')
plt.close()
print("Volcano plots saved.")

# ── PCA ───────────────────────────────────────────────────────────────────────
if os.path.exists(PCA_FILE):
    pca = pd.read_csv(PCA_FILE, sep='\t')
    pc1_col = [c for c in pca.columns if 'PC1' in c][0]
    pc2_col = [c for c in pca.columns if 'PC2' in c][0]

    pca['condition'] = pca.iloc[:, 0].apply(
        lambda s: 'NT_Control'    if 'NT_Control'    in str(s) else
                  'NT_HQ461'      if 'NT_HQ461'      in str(s) else
                  'I733V_Control' if 'I733V_Control' in str(s) else
                  'I733V_HQ461'   if 'I733V_HQ461'  in str(s) else 'other'
    )

    colors_map   = {'NT_Control': '#2166ac', 'NT_HQ461': '#92c5de',
                    'I733V_Control': '#d6604d', 'I733V_HQ461': '#f4a582'}
    cond_labels  = {'NT_Control': 'NT Control', 'NT_HQ461': 'NT + HQ461',
                    'I733V_Control': 'I733V Control', 'I733V_HQ461': 'I733V + HQ461'}

    fig, ax = plt.subplots(figsize=(5, 4.5))
    for cond, grp in pca.groupby('condition'):
        ax.scatter(grp[pc1_col], grp[pc2_col],
                   label=cond_labels.get(cond, cond),
                   color=colors_map.get(cond, 'gray'), s=80, zorder=3)

    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel(pc1_col, fontsize=11)
    ax.set_ylabel(pc2_col, fontsize=11)
    ax.set_title('PCA — LNCaP CDK12 RNA-seq', fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/figures/PCA.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{OUTDIR}/figures/PCA.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("PCA saved.")

# ── GSEA ──────────────────────────────────────────────────────────────────────
import gseapy as gs

# Published curated DDR gene sets
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

# Custom CDK12/13 signatures
CDK12_LOSS_DN       = ["BRCA1","BRCA2","ATR","ATRIP","TOPBP1","ETAA1","FANCA","FANCC","FANCD2","FANCI","FANCL","FANCM","RAD51","RAD51C","RAD51D","PALB2","BARD1","NBN","MRE11","RAD50","CHEK1","RPA1","RPA2","RPA3","RFC1","RFC2","RFC3","RFC4","RFC5","PCNA","POLD1","POLD2","POLD3","POLE","POLE2","LIG1","LIG3","XRCC1","XRCC2","XRCC3","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC45","CDT1","ORC1","EME1","MUS81","SLX4","GEN1","BLM","WRN","EXO1","DNA2"]
CDK12_LOSS_UP       = ["CGAS","STING1","IRF3","IRF7","IFNB1","IFNA1","IFNA2","MX1","MX2","OAS1","OAS2","OAS3","OASL","ISG15","ISG20","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","STAT1","STAT2","IRF9","DDX58","IFIH1","TRIM56","TRIM32","RSAD2","HERC5","XAF1"]
CDK13_LOSS_DN       = ["MCL1","CCND1","CDK4","CDK6","MYC","MYCN","SRSF1","SRSF2","SRSF3","SRSF7","U2AF1","U2AF2","SF3B1","ELL","ELL2","AFF1","AFF4","CDC7","TIPIN","CLASPIN","E2F1","E2F2","E2F3"]
CDK12_CDK13_SHARED  = ["CCNK","BRCA1","BRCA2","ATR","FANCD2","FANCI","CHEK1","MCL1","CCND1","MYC","WEE1","CDC25A","POLR2A","SUPT5H","SUPT4H1","PAF1","CTR9","LEO1","CDC73","SPT6"]
HQ461_DN            = sorted(set(CDK12_LOSS_DN + CDK13_LOSS_DN + CDK12_CDK13_SHARED))

CUSTOM_SETS = {
    "CDK12_LOSS_DN":        CDK12_LOSS_DN,
    "CDK12_LOSS_UP":        CDK12_LOSS_UP,
    "CDK13_LOSS_DN":        CDK13_LOSS_DN,
    "CDK12_CDK13_SHARED_DN": CDK12_CDK13_SHARED,
    "HQ461_CDK12_CDK13_DN": HQ461_DN,
}

def run_gsea(results_dict, gene_sets, label):
    """Run preranked GSEA for all comparisons and return NES/FDR DataFrames."""
    gsea_res = {}
    for name, (res, title) in results_dict.items():
        rc = res.dropna(subset=['stat','gene_name']).drop_duplicates(subset='gene_name')
        ranking = rc[['gene_name','stat']].set_index('gene_name').sort_values('stat', ascending=False)
        pre = gs.prerank(
            rnk=ranking, gene_sets=gene_sets,
            permutation_num=1000, seed=42,
            min_size=5, max_size=500, verbose=False, outdir=None,
        )
        gsea_res[title] = pre.res2d.set_index('Term')
        print(f"  {label} GSEA done: {title}")
    return gsea_res

def plot_gsea_heatmap(gsea_res, sig_order, outpath, title):
    comp_labels = list(gsea_res.keys())
    nes_df = pd.DataFrame({l: gsea_res[l]['NES']       for l in comp_labels}).loc[sig_order]
    fdr_df = pd.DataFrame({l: gsea_res[l]['FDR q-val'] for l in comp_labels}).loc[sig_order]
    nom_df = pd.DataFrame({l: gsea_res[l]['NOM p-val'] for l in comp_labels}).loc[sig_order]

    vmax = max(abs(nes_df.values.max()), abs(nes_df.values.min()))
    fig, ax = plt.subplots(figsize=(max(7, len(sig_order) * 1.7), 3.8))
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im   = ax.imshow(nes_df.values.T, cmap='RdBu_r', norm=norm, aspect='auto')

    ax.set_xticks(range(len(sig_order)))
    ax.set_xticklabels(sig_order, rotation=30, ha='right', fontsize=9)
    ax.set_yticks(range(len(comp_labels)))
    ax.set_yticklabels(comp_labels, fontsize=9)

    for i in range(len(sig_order)):
        for j in range(len(comp_labels)):
            nes = nes_df.values[i, j]
            fdr = fdr_df.values[i, j]
            nom = nom_df.values[i, j]
            stars = '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else ('+' if nom < 0.05 else '')
            ax.text(i, j, f'{nes:.2f}{stars}', ha='center', va='center', fontsize=8.5,
                    color='white' if abs(nes) > vmax * 0.6 else 'black')

    plt.colorbar(im, ax=ax, shrink=0.7, label='NES')
    ax.set_title(title, fontsize=11)
    plt.tight_layout(pad=2.0)
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

print("\nRunning GSEA — DDR (published gene sets)...")
ddr_gsea = run_gsea(results, DDR_SETS, 'DDR')
plot_gsea_heatmap(
    ddr_gsea,
    sig_order=['GOBP_DSB_REPAIR_VIA_HR','KEGG_FANCONI_ANEMIA_PATHWAY',
               'GOBP_DOUBLE_STRAND_BREAK_REPAIR','GOBP_DNA_REPAIR'],
    outpath=f'{OUTDIR}/figures/GSEA_DDR_heatmap',
    title='DDR Gene Set GSEA — LNCaP CDK12 HQ461\n(negative NES = enriched in second condition)',
)

print("\nRunning GSEA — CDK12/13 custom signatures...")
custom_gsea = run_gsea(results, CUSTOM_SETS, 'CDK12/13')
plot_gsea_heatmap(
    custom_gsea,
    sig_order=["CDK12_LOSS_DN","CDK13_LOSS_DN","CDK12_CDK13_SHARED_DN",
               "HQ461_CDK12_CDK13_DN","CDK12_LOSS_UP"],
    outpath=f'{OUTDIR}/figures/GSEA_CDK12_13_heatmap',
    title='CDK12/13 Signature GSEA — LNCaP CDK12 HQ461\n(negative NES = enriched in second condition)',
)

# ── Hallmark GSEA dot plot — NT: Control vs HQ461 ────────────────────────────
print("\nRunning GSEA — MSigDB Hallmarks (NT: Control vs HQ461)...")

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

# Rank NT HQ461 vs Control
nt_res = results['NT_HQ461_vs_NT_Control'][0]
nt_rc  = nt_res.dropna(subset=['stat','gene_name']).drop_duplicates(subset='gene_name')
nt_rank = nt_rc[['gene_name','stat']].set_index('gene_name').sort_values('stat', ascending=False)

hallmark_pre = gs.prerank(
    rnk=nt_rank, gene_sets=HALLMARKS,
    permutation_num=1000, seed=42,
    min_size=5, max_size=500, verbose=False, outdir=None,
)
hdf = hallmark_pre.res2d.copy()
hdf['Term_clean'] = hdf['Term'].str.replace('HALLMARK_','').str.replace('_',' ').str.title()
hdf['neg_log10_fdr'] = -np.log10(hdf['FDR q-val'].clip(1e-3))
hdf_plot = hdf.sort_values('NES').reset_index(drop=True)

max_y  = hdf_plot['neg_log10_fdr'].max()
sizes  = 40 + 310 * (hdf_plot['neg_log10_fdr'] / max_y)

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
ax.set_title('WT LNCaP: Control vs HQ461\nMSigDB Hallmark GSEA', fontsize=13)
ax.tick_params(axis='x', labelsize=12)
ax.set_xlim(-3, 3)

for fdr_val, label in [(0.3, 'FDR=0.30'), (0.1, 'FDR=0.10'), (0.01, 'FDR=0.01')]:
    s = 40 + 310 * (-np.log10(fdr_val) / max_y)
    ax.scatter([], [], s=s, c='grey', alpha=0.7, edgecolors='black', linewidths=0.4, label=label)
ax.legend(title='Dot size\n= -log₁₀(FDR)', fontsize=9, title_fontsize=9,
          frameon=True, framealpha=0.9, edgecolor='#CCCCCC',
          loc='lower right', handletextpad=1.0)

plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/hallmark_gsea_dotplot.pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{OUTDIR}/figures/hallmark_gsea_dotplot.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {OUTDIR}/figures/hallmark_gsea_dotplot")

print(f"\n✓ All done. Results in: {OUTDIR}")
