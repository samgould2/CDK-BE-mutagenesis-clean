"""
make_gsea_svgs.py
─────────────────
Regenerates SVG (+ PDF + PNG) for the Hallmark and Reactome GSEA dot plots,
and the CDK12/13 custom-signature heatmap, from the existing DESeq2 CSVs.

Run from the project root:
    python3 scripts/make_gsea_svgs.py 2>&1 | grep -v findfont

Outputs (overwrite in place):
    figures/04_gsea_cdk12/GSEA_CDK12_13_heatmap.{pdf,png,svg}
    figures/05_hallmark_gsea/hallmark_gsea_*.{pdf,png,svg}
    figures/06_reactome/reactome_*.{pdf,png,svg}
"""

import os, warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import gseapy as gs

plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'Helvetica Neue',
                                    'Helvetica', 'DejaVu Sans', 'sans-serif']
plt.rcParams['pdf.fonttype']    = 42
plt.rcParams['ps.fonttype']      = 42
plt.rcParams['svg.fonttype']     = 'none'

_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(_ROOT, 'pc3_rnaseq_results')
RESDIR = os.path.join(OUTDIR, 'deseq2_results')

def find_csv(fname, base_dir):
    direct = os.path.join(base_dir, f'{fname}.csv')
    if os.path.exists(direct):
        return direct
    for entry in os.scandir(base_dir):
        if entry.is_dir():
            candidate = os.path.join(entry.path, f'{fname}.csv')
            if os.path.exists(candidate):
                return candidate
    raise FileNotFoundError(f"Cannot find {fname}.csv under {base_dir}")

COMP_KEYS = [
    'NT_HQ461_vs_NT_Control', 'NT_ZLY025_vs_NT_Control',
    'NT_vs_I733V_Control', 'NT_vs_G822E_Control', 'NT_vs_G731K_Control',
    'NT_vs_I733V_HQ461', 'NT_vs_G822E_HQ461', 'NT_vs_G731K_HQ461',
    'NT_vs_I733V_ZLY025', 'NT_vs_G822E_ZLY025', 'NT_vs_G731K_ZLY025',
    'I733V_HQ461_vs_I733V_Control', 'G822E_HQ461_vs_G822E_Control', 'G731K_HQ461_vs_G731K_Control',
    'I733V_ZLY025_vs_I733V_Control', 'G822E_ZLY025_vs_G822E_Control', 'G731K_ZLY025_vs_G731K_Control',
]

COMP_LABELS = {
    'NT_HQ461_vs_NT_Control':         'NT: HQ461 vs Control',
    'NT_ZLY025_vs_NT_Control':        'NT: ZLY-025 vs Control',
    'NT_vs_I733V_Control':            'NT vs I733V — no drug',
    'NT_vs_G822E_Control':            'NT vs G822E — no drug',
    'NT_vs_G731K_Control':            'NT vs G731K — no drug',
    'NT_vs_I733V_HQ461':              'HQ461: NT vs I733V',
    'NT_vs_G822E_HQ461':              'HQ461: NT vs G822E',
    'NT_vs_G731K_HQ461':              'HQ461: NT vs G731K',
    'NT_vs_I733V_ZLY025':             'ZLY-025: NT vs I733V',
    'NT_vs_G822E_ZLY025':             'ZLY-025: NT vs G822E',
    'NT_vs_G731K_ZLY025':             'ZLY-025: NT vs G731K',
    'I733V_HQ461_vs_I733V_Control':   'I733V: HQ461 vs Control',
    'G822E_HQ461_vs_G822E_Control':   'G822E: HQ461 vs Control',
    'G731K_HQ461_vs_G731K_Control':   'G731K: HQ461 vs Control',
    'I733V_ZLY025_vs_I733V_Control':  'I733V: ZLY-025 vs Control',
    'G822E_ZLY025_vs_G822E_Control':  'G822E: ZLY-025 vs Control',
    'G731K_ZLY025_vs_G731K_Control':  'G731K: ZLY-025 vs Control',
}

print("Loading DESeq2 results...")
results = {}
for fname in COMP_KEYS:
    path = find_csv(fname, RESDIR)
    results[fname] = (pd.read_csv(path, index_col=0), COMP_LABELS[fname])
print(f"  Loaded {len(results)} comparisons.\n")

# ── CDK12/13 custom sets ──────────────────────────────────────────────────────
CDK12_LOSS_DN      = ["BRCA1","BRCA2","ATR","ATRIP","TOPBP1","ETAA1","FANCA","FANCC","FANCD2","FANCI","FANCL","FANCM","RAD51","RAD51C","RAD51D","PALB2","BARD1","NBN","MRE11","RAD50","CHEK1","RPA1","RPA2","RPA3","RFC1","RFC2","RFC3","RFC4","RFC5","PCNA","POLD1","POLD2","POLD3","POLE","POLE2","LIG1","LIG3","XRCC1","XRCC2","XRCC3","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC45","CDT1","ORC1","EME1","MUS81","SLX4","GEN1","BLM","WRN","EXO1","DNA2"]
CDK12_LOSS_UP      = ["CGAS","STING1","IRF3","IRF7","IFNB1","IFNA1","IFNA2","MX1","MX2","OAS1","OAS2","OAS3","OASL","ISG15","ISG20","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","STAT1","STAT2","IRF9","DDX58","IFIH1","TRIM56","TRIM32","RSAD2","HERC5","XAF1"]
CDK13_LOSS_DN      = ["MCL1","CCND1","CDK4","CDK6","MYC","MYCN","SRSF1","SRSF2","SRSF3","SRSF7","U2AF1","U2AF2","SF3B1","ELL","ELL2","AFF1","AFF4","CDC7","TIPIN","CLASPIN","E2F1","E2F2","E2F3"]
CDK12_CDK13_SHARED = ["CCNK","BRCA1","BRCA2","ATR","FANCD2","FANCI","CHEK1","MCL1","CCND1","MYC","WEE1","CDC25A","POLR2A","SUPT5H","SUPT4H1","PAF1","CTR9","LEO1","CDC73","SPT6"]
HQ461_DN           = sorted(set(CDK12_LOSS_DN + CDK13_LOSS_DN + CDK12_CDK13_SHARED))

CUSTOM_SETS = {
    "CDK12_LOSS_DN":         CDK12_LOSS_DN,
    "CDK12_LOSS_UP":         CDK12_LOSS_UP,
    "CDK13_LOSS_DN":         CDK13_LOSS_DN,
    "CDK12_CDK13_SHARED_DN": CDK12_CDK13_SHARED,
    "HQ461_CDK12_CDK13_DN":  HQ461_DN,
}

def run_gsea(results_dict, gene_sets, label, max_size=500):
    gsea_res = {}
    for name, (res, title) in results_dict.items():
        rc = res.dropna(subset=['stat', 'gene_name']).drop_duplicates(subset='gene_name')
        ranking = rc[['gene_name', 'stat']].set_index('gene_name').sort_values('stat', ascending=False)
        pre = gs.prerank(rnk=ranking, gene_sets=gene_sets,
                         permutation_num=1000, seed=42,
                         min_size=5, max_size=max_size, verbose=False, outdir=None)
        gsea_res[title] = pre.res2d.set_index('Term')
        print(f"  {label} done: {title}")
    return gsea_res

def plot_gsea_heatmap(gsea_res, sig_order, outpath, title, group_separators=None):
    comp_labels = list(gsea_res.keys())
    nes_df = pd.DataFrame({l: gsea_res[l]['NES']       for l in comp_labels}).loc[sig_order]
    fdr_df = pd.DataFrame({l: gsea_res[l]['FDR q-val'] for l in comp_labels}).loc[sig_order]
    nom_df = pd.DataFrame({l: gsea_res[l]['NOM p-val'] for l in comp_labels}).loc[sig_order]
    vmax = max(abs(nes_df.values.max()), abs(nes_df.values.min()))
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu_r
    n_rows, n_cols = len(comp_labels), len(sig_order)
    cell_w, cell_h = 1.7, 0.7
    fig, ax = plt.subplots(figsize=(n_cols * cell_w + 2, n_rows * cell_h + 1.5))
    ax.set_xlim(0, n_cols); ax.set_ylim(0, n_rows)
    ax.set_aspect('equal'); ax.axis('off')
    for i, gene_set in enumerate(sig_order):
        for j, comp in enumerate(comp_labels):
            nes = nes_df.loc[gene_set, comp]
            fdr = fdr_df.loc[gene_set, comp]
            nom = nom_df.loc[gene_set, comp]
            color = cmap(norm(nes))
            ax.add_patch(plt.Rectangle((i, n_rows - j - 1), 1, 1, color=color, zorder=1))
            stars = '***' if fdr < 0.001 else '**' if fdr < 0.01 else '*' if fdr < 0.05 else ('+' if nom < 0.05 else '')
            tc = 'white' if abs(nes) > vmax * 0.6 else 'black'
            ax.text(i + 0.5, n_rows - j - 0.5, f'{nes:.2f}{stars}',
                    ha='center', va='center', fontsize=8.5, color=tc, zorder=2)
    for i, gs_name in enumerate(sig_order):
        ax.text(i + 0.5, -0.15, gs_name, ha='right', va='top', fontsize=9, rotation=30)
    for j, comp in enumerate(comp_labels):
        ax.text(-0.05, n_rows - j - 0.5, comp, ha='right', va='center', fontsize=9)
    if group_separators:
        for sep in group_separators:
            y = n_rows - sep
            ax.axhline(y=y, xmin=0, xmax=1, color='white', linewidth=2, zorder=3)
            ax.axhline(y=y, xmin=0, xmax=1, color='black', linewidth=0.8, linestyle='--', zorder=4)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02, aspect=15); cbar.set_label('NES', fontsize=9)
    ax.set_title(title, fontsize=11, pad=20)
    plt.tight_layout(pad=2.0)
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.{{pdf,png,svg}}")


def plot_gsea_heatmap_T(gsea_res, sig_order, outpath, title, col_label_sep='|'):
    """Transposed: gene sets = rows, comparisons = columns."""
    comp_labels = list(gsea_res.keys())
    nes_df = pd.DataFrame({l: gsea_res[l]['NES']       for l in comp_labels}).loc[sig_order]
    fdr_df = pd.DataFrame({l: gsea_res[l]['FDR q-val'] for l in comp_labels}).loc[sig_order]
    nom_df = pd.DataFrame({l: gsea_res[l]['NOM p-val'] for l in comp_labels}).loc[sig_order]
    vmax = max(abs(nes_df.values.max()), abs(nes_df.values.min()))
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu_r
    n_rows, n_cols = len(sig_order), len(comp_labels)
    cell_h, cell_w = 0.65, 1.8
    fig, ax = plt.subplots(figsize=(n_cols * cell_w + 4.5, n_rows * cell_h + 2.2))
    ax.set_xlim(0, n_cols); ax.set_ylim(0, n_rows)
    ax.set_aspect('equal'); ax.axis('off')
    clean = lambda s: s.replace('HALLMARK_','').replace('_',' ').title()
    for i, gs_name in enumerate(sig_order):
        row_idx = n_rows - i - 1
        for j, comp in enumerate(comp_labels):
            nes = nes_df.loc[gs_name, comp]; fdr = fdr_df.loc[gs_name, comp]; nom = nom_df.loc[gs_name, comp]
            ax.add_patch(plt.Rectangle((j, row_idx), 1, 1, color=cmap(norm(nes)), zorder=1))
            stars = ('***' if fdr<0.001 else '**' if fdr<0.01 else '*' if fdr<0.05 else ('+' if nom<0.05 else ''))
            tc = 'white' if abs(nes) > vmax * 0.6 else 'black'
            ax.text(j+0.5, row_idx+0.5, f'{nes:.2f}{stars}', ha='center', va='center', fontsize=9, color=tc, zorder=2)
    for i, gs_name in enumerate(sig_order):
        label = clean(gs_name) if gs_name.startswith('HALLMARK_') else gs_name
        ax.text(-0.08, n_rows-i-0.5, label, ha='right', va='center', fontsize=10)
    for j, comp in enumerate(comp_labels):
        parts = comp.split(col_label_sep, 1)
        if len(parts) == 2:
            ax.text(j+0.5, n_rows+0.45, parts[0], ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(j+0.5, n_rows+0.05, parts[1], ha='center', va='bottom', fontsize=10, fontweight='bold')
        else:
            ax.text(j+0.5, n_rows+0.25, comp, ha='center', va='bottom', fontsize=10, fontweight='bold')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.02, aspect=18); cbar.set_label('NES', fontsize=10)
    ax.set_title(title, fontsize=11, pad=30)
    plt.tight_layout(pad=1.5)
    plt.savefig(f'{outpath}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.{{pdf,png,svg}}")

# ── CDK12/13 heatmap ─────────────────────────────────────────────────────────
print("Running CDK12/13 custom-set GSEA...")
FIGDIR_04 = os.path.join(OUTDIR, 'figures', '04_gsea_cdk12')
os.makedirs(FIGDIR_04, exist_ok=True)
custom_gsea = run_gsea(results, CUSTOM_SETS, 'CDK12/13')
plot_gsea_heatmap(
    custom_gsea,
    sig_order=["CDK12_LOSS_DN","CDK13_LOSS_DN","CDK12_CDK13_SHARED_DN",
               "HQ461_CDK12_CDK13_DN","CDK12_LOSS_UP"],
    outpath=os.path.join(FIGDIR_04, 'GSEA_CDK12_13_heatmap'),
    title='CDK12/13 Signature GSEA — PC3 CDK12 v2\n'
          '(negative NES = enriched in second/denominator condition)',
)

# ── Hallmark gene sets ────────────────────────────────────────────────────────
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

FIGDIR_05 = os.path.join(OUTDIR, 'figures', '05_hallmark_gsea')
os.makedirs(FIGDIR_05, exist_ok=True)

def hallmark_dotplot(res_key, res_label, outpath):
    res_r = results[res_key][0]
    rc = res_r.dropna(subset=['stat', 'gene_name']).drop_duplicates(subset='gene_name')
    rank = rc[['gene_name', 'stat']].set_index('gene_name').sort_values('stat', ascending=False)
    pre = gs.prerank(rnk=rank, gene_sets=HALLMARKS, permutation_num=1000, seed=42,
                     min_size=5, max_size=500, verbose=False, outdir=None)
    hdf = pre.res2d.copy()
    hdf['Term_clean']    = hdf['Term'].str.replace('HALLMARK_', '').str.replace('_', ' ').str.title()
    hdf['neg_log10_fdr'] = -np.log10(hdf['FDR q-val'].clip(1e-3))
    hdf_plot = hdf[hdf['FDR q-val'] < 0.25].sort_values('NES').reset_index(drop=True)
    if hdf_plot.empty:
        print(f"  No significant gene sets (FDR<0.25) for {res_label} — skipping")
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
    plt.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.{{pdf,png,svg}}")

print("\nRunning Hallmark GSEA...")
hallmark_dotplot('NT_HQ461_vs_NT_Control',  'NT: HQ461 vs Control',   f'{FIGDIR_05}/hallmark_gsea_NT_HQ461')
hallmark_dotplot('NT_ZLY025_vs_NT_Control', 'NT: ZLY-025 vs Control', f'{FIGDIR_05}/hallmark_gsea_NT_ZLY025')
hallmark_dotplot('NT_vs_I733V_Control', 'NT vs I733V — no drug', f'{FIGDIR_05}/hallmark_gsea_NT_vs_I733V_Control')
hallmark_dotplot('NT_vs_G822E_Control', 'NT vs G822E — no drug', f'{FIGDIR_05}/hallmark_gsea_NT_vs_G822E_Control')
hallmark_dotplot('NT_vs_G731K_Control', 'NT vs G731K — no drug', f'{FIGDIR_05}/hallmark_gsea_NT_vs_G731K_Control')

# ── Reactome ──────────────────────────────────────────────────────────────────
FIGDIR_06 = os.path.join(OUTDIR, 'figures', '06_reactome')
os.makedirs(FIGDIR_06, exist_ok=True)

def reactome_dotplot(res_key, res_label, outpath, gene_sets, top_n=30):
    res_r = results[res_key][0]
    rc = res_r.dropna(subset=['stat', 'gene_name']).drop_duplicates(subset='gene_name')
    rank = rc[['gene_name', 'stat']].set_index('gene_name').sort_values('stat', ascending=False)
    pre = gs.prerank(rnk=rank, gene_sets=gene_sets, permutation_num=1000, seed=42,
                     min_size=5, max_size=500, verbose=False, outdir=None)
    rdf = pre.res2d.copy()
    rdf['Term_clean'] = (rdf['Term'].str.replace(r'^REACTOME[_ ]', '', regex=True)
                         .str.replace('_', ' ').str.title())
    rdf['neg_log10_fdr'] = -np.log10(rdf['FDR q-val'].clip(1e-4))
    sig = rdf[rdf['FDR q-val'] < 0.25].copy()
    if sig.empty:
        print(f"  No significant Reactome terms (FDR<0.25) for {res_label} — skipping")
        return
    sig = sig.reindex(sig['NES'].abs().sort_values(ascending=False).index)
    sig = sig.head(top_n).sort_values('NES').reset_index(drop=True)
    max_y = sig['neg_log10_fdr'].max() or 1
    sizes = 40 + 310 * (sig['neg_log10_fdr'] / max_y)
    fig_h = max(6, len(sig) * 0.38)
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
    plt.savefig(f'{outpath}.svg', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}.{{pdf,png,svg}}")

# ── Hallmark heatmap — NT vs variants, Control (no drug) ──────────────────────
print("\nRunning Hallmark GSEA for NT vs variant Control comparisons...")
_ctrl_comps = {
    'NT_vs_G731K_Control': results['NT_vs_G731K_Control'],
    'NT_vs_G822E_Control': results['NT_vs_G822E_Control'],
    'NT_vs_I733V_Control': results['NT_vs_I733V_Control'],
}
_hallmark_ctrl_gsea = run_gsea(_ctrl_comps, HALLMARKS, 'Hallmark-Control')

_HALLMARK_ORDER = [
    'HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT',
    'HALLMARK_MYC_TARGETS_V1', 'HALLMARK_MYC_TARGETS_V2', 'HALLMARK_CELL_CYCLE',
    'HALLMARK_DNA_REPAIR', 'HALLMARK_P53_PATHWAY', 'HALLMARK_APOPTOSIS',
    'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
    'HALLMARK_HYPOXIA', 'HALLMARK_ANDROGEN_RESPONSE',
    'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
]

_hallmark_ctrl_clean = {
    'NT vs|G731K': _hallmark_ctrl_gsea['NT vs G731K — no drug'],
    'NT vs|G822E': _hallmark_ctrl_gsea['NT vs G822E — no drug'],
    'NT vs|I733V': _hallmark_ctrl_gsea['NT vs I733V — no drug'],
}

plot_gsea_heatmap_T(
    _hallmark_ctrl_clean,
    sig_order=_HALLMARK_ORDER,
    outpath=os.path.join(FIGDIR_05, 'hallmark_gsea_NT_vs_variants_Control'),
    title='Hallmark GSEA — NT vs CDK12 Variants (no drug)\n'
          '(negative NES = enriched in variant)',
    col_label_sep='|',
)


print("\nRunning Reactome GSEA...")
try:
    REACTOME = gs.get_library('Reactome_2022')
    print(f"  Loaded Reactome_2022: {len(REACTOME)} gene sets")
    reactome_dotplot('NT_vs_G822E_Control', 'NT vs G822E — no drug', f'{FIGDIR_06}/reactome_NT_vs_G822E', REACTOME)
    reactome_dotplot('NT_vs_G731K_Control', 'NT vs G731K — no drug', f'{FIGDIR_06}/reactome_NT_vs_G731K', REACTOME)
    reactome_dotplot('NT_vs_I733V_Control', 'NT vs I733V — no drug', f'{FIGDIR_06}/reactome_NT_vs_I733V', REACTOME)
except Exception as e:
    print(f"  Reactome download failed ({e}) — skipping Reactome GSEA")

print("\nDone.")
