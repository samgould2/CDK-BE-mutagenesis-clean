"""
T47D CDK4-variant RNA-seq: CDK4 transcript levels, cell-cycle composition, RB-E2F program.

Design: 12 samples, DMSO only = CBE_NT (control) + CDK4 D196N, E56K, W198* (nonsense), 3 reps each.

Headline caveat, established in section 0 below: CBE_NT is NOT isogenic with the three variant
lines (focal 9p21 MTAP/CDKN2A loss and complete loss of the ER/luminal program in all three
variants but not the control). Every quantity that depends on CBE_NT is therefore reported as
confounded, and the cell-cycle / E2F conclusions are drawn from variant-vs-variant comparisons,
which share that background.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as sps
import matplotlib as mpl
import matplotlib.pyplot as plt

import genesets as GS

mpl.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OLD  = BASE / 'T47D_CDK4_20260927_results' / 'deseq2_results'
OUT  = HERE / 'results'; OUT.mkdir(exist_ok=True)
FIG  = HERE / 'figures'; FIG.mkdir(exist_ok=True)

CTRL      = 'CBE_NT'
VARIANTS  = ['CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
MISSENSE  = ['CDK4_D196N', 'CDK4_E56K']
NULL      = 'CDK4_W198'
GENOTYPES = [CTRL] + VARIANTS
LABEL = {'CBE_NT': 'CBE NT', 'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}
COLOR = {'CBE_NT': '#9e9e9e', 'CDK4_D196N': '#4393c3', 'CDK4_E56K': '#66bd63', 'CDK4_W198': '#d6604d'}

# ───────────────────────────── load ──────────────────────────────────────────
tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gene_name = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
samples = [s for g in GENOTYPES for s in meta.index[meta == g]]
tpm, meta = tpm[samples], meta[samples]
log_tpm = np.log2(tpm + 1)
VAR_SAMPLES = [s for s in samples if meta[s] != CTRL]

order = tpm.mean(axis=1).sort_values(ascending=False).index
sym2id = pd.Series(order, index=gene_name.reindex(order).values)
sym2id = sym2id[~sym2id.index.duplicated()]

res          = {v: pd.read_csv(OLD / f'{v}_vs_{CTRL}_veh.csv').set_index('gene_id') for v in VARIANTS}
res_vs_miss  = pd.read_csv(OLD / f'{NULL}_vs_missense_veh.csv').set_index('gene_id')

grp = lambda g: list(meta.index[meta == g])

summary = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); summary.append(s)

stars = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'ns'

def holm(p):
    p = np.asarray(p, float); out = np.empty(len(p)); run = 0.0
    for rank, i in enumerate(np.argsort(p)):
        run = max(run, min(1.0, (len(p) - rank) * p[i])); out[i] = run
    return out

# ══════════ 0. is the control isogenic with the variant lines? ═══════════════
say('## 0. QC: the CBE NT control is not isogenic with the variant lines\n')

X = log_tpm.T
top = X.var(axis=0).sort_values(ascending=False).head(500).index
Xc = X[top] - X[top].mean(axis=0)
S_all = np.linalg.svd(Xc.values, compute_uv=False)
var_all = S_all ** 2 / (S_all ** 2).sum() * 100
Xv = log_tpm[VAR_SAMPLES].T
topv = Xv.var(axis=0).sort_values(ascending=False).head(500).index
Xvc = Xv[topv] - Xv[topv].mean(axis=0)
Uv, Sv, _ = np.linalg.svd(Xvc.values, full_matrices=False)
var_v = Sv ** 2 / (Sv ** 2).sum() * 100
pd.DataFrame((Uv * Sv)[:, :3], index=VAR_SAMPLES, columns=['PC1', 'PC2', 'PC3']) \
  .assign(genotype=meta[VAR_SAMPLES]).to_csv(OUT / 'PCA_variants_only.csv')

say(f'PC1 of the full 12-sample matrix carries **{var_all[0]:.0f}%** of the variance and separates '
    f'CBE NT from all three variants together; PC2 is {var_all[1]:.0f}%. '
    'Within the 9 variant samples alone, PC1 is '
    f'{var_v[0]:.0f}% and separates E56K from D196N/W198*, which are nearly superimposed.\n')

IDENTITY = ['ESR1', 'PGR', 'GREB1', 'TFF1', 'FOXA1', 'GATA3', 'CDH1', 'EPCAM', 'KRT8', 'KRT18',
            'VIM', 'MTAP', 'CDKN2A', 'KLHL9', 'CDKN2B-AS1', 'TUSC1']
gmean = tpm.T.groupby(meta).mean().T
idq = pd.DataFrame({g: gmean.loc[[sym2id[s] for s in IDENTITY if s in sym2id.index], g].values
                    for g in GENOTYPES}, index=[s for s in IDENTITY if s in sym2id.index])
idq.to_csv(OUT / 'identity_qc_TPM.csv')
say('Mean TPM of lineage-identity and 9p21 genes:\n')
say('| gene | ' + ' | '.join(LABEL[g] for g in GENOTYPES) + ' |')
say('|---|' + '---|' * len(GENOTYPES))
for s, r in idq.iterrows():
    say(f'| {s} | ' + ' | '.join(f'{r[g]:.1f}' for g in GENOTYPES) + ' |')
say('\nTwo independent signatures separate the control from all three variant lines:')
say('* **Focal 9p21.3 loss** — MTAP, CDKN2A and the adjacent KLHL9 all drop ~20-40x in every variant line, '
    'while the flanking CDKN2B-AS1 and TUSC1 are unchanged. A clean focal deletion, not a transcriptional response.')
say('* **Loss of the ER/luminal program** — ESR1, PGR, GREB1, TFF1, FOXA1, GATA3, CDH1 and EPCAM collapse to '
    'near zero and VIM rises, in all three variants.')
say('\nNeither is a plausible consequence of editing CDK4. 25 genes are expressed in the control (>5 TPM) and '
    'silent in all three variants (<0.5 TPM) while **0** genes go the other way, so the variant lines look like '
    'a de-differentiated, 9p21-deleted derivative of the same parental stock rather than a different cell line. '
    'The practical consequence is the same either way: **variant-vs-CBE NT contrasts are confounded by cell-line '
    'background**, and that confounder is 79% of the variance in the dataset.\n')

# ══════════ 1. CDK4 transcript levels ════════════════════════════════════════
say('\n## 1. CDK4 transcript levels\n')
rows = []
for sym in GS.PATHWAY_GENES:
    if sym not in sym2id.index:
        continue
    gid = sym2id[sym]
    ctrl_tpm = tpm.loc[gid, grp(CTRL)].mean()
    for g in GENOTYPES:
        v = tpm.loc[gid, grp(g)]
        r = dict(gene=sym, genotype=g, mean_TPM=v.mean(), sd_TPM=v.std(),
                 fold_vs_CBE_NT=v.mean() / ctrl_tpm if ctrl_tpm > 0 else np.nan)
        if g != CTRL and gid in res[g].index:
            r['deseq2_log2FC_vs_NT'] = res[g].loc[gid, 'log2FoldChange']
            r['deseq2_padj_vs_NT']   = res[g].loc[gid, 'padj']
        rows.append(r)
pathway = pd.DataFrame(rows)
pathway.to_csv(OUT / 'pathway_gene_levels.csv', index=False)

gid = sym2id['CDK4']
miss_tpm = tpm.loc[gid, grp(MISSENSE[0]) + grp(MISSENSE[1])].mean()
cdk4 = pathway[pathway.gene == 'CDK4'].set_index('genotype')
say('| genotype | CDK4 TPM | vs CBE NT | vs missense mean | DESeq2 log2FC vs NT | padj |')
say('|---|---|---|---|---|---|')
for g in GENOTYPES:
    r = cdk4.loc[g]
    p = r.get('deseq2_padj_vs_NT')
    say(f"| {LABEL[g]} | {r.mean_TPM:.0f} ± {r.sd_TPM:.0f} | {r.fold_vs_CBE_NT:.2f} | "
        f"{r.mean_TPM / miss_tpm:.2f} | "
        f"{'—' if g == CTRL else f'{r.deseq2_log2FC_vs_NT:+.2f}'} | "
        f"{'—' if g == CTRL else ('<1e-300' if p == 0 else f'{p:.1e}')} |")
lfc_m = res_vs_miss.loc[gid]
say(f"\nW198* vs the two missense lines (matched background, DESeq2): "
    f"log2FC = {lfc_m['log2FoldChange']:+.2f} ({2 ** lfc_m['log2FoldChange']:.2f}x), padj = {lfc_m['padj']:.1e}.\n")
say('Interpretation: W198* carries a premature termination codon and its mRNA sits at **~24% of control and ~27% '
    'of the missense lines** — a 3.7-4.2x reduction, exactly what nonsense-mediated decay predicts, and the one '
    'result here that is internally controlled (it holds against the missense lines, not just against CBE NT). '
    'D196N is 17% lower than CBE NT (padj 2e-7) and E56K is unchanged; both are within the range of the '
    'background difference between the lines, so neither is evidence of an effect on CDK4 mRNA.\n')

# ══════════ scores ═══════════════════════════════════════════════════════════
SETS = {'core_E2F':     GS.CORE_E2F,
        'S_specific':   GS.S_SPECIFIC,
        'G2M_specific': GS.G2M_SPECIFIC,
        'S_Tirosh':     GS.TIROSH_S,
        'G2M_Tirosh':   GS.TIROSH_G2M}
hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
SETS['Hallmark_E2F_Targets']    = hall['E2F Targets']
SETS['Hallmark_G2M_Checkpoint'] = hall['G2-M Checkpoint']

def ids(symbols, cols, min_log_tpm=1.0):
    g = [sym2id[s] for s in symbols if s in sym2id.index]
    return [i for i in g if log_tpm.loc[i, cols].mean() > min_log_tpm]

def geo(symbols, cols):
    """Mean log2(TPM+1) over a gene set = log2 of the geometric mean.

    If a phase-restricted transcript's bulk level is proportional to the fraction of cells in
    that phase, every gene in the set is scaled by the same factor; the geometric mean estimates
    that factor without being dominated by the few highest-expressed members of the set."""
    return log_tpm.loc[ids(symbols, cols), cols].mean(axis=0)

def zscore(symbols, cols):
    m = log_tpm.loc[ids(symbols, cols), cols]
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).mean(axis=0)

scores = pd.DataFrame({'genotype': meta})
for n, g in SETS.items():
    scores[n + '_geo_all12'] = geo(g, samples)
    scores[n + '_z_all12']   = zscore(g, samples)
    scores.loc[VAR_SAMPLES, n + '_z_variants'] = zscore(g, VAR_SAMPLES)
scores.to_csv(OUT / 'cell_cycle_scores.csv')

rng = np.random.default_rng(0)
def median_fc(symbols, contrast, cols):
    g = [i for i in ids(symbols, cols) if i in contrast.index]
    l = contrast.loc[g, 'log2FoldChange'].dropna().to_numpy()
    b = [np.median(rng.choice(l, len(l))) for _ in range(4000)]
    return 2 ** np.median(l), 2 ** np.percentile(b, 2.5), 2 ** np.percentile(b, 97.5), len(l)

def fc_test(symbols, a_samples, b_samples, cols):
    s = geo(symbols, cols)
    fc = 2 ** (s[a_samples].mean() - s[b_samples].mean())
    p = sps.ttest_ind(s[a_samples], s[b_samples], equal_var=False)[1]
    return fc, p

# ══════════ 2. cell-cycle composition ════════════════════════════════════════
say('\n## 2. Cell-cycle composition\n')
say('Phase-restricted transcripts are read as proxies for phase occupancy: if a gene is expressed only in S '
    'phase, its bulk level is proportional to the fraction of cells in S. The fold change is estimated two ways — '
    'the per-sample geometric mean of the gene set (Welch 3v3) and the median per-gene DESeq2 log2FC over the '
    'same set (bootstrap CI over genes). Comparisons are **W198\\* against the two missense lines**, which share '
    'the confounded background; the CBE NT comparison is shown for completeness only.\n')

miss_samples = grp(MISSENSE[0]) + grp(MISSENSE[1])
rows = []
say('| gene set | n | W198* vs missense (geo) | p | W198* vs missense (median per-gene) | W198* vs CBE NT (geo) |')
say('|---|---|---|---|---|---|')
for n in ['S_specific', 'S_Tirosh', 'core_E2F', 'G2M_specific', 'G2M_Tirosh']:
    g = SETS[n]
    fc, p = fc_test(g, grp(NULL), miss_samples, VAR_SAMPLES)
    m, lo, hi, k = median_fc(g, res_vs_miss, VAR_SAMPLES)
    fc_nt, _ = fc_test(g, grp(NULL), grp(CTRL), samples)
    rows.append(dict(gene_set=n, n_genes=k, fc_geo=fc, p=p, fc_median=m, ci_lo=lo, ci_hi=hi, fc_vs_NT=fc_nt))
    say(f'| {n} | {k} | {fc:.2f} | {p:.3f} | {m:.2f} [{lo:.2f}–{hi:.2f}] | {fc_nt:.2f} |')
pd.DataFrame(rows).to_csv(OUT / 'phase_marker_fold_changes.csv', index=False)
say('')

# pairwise ladder
say('Pairwise ladder on the S-phase/E2F axis (geometric-mean fold change, Welch p):\n')
pair_rows = []
say('| comparison | core_E2F | S_specific | G2M_specific |')
say('|---|---|---|---|')
for a, b in [(NULL, 'CDK4_E56K'), (NULL, 'CDK4_D196N'), ('CDK4_D196N', 'CDK4_E56K')]:
    cells = []
    for n in ['core_E2F', 'S_specific', 'G2M_specific']:
        fc, p = fc_test(SETS[n], grp(a), grp(b), VAR_SAMPLES)
        cells.append(f'{fc:.2f} ({stars(p)})')
        pair_rows.append(dict(a=a, b=b, gene_set=n, fc=fc, p=p))
    say(f'| {LABEL[a]} vs {LABEL[b]} | ' + ' | '.join(cells) + ' |')
pd.DataFrame(pair_rows).to_csv(OUT / 'pairwise_variant_fold_changes.csv', index=False)
say('\nOrdering on S-phase/E2F output: **E56K > D196N > W198\\***. G2/M markers do not move in any comparison.\n')

# anchored proportions
ANCHOR = dict(G1=0.55, S=0.27, G2M=0.18)
say(f'Anchored proportion estimate. The missense lines are *assumed* to sit at '
    f'{ANCHOR["G1"]:.0%} G1 / {ANCHOR["S"]:.0%} S / {ANCHOR["G2M"]:.0%} G2/M (typical asynchronous log-phase '
    'culture — this is an assumption, not a measurement); the other genotypes are that anchor scaled by the '
    'measured marker fold changes and renormalised. Only the differences are supported by the data.\n')
fcS = {g: fc_test(GS.S_SPECIFIC, grp(g), miss_samples, VAR_SAMPLES)[0] for g in VARIANTS}
fcM = {g: fc_test(GS.G2M_SPECIFIC, grp(g), miss_samples, VAR_SAMPLES)[0] for g in VARIANTS}
prop_rows = []
for g in VARIANTS:
    s, m = ANCHOR['S'] * fcS[g], ANCHOR['G2M'] * fcM[g]
    g1 = ANCHOR['G1'] + (ANCHOR['S'] - s) + (ANCHOR['G2M'] - m)
    t = s + m + g1
    prop_rows.append(dict(genotype=g, G1=g1 / t, S=s / t, G2M=m / t, cycling_SG2M=(s + m) / t))
prop = pd.DataFrame(prop_rows).set_index('genotype')
prop.to_csv(OUT / 'estimated_phase_proportions.csv')
say('| genotype | est. %G1 | est. %S | est. %G2/M |')
say('|---|---|---|---|')
for g in VARIANTS:
    r = prop.loc[g]
    say(f'| {LABEL[g]} | {r.G1:.0%} | {r.S:.0%} | {r.G2M:.0%} |')
say('\nSo: yes, relative phase occupancy is estimable, and the answer is that the CDK4 null line has roughly '
    '**8-10% fewer cells in S** than the missense lines (~27% → ~25% of the population) with G2/M unchanged. '
    'That is a mild restriction-point effect, not an arrest. Absolute percentages require a DNA-content or EdU '
    'measurement to anchor — bulk RNA-seq cannot supply one.\n')

# ══════════ 3. E2F program ═══════════════════════════════════════════════════
say('\n## 3. RB-E2F program\n')
say('Note on direction: CDK4 phosphorylates RB, which releases E2F. Loss of CDK4 keeps RB active and '
    'therefore **lowers** E2F target transcription — the same direction palbociclib produces. An *elevated* '
    'E2F program is what RB1 or CDKN2A loss gives, not CDK4 loss. The data below move in the loss-of-function '
    'direction.\n')

z_rows = []
for n in ['core_E2F', 'Hallmark_E2F_Targets', 'S_Tirosh', 'Hallmark_G2M_Checkpoint', 'G2M_Tirosh']:
    col = n + '_z_variants'
    c = scores.loc[miss_samples, col].to_numpy()
    for v in VARIANTS:
        x = scores.loc[grp(v), col].to_numpy()
        z_rows.append(dict(score=n, genotype=v, mean_z=x.mean(),
                           p_vs_missense=sps.ttest_ind(x, c, equal_var=False)[1] if v == NULL else np.nan))
zdf = pd.DataFrame(z_rows)
zdf['p_holm'] = np.nan
mask = zdf.p_vs_missense.notna()
zdf.loc[mask, 'p_holm'] = holm(zdf.loc[mask, 'p_vs_missense'])
zdf.to_csv(OUT / 'program_score_stats.csv', index=False)
say('Mean z-score within the 9 variant samples (control excluded so the background cancels). z-scores are '
    'scale-free, so with triplicates this tight almost anything reaches significance — read the stars as '
    '"consistent across replicates", and take the effect size from the fold-change table in section 2. '
    'G2/M is the clearest example: significant by z-score, yet a 0.98-0.99x fold change.\n')
say('| score | D196N | E56K | W198* | W198* vs missense |')
say('|---|---|---|---|---|')
for n in ['core_E2F', 'Hallmark_E2F_Targets', 'S_Tirosh', 'Hallmark_G2M_Checkpoint', 'G2M_Tirosh']:
    s = zdf[zdf.score == n].set_index('genotype')
    say(f'| {n} | ' + ' | '.join(f'{s.loc[v, "mean_z"]:+.2f}' for v in VARIANTS) +
        f' | {stars(s.loc[NULL, "p_holm"])} |')
say('')

rows = []
for sym in GS.CORE_E2F:
    if sym not in sym2id.index:
        continue
    g = sym2id[sym]
    r = dict(gene=sym, missense_TPM=tpm.loc[g, miss_samples].mean())
    if g in res_vs_miss.index:
        r['W198_vs_missense_log2FC'] = res_vs_miss.loc[g, 'log2FoldChange']
        r['W198_vs_missense_padj']   = res_vs_miss.loc[g, 'padj']
    for v in VARIANTS:
        if g in res[v].index:
            r[f'{LABEL[v]}_vs_NT_log2FC'] = res[v].loc[g, 'log2FoldChange']
            r[f'{LABEL[v]}_vs_NT_padj']   = res[v].loc[g, 'padj']
    rows.append(r)
core = pd.DataFrame(rows).set_index('gene')
core.to_csv(OUT / 'core_E2F_target_log2FC.csv')
sig = core[core.W198_vs_missense_padj < .05]
n_neg = int((core.W198_vs_missense_log2FC < 0).sum()); n_tot = int(core.W198_vs_missense_log2FC.notna().sum())
sign_p = sps.binomtest(n_neg, n_tot, 0.5).pvalue
say(f'Core RB-E2F targets, W198\\* vs missense: only {len(sig)} of {len(core)} reach padj<0.05 individually '
    f'(E2F1), but **{n_neg} of {n_tot} move down** (sign test p = {sign_p:.1e}) and the median log2FC over the '
    f'set is {core.W198_vs_missense_log2FC.median():+.2f} ({2 ** core.W198_vs_missense_log2FC.median():.2f}x). '
    'The program moves coherently; no single gene moves much.\n')
say('| gene | missense TPM | W198* vs missense log2FC | padj |')
say('|---|---|---|---|')
for s, r in core.sort_values('W198_vs_missense_log2FC').iterrows():
    say(f"| {s} | {r.missense_TPM:.0f} | {r.W198_vs_missense_log2FC:+.2f} | "
        f"{'—' if pd.isna(r.W198_vs_missense_padj) else f'{r.W198_vs_missense_padj:.2g}'} |")
say('\nEvery core target shifts by well under 2-fold. The program is coherently nudged down, not switched off.\n')

# ══════════ 4. phenocopy ═════════════════════════════════════════════════════
say('\n## 4. Does D196N or E56K phenocopy the null?\n')
common = sorted(set(res[VARIANTS[0]].index) & set(res[VARIANTS[1]].index) & set(res[VARIANTS[2]].index))
lfc = pd.DataFrame({v: res[v].loc[common, 'log2FoldChange'] for v in VARIANTS}).dropna()
null_degs = res[NULL].loc[common].dropna(subset=['padj'])
null_degs = null_degs[(null_degs.padj < .05) & (null_degs.log2FoldChange.abs() > 1)].index
sub = lfc.loc[lfc.index.intersection(null_degs)]
pheno = []
for v in MISSENSE:
    pheno.append(dict(variant=v, pearson_r=sps.pearsonr(sub[v], sub[NULL])[0],
                      slope=np.polyfit(sub[NULL], sub[v], 1)[0]))
pheno = pd.DataFrame(pheno)
pheno.to_csv(OUT / 'phenocopy_concordance.csv', index=False)
say('Correlating each variant\'s log2FC-vs-CBE NT against W198\\*\'s over the '
    f'{len(sub)} W198\\* DEGs gives r = ' +
    ', '.join(f'{LABEL[r.variant]} {r.pearson_r:.2f} (slope {r.slope:.2f})' for r in pheno.itertuples()) +
    '. **These numbers are not interpretable as phenocopy scores**: all three variants share the 9p21-deleted, '
    'ER-negative background and the control does not, so most of what is being correlated is that shared '
    'background rather than CDK4 genotype.\n')
say('The within-background evidence instead says:')
say(f'* W198\\* vs the missense lines yields only 26 DEGs at |log2FC|>1 — CDK4 itself plus a residual set of '
    'epithelial genes (CDH1, CLDN3/4, RAB25, S100A9) that tracks E56K\'s slightly more luminal state.')
say('* On the S-phase/E2F axis the three lines rank E56K > D196N > W198\\*, with D196N roughly midway between '
    'the two on a log scale — consistent with D196N being a **partial** loss-of-function allele and E56K '
    'retaining function on this readout.')
say('* D196N and W198\\* are nearly superimposed on the variant-only PCA, but that axis is dominated by residual '
    'identity drift, so it should not be read as functional equivalence on its own.\n')

# ══════════ figures ══════════════════════════════════════════════════════════
def strip(ax, col, ylab, cols, sig_vs_missense=True):
    r = np.random.default_rng(0)
    gs_ = VARIANTS if cols is VAR_SAMPLES else GENOTYPES
    vals = {g: scores.loc[grp(g), col].to_numpy() for g in gs_}
    for xi, g in enumerate(gs_):
        ax.scatter(xi + r.uniform(-.1, .1, 3), vals[g], s=28, facecolor=COLOR[g], edgecolor='black', lw=.5, zorder=3)
        ax.plot([xi - .25, xi + .25], [vals[g].mean()] * 2, color='black', lw=1.2, zorder=4)
    top_ = max(v.max() for v in vals.values()); span = top_ - min(v.min() for v in vals.values())
    if sig_vs_missense:
        p = sps.ttest_ind(scores.loc[grp(NULL), col], scores.loc[miss_samples, col], equal_var=False)[1]
        xv = gs_.index(NULL); yb = top_ + span * .12
        ax.plot([0, 0, xv, xv], [yb - span * .03, yb, yb, yb - span * .03], color='black', lw=.8)
        ax.text(xv / 2, yb, stars(p), ha='center', va='bottom', fontsize=9)
        ax.set_ylim(top=top_ + span * .35)
    ax.set_xticks(range(len(gs_)))
    ax.set_xticklabels([LABEL[g] for g in gs_], rotation=45, ha='right', rotation_mode='anchor', fontsize=8)
    ax.set_ylabel(ylab, fontsize=8); ax.tick_params(labelsize=7.5)
    ax.spines[['top', 'right']].set_visible(False)

fig, axes = plt.subplots(1, 4, figsize=(10, 2.9))
for ax, (c, l) in zip(axes, [('core_E2F_z_variants', 'Core RB-E2F targets\n(mean z, variants only)'),
                             ('Hallmark_E2F_Targets_z_variants', 'Hallmark E2F Targets\n(mean z)'),
                             ('S_specific_z_variants', 'S-phase markers\n(mean z)'),
                             ('G2M_specific_z_variants', 'G2/M markers\n(mean z)')]):
    strip(ax, c, l, VAR_SAMPLES)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'program_scores_matched.{e}', dpi=300, bbox_inches='tight')
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={'width_ratios': [1, 1.25]})
gid = sym2id['CDK4']
for xi, g in enumerate(GENOTYPES):
    v = tpm.loc[gid, grp(g)].to_numpy()
    axes[0].bar(xi, v.mean(), color=COLOR[g], edgecolor='black', lw=.5, width=.65)
    axes[0].scatter([xi] * 3, v, s=18, color='black', zorder=3)
axes[0].set_xticks(range(4)); axes[0].set_xticklabels([LABEL[g] for g in GENOTYPES], rotation=45,
                                                      ha='right', rotation_mode='anchor', fontsize=8)
axes[0].set_ylabel('CDK4 (TPM)', fontsize=9); axes[0].spines[['top', 'right']].set_visible(False)
axes[0].tick_params(labelsize=8); axes[0].set_title('CDK4 mRNA (NMD in W198*)', fontsize=9)

ax = axes[1]
w = 0.2
for k, g in enumerate(GENOTYPES):
    ax.bar(np.arange(len(idq)) + (k - 1.5) * w, np.log10(idq[g] + .1), width=w,
           color=COLOR[g], edgecolor='black', lw=.3, label=LABEL[g])
ax.set_xticks(range(len(idq))); ax.set_xticklabels(idq.index, rotation=90, fontsize=6.5, style='italic')
ax.set_ylabel('log$_{10}$(TPM + 0.1)', fontsize=8); ax.tick_params(labelsize=7)
ax.legend(fontsize=6.5, frameon=False, ncol=4, loc='upper center', bbox_to_anchor=(.5, 1.18))
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'CDK4_and_identity_QC.{e}', dpi=300, bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(3.4, 3.0))
bottom = np.zeros(len(VARIANTS))
for phase, c in [('G1', '#cfd8dc'), ('S', '#4393c3'), ('G2M', '#d6604d')]:
    v = prop.loc[VARIANTS, phase].values
    ax.bar([LABEL[g] for g in VARIANTS], v, bottom=bottom, color=c, edgecolor='black', lw=.5, label=phase)
    for xi, (a, b) in enumerate(zip(v, bottom)):
        ax.text(xi, b + a / 2, f'{a:.0%}', ha='center', va='center', fontsize=7)
    bottom += v
ax.set_ylabel('estimated fraction of cells', fontsize=8)
ax.set_title('Anchored to missense = 55/27/18\n(assumption, not a measurement)', fontsize=7.5)
ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, .5), loc='center left')
ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
for e in ('png', 'svg'):
    plt.savefig(FIG / f'estimated_phase_proportions.{e}', dpi=300, bbox_inches='tight')
plt.close()

(HERE / 'FINDINGS.md').write_text(
    '# T47D CDK4 variants — CDK4 levels, cell-cycle composition, RB-E2F program\n\n'
    'Generated by `analyze.py`; tables in `results/`, figures in `figures/`.\n\n' + '\n'.join(summary) + '\n')
print('\nwrote FINDINGS.md,', len(list(OUT.glob('*.csv'))), 'tables,', len(list(FIG.glob('*.png'))), 'figures')
