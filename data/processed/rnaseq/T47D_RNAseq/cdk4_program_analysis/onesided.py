"""
One-sided vs two-sided testing of the cell-cycle / E2F scores.

Pre-specified direction: loss of CDK4 reduces RB phosphorylation, so E2F and S-phase
scores are tested in the LOWER tail (variant < reference). Reported for both reference
schemes (CBE NT, and the missense lines) so the two can be compared directly.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats as sps

import genesets as GS

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OLD  = BASE / 'T47D_CDK4_20260927_results' / 'deseq2_results'
OUT  = HERE / 'results'

CTRL     = 'CBE_NT'
VARIANTS = ['CDK4_D196N', 'CDK4_E56K', 'CDK4_W198']
MISSENSE = ['CDK4_D196N', 'CDK4_E56K']
NULL     = 'CDK4_W198'
LABEL = {'CBE_NT': 'CBE NT', 'CDK4_D196N': 'D196N', 'CDK4_E56K': 'E56K', 'CDK4_W198': 'W198*'}

tpm = pd.read_csv(BASE / 'salmon.merged.gene_tpm.tsv', sep='\t', index_col=0)
gene_name = tpm.pop('gene_name')
meta = pd.Series(tpm.columns.str.extract(r'^T47D_(.+)_veh_\d+$')[0].values, index=tpm.columns)
log_tpm = np.log2(tpm + 1)
order = tpm.mean(axis=1).sort_values(ascending=False).index
sym2id = pd.Series(order, index=gene_name.reindex(order).values)
sym2id = sym2id[~sym2id.index.duplicated()]
grp = lambda g: list(meta.index[meta == g])
VAR_SAMPLES = [s for s in tpm.columns if meta[s] != CTRL]
MISS_SAMPLES = grp(MISSENSE[0]) + grp(MISSENSE[1])

res         = {v: pd.read_csv(OLD / f'{v}_vs_{CTRL}_veh.csv').set_index('gene_id') for v in VARIANTS}
res_vs_miss = pd.read_csv(OLD / f'{NULL}_vs_missense_veh.csv').set_index('gene_id')

hall = json.loads((BASE / 'T47D_CDK4_20260927_results' / 'hallmark_gene_sets_cache.json').read_text())
SETS = {'core_E2F': GS.CORE_E2F, 'Hallmark_E2F_Targets': hall['E2F Targets'],
        'S_specific': GS.S_SPECIFIC, 'S_Tirosh': GS.TIROSH_S,
        'G2M_specific': GS.G2M_SPECIFIC, 'G2M_Tirosh': GS.TIROSH_G2M,
        'Hallmark_G2M_Checkpoint': hall['G2-M Checkpoint']}

def holm(p):
    p = np.asarray(p, float); out = np.empty(len(p)); run = 0.0
    for rank, i in enumerate(np.argsort(p)):
        run = max(run, min(1.0, (len(p) - rank) * p[i])); out[i] = run
    return out

stars = lambda p: '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < .05 else 'ns'

def ids(symbols, cols):
    g = [sym2id[s] for s in symbols if s in sym2id.index]
    return [i for i in g if log_tpm.loc[i, cols].mean() > 1.0]

def zscore(symbols, cols):
    m = log_tpm.loc[ids(symbols, cols), cols]
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0).mean(axis=0)

def geo(symbols, cols):
    return log_tpm.loc[ids(symbols, cols), cols].mean(axis=0)

rng = np.random.default_rng(0)
def median_fc(symbols, contrast, cols):
    g = [i for i in ids(symbols, cols) if i in contrast.index]
    l = contrast.loc[g, 'log2FoldChange'].dropna().to_numpy()
    b = np.array([np.median(rng.choice(l, len(l))) for _ in range(10000)])
    # one-sided bootstrap p for "median log2FC < 0"
    return 2 ** np.median(l), 2 ** np.percentile(b, 2.5), 2 ** np.percentile(b, 97.5), (b >= 0).mean(), len(l)

rows = []
for ref_name, cols, ref_samples_fn in [('vs CBE NT', list(tpm.columns), lambda v: grp(CTRL)),
                                       ('vs missense', VAR_SAMPLES, lambda v: MISS_SAMPLES)]:
    for sname, genes in SETS.items():
        z = zscore(genes, cols); g = geo(genes, cols)
        variants = VARIANTS if ref_name == 'vs CBE NT' else [NULL]
        p2, p1 = [], []
        for v in variants:
            x, c = z[grp(v)].to_numpy(), z[ref_samples_fn(v)].to_numpy()
            t, pp = sps.ttest_ind(x, c, equal_var=False)
            p2.append(pp)
            p1.append(pp / 2 if t < 0 else 1 - pp / 2)      # lower-tail, pre-specified
            rows.append(dict(reference=ref_name, score=sname, variant=v,
                             fold_change=2 ** (g[grp(v)].mean() - g[ref_samples_fn(v)].mean()),
                             delta_z=x.mean() - c.mean(), p_two=pp, p_one=p1[-1]))
        n = len(variants)
        for r, a2, a1 in zip(rows[-n:], holm(p2), holm(p1)):
            r['p_two_holm'], r['p_one_holm'] = a2, a1
df = pd.DataFrame(rows)
df.to_csv(OUT / 'onesided_vs_twosided_stats.csv', index=False)

for ref in ['vs CBE NT', 'vs missense']:
    d = df[df.reference == ref]
    print(f'\n═══ {ref} ═══')
    print(f"{'score':24s} {'variant':8s} {'FC':>6s} {'p2':>8s} {'p2_holm':>9s} {'p1':>8s} {'p1_holm':>9s}   verdict")
    for r in d.itertuples():
        flip = '  <-- ns -> sig' if (r.p_two_holm >= .05 and r.p_one_holm < .05) else ''
        print(f'{r.score:24s} {LABEL[r.variant]:8s} {r.fold_change:6.2f} {r.p_two:8.4f} {r.p_two_holm:9.4f} '
              f'{r.p_one:8.4f} {r.p_one_holm:9.4f}   {stars(r.p_two_holm):>3s} -> {stars(r.p_one_holm):>3s}{flip}')

print('\n═══ gene-level estimator (median per-gene log2FC, bootstrap over genes) ═══')
print(f"{'reference':12s} {'score':24s} {'variant':8s} {'medFC':>6s} {'95% CI':>14s} {'1-sided p':>10s}  n")
for sname in ['core_E2F', 'S_specific', 'S_Tirosh', 'G2M_specific', 'G2M_Tirosh']:
    for v in VARIANTS:
        m, lo, hi, p, n = median_fc(SETS[sname], res[v], list(tpm.columns))
        print(f"{'vs CBE NT':12s} {sname:24s} {LABEL[v]:8s} {m:6.2f} {f'[{lo:.2f}-{hi:.2f}]':>14s} {p:10.4f}  {n}")
    m, lo, hi, p, n = median_fc(SETS[sname], res_vs_miss, VAR_SAMPLES)
    print(f"{'vs missense':12s} {sname:24s} {LABEL[NULL]:8s} {m:6.2f} {f'[{lo:.2f}-{hi:.2f}]':>14s} {p:10.4f}  {n}")
