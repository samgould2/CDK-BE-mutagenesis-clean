"""Controls for the gene-length elongation analysis.

Gene length is correlated with expression level and with intron number, so a raw
length-vs-log2FC trend could in principle be driven by either. This script quantifies
how much of the length effect survives those controls.

1. Expression stratification — length tertiles recomputed INSIDE each expression tertile,
   so the comparison is between genes of similar abundance. Note the global
   top-decile-vs-median "gap" used in the figures is NOT valid within a stratum: length
   and expression are anticorrelated, so within the low-expression stratum the median gene
   is already long and the gap collapses for the wrong reason. Within-stratum
   top-vs-bottom tertile is used here instead.
2. Multiple regression of log2FC on log10(span), log10(introns+1) and log10(baseMean),
   to see how the length coefficient changes as covariates are added.

Outputs: results/elongation_controls_stratified.csv
         results/elongation_controls_regression.csv

Run:  cd reanalysis/scripts && python 08_elongation_controls.py
"""
import os
import numpy as np
import pandas as pd
from _style import RES, GENOS, DRUGS, drug_con, expressed_pc

pc = expressed_pc()
ref = drug_con('NT', 'HQ461').reindex(pc.index)
expr_tertile = pd.qcut(np.log10(ref.baseMean), 3, labels=['low', 'mid', 'high'])

print(f'{len(pc):,} genes; corr log10(span) vs log10(baseMean) = '
      f'{np.corrcoef(np.log10(pc.span_bp), np.log10(ref.baseMean))[0, 1]:.3f}')
print(f'corr log10(span) vs log10(introns+1) = '
      f'{np.corrcoef(np.log10(pc.span_bp), np.log10(pc.n_introns + 1))[0, 1]:.3f}\n')


def tertile_effect(lfc, genes):
    """median log2FC(longest third) - median(shortest third), tertiles over `genes`."""
    t = pd.qcut(genes.span_bp.rank(method='first'), 3, labels=False)
    return np.median(lfc[t == 2].dropna()) - np.median(lfc[t == 0].dropna())


rows = []
for d in DRUGS:
    for g in GENOS:
        lfc = drug_con(g, d).log2FC_shrunk.reindex(pc.index)
        rec = {'drug': d, 'genotype': g, 'global': tertile_effect(lfc, pc)}
        for t in ['low', 'mid', 'high']:
            m = (expr_tertile == t).values
            rec[f'expr_{t}'] = tertile_effect(lfc[m], pc[m])
        rows.append(rec)
strat = pd.DataFrame(rows)
strat.to_csv(os.path.join(RES, 'elongation_controls_stratified.csv'), index=False)
print('Long-minus-short tertile log2FC, globally and within expression strata:')
print(strat.round(3).to_string(index=False))

X = pd.DataFrame({'logspan': np.log10(pc.span_bp),
                  'logintron': np.log10(pc.n_introns + 1),
                  'logexpr': np.log10(ref.baseMean)})
reg_rows = []
for d in DRUGS:
    for g in GENOS:
        y = drug_con(g, d).log2FC_shrunk.reindex(pc.index)
        for model in (['logspan'], ['logspan', 'logexpr'],
                      ['logspan', 'logintron'], ['logspan', 'logintron', 'logexpr']):
            df = X[model].join(y.rename('y')).dropna()
            A = np.column_stack([np.ones(len(df))] + [df[c].values for c in model])
            coef = np.linalg.lstsq(A, df.y.values, rcond=None)[0]
            reg_rows.append({'drug': d, 'genotype': g, 'model': ' + '.join(model), 'n': len(df),
                             **{f'beta_{c}': v for c, v in zip(model, coef[1:])}})
reg = pd.DataFrame(reg_rows)
reg.to_csv(os.path.join(RES, 'elongation_controls_regression.csv'), index=False)
print('\nlog2FC regression coefficients (NT only shown):')
print(reg[reg.genotype == 'NT'].round(3).to_string(index=False))
