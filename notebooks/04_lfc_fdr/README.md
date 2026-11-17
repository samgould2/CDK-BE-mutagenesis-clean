# 04 — log fold change and FDR

`lfc_fdr_tables.ipynb` turns count tables into the per-guide, per-condition
log fold changes and FDRs under `data/processed/lfc_fdr_tables/`, and the
master tables the lollipops and heatmaps read.

`cdkbe.lfc_fdr` holds the statistics and `docs/METHODS.md` §1.2 describes
them, including the two baseline artefacts that have to be avoided when
choosing a contrast.
