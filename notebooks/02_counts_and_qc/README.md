# 02 — counts and QC

`count_tables.ipynb` aggregates the per-sample guide counts in
`data/raw/screen_counts/` into the barcode, protospacer and matched count
tables under `data/processed/count_tables/`, and the reads-per-million tables
beside them.

Screen QC panels — cumulative population doublings, library representation,
LFC correlation between arms — are drawn from those tables by
`cdkbe.source_data.cpd` and `cdkbe.lfc_correlation`, and appear in
`09_supplementary/`.
