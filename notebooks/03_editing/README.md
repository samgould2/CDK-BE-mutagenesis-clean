# 03 — editing

`editing_tables.ipynb` turns the per-sample sensor-editing calls in
`data/raw/screen_editing/` into the per-gRNA tables under
`data/processed/editing/`: the HGVSp quantification, the single-amino-acid
variant calls, and their canonical variants, which keep only the edit the
editor itself makes.

Canonical versus non-canonical matters: several published panels were made
before canonical filtering existed. `docs/DISCREPANCIES.md` records which.
