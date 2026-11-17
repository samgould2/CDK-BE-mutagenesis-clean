# Reproducing the analysis

## Setup

```bash
pip install -e .
python3 -m pytest tests/ -q
```

The tests read real data, so a pass also confirms `data/` is intact.

## Order

Notebooks are numbered in dependency order. Within a directory they are
independent of each other unless a name says otherwise.

| directory | reads | writes |
|---|---|---|
| `01_library_design` | `data/reference/library_generation` | the guide library |
| `02_counts_and_qc` | `data/raw/screen_counts` | `data/processed/count_tables` |
| `03_editing` | `data/raw/screen_editing` | `data/processed/editing` |
| `04_lfc_fdr` | count tables + editing | `data/processed/lfc_fdr_tables` |
| `05_rnaseq` | `data/processed/rnaseq` | DESeq2 results, GSEA |
| `06_structure_and_conservation` | `data/reference/KLIFS`, `conservation` | per-position tables |
| `07_other_analyses` | various | splicing, off-target, clinical |
| `08_figures` | processed tables only | `figures/`, source-data workbooks |

`08_figures` never recomputes anything. If a number in a figure looks wrong,
the fix belongs upstream, not in the plotting notebook.

## Rerunning one figure

Each notebook in `08_figures` covers one manuscript figure and writes both the
panels and that figure's source-data sheet. Nothing in it reaches outside
`data/processed/` and `data/reference/`, so it runs without the raw tiers
present.

To write the panels without opening a notebook, and under their manuscript
names:

```bash
python scripts/render_figures.py 5          # one figure
python scripts/render_figures.py 5j --png   # one panel, with a PNG beside it
python scripts/render_figures.py            # everything
```

## Rewriting the source-data workbooks

`source_data_tables/` is the same panels as numbers: one workbook per figure,
one sheet per panel, with an INDEX sheet saying what each sheet is.

```bash
python scripts/build_source_data.py          # all of them
python scripts/audit_source_data.py          # what each panel prints vs its sheet
```

The audit is a prompt to look, not a test: it reports every word a panel
prints that is not in its sheets, and most of those are axis titles rather
than data. Its own docstring says how to read the output.

## The ChimeraX panels

```bash
python3 scripts/chimerax/render.py --list    # what would render, and for which figure
python3 scripts/chimerax/render.py           # all of them
python3 scripts/chimerax/render.py CDK9_fig2b_labeled
```

Rendering needs the ChimeraX application, because saving an image needs an
OpenGL context. The `.cxc` command files are generated: change
`scripts/chimerax/generate_*.py`, not the files under `cxc/`.

## The Sanger index

```bash
python3 scripts/make_sanger_index.py         # INDEX.csv, the FASTA, the README, the supplement
python3 scripts/make_sanger_index.py --check # compare with what is on disk
```

## Rerunning from raw

```bash
# counts -> LFC/FDR
jupyter nbconvert --execute --inplace notebooks/02_counts_and_qc/*.ipynb
jupyter nbconvert --execute --inplace notebooks/04_lfc_fdr/*.ipynb

# raw editing -> processed editing  (slow: reads ~300 MB of gzipped per-read calls)
jupyter nbconvert --execute --inplace notebooks/03_editing/*.ipynb
```

Both write into `data/processed/`, overwriting what shipped. `data/MANIFEST.csv`
has the SHA-256 of every shipped file, so you can confirm a rerun reproduced
them; `python3 scripts/verify_data.py` does the comparison, `--refresh` records
a rebuild, and `--prune` drops the row of a file you removed on purpose.

## The figure notebooks

They ship executed, so every panel is visible without running anything. To
redo that after changing a drawer:

```bash
python3 scripts/run_notebooks.py --inplace   # the figure notebooks, output kept
python3 scripts/run_notebooks.py --all       # the processing notebooks too
```

Without `--inplace` the executed copies go to a temporary directory, which is
what you want when you only mean to check that they still run.

## A note on the original notebooks

This analysis was first written as a set of exploratory notebooks, and the
modules here were extracted from them. Those originals are **not** in this
repository: several cannot be run top to bottom in any case --
`FACS_dose_response_2024_8_22.ipynb` redefines `log_logistic_4_param` eleven
times with different normalisation, and at least one cell reads a `df` that no
preceding cell in file order produces.

They are still cited. A docstring that says *"from
`09_analysis_and_viz/analysis.ipynb` cell 137"* is naming the original
notebook and cell a function came from, so that a number in a figure can be
traced back to where it was first computed. The path is a citation, not a file
in this checkout.
