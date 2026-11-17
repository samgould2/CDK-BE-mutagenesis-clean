# CDK base-editing mutagenesis screens

Code and data for the CDK base-editing mutagenesis manuscript: pooled
base-editor tiling screens across nine CDKs, the drug-resistance variants they
recover, and the structural, transcriptional and single-cell analyses built on
them.

The repository is **self-contained**. Clone it, `pip install -e .`, and every
panel this analysis produces can be rebuilt from the data in `data/` — no
external drives, no lab file shares, no paths outside the checkout. A test
asserts that mechanically, because it is the kind of property that breaks
silently.

| | |
|---|---|
| published panels | 162 — 57 main, 105 supplementary |
| redrawn from data here | **118** (44 main, 74 supplementary) |
| source-data tables | 119 panels, in 22 workbooks — one per figure |
| data shipped | ~1.7 GB, every file checksummed in `data/MANIFEST.csv` |
| tests | 769, run against the real data |

Every panel is listed in [`docs/COVERAGE.md`](docs/COVERAGE.md) with its
status. The ones not redrawn here are western blots, micrographs, Sanger
traces and schematics, and ten panels that are included as the images that
were published.

---

## Quick start

```bash
pip install -e .
```

```bash
python3 -m pytest tests/ -q
```

The tests read the real data, so a pass also confirms `data/` is intact.

```bash
python3 scripts/render_figures.py
```

That writes every panel the repository can draw into `figures/manuscript/`,
named for its place in the manuscript — `main_fig4a_CDK2.pdf`,
`supp_fig12c.pdf` — with `figures/manuscript/INDEX.csv` recording which
function drew each file.

```bash
python3 scripts/build_source_data.py
```

That writes the same panels as numbers: `source_data_tables/Figure_4.xlsx` and
its siblings, one sheet per panel, each workbook opening with an INDEX sheet
saying what every sheet is and where it came from.

```bash
python3 scripts/run_notebooks.py --inplace
```

That re-executes the figure notebooks and keeps their output, which is how the
rendered panels in `notebooks/08_figures/` and `notebooks/09_supplementary/`
got there. `tests/test_self_contained.py` checks that what they stored names no
path outside the repository and holds no error.

Python 3.9 or newer. `pip install -e .` is enough for every figure; the extras
in `pyproject.toml` cover the alignment panels (`msa`), re-fetching reference
data (`fetch`), re-deriving the single-cell statistics from objects that are
not shipped (`singlecell`), the read-level pipeline (`pipeline`) and the test
suite (`dev`).

---

## How reproduction works

Four tiers, each derived from the one above it, and each shipped so you can
start anywhere.

```
data/raw/            per-sample guide counts, sensor-editing calls, flow
                     exports, bench plates, Sanger traces
      |   notebooks/02_counts_and_qc, 03_editing
      v
data/processed/      aggregated counts, per-gRNA editing, LFC + FDR per guide
                     per condition, DESeq2 results, single-cell statistics
      |   src/cdkbe/source_data/   <- one module per group of panels
      v
figures/manuscript/  the panels, named for the manuscript
source_data_tables/  the same panels as numbers, one workbook per figure
```

`notebooks/` walks the same path in dependency order, and
[`docs/REPRODUCE.md`](docs/REPRODUCE.md) gives the run order and the commands
for rebuilding each stage.

Two rules the layout enforces:

* **`08_figures` never recomputes a screen statistic.** Figure notebooks read
  `data/processed/` and plot. If a number in a figure looks wrong, the fix
  belongs upstream.
* **Every panel is drawn by a function, not by a notebook.** The notebooks
  call `cdkbe.source_data`; so does `scripts/render_figures.py`, and so does
  the workbook builder. A panel therefore cannot differ between the notebook
  that shows it, the file that ships and the table behind it.

---

## What is in here

| path | what it is |
|---|---|
| `src/cdkbe/` | the package: 91 modules. `source_data/` has one module per group of panels, holding both the table and the drawing |
| `notebooks/` | 26 notebooks, numbered in dependency order. `08_figures/` and `09_supplementary/` ship executed, so every panel is visible on GitHub without running anything |
| `scripts/` | command-line entry points — render figures, build the workbooks, verify the data, drive ChimeraX |
| `data/` | the four tiers below, ~1.7 GB, with `MANIFEST.csv` recording all 3,044 files' sizes and SHA-256s, and each file's original source path where it has one |
| `source_data_tables/` | 22 workbooks, one per figure, one sheet per panel, plus the Sanger supplement |
| `tests/` | 769 checks against the real data, including byte-level comparisons with published artwork |
| `pipeline/` | the read-level pipeline upstream of `data/raw/` — sensor extraction, CRISPResso, HGVSp quantification — with the cluster job files it ran as |
| `docs/` | coverage, methods, data description, known differences, per-figure and per-module indexes |
| `figures/` | render output. Not committed: a full ChimeraX render alone is 119 MB |

### The package

`cdkbe` is organised so that a panel's numbers and its drawing live together.

* `cdkbe.io` — compression-agnostic read/write, so nothing has to know whether
  a table is gzipped
* `cdkbe.paths` — every location resolved from the repository root, so a
  notebook works wherever it is launched from
* `cdkbe.counts`, `.lfc`, `.lfc_fdr`, `.editing`, `.master_table` — the screen
  statistics
* `cdkbe.source_data.*` — one module per group of panels; each exposes the
  table (`PANELS`) and the drawer (`FIGURES`), which the figure notebooks, the
  renderer and the workbook builder all share
* `cdkbe.beatseq`, `.essentiality`, `.conservation` — the single-cell, early
  essentiality and KLIFS/MSA analyses

[`docs/MODULES.md`](docs/MODULES.md) is a one-line index of all 91.

---

## The data

[`docs/DATA.md`](docs/DATA.md) describes each tier in detail; the short
version:

| tier | holds |
|---|---|
| `raw/` | as close to the sequencer and the bench as this repository goes: per-sample guide counts, sensor-editing calls, flow exports, MTT plates, xenograft measurements, 288 Sanger traces, and the finished panels collaborators delivered |
| `processed/` | everything derived, and regenerable by the notebooks here: count tables, per-gRNA editing, LFC + FDR, DESeq2 results, single-cell statistics |
| `reference/` | external annotation not produced here: the KLIFS pocket alignment, conservation, AlphaMissense, gnomAD, the guide library, OncoKB/CIViC exports, and the ChimeraX renders the structure figures use |
| `external/` | not shipped. `NOT_INCLUDED.csv` lists 4.5 GB of per-sample intermediates and archival material with its size and how to regenerate or obtain it. Nothing a main figure needs is in that list |

**Raw sequencing reads** are deposited at NCBI under BioProject
[PRJNA1390412](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1390412). `data/raw/`
here begins one step downstream of them — per-sample guide counts and per-sample
sensor editing calls — which is where `pipeline/` leaves off. To start from the
reads, fetch them from the BioProject and run `pipeline/`.

```bash
python3 scripts/verify_data.py             # check data/ against the manifest
```

```bash
python3 scripts/verify_data.py --refresh   # record the current state
```

Rerunning a notebook rewrites what it produces, so drift under
`data/processed/` after a rebuild is expected. Drift under `data/raw/` is not —
nothing should ever write there.

### Sanger sequencing

All 288 traces ship, in three groups: the edited locus amplified out of
treated cells (187) and the library cloning colonies (45), both under
`data/raw/sanger/`, and the validation guide plasmids (56), which sit with the
guide designs under `data/reference/validation_guides/`.
`data/raw/sanger/INDEX.csv` indexes all three, one row per `.ab1` — the sample,
primer, gene, variant, treatment and read direction parsed out of the
filename, plus the run date, instrument, well, dye, read length and median
Phred read out of the ABIF header.

`basecalls.fasta.gz` makes the set searchable without trace software, but a
base-edited locus is a mixed peak that a basecall flattens, so quantification
goes back to the `.ab1`.

### Structure panels

The ChimeraX renders behind figures 1f, 2c, 2j, 3c, 4c–4e, 6 and
supplementary 15e cannot be drawn from a table, and ChimeraX is a GUI
application. So the renders those panels use ship under
`data/reference/structure_renders/` and the code reads them from there — a
clone with no ChimeraX still assembles every figure.

To re-render them you need the application:

```bash
python3 scripts/chimerax/render.py --list
```

```bash
python3 scripts/chimerax/render.py
```

```bash
python3 scripts/sync_structure_renders.py
```

The `.cxc` command files are generated — edit `scripts/chimerax/generate_*.py`,
not the files under `cxc/`.

---

## The tests

The suite is not a formality. It pins the things that would otherwise drift
without anyone noticing:

* **panels against published artwork.** Where a published SVG was available,
  the redraw is compared to it character for character and fill for fill — the
  competition-assay heatmaps and the T47D cell-cycle panels reproduce exactly.
* **statistics against the values printed on the figures.** Every p-value on
  the in vivo panels, the DEG counts on 3f, the bar labels on 3g, the eighteen
  counts on 4b.
* **self-containment.** No file may resolve a path outside the repository —
  and no notebook may have stored one in its output, which is how an absolute
  path would reach a reader who never ran anything.
* **the shipped notebooks.** None holds an error output: what you see rendered
  is what ran.
* **declared dependencies.** Every third-party import must be in
  `pyproject.toml`; the check walks the import graph rather than a list.
* **the shipped workbooks.** Every panel with source data has a sheet, every
  sheet has row labels, every INDEX line says what its sheet is.
* **the data.** The tests read `data/`, so a pass is also a statement that the
  tables are the ones the figures were made from.

---

## Documentation

| file | what it is |
|---|---|
| [`docs/REPRODUCE.md`](docs/REPRODUCE.md) | the run order, and how to rebuild each stage |
| [`docs/METHODS.md`](docs/METHODS.md) | what each stage of the analysis does, and the conventions that bite |
| [`docs/COVERAGE.md`](docs/COVERAGE.md) | every published panel and its status |
| [`docs/FIGURES.md`](docs/FIGURES.md) | every panel: what it shows, what draws it, what the file is called |
| [`docs/MODULES.md`](docs/MODULES.md) | one-line index of every module |
| [`docs/DATA.md`](docs/DATA.md) | the data tree, tier by tier |
| [`docs/DISCREPANCIES.md`](docs/DISCREPANCIES.md) | every place the published figure, the mapping and the code disagree |

`FIGURES.md`, `MODULES.md` and `COVERAGE.md` are generated:

```bash
python3 scripts/make_figure_docs.py
```

```bash
python3 scripts/make_coverage_doc.py
```

### Known differences

[`docs/DISCREPANCIES.md`](docs/DISCREPANCIES.md) is worth reading before
comparing a regenerated panel with the published one. It records every place
the two differ, what was checked, and what this repository does about it — a
colour scale centred differently, a filter applied before canonical editing
calls existed, a selection rule reconstructed from the figure it produced.
None of them blocks a rebuild, and several are choices with two defensible
answers rather than errors.

Ten panels are marked **published image** in `docs/COVERAGE.md`: they are
included as the figure that was published rather than redrawn from a table.
Six others in the same group are redrawn — four from the flow-cytometry export
behind them, and two from values read back out of the published figure, which
is exact for a plot drawn against a labelled linear axis (`cdkbe.svg_values`).

---

## Citing

If you use this code or data, please cite the manuscript. (Citation to be
added on publication.) The raw sequencing reads are NCBI BioProject
[PRJNA1390412](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1390412).

## Licence

MIT — see [`LICENSE`](LICENSE). The copyright line names the repository's
author; change it if it should name the laboratory or the institution instead.
