# Data

Three tiers, each regenerable from the one above it.

```
raw/        what came off the sequencer and the bench, per sample
processed/  everything the notebooks derive from raw/
reference/  external annotation, not produced by this project
external/   what is NOT shipped, with its size and how to get it back
```

`data/MANIFEST.csv` lists every shipped file -- all 3,044 of them -- with its
size before and after compression and a SHA-256. A file that came through the
original ingest also carries its source path in the tree it came from; one
written later, here, carries an empty source and says so in its note.

```bash
python3 scripts/verify_data.py           # compare the tree with the manifest
python3 scripts/verify_data.py --refresh # record a rebuild's new checksums
```

`--prune` drops the row of a file removed on purpose and `--adopt` adds a row
for one that has none, so neither kind of change has to be made by hand.

## Compression

Large text tables ship gzipped — the sensor-editing CSVs go from 1.8 GB to
298 MB, about 7x. **No notebook mentions compression.** `cdkbe.io.read` takes a
logical path and finds whichever of `.csv`, `.csv.gz` or `.zip` is on disk:

```python
from cdkbe import io, paths
df = io.read(paths.RAW_EDITING / "CBE_subpool1" / "some_sample.csv")   # finds .csv.gz
```

So a file can be re-compressed, or decompressed, without touching any code.

## The chains you can rerun

**counts -> LFC/FDR.** `data/raw/screen_counts/` holds per-sample guide counts
exactly as produced by the demultiplexing step. `notebooks/02_counts_and_qc`
aggregates them into `data/processed/count_tables/`, `notebooks/04_lfc_fdr`
turns those into `data/processed/rpm_tables/` and
`data/processed/lfc_fdr_tables/`. Every figure that shows an LFC reads the
last of those.

**raw editing -> processed editing.** `data/raw/screen_editing/` holds the
per-sample, per-read sensor calls. `notebooks/03_editing` reduces them to the
per-gRNA tables in `data/processed/editing/` — `*_full_canonical`, `*_HGVSp`,
and the single-amino-acid-variant tables in `SAVs_canonical/`. The per-sample
intermediates between those two points are the single largest tree in the
original repository (1.8 GB) and are **not** shipped; they are rebuilt by
rerunning that notebook.

## Experiments contributed from elsewhere

Two trees under `raw/` were produced outside this project's own sequencing and
bench work.

`raw/competition_assays/` holds the flow competition-assay fold changes behind
figure 2f and supplementary 7c, 7e, 10f, 10g, 11b and 11i: three replicates per
condition, with the contributors' own `README.txt` explaining each file. Those
replicate values are the form the experiment is held in here, and every panel
built on them is drawn from them.

Two delivered RNA-seq trees sit under `processed/rnaseq/`, as `PC3_RNAseq/`
and `LNCaP_RNAseq/`. Both were renamed on arrival -- the delivered names
carried spaces, which makes every path that quotes them fragile -- so their
own READMEs and scripts, left as delivered, still use the old names. Nothing
else changed in either tree, and `data/MANIFEST.csv` records each file's
delivered path.

`raw/collaborator_panels/` holds the LNCaP and PC3 CDK12 packages
(supplementary 11c, 11D-11K) and the T47D and MCF7 panels (figure 4g, 4i-4m
and supplementary 13a-13j). The first group carries its values, its code and
its published image, one folder per panel, and redraws here. Of the sixteen
T47D/MCF7 panels, six redraw: the four cell-cycle panels from their FlowJo
export, and figure 4g and supplementary 13a from values read back out of the
published figure (`cdkbe.svg_values`). The other ten are included as the
images that were published; `docs/COVERAGE.md` marks them **published
image**.

## Sanger sequencing

`raw/sanger/` holds every Sanger trace this repository ships, with one table
over all of them. `raw/sanger/INDEX.csv` has a row per `.ab1`: the sample
label and sequencing primer the vendor recorded, the gene, variant,
treatment, selection state and read direction parsed out of that label, and
the run date, instrument, well, order, dye, read length and median Phred read
out of the ABIF header. `scripts/make_sanger_index.py` writes it and
documents every parsing rule.

Three groups. `endogenous_editing/` is the edited locus amplified out of
treated cells -- the validation and post-selection genotyping behind figures
2, 3 and 4 -- with the sequencing primer tables and the SnapGene transcript
maps beside it. `library_colonies/` is the library cloning QC, read at the U6
promoter, with the protospacer and barcode called per colony. The validation
guide plasmids are the third: they stay under `reference/validation_guides/`
beside the designs they confirm, and are indexed with the rest.

`raw/sanger/basecalls.fasta.gz` carries the vendor basecall for every trace,
so the set can be searched without trace software -- but a base-edited locus
is a mixed peak, which a basecall flattens, so quantification goes back to
the `.ab1`. `source_data_tables/Supplementary_Data_Sanger_sequencing.xlsx` is
the same index as a two-sheet workbook, for the supplement.

## What is not here

`data/external/NOT_INCLUDED.csv` lists 4.5 GB of material left out, each row
with its size, what it contains and how to obtain or regenerate it. The
largest entries are the per-sample editing intermediates (1.8 GB, regenerable),
every figure ever rendered during the project (865 MB), and the off-target
search output (689 MB).

Nothing in that list is needed to reproduce a main figure.

### Raw sequencing reads

The raw reads are not in this repository. They are deposited at NCBI under
**BioProject [PRJNA1390412](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1390412)**.

`data/raw/` begins one step downstream of them: per-sample guide counts and
per-sample sensor editing calls, which is where `pipeline/` leaves off. To go
back further -- demultiplexing, sensor extraction, CRISPResso, HGVSp
quantification -- fetch the reads from the BioProject and run `pipeline/`;
`pipeline/README.md` gives the steps and the cluster job files they ran as.

## Self-containment

Every path any notebook, module or script resolves lies inside this repository.
`tests/test_self_contained.py` enforces that mechanically: it parses each `.py`,
`.ipynb`, `.cxc` and `.sh` under `src/`, `tests/`, `scripts/`, `notebooks/` and
`pipeline/` and fails on a string literal naming a home directory or one of the
trees this repository was assembled from. It also checks that every path
constant in `cdkbe.paths` and `cdkbe.beatseq.paths` points at a file that is
actually here.

One file is exempt, for a stated reason: `cdkbe.beatseq.paths._ABSENT`
documents the ~34 GB of single-cell `.h5ad` objects that are deliberately not
vendored -- no figure reads them, every BEAT-seq panel comes from the derived
tables under `data/processed/beatseq`.

The script that assembled `data/` from the original trees is not part of this
repository: it could only run with those trees on disk. What it recorded does
ship -- `data/MANIFEST.csv` gives every file's original source path, its size
before and after compression, and its SHA-256.

The ChimeraX command files under `scripts/chimerax/cxc*/` are generated and
gitignored: ChimeraX silently aborts a render given a relative path, so they
must carry an absolute `cd` to wherever the repository sits. Regenerate them
with the `scripts/chimerax/generate_*.py` scripts; the structures they open
ship in `scripts/chimerax/structures/`.

## Tiers in detail

### `data/raw/`

| directory | size | files | contents |
|---|---:|---:|---|
| `screen_editing` | 312 MB | 267 | per-sample sensor editing calls |
| `screen_counts` | 211 MB | 305 | per-sample guide counts |
| `flow` | 185 MB | 125 | flow cytometry exports |
| `collaborator_panels` | 53 MB | 73 | panels contributed from elsewhere, with whatever data and code came with them |
| `sanger` | 46 MB | 260 | 232 `.ab1` traces, their basecalls and `INDEX.csv`, which indexes all 288 |
| `cdk9_boltz` | 8 MB | 252 | Boltz-2 predictions of CDK9 bound to four ligands |
| `mtt_dose_response` | 2 MB | 67 | MTT/FACS dose-response plates |
| `screen_classification` | 0 MB | 122 | guide classification |
| `competition_assays` | 0 MB | 15 | competition time courses, by cell line |
| `invivo_xenograft` | 0 MB | 5 | tumour volumes, Ki67, necrosis |
| `sample_key.xlsx` | 0 MB | 1 | sample -> condition key |

### `data/processed/`

| directory | size | files | contents |
|---|---:|---:|---|
| `beatseq` | 194 MB | 323 | BEAT-seq derived tables: matrices, plus the delivered `outputs/` and `metadata/` trees |
| `rnaseq` | 184 MB | 489 | DESeq2 results, TPM matrices, GSEA |
| `lfc_fdr_tables` | 147 MB | 189 | LFC + FDR per guide per condition, hit tables, T0-vs-plasmid essentiality |
| `editing` | 146 MB | 45 | per-gRNA editing and single-aa variant calls |
| `splicing` | 27 MB | 11 | SpliceAI scores |
| `rpm_tables` | 20 MB | 27 | reads per million |
| `count_tables` | 12 MB | 27 | aggregated counts (barcode/proto/matched) |
| `invivo` | 0 MB | 12 | xenograft tumour volume and IHC summaries |
| `resistance_patterns` | 0 MB | 3 | the KLIFS alignments the cross-CDK panels read |
| `validation` | 0 MB | 3 | allele frequencies for figures 4f, 4g and supplementary 13a |
| `cdk9_boltz` | 0 MB | 1 | the summarised interaction fingerprints |

### `data/reference/`

| directory | size | files | contents |
|---|---:|---:|---|
| `structure_renders` | 93 MB | 67 | ChimeraX renders placed in the figures (not regenerable without ChimeraX) |
| `conservation` | 81 MB | 126 | ortholog + kinome conservation |
| `KLIFS` | 67 MB | 28 | KLIFS pocket alignment, OncoKB resistance |
| `validation_guides` | 22 MB | 143 | validation guide designs, and the 56 guide-plasmid `.ab1` traces |
| `library_generation` | 11 MB | 32 | CDK sequences, edit libraries |
| `alphamissense` | 6 MB | 8 | AlphaMissense |
| `gnomAD` | 4 MB | 9 | population variants |
| `resistance_databases` | 1 MB | 2 | CIViC/OncoKB exports |
| `library` | 1 MB | 2 | final gRNA library |
| `beatseq_panels` | 0 MB | 5 | BEAT-seq panel PDFs referenced by supplementary figures |
| `cell_counts` | 0 MB | 1 | cells per screen arm |
| `disorder` | 0 MB | 9 | per-residue disorder |
