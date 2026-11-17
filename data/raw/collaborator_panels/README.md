# Panels delivered by collaborators

These are the figure panels that were produced outside this repository and
handed over as finished work. They are vendored here so the repository is
self-contained: nothing under `src/` or `notebooks/` reaches outside the repo
for them.

Two kinds of thing sit here, and the difference matters for what can be
reproduced:

* **`LNCaP/` and `PC3/` carry their own data and code.** Each panel folder
  holds the values that were plotted, the script that plotted them, and the
  published image. `src/cdkbe/source_data/lncap_pc3.py` reads the values back
  and redraws the panels, and the published image is what the tests compare
  against.
* **`T47D_MCF7/` holds published figures.** Sixteen `.svg`/`.png` files. Six
  of them redraw here: the four cell-cycle panels from the FlowJo export in
  `data/raw/collaborator_panels/T47D_cell_cycle/`, and figure 4g and
  supplementary 13a from values read back out of the figure itself
  (`cdkbe.svg_values`). The other ten are included as the images that were
  published; `docs/COVERAGE.md` marks them **published image**.

## Provenance

Copied on 2026-11-16 from the lab Dropbox, at the paths the manuscript's
source-data mapping gives. Each set was delivered from a collaborator's own
folder; what matters here is which delivery a panel belongs to, so that is
what is recorded rather than the folder it sat in.

| here | delivery |
|---|---|
| `T47D_MCF7/` | the T47D and MCF7 panels: competition, cell cycle, Sanger logos and RNA-seq |
| `LNCaP/` | the LNCaP CDK12 I733V package: DESeq2, qPCR, GSEA, logo plots |
| `PC3/PC3_panel_k/` | the PC3 CDK12 volcano package |
| `PC3/Cdk12_logo_plots/` | the PC3 CDK12 amplicon logo plots |

`data/MANIFEST.csv` records each file's source and SHA-256.

The one FlowJo workspace here, `T47D_cell_cycle/20260908_T47D_cell_cycle_7dpt.wsp`,
had the absolute path rewritten out of each of its sample URIs: they began with
a home directory and a personal Dropbox folder, and now start at the
acquisition directory. The gating is unchanged, and no code reads the
workspace -- the panels come from the exported table beside it.

`LNCaP/README.txt`, `PC3/PC3_panel_k/README.txt` and
`PC3/Cdk12_logo_plots/README.md` are the authors' own notes, left exactly as
delivered. They are the authority on how each panel was computed — including
the two corrections they record, a `lfc_shrink` reference-level bug in the PC3
DESeq2 contrasts and a reverse-strand CIGAR bug in the logo-plot pileup.

The competition assays that go with these figures are not here; they live in
`data/raw/competition_assays/`, which has its own README.

## Which panel is which

| panel | folder |
|---|---|
| supplementary 11c | `LNCaP/qPCR_heatmap/` |
| supplementary 11D | `LNCaP/RNAseq_panel_d/` |
| supplementary 11E | `LNCaP/RNAseq_panel_e/` |
| supplementary 11F | `LNCaP/RNAseq_panel_f/` |
| supplementary 11G | `LNCaP/RNAseq_panel_g/` |
| supplementary 11H | `PC3/Cdk12_logo_plots/` |
| supplementary 11K | `PC3/PC3_panel_k/` |
| figure 4g, 4i–4m; supplementary 13a–13j | `T47D_MCF7/` |

The supplementary 11 letters are not guesswork: the delivered folders are named
`RNAseq_panel_d` … `RNAseq_panel_g` after the panels they make, and the
deposition README names the qPCR heatmap as the one before them. The T47D and
MCF7 files are matched to panels by the mapping workbook alone, one file per
row.
