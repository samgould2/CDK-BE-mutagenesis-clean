# Methods

What each stage of the analysis does, in the order the data moves through it,
and which module to read for the detail. Every statement here is about code that
is in this repository; where a step's inputs are not, that is said explicitly.

`docs/FIGURES.md` lists the panels themselves and what draws each one.
`docs/DATA.md` describes the data tree. `docs/REPRODUCE.md` is the run order.

---

## 1. Pooled base-editing screens

### 1.1 Counts

`cdkbe.counts` turns per-sample sequencing counts into per-screen count and RPM
tables. Each screen yields three count flavours, differing in how strictly a read
had to match:

| flavour | reads assigned by |
|---|---|
| `barcode` | the sensor barcode |
| `proto` | the protospacer |
| `matched` | barcode and protospacer agreeing |

Most analyses use `barcode`. The compound-mutant screen uses `proto`, because its
guides share barcodes — a detail that silently changes which hits appear if it is
got wrong.

### 1.2 Log fold change and FDR

`cdkbe.lfc_fdr` takes counts and RPM to a log2 fold change, an empirical p-value
and a Benjamini–Hochberg FDR per guide. Three properties of the original are
preserved deliberately and documented in the module:

* **The bootstrap was unseeded.** The null is drawn with `np.random.choice` and no
  seed was ever set, so two runs of the original give different p-values. `seed`
  is now an explicit argument. A rerun is therefore reproducible, but recomputed
  FDRs will not match the shipped tables to the last digit.
* **`one_sided_FDR_v2` is two-sided** despite its name: it computes a high-tail
  and a low-tail p-value and keeps the smaller. It is called `empirical_fdr` here.
* **The function was defined twice** in the source notebook, the second version
  incompatible with the first. Only the first is ported.

The shipped LFC/FDR tables remain the canonical record; `cdkbe.master_table`
stacks every screen and baseline into one long-format table per editor.

### 1.3 Sensor editing

`cdkbe.editing_pipeline` reduces per-sample sensor calls to four tables per
screen: every distinct edit with its DNA change, the same restricted to canonical
edits, and reads per (gRNA, protein consequence) in both flavours. **Canonical**
means the edit is the one the editor makes (`A>G` for ABE, `C>T` for CBE).
`Canonical_window` records whether it also landed in the expected window but is
*not* part of the filter.

`cdkbe.editing` collapses those to one variant call per gRNA. It has one trap
worth knowing:

> **Two SAV conventions exist.** The shipped `SAVs_canonical/*.csv` are built
> *without* first dropping reads that carry a substitution the editor cannot
> make; `logo_plots_KB_compound_mut.ipynb` applies that filter and gets a
> different table (16,026 rows against 9,850 for CBE subpool 1).
> `editing.savs_from_full_canonical(..., own_edits_only=)` reproduces either.
> The default reproduces the shipped tables. Picking the wrong one changes which
> variant a guide is attributed to, without any error.

### 1.4 Figures from the screens

* `cdkbe.lollipop` — median LFC per guide against codon position, one panel per
  gene. Marker area scales with −log10(FDR); guides clearing both the FDR and LFC
  cutoffs are opaque and labelled. Shaded bands are KLIFS pocket positions.
  The gRNA denominator printed in each title was a hard-coded literal in the
  source notebook; see `docs/DISCREPANCIES.md` item 12 for where it disagrees
  with this repository's own counts.
* `cdkbe.logo_panel` — per-hit LFC bars beside a sequence logo of what that guide
  actually edits. Each logo column sums to 100: the observed edits at their
  measured %Reads, the reference residue taking the remainder, so glyph height is
  the fraction of sensor reads carrying that residue.
* `cdkbe.heatmap` — the cross-drug resistance heatmap and hits-per-drug bars.
  Row linkage is computed on the drug columns *before* the T0 column is appended,
  so baseline fitness never drives the drug clustering.
* `cdkbe.lfc_correlation` — LFC correlation clustermaps.

---

## 2. Validation assays

### 2.1 MTT dose response

`cdkbe.dose_response` fits four-parameter log-logistic curves. The source
notebook redefined the fitting function eleven times with different bounds and
normalisations; the definition ported here is its **cell 125**, the one live for
every panel drawn from it.

Viability is `(measurement − background) / no_drug_avg`, where `no_drug_avg` is
that cell line's own DMSO wells minus the same background.

> **Three background conventions are in use**, and they are not cosmetic — each
> moves every viability value on the panel:
>
> | convention | background taken from |
> |---|---|
> | `none` | nothing subtracted |
> | `cell_line_no_cells` | wells whose *cell line* is "No Cells" |
> | `drug_no_cells` | wells whose *drug* is "No Cells" |
>
> Supplementary 2a uses `none`, 2b and 2c use `cell_line_no_cells`, and figure 3d
> and supplementary 9d–9g use `drug_no_cells`. The convention is named per panel
> in each module's `SPEC` rather than assumed.

LC50 equality between two curves is tested with an extra-sum-of-squares F test,
Holm-corrected across the comparisons on a panel.

### 2.2 In vivo

`cdkbe.invivo_panel` draws the treatment-grouped endpoint panel for the xenograft
cohorts. x is the treatment condition, so the WT and mutant tumours of the same
mouse sit together and grey lines can join matched flanks; light fill is WT,
solid is mutant. `cdkbe.source_data.leukemia_burden` reads the B-ALL
transplantation flow exports directly.

### 2.3 Competition assays

A competition assay mixes edited (fluorescent) cells with unedited ones and
reads the edited fraction by flow at T0 and at the end of treatment. The number
every panel shows is the ratio of the two, so a fold change above 1 means the
variant was enriched under drug.

`data/raw/competition_assays/` holds the delivered fold changes, three
replicates per condition, with the collaborator's own file-by-file guide in
`README.txt`. `cdkbe.source_data.competition_assays` computes everything else
the panel prints: the mean, the SEM over the three replicates, and a two-sided
unpaired t-test of each variant against the panel's non-targeting arm within
the same treatment, starred `****` below 1e-4, `***` below 1e-3, `**` below
0.01, `*` below 0.05 and `ns` otherwise. That is the test the one delivered
file that carries its own statistics names in its header, and it reproduces
every star on all five published heatmaps.

The heatmaps share one colour scale, `RdBu_r` from 0 to 5.5 centred at 1 --
seaborn's `center` truncates the colormap rather than bending the norm, which
is why the colourbar's ticks stay evenly spaced. Two panels orient differently
and the module follows each: figure 2f puts variants across and treatments
down, the CDK12 family does the opposite and rules a thick line between gene
groups.

Supplementary 7c and 7e are the mouse B-ALL assays, drawn as bars with their
replicates on top. This repository holds their values but not the published
figure, so that style is its own rather than a reproduction; see
`docs/DISCREPANCIES.md` item 24.

### 2.4 The LNCaP and PC3 CDK12 panels

Supplementary 11c and 11D-11K were produced outside this repository and came
as complete packages: the values that were plotted, the script that plotted
them and the published image, one folder per panel, all under
`data/raw/collaborator_panels/`. `cdkbe.source_data.lncap_pc3` reads the values
back and redraws each panel with the parameters those scripts used.

**No statistic is recomputed there.** The DESeq2 fits (pyDESeq2 0.6.0, genes
with at least 10 counts in at least 3 samples, `~condition`, Cook's refit,
Wald test, unshrunk MLE log2FC), the preranked GSEA runs (gseapy 1.3.2, 1000
permutations, seed 42) and the qPCR 2^-ddCq with its Bonferroni correction are
the contributing authors' own. What this repository checks is that the redraw
still says what the packages say -- every DE count their READMEs print comes
back exactly. Supplementary 11D's PCA is from nf-core/rnaseq 3.25.0's DESeq2 QC
step, not from pyDESeq2, which is why its axes are labelled 31% and 15%.


---

## 3. Structure and conservation

* `cdkbe.conservation` — per-residue evolutionary conservation for the nine
  screened CDKs from a deep, taxon-stratified ortholog search, plus the
  kinome-wide version at the 85 KLIFS pocket positions, and the tests asking
  whether conservation predicts what a base edit does in the screen.
* `cdkbe.essentiality` — the early (T0 vs plasmid) residue-essentiality score,
  its corrections and the orthogonal checks on them, peak finding along the
  kinase domain, the KLIFS-position detail figures, and the dbNSFP comparison.
* `scripts/chimerax/` — structure renders, behind figures 1f, 2j, 3c, 4c, 4d,
  4e and 6. The `.cxc` command files are **generated** rather than committed:
  ChimeraX silently aborts a render when a path is relative, so each file needs
  an absolute `cd` into this checkout.

  ```bash
  python3 scripts/chimerax/render.py          # every panel a figure needs
  python3 scripts/chimerax/render.py --list   # what would run, and why
  ```

  > **This has to use the application binary, not `bin/ChimeraX`.** Saving an
  > image needs an OpenGL context, and `--nogui` has none — it fails with
  > `LimitationError: Unable to save images because OpenGL rendering is not
  > available`. `ChimeraX.app/Contents/MacOS/ChimeraX` does have one and runs
  > the same commands to completion. `--nogui` is still right for anything that
  > does not save an image, such as `validate_structures.py`.

  The renders take about ten minutes in six jobs, chosen to cover every needed
  image in as few ChimeraX launches as possible. `scripts/render_figures.py`
  then copies them under their manuscript names — `main_fig2j.png`,
  `main_fig4c_CDK4_active_site_zoom_labeled.png` — rather than re-encoding a
  3,000-pixel render into a PDF.

  ### The leader-line panels (figure 6 and supplementary 15e)

  Two panels are not command files. `label_consensus.py` opens the structure,
  projects each marked residue through the **live camera** into frame
  coordinates, parks its label in the left or right gutter, and joins the two
  with a leader line; then it writes the `2dlabels` commands, the rendered PNG,
  and a layered SVG whose structure is a raster and whose labels and leaders are
  real `<text>` and `<path>`, so they stay editable in Illustrator.

  ```bash
  # figure 6: one line per site
  ChimeraX --cmd "runscript label_consensus.py core7 variants \
      consensus_strong recurrent,strong_singleton"
  # supplementary 15e: a block per peak, and the camera pulled back for them
  ChimeraX --cmd "runscript label_consensus.py core7 full peaks recurrent \
      panel=peaks.cxc sites=MSA_peak_sites_labels.csv \
      stem=peaks_leader zoom=0.66"
  ```

  > **The window's shape is part of the figure.** Every label position is a
  > fraction of the ChimeraX frame, so the published 1400×796 window is set in
  > the script rather than left to the display — ChimeraX honours a taller
  > request on a taller screen, and a square frame stretches every vertical gap
  > by 1.76 and distorts the structure underneath. `window=WxH` overrides it.

  `render.py` runs both as part of its normal pass; it reads
  `docs/source_data_mapping_supp.csv` as well as the main one, so a
  supplementary panel that uses a render is picked up the same way.

---

## 4. Single-cell CROP-seq (figure 5)

164 guides (92 ABE + 72 CBE) against CDK2/4/6/7/8/9/12/13/19 plus non-targeting
controls, in A549 and T47D, untreated, 10x 3′ SC3Pv4. The design, including the
hazards below, is recorded in `cdkbe.beatseq.paths`.

### 4.1 What is and is not in this repository

The delivered bundle's **derived tables are all copied in**, under
`data/processed/beatseq/` — 248 of them. The single-cell objects themselves are
not: roughly 34 GB of `.h5ad`, the last row of
`data/external/NOT_INCLUDED.csv`.
Nothing that draws a panel or recomputes a shipped statistic opens one.
`cdkbe.beatseq.singlet_io` is the only module that would, and it is named
`singlet_io` rather than `io` so it cannot be confused with `cdkbe.io`, the
compression-agnostic table reader.

### 4.2 Normalisation and the perturbation statistic

`cdkbe.beatseq.norm` differs in three ways from the pseudobulk that came with
the data, each of which changes the result:

1. **Leave-one-out bias.** The shipped `control_norm` layer z-scores against the
   whole clean-control pool, and each NT guide's own cells are inside that pool.
   Its pseudobulk is shrunk toward zero while targeting guides are not, which
   manufactures a ~14% "targeting > control" difference out of nothing. The
   control pool for guide *g* here always excludes *g*'s own cells.
2. **Incompatible gene universes.** The per-arm sets shipped with the data span
   5,563–10,331 genes. One universe is used: genes detected in ≥5% of cells in
   every arm.
3. **An n-unbiased perturbation statistic**, *D*, in place of one whose value
   depends on how many cells a guide happens to have.

### 4.3 Severity

`cdkbe.beatseq.nullcal` calibrates *D* against a **size-matched** pseudo-guide
null:

> severity = (D − μ(n_g)) / σ(n_g), with μ and σ interpolated in log n

The earlier pooled null was justified on the grounds that *D* is n-unbiased.
Its *mean* is; its *standard deviation* is not, which is why the null has to be
size-matched. The original's hard-coded draw sizes also failed to span the real
guide range — 54%, 69%, 24% and 0% of guides in the four arms fell outside any
size the null was drawn at — so sizes are now drawn from the observed cell
counts.

`cdkbe.beatseq.severity` assembles this and runs the pre-registered tests.
**It recomputes the shipped `severity.csv` from the shipped tables to within
1e-12**, which `tests/test_figure5.py` pins: the statistic figure 5j and 5k plot
is derived here, not merely read from a CSV.

### 4.4 Known data hazards

* **Four targeting protospacers appear in both editor arms installing different
  alleles** (including CDK9 L156F / I157V). Same spacer, same integration
  biology, different amino-acid change — the cleanest controlled contrasts in the
  experiment. Any cross-arm join must key on `(gRNA_id, editor)`; joining on the
  spacer alone silently collapses these pairs.
* **Two NT controls share a protospacer**, so one is never recovered and the
  effective number of independent NT units is 9, not 10. That is the hard ceiling
  on any guide-level permutation null.
* **gem_group is perfectly nested inside arm.** Both T47D arms have exactly one,
  so T47D cannot support any statistic relying on technical replication.
* Ribosomal, mitochondrial and immediate-early genes are excluded from global
  fits: the first two are short, highly expressed and depth-sensitive and will
  manufacture a gene-length slope on their own; IEGs are simultaneously the best
  short-half-life markers and the worst dissociation-stress artifact.

### 4.5 Panels

`cdkbe.beatseq.fig_panels` draws all 46 source panels, each as its own file for
assembly in Illustrator — PDF and SVG carry live text, point clouds are
rasterised at 400 dpi and nothing else is. Its blocks are lettered independently
of the manuscript, so `panel_f` is **not** figure 5f.
`cdkbe.source_data.figure5.SPEC` is the mapping between the two, with the
evidence for each entry; `tests/test_figure5.py` re-runs that evidence.

---

## 4b. Predicted CDK9-ligand complexes (supplementary 6)

`cdkbe.boltz` analyses Boltz-2 predictions of CDK9 bound to four ligands
(atuveciclib, KB-0742, the PROTAC KI-CDK9d-32, voruciclib) in three backgrounds
(WT, L156F, I157V), 3 seeds x 3 structural samples = **108 predicted
structures**. Confidence is recorded per structure, affinity per seed, so there
are 36 affinity records, not 108.

Each protein residue gets a continuous contact score per interaction channel,

> S = 1 / (1 + exp[(d_min - d_mid) / tau])

over five channels (side-chain carbon, backbone polar, side-chain polar,
charged, aromatic carbon) plus a steric-clash descriptor taken from the single
most overlapping atom pair. Features whose maximum across all 108 structures
falls below 0.15 are dropped, leaving **132**; restricted to the WT ATP pocket,
**58**. Both counts are stated independently in the analysis README and both
come out of this port exactly.

The PCA is weighted: each feature's within-ligand standard deviation is averaged
over ligands and the feature scaled by 1 / (0.15 + that), so contacts that
fluctuate between predictions of the same ligand do not dominate the leading
components.

> **Two things here will catch you out.** The notebook defines the ATP pocket
> twice -- a labelled region map covering 32 residues, and a flat set of 37 that
> every PCA actually uses. And an SVD's axis signs are arbitrary, so the same
> data can come back mirrored on a different machine; `ifp.SIGN_ANCHORS` pins
> each component to the orientation that was published.

These scores measure atom-class proximity. They are not validated hydrogen
bonds, electrostatics or stacking geometries, and the binding probabilities are
model output rather than measurements.

---

## 5. Bulk RNA-seq

* `cdkbe.source_data.rnaseq` / `rnaseq_figures` — the PC3 series behind figures
  3f–3h, and the DESeq2 volcanoes in supplementary 5e.
* `cdkbe.source_data.rnaseq_t47d` / `rnaseq_t47d_figures` — the T47D CDK4 allelic
  series behind figures 5f–5i. Both volcanoes use contrasts refit with DESeq2
  independent filtering disabled, so the two panels share an identical gene
  universe (11,489 genes) and their DEG counts are comparable. The E2F score is
  the mean z-score across the 198 expressed Hallmark E2F target genes, tested
  one-sided against the control and Holm-corrected over the three variants.

---

## 6. Conventions that bite

**`pdf.fonttype` is read when a figure is *written*, not when it is drawn.**
A figure drawn under the right rcParams and saved outside that context comes out
with Type 3 text instead of TrueType — the same picture, a different glyph
encoding, and no longer matching the published file. Everything that writes a
figure here forces it at save time; `scripts/render_figures.py` does it centrally.

**Compression is not part of a table's identity.** `cdkbe.io.resolve` accepts
`foo.csv` and `foo.csv.gz` interchangeably, and `cdkbe.io.write` gzips anything
over 1 MB. Code that calls `pd.read_csv` on a literal `.csv` path will miss a
gzipped twin — this is how `pseudo_guide_null.csv` and
`cdk4_program_analysis/pub_fig_square.py` were both broken until they were routed
through `io`.

**Figure numbering is not module numbering.** Panel letters in `fig_panels`, in
the 2×2 RNA-seq grid, and in the manuscript are three different schemes. Treat
any coincidence between them as a coincidence.

### 2.5 T47D cell cycle

T47D lines carrying a base-edited CDK variant were treated for seven days with
DMSO, palbociclib or serum starvation, stained, and gated into G0, G1, S and
G2/M. `data/raw/collaborator_panels/T47D_cell_cycle/` holds the FlowJo export
of those four gate frequencies -- six cell lines by three treatments by three
replicates -- with the delivered notebook, its style module and the FlowJo
workspace beside it. `cdkbe.source_data.t47d_cell_cycle` draws figure 4k, 4l
and supplementary 13c, 13d from it.

Two points of arithmetic. The four gates are siblings on one parent
population, so they do not sum to 100% (totals run 91 to 99%); the stacked
panels rescale each condition to 100, the %S-phase panels do not. And the
significance brackets are Welch's two-sided t-test of each variant against its
own non-targeting arm, Bonferroni-corrected across **every comparison the panel
draws** -- four per panel, not two per treatment.

Serum starvation is in the export and is not drawn on any of the four
published panels.

### 2.6 Allele frequency across selection

The same measurement appears three times: a base-edited population is
amplified at the edited locus, Sanger sequenced before and after drug
selection, and the variant allele's share of the trace is read off. Figure 4f
is the A549 CDK2/4/6 panel, figure 4g the T47D panel under palbociclib, and
supplementary 13a the MCF7 panel under palbociclib.

**None of the three is recomputed here; all three are held as tables.** 4f's
numbers were literals in the original notebook, transcribed to
`data/processed/validation/allele_frequency_4f.csv` (`docs/DISCREPANCIES.md`
item 20). 4g and 13a were read back out of the published figures: both are
matplotlib line plots against a labelled linear y axis, so `cdkbe.svg_values`
inverts the axis and reads each line's two vertices back out, exactly.

That inversion is validated before it is used.
`scripts/recover_selection_tables.py` first runs it over the two T47D
cell-cycle panels, whose published SVGs *and* underlying FlowJo export are
both in this repository, and refuses to write anything unless the recovered
bar heights equal the measured gate frequencies. They agree to 2e-7
percentage points, which is the six-decimal precision of the SVG coordinates.
Fifteen of the sixteen recovered values are whole percentages, as 4f's
hand-quantified numbers are; the one fraction, MCF7 CDK2 S261G after
selection, is 4.264871%.

Nothing in the repository reads colour out of a published image. A value
encoded in a fill is quantised to eight bits per channel, and for figure 4m it
is a within-arm z-score rather than the expression matrix behind it, so that
panel is included as the image that was published rather than recovered
(`docs/DISCREPANCIES.md` item 32).

### 2.7 Sanger sequencing

Three sets of traces ship, 288 `.ab1` files in all, as the sequencing vendors
delivered them: the edited locus amplified out of treated cells
(`data/raw/sanger/endogenous_editing/`), the library cloning colonies read at
the U6 promoter (`data/raw/sanger/library_colonies/`), and the individual
validation guide plasmids, also at U6, which sit under
`data/reference/validation_guides/` beside their designs.

The amplicons were read in both directions with the primers in
`endogenous_editing/cdk_seq_primers.xlsx` and `cdk246_seq_primers.xlsx`; the
SnapGene transcript maps beside them show where each primer sits. Capillary
runs are BigDye v3 on a 3730xl; four deliveries are Plasmidsaurus and carry no
run date in their headers.

`data/raw/sanger/INDEX.csv` is one row per trace, written by
`scripts/make_sanger_index.py`, which parses the sample labels and reads the
run metadata out of each ABIF header. Where a parsed claim can be checked it
is: 55 of the 56 guide ids read out of filenames exist in the screening
library, and the one that does not is flagged rather than corrected
(`docs/DISCREPANCIES.md` item 33).

Quantification goes back to the chromatogram. `basecalls.fasta.gz` is there to
make the set searchable, but a base-edited locus is a mixed peak and a
basecall flattens it.
