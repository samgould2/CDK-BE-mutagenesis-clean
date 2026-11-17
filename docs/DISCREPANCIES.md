# Reproduction notes

Every place where the published figure, the source-data mapping and the code
that regenerates it differ, what was checked, and what this repository does
about each one.

Nothing here blocks a rebuild. Most are choices with two defensible answers --
which of two layouts to publish, whether to reproduce a figure or correct it --
and they are written down so that a reader comparing a regenerated panel with
the published one can see immediately why the two differ. An item marked
RESOLVED is kept because the reasoning is the useful part.

Ordered by how much they matter.

---

## 1. Figure 3d: which plate and which fit the panel uses

**Status: reproduces the published panel exactly.**

The original cell reads `allele_specific = df`, so which plate it draws depends
on which assignment is live. It is cell 136, `2025_8_31_CDK12_13.csv`, carrying
exactly the eight cell lines the panel asks for (`CDK12 I733V (3)` through
`NT CBE`) and the three CDK12/13 inhibitors. The fit is the
`log_logistic_4_param` defined in cell 125, the last before it, which is also
the one live for supplementary 9d-9g.

Worth stating because the notebook defines that function eleven times with
different normalisation, so the fit a panel uses is a question with a real
answer rather than a formality.

Ported to `cdkbe.dose_response` and `cdkbe.source_data.validation_mtt`, the
panel reproduces the published figure exactly — all eighteen sub-panels, and
every LC50 fold shift and significance mark:

| | BSJ-4-116 | CDK12-IN-2 | HQ461 |
|---|---|---|---|
| CDK12 I733V | 4.0x **** | 2.5x *** | 2.2x ** |
| CDK12 G731K | 2.4x **** | 2.9x ** | 3.4x ** |
| CDK12 M816I_D817N | 5.5x **** | 3.5x **** | 2.8x **** |
| CDK12 G822K | 0.6x ** | 1.2x ns | 2.8x **** |
| CDK12 G822E | 0.8x ns | 1.6x ns | 2.5x **** |
| CDK13 I711V | 3.4x **** | 2.1x **** | 1.7x *** |

`tests/test_supplementary.py::test_figure_3d_statistics_match_the_published_panel`
pins all eighteen. **Nothing is outstanding.**

---

## 2. The editing tables and the raw editing calls are from different runs

`data/raw/screen_editing/` holds a **redo** run; `data/processed/editing/` was
built from the earlier one, whose per-sample outputs are not in the repository.
Rebuilding gives 2.4% fewer reads (345.4M vs 351.9M on CBE_subpool1), and no
subset of the available files reproduces the shipped total.

**It changes no conclusion**: 2,439 of 2,439 top variant calls are identical
and the WT fraction correlates 0.9999. Depth differs, biology does not.

**Decision:** leave as is, or re-run the pipeline so raw and processed agree.

---

## 3. The shipped FDR values come from an unseeded bootstrap

`np.random.choice` draws the null distribution, and the run that produced the
shipped tables set no seed, so those exact FDR values cannot be reproduced
draw for draw. How much that matters was measured rather than assumed.

Measured on ABE_subpool1 / barcode / DMSO: LFC columns are exact, FDR
correlates 0.992–0.999, the largest single difference reaches 0.22, and
**98.7–100% of hit calls agree** at FDR < 0.1 — about one guide in 80 changes
status between runs.

**Decided:** the shipped tables stay the tables of record.
`notebooks/04_lfc_fdr/` rebuilds equivalents with an explicit seed and refuses
to overwrite them.

---

## 4. Figure 1f has two layouts, and the published one is the older

The cell that names `res_ABE_CBE_combined_length_scaled.pdf` offsets each panel
so the kinase domains line up across genes. Measured on the shipped PDF, all
nine panels share one left edge (68.4 pt for every gene; the code produces a
260 pt spread). The figure predates the change.

**Currently:** `figure_1f(align_domains=True)` by default, following the code.
`align_domains=False` reproduces the published panel. Both are tested.

**Open:** which one goes in the manuscript.

---

## 5. Figure 2b used the non-canonical editing tables

`logo_plots.ipynb` predates canonical filtering and reads `SAVs/` and
`*_HGVSp`, where a bystander C>T on an ABE guide still counts. The canonical
tables keep only the editor's own edit.

**Effect:** the hit set and every label are identical; the logo windows are
slightly wider in the published panel (`VFITRDGV` vs `VFITRD` on one guide,
`VLVTR` vs `LVTR` on another).

**Decided, 11 November 2026:** canonical everywhere.
`draw_2b(canonical=False)` still reproduces the published panel.

---

## 5b. Figure 4p: closed -- the selection rule was a significance filter

*(This panel was figure 4n until the November relabelling; see item 28.)*

Both `clinical_vars_cdk6.pdf` and `clinical_vars_cdk6_new.pdf` on disk show
**five** rows -- V47A, R46W, R46W_Q48*, R168H, A167T -- and for a long time
nothing in the notebooks produced that set. `clinical.ipynb` cells 19-25 take
CDK6 guides editing the clinically reported codons (CBE 46/168, ABE 46/47/168)
at >=10% sensor editing and then the top 15 by palbociclib LFC, which yields 14
guides with the canonical tables and 18 without -- several guides per variant,
plus R46Q and R168C, none of which are in the published panel. The published
five are not the top five by LFC either: Q48* outranks R168H and A167T.

**The missing step was a significance filter, and it reproduces the panel
exactly.** Of the 14 guides editing a clinically reported codon, exactly five
have FDR < 0.1 in at least one of the screen's six drug arms:

| guide | variant | best FDR | LFC palbociclib |
|---|---|---:|---:|
| `gRNA_CDK6_targ_8184_ABE` | V47A | 0.0566 | +2.38 |
| `gRNA_CDK6_targ_8422_CBE` | R46W | 0.0534 | +0.38 |
| `gRNA_CDK6_targ_8419_CBE` | R46W_Q48* | 0.0939 | +0.06 |
| `gRNA_CDK6_targ_8100_CBE` | R168H | 0.0970 | -0.99 |
| `gRNA_CDK6_targ_8101_CBE` | A167T | 0.0973 | -1.27 |

One guide per variant, which is why the panel has no duplicates, and ordering
what is left by palbociclib LFC gives the published order. The sixth-best FDR
is R168C at 0.137, so the set is not sensitive to where between 0.1 and 0.13
the line is drawn.

Two details matter and both were worth checking rather than assuming:

* **The threshold is 0.1, not 0.05.** At 0.05 nothing clears -- the smallest
  FDR among all 14 candidates is 0.0534. 0.1 is this repository's threshold for
  a resistance hit throughout; it is printed on the face of supplementary 8c
  and 9b.
* **All six drug arms count, not only the three CDK4/6 inhibitors the panel
  draws.** R46W's significant arm is one of the other three: across
  palbociclib, ribociclib and abemaciclib its best FDR is 0.909.

`_clinical_hits` now applies the filter (`CLINICAL_MAX_FDR`,
`CLINICAL_CONDITIONS`) and returns those five rows in that order. The `n=15`
cap the notebook used is kept but no longer binds.

---

## 28. Figure 4's last three panels were relabelled

The November pass moved them round by one, which is a cycle rather than a swap:

| was | is now | panel |
|---|---|---|
| 4n | **4p** | CDK6 clinical resistance variants |
| 4o | **4n** | xenograft tumour volume, CDK6 R46W |
| 4p | **4o** | xenograft Ki67, CDK6 R46W |

The mapping rows, the registry keys, the drawer names and the files the
drawers write all follow the new letters, so `main_fig4n.pdf` is now the
tumour-volume panel. Anything citing the old letters -- including earlier
entries in this file -- means the panel, not the letter.

---

## 6. Three figure filenames in the mapping do not match the code

| mapping says | the code writes | note |
|---|---|---|
| `cdk9_top_vars_v2.pdf` (2b) | `cdk9_top_vars_v3.pdf` | no code writes `_v2` |
| `clinical_vars_cdk6_new.pdf` (4p, was 4n) | `clinical_vars_cdk6.pdf` | and that line is commented out |
| `Fig2_DEG_counts.pdf` (3f) | — | the panel contains the baseline sub-panel, which only exists in `Fig2b_DEG_counts_with_baseline` |

3f is matched to `Fig2b`, whose counts (21, 62, 244, 408, 230, 360) appear in
the panel and exist nowhere else; the other two are matched to the code.

**Open:** confirmation, particularly for 2b — if the manuscript uses
the `_v2` render, it was made by a version of that cell that no longer exists.

---

## 7. Figure 2b's published PDF is wider than its figsize

The cell says `figsize=(7, 8)`; the published PDF measures about 8 inches wide.
That is the signature of `bbox_inches="tight"`, which expands the canvas so the
long compound labels (`D104N_C106Y`) are not clipped. Saved tight here, which
reproduces it.

---

## 8. Why the analysis was re-implemented rather than re-executed

The original exploratory notebooks are not part of this repository, and the
modules here were written from them rather than around them. The reason is
worth recording, because it is also why each module's docstring cites a
notebook and a cell: in several of those notebooks the state a cell depends on
is not determined by file order.

- a function is redefined through the file, so which definition a cell uses
  depends on execution order rather than position (item 1)
- `one_sided_FDR_v2` is defined twice with incompatible signatures, and the
  file as last saved narrowed its sample list to two of the nine screens
- in the sensor pipeline, the lines that write the output tables are commented
  out
- the figure 4f cell has its `savefig` commented out, though its data is
  written literally in the cell, so the panel itself is unambiguous (item 20)

Each of those is settled in the port, and the settlement is what the tests
pin.

---

## 9. `one_sided_FDR_v2` computes a two-sided test

It draws both tails and keeps the smaller. The name is the only thing wrong
with it, so it is ported as `empirical_fdr`. No numbers change.

---

## 11. Figure 5: which BEAT-seq panel drew each published panel

**Status: RESOLVED. All eleven panels are mapped, redrawn and pinned by tests.**

This item previously said the mapping was only partly established. It no
longer needs confirming: each panel was settled against
something only the right candidate satisfies, and
`tests/test_figure5.py` re-runs every one of those checks. The mapping lives in
`cdkbe.source_data.figure5.SPEC`, with the evidence recorded beside each entry.

| panel | source | how it was settled |
|---|---|---|
| 5a | `PANEL_i2_library_donut_CBE` + `PANEL_i1_library_donut_ABE` | all eight published counts (72/92 gRNAs, 17/16/29/10 and 11/10/61/10) are in these two panels |
| 5b | `PANEL_a1_umap_target_A549` + `PANEL_a3_umap_phase_A549` | pixel-identical embedded rasters; published header "n=18,000 cells" |
| 5c | `PANEL_a2_umap_target_T47D` + `PANEL_a4_umap_phase_T47D` | pixel-identical rasters; "n=17,339 cells" |
| 5d | `PANEL_c1s_editing_by_variant_class_combined` | one axes carrying both cell lines and both brackets (p = 0.018, 0.035); `PANEL_c1` draws the same numbers as two axes |
| 5e | `PANEL_c3_editing_vs_cellcycle_CDK4_T47D` | "percentage points" and "phase-adjusted" appear in the `.ai` and in no other candidate |
| 5f-5i | `pub_fig_square_4panel_mid.pdf`, panels a-d | of the three layouts, only the mid one puts 80/160/240 on the volcano axes; the `.ai` has those and not 60/120/180 |
| 5j | `PANEL_d12w_severity_by_CDK_stacked` | stacked, so not d1/d1s/d2/d2s; "severity (SD of the size-matched null)" is absent from the `.ai`; and the published plot area is 1.73:1, which is d12w (1.734) not d12 (1.164) |
| 5k | `PANEL_e2n_...T47D` + `PANEL_e1n_...A549` | the `.ai` says "max drug log2FC"; e1/e2 say "maximum drug log2FC in the parent screen" |

An automated text-similarity pass over the candidates agreed on 5d and on 5k's
A549 half and disagreed elsewhere -- it put 5e on `PANEL_d12w` and 5j on
`PANEL_d1s`, and it had no way to see that 5a, 5b, 5c and 5k are each built
from two source panels. The mapping above rests on the per-panel evidence in
the table, not on that pass.

Every redrawn part is text-identical to the published panel file it came from,
checked against extracts vendored in `tests/reference/figure5/`.

**What is still not in the repository:** figure 5a's experimental schematic --
the cells, the arrows, the "+7 days" -- was drawn by hand in Illustrator. Only
the two donuts beside it come from code.

---

## 12. Lollipop titles: where the gRNA denominator comes from

**Status: reproduced as published; whether to change it is open.**

Every lollipop title ends "ABE: n/N gRNAs, CBE: n/N gRNAs". In
`lollipop_plots_canonical.ipynb` the numerator is computed but **N is a literal
typed into the f-string** (cells 8, 9, 10 -- one per subpool):

| subpool | published N | targeting guides | whole pool |
|---|---:|---:|---:|
| subpool1 (F1-R1) | 2,439 | 2,439 | 2,730 |
| subpool2 (F2-R2) | 4,533 | 4,533 | 4,868 |
| subpool3 (F3-R3) | **1,655** | **1,505** | **1,655** |

subpool1 and subpool2 print the targeting-guide count. subpool3 prints the
whole pool, which also counts its 79 non-targeting and 71 intron controls --
guides that cannot appear in the numerator, since that counts targeting guides
with sensor information. So figure 12a-f and the other subpool-3 lollipops
quote a denominator 150 larger than the set the numerator is drawn from.

The numerators are unaffected, and so is every point, hit and label on every
panel. This is the title text only.

`cdkbe.source_data.screen_figures.PUBLISHED_N_LIBRARY` holds the three literals
so the panels reproduce what was published. Passing `n_library=None` computes
the consistent value instead (1,505 for subpool3); switching the default would
correct the figures rather than reproduce them.

---

## 13. Supplementary 14 and 15c: panels held as tables but not as code

**Status: ten of the eleven supplementary-14 panels are now redrawn; 15c is
reconstructed. Only 14k remains undrawn.**

### Supplementary 14 (eleven panels)

These come from five composite PDFs in the BEAT-seq delivery's
`plotting/figures/` tree. That tree shipped **outputs only** -- no plotting
code came with the bundle -- so the panels are redrawn here from the tables the
bundle saved, against the numbers the published PDFs print on their own faces.

**The sub-panel letters are settled, and this document previously said they
were not.** The composites carry their own letters, printed in the panel
titles. `experiment_summary.pdf` is lettered a-d, and those four are
supplementary 14a-14d in that order; `de_summary.pdf` is lettered a-d too, and
those are 14g-14j:

| panel | the composite's own title |
|---|---|
| 14a | a  Cells recovered (138,363 post-QC) |
| 14b | b  Guide assignment (101,805 singlets) |
| 14c | c  Library coverage (recovered / designed) |
| 14d | d  Representation depth |
| 14e | A549/ABE (n=91) ... -- `cells_per_guide_per_arm.pdf` |
| 14f | Designed consequence -- `de_by_variant_class.pdf` |
| 14g | a  DE calls (8/288 active) |
| 14h | b  Not power-limited (rho +0.06, -0.03, -0.05, +0.05) |
| 14i | c  Same guide, both lines (rho = 0.26) |
| 14j | d  Active guide-arms |
| 14k | All arms -- pseudobulk transcriptional correlation |

Nothing was guessed to get there, and nothing had to be asked. Every figure in
every one of those titles regenerates from the delivered tables:

| check | panel prints | the port gives |
|---|---|---|
| cells passing QC, four arms | 93 / 95 / 94 / 93 % | the same |
| post-QC cells in 14a's title | 138,363 | the same |
| singlets, four arms | 73 / 70 / 86 / 88 % | the same |
| singlets in 14b's title | 101,805 | the same |
| guides recovered | 91/92, 71/72, 91/92, 71/72 | the same |
| cells per guide, medians | 509 / 537 / 138 / 71 | the same |
| ... means | 500 / 549 / 136 / 70 | the same |
| active by consequence | 4/232, 2/34, 2/22 | the same |
| DE active per arm | 2/82, 0/62, 3/82, 3/62 | the same |
| cells against DE genes, Spearman | +0.06, -0.03, -0.05, +0.05 | the same |
| the two cell lines, Spearman | 0.26 | the same |
| 14j's eight bars, cells and edit rates | 122 cells / 66% ... | the same |

Two tables had to be vendored for this beyond what the figures directory
shipped: `pseudobulk_guide_annotation.csv` (324 guide-arms with their DE call,
variant class, cells and edit rate) and `guide_de_summary.csv`, both from the
delivery's `05_primary_analysis/outputs/tables/de/`.

14e's histogram bins are the one thing no PDF prints, and they were recovered
rather than chosen: **29 logarithmic bins across the pooled range reproduce all
four arms' bar heights exactly**, and 28 and 30 both miss. `tests/
test_supplementary.py` pins all four arms' bar sequences.

**14k is still not redrawn**, and this one is a real gap rather than an
unanswered question. The published clustermap is over "2,607 shared genes"; the
delivered pseudobulk matrix has 324 guide-arms by **4,779** genes with no
missing values, so no "shared across arms" rule cuts it to 2,607, and the
gene filter is in none of the code here. Guessing a gene set would change
every correlation on the page, so `beatseq_qc.source_14k()` returns the
annotation tracks and the null bands and stops there.

**What would settle 14k:** the gene selection used for the pooled clustermap,
or the plotting code if it exists outside the delivered bundle.

### Supplementary 15c

**Status: RESOLVED by reconstruction, and pinned against the published panel.**

The published `KLIFS_resistance_MSA.pdf` is a **two-panel** figure: the CDK
base-editing screen on the 85 KLIFS positions on top, the ten OncoKB clinical
kinases on the same positions underneath, a shared consensus track between them,
two legends and a footnote.

`CDK_resistance_patterns.ipynb` cell 14 -- the only cell that writes that
filename -- draws the **CDK block alone**, under a different title. The
two-panel version is not in the notebook, not in
`canonical_resistance_figures.py`, and not in the IPython history. That code
still appears to be lost, the same way the `filtered_top_only` tables' code is.

Everything the lower block needs was, however, already in this repository:

* the ten kinases' 85-position pocket sequences, from each one's own
  `KLIFS2UniProtIdx` map in `data/reference/KLIFS/oncokb_klifs/`;
* the per-position drug counts, from `klifs_msa.oncokb_by_position()`, which
  reproduced the published footnote (30 of 104 mutations outside the numbered
  positions) before anything was drawn.

So `klifs_msa.supp_15c` now draws the whole panel. Both blocks go through
pymsaviz's own renderer onto axes of one figure, so they cannot drift apart, and
`two_panel=False` still gives cell 14's single block verbatim.

**This is a reconstruction of a lost figure, not a recovery of lost code**, and
it is worth checking as such. The evidence that it is right: the rebuilt PDF and
the published one contain **the same 2,189 characters**, with no character in
one and not the other -- every residue of all seventeen rows, both position
axes, the consensus sequence, both legends, the title and the footnote. Only
`pdftotext`'s reading order differs, because the row labels and the subdomain
annotations sit at the same height. `tests/test_supp15c.py` pins that for both
the `core7` and `all9` sets, and pins the drug counts behind the shading, which
text cannot cover.

**What is still worth checking:** the colour scales and the diamond markers were
read off the published figure rather than from code. The diamonds mark positions
carrying resistance in *both* blocks; if they meant something narrower in the
original, say so and the rule is one line to change.

---

## 14. The SAV tables use two different filters, depending on the screen

**Status: both reproduce; nothing outstanding.**

A single amino-acid variant table splits compound sensor outcomes (`R46W_K48E`)
into their parts and shares the reads out to each. Whether reads carrying a
substitution the editor *cannot* make are dropped first differs between screens:

| table | filter | rows (CBE subpool 1) |
|---|---|---:|
| shipped `SAVs_canonical/*.csv` | none | 16,026 |
| `logo_plots_KB_compound_mut.ipynb` | editor's own edits only | fewer |

So the main-figure panels (2b, 2i, 4n) and supplementary 8c are built on the
unfiltered convention, and supplementary 9b on the filtered one -- because that
is what each published panel used.

`cdkbe.editing.savs_from_full_canonical` takes `own_edits_only` for this, and
with the default (`False`) reproduces every shipped `SAVs_canonical` table
exactly. Both conventions are pinned in `tests/test_supplementary.py`.

This is worth knowing when building a new panel: picking the wrong one changes
which variant a guide is attributed to, silently.

---

## 15. Supplementary 2c: three of six LC50s do not reproduce

**Status: open, but nothing else depends on it.**

The three dose-response grids each use a **different background convention**,
because `FACS_dose_response_2024_8_22.ipynb` redefines the fit between them:

| panel | plate | background taken from |
|---|---|---|
| 2a (cell 158) | `2024_09_27.csv` | nothing -- `no_cell_avg = 0` |
| 2b, 2c (cell 162) | `2025_4_9_csv.csv`, `2025_9_9_CDK246.csv` | wells whose **Cell Line** is "NO CELLS" |
| 3d, 9d-g (cell 125) | `2025_8_31_CDK12_13.csv`, `2025_8_10.csv` | wells whose **Drug** is "No Cells" |

With those conventions applied, 2a reproduces all six published LC50s exactly
and 2b all three. 2c reproduces three of six:

| drug | published | reproduced |
|---|---:|---:|
| Ribociclib | 2,754 | 2,754 |
| Atirmociclib | 5,542 | 5,542 |
| INX-315 | 1,687 | 1,687 |
| Palbociclib | 6,298 | **4,951** |
| Abemaciclib | 3,354 | **2,845** |
| Tagtociclib | 3,653 | **2,972** |

The three that disagree share the plate, the dose grid (8 doses x 3 wells) and
the fit with the three that agree, and `2025_9_9_CDK246.csv` has no exclusion
column. No MTT plate in the repository reproduces the published values for
those three. The most likely explanation is that the published panel was drawn
from an earlier version of that file, the same situation as item 2.

One detail that has to be right: the background rows are matched as
`"NO CELLS"` in the original and spelled `"No cells"` in that file, so the
match is made case-insensitively here. Without that the background is `NaN`
and so is every viability on the panel.

**What would settle it:** the version of `2025_9_9_CDK246.csv` that was live
when the figure was made, if a copy survives.

---

## 16. Two supplementary-5 panels are keyed in upper case

`supp_figures.SUPP_SPEC` keys supplementary 5's two lollipops as `S5A` and
`S5B`. Every other key in both registries is lower case, and so is every panel
letter in the manuscript's own figures, so `scripts/render_figures.py`
lower-cases the letter when it names a file: these come out as `supp_fig5a.pdf`
and `supp_fig5b.pdf`, beside `supp_fig5e.pdf` from the same figure. The registry
keys are left as they are, because tests pin them.

**What would settle it:** confirm that supplementary 5's panels are lettered a,
b, ... e in the manuscript and not A, B, ... e. If they really are upper case,
the normalisation in `render_figures.stem()` should come out and the other keys
in that figure should change instead. Nothing else depends on it.

---

## 17. Compression makes a literal `.csv` path fail

`data/processed/rnaseq/T47D_RNAseq/cdk4_program_analysis/pub_fig_square.py` is
the script that drew figures 5f-5i. It reads its inputs with `pd.read_csv` on
literal `.csv` paths, and those tables are gzipped in this repository, so it
cannot run here as written.

It is kept unchanged as the record of what was run. The version that runs is
`cdkbe.source_data.rnaseq_t47d_figures`, which is the same code reading through
`cdkbe.io` and with the three-layout loop behind `__main__`.
`tests/test_figure5.py` pins it against the shipped PDF, so the two cannot drift
apart silently.

The same thing had happened in `cdkbe.beatseq.fig_panels`, where
`pseudo_guide_null.csv` is gzipped on disk and `panel_j` raised
`FileNotFoundError` until its table paths were routed through `io.resolve`.
**Read every table through `cdkbe.io`**, which resolves the compression, rather
than through a literal `.csv` path -- see `docs/METHODS.md` section 6.

---

## 18. The manuscript's own mapping names a different panel for figure 5k

`docs/source_data_mapping_main.csv` -- the panel-to-source record kept beside the
manuscript, and vendored here -- agrees with `cdkbe.source_data.figure5.SPEC` on
ten of figure 5's eleven panels, including every two-part one. That is an
independent confirmation of the mapping, written without reference to any of
this repository's code.

On 5k it names `PANEL_e1_severity_vs_drug_lfc_A549` for the A549 half, paired
with `PANEL_e2n_severity_vs_drug_lfc_T47D` for the T47D half. The published
panel is `PANEL_e1n`:

* `PANEL_e1` labels its axes "maximum drug log2FC in the parent screen" and
  "severity (SD of the size-matched null)". **Neither string appears anywhere in
  `figure5_v2.ai`.**
* `PANEL_e1n` labels them "max drug log2FC" and "severity". "max drug log2 FC"
  appears in the figure **twice** -- once for each half of 5k.
* The mapping's pairing would also put two different axis-label conventions on
  the two halves of one panel.

So this is a one-character slip in the mapping file (`e1` for `e1n`), not an
error here. `tests/test_figure5.py` pins both the agreement and this single
exception, so if the mapping file is corrected the test will say so.

**Nothing needs to change in the code.** The mapping file is worth correcting if
it is going to be submitted with the manuscript.

---

## 19. The supplementary figure numbering in the manuscript draft is not the current one

`Gould_et_al_CDK_manuscript.rtfd` (dated 15 August) carries legends for five
supplementary figures. `source_data_mapping_supp.xlsx` (November) covers
sixteen, and they do not line up -- the draft's Supplementary Figure 3 is
"additional information about subpool 1", whose panels are spread across the
workbook's S3, S5, S8 and S9.

The workbook is the later document and is what this repository follows. The
draft is noted here only so that nobody re-maps panels from it by mistake.

**There is no supplementary figure artwork in any of the three repositories** --
only the six main-figure `.ai` files. That is why the panels in item 13, and the
"Missing" rows in the workbook, cannot be identified from the files on disk.

---

## 20. Figure 4f's numbers were literals in the notebook

`analysis.ipynb` cell 137 draws figure 4f from two Python lists typed into the
cell -- eleven variants with their pre- and post-selection allele frequencies.
No upstream table produces them, here or in the original tree, and nothing
recomputes them from the Sanger traces under
`data/reference/validation_guides/2025_1_15_CDK_validation_guide_seq/`.

They are transcribed verbatim into
`data/processed/validation/allele_frequency_4f.csv`, which matters because
the original notebooks are not shipped and those 22 numbers exist
nowhere else.
`cdkbe.source_data.allele_frequency` is therefore a transcription, not a
derivation. If the quantification is ever redone from the traces, that table is
what the result has to be checked against.

---

## 21. The CDK9 Boltz-2 analysis has two ATP-pocket definitions and an arbitrary axis sign

Both are in `cdk9_ifp.ipynb`, and both are easy to get wrong when porting.

**Two pockets.** Cell 15 builds `CDK9_ATP_POCKET` from a labelled `CDK9_REGIONS`
map: 32 residues. Cell 19 then defines its own flat `ATP_POCKET`: 37 residues,
adding 47, 49, 65, 67 and 150, which no region names. **Every PCA in the
notebook uses cell 19's set**, and so does the README's methods description, so
that is the one `cdkbe.boltz.ifp.ATP_POCKET` carries. Using the other gives 53
pocket features instead of the published 58.

**Arbitrary signs.** `np.linalg.svd` may return a component or its negative;
both describe the same fit. On this machine PC2 came back mirrored against the
published panels, which flips supplementary 6h top-to-bottom and reverses 6j's
ordering and colours. `ifp.SIGN_ANCHORS` pins each component by the sign of one
named loading in the published figure -- `167:charged_polar` for PC1,
`154:polar_sc` for PC2 -- so the panels reproduce rather than mirror.

Nothing here is outstanding; it is recorded because a future re-run on another
machine would otherwise look like a real disagreement.

---

## 22. The leader-line panels: what was fixed, and one difference left

**Status: supplementary 15e now reproduces the published panel exactly.
Figure 6 reproduces it and adds one site.**

`scripts/chimerax/label_consensus.py` draws the two panels whose labels sit in
the frame's gutters joined to their residue by a leader line: figure 6 (the
consensus panel) and supplementary 15e (the peaks panel). Three corrections
were needed before it would run and reproduce them, and none of the three
announced itself:

1. **`p["h"] = LINE_H`, and `LINE_H` is defined nowhere.** A `NameError` that
   stopped the script dead, so neither panel could be regenerated at all. The
   file defines `LINE_H_PER_PT` next to `CHAR_W_PER_PT`, and the line above
   sets `char_w = CHAR_W_PER_PT * font`, so the height half of that pair is
   what belongs there. The published panels' own geometry confirms it: the
   gaps between the peak panel's blocks are 70.9 and 102.7 px for a 4-line and
   a 6-line block, and `(n_lines * LINE_H_PER_PT + BLOCK_GAP_PER_PT) * font`
   gives 70.84 and 102.68.

2. **The SVG writer always emitted a square viewBox.** Label positions are
   fractions of the ChimeraX frame, and every published panel was composed in a
   1400x796 window, so writing them into a 1400x1400 box stretched every
   vertical gap by 1.76 and distorted the structure image underneath. This was
   the "looser spacing" this document used to record as unexplained. `WINDOW`
   is now set to the frame the published panels were made in rather than left
   to whatever the display allows, because ChimeraX honours a taller request on
   a taller screen and that moves every label.

3. **`full` mode was missing from `main()`.** Supplementary 15e's labels are
   multi-line blocks -- a header, then one line per contributing gene, wrapped
   at `WRAP` -- and `variant_block`, `MAX_W_FULL`, `BLOCK_GAP_PER_PT` and
   `spread`'s `key_h` were all still in the file with nothing calling them.
   The mode is restored, and the panel's `zoom` (0.66, which pulls the
   structure back to clear the gutters for the wide blocks) was recovered by
   fitting the leader-line endpoints to the published ones.

With all three fixed, supplementary 15e comes back **label for label identical
to the published SVG**, with every coordinate inside half a pixel of a
1400-wide frame and the same seven ramp colours. The file is even the same
size. `tests/test_leader_panels.py` pins it against the published artwork,
vendored under `tests/reference/leader_panels/`.

**Figure 6 is right except for one site.** The regenerated panel carries all 32
of the published panel's labels and one more: `CDK7 F91L`, a `strong_singleton`
at the gatekeeper (KLIFS GK:45, CDK2 F80) in the current
`consensus_sites.csv`. It is the weakest of the eight strong singletons by FDR
(4.7e-03, on one guide), so the likeliest explanation is that the published
figure predates the site clearing the threshold. The extra label pushes the
rest of the stack, so the other 32 sit in different places; their text is
unchanged.

**What would settle it:** whether figure 6 should show the CDK7 F91L
gatekeeper site. If it should not, the tier rule that excludes it is one line.

---

## 23. The competition-assay heatmaps reproduce, except for one panel's colour scale

Seven panels arrived in the November mapping workbook with sources in the lab
Dropbox: figure 2f and supplementary 7c, 7e, 10f, 10g, 11b and 11i. Five of
them are fold-change heatmaps, and for all five the **published SVG was
delivered beside the data**, written by matplotlib with `svg.fonttype: none`,
so every cell's fill is a hex string in the file and every annotation is real
text. That makes the comparison exact rather than visual, and the SVGs are
vendored under `tests/reference/competition_assays/`.

Nothing about the colour scale was written down, so it was read back out of
those files. The colourbar runs 0 to 5.5 with evenly spaced integer ticks -- a
plain linear norm, not a two-slope one -- and the ramp it draws is `RdBu_r`
truncated to 0.389-1.0, which is exactly what seaborn's `center=1` does to that
colormap at that vmin and vmax. So the family is
`cmap="RdBu_r", vmin=0, vmax=5.5, center=1`.

With that, **four of the five come back byte-identical** in every cell colour
and every character of text: figure 2f, supplementary 10f, 10g and 11i.

**Supplementary 11b is the exception, and only in its four fill colours.** Its
published SVG is truncated to 0.3325-1.0 instead, which at this vmin and vmax
means a `center` near 1.38 rather than 1.0 -- close enough to that panel's own
non-targeting control mean (1.3700) to look deliberate, but nothing says so.
The numbers, the SEMs, the stars, the labels and the layout are identical; the
cells are a shade warmer in the published panel than they are here.
`tests/test_competition_assays.py` pins the difference so it cannot be mistaken
for drift later.

**What would settle it:** whether supplementary 11b should share the family's
colour scale (centred at 1, as it is drawn here) or keep the published panel's
warmer one.

Two smaller things settled themselves on the way:

* **The CDK13 files repeat the CDK12 files' control arm.** `Cdk13 FC.txt` and
  `Cdk12 FC.txt` carry the same three NT replicates, as do the two ZLY files
  for their ABE NT arm, which is why supplementary 10f and 10g draw one
  non-targeting row rather than two. Pinned by a test.
* **Supplementary 10g drops a row its file carries.** The ZLY file has both an
  ABE and a CBE non-targeting arm; the published panel shows the ABE arm, as a
  row labelled simply `NT`, and measures every star against it. Keeping the CBE
  arm would add a row the manuscript does not have.

The significance convention also differs between the two orientations, and both
are reproduced as published: the CDK12-family panels mark every variant cell in
every treatment, including `ns`, while figure 2f marks only the treated row and
marks neither non-targeting arm.

## 24. Supplementary 7c and 7e: numbers here, published figure not

The two mouse B-ALL competition assays are held as a summary of three
replicates each -- `mCdk9 Summary Plot .txt` and `mCdk7 Comp Assay Summary.txt`
-- without the published panel or the code that drew it. The
numbers are exact and both arms separate cleanly (p = 0.0075 and p = 0.0085,
two-sided unpaired t-test), but **the panel style here is this repository's own
convention**, the bar-and-points style the other in vivo panels use, not a
reproduction of something published.

The FlowJo exports those three replicates were computed from exist alongside
the summaries, each an `.xls` with the gating tree and its fold-change columns.
They are not included here: the per-replicate fold changes in the exports do
not match the summary's three values, so a normalisation step sits between
them that this repository cannot reconstruct. The summary values are what the
panels use.

**What would settle it:** the published supplementary 7c and 7e panels, or
whoever drew them saying what they look like.

## 25. The LNCaP and PC3 packages carry two corrections of their own

Supplementary 11c and 11D-11K arrived as self-describing packages -- the values
that were plotted, the script that plotted them, and the published image, one
folder per panel. `src/cdkbe/source_data/lncap_pc3.py` redraws each panel from
the delivered values, and every DE count its READMEs print comes back exactly
(780, 39, 286 and 1 differentially expressed genes for the four LNCaP
contrasts; 141/52, 764/429 and 696/363 up and down for the three PC3 ones).

Two corrections are recorded by the contributing authors and are carried with
the data rather than re-made here. Both matter for reading the figures:

1. **The PC3 contrasts were refitted without `lfc_shrink`.** The original
   `pc3_deseq2_gsea.py` called `lfc_shrink(coeff='condition[T.<test>]')`, which
   shrinks against the *design* reference (`G731K_Control`), not the reference
   named in the contrast. `log2FoldChange` in `NT_vs_I733V_Control` and
   `NT_vs_G822E_Control` was therefore effectively NT versus G731K. `padj` and
   `stat` were never affected, and the significance calls agree for more than
   99.9% of genes. The shipped tables are the unshrunk MLE from each Wald
   contrast.
2. **The logo-plot pileup was corrected for reverse-strand reads.** Roughly
   40-45% of the Nanopore reads align to the reverse strand, and `mappy`
   reports query positions against the reverse complement, so the CIGAR walk
   has to reverse-complement the read (sequence *and* quality string) first.
   Without it, about half the base calls come out miscalled -- including in
   the non-targeting controls, which should read 100% wild type.

One artifact in the logo plots is deliberately left alone: about 5.4% of reads
at residue 822 in the G822E/K sample call a stop codon, which survives a
Q20-to-Q30 filter, is absent from the non-targeting control on the identical
sequence context (0.04%), and is not explained by a nearby indel. The delivered
analysis leaves those reads in their own bucket rather than folding them back
into wild type or G822E, on the grounds that a genuine CBE bystander
transversion would be signal and not noise. A test pins that they are still
there.

Finally, supplementary 11D's PCA is **not** from pyDESeq2: it is nf-core/rnaseq
3.25.0's own DESeq2 QC step, on variance-stabilised counts over the 500 most
variable genes, which is why its axis labels say 31% and 15% and no code here
recomputes them.

## 26. Sixteen T47D and MCF7 panels are held as published figures

Figure 4g, 4i-4m and supplementary 13a-13j were contributed as finished
figures: fourteen `.svg` and two `.png`, under
`data/raw/collaborator_panels/T47D_MCF7/`. Six of them are redrawn here --
items 27 and 31 -- and `docs/COVERAGE.md` marks the other ten **published
image**.

For those ten, this repository has the figure and not the table. Some of the
SVGs carry their labels as real text -- enough to say that
`PCA.svg` is over PC1 52% and PC2 25%, that `volcano_variant_vs_control.svg`
contrasts NT against CDK2 S261G and three more variants, that the stacked
`T47D_cell_cycle_*` panels are percent S phase with their own significance
stars, and that `GSEA_variant_vs_control_palbo_dotplot.svg` is five Hallmark
sets across four variants -- but the plotted values are geometry, not text.
Where the plotted values are geometry against a labelled linear axis they can
be read back exactly, which is what item 31 does for two of them. Where they
are encoded some other way -- a heatmap's fill, a dot plot's area -- they
cannot, and the panel stays the published figure.

One classification rule is worth recording, because the distinction is easy
to lose: `CDK_sanger_aa_logo_plots_all.svg` (supplementary 13b) was counted as
"no data" while the coverage script looked for the word *Sanger* anywhere in
the source cell. A Sanger **trace** has no underlying table; a logo plot drawn
**from** Sanger reads does. The rule now reads that word only when the cell
begins with "None".

## 27. Four of the sixteen came with their measurements

Figure 4k, 4l and supplementary 13c, 13d are the T47D cell-cycle panels, and
the experiment behind them is complete: the FlowJo export the panels plot (54
rows -- six cell lines by three treatments by three replicates), the notebook
that drew them, the style module it imports, the FlowJo workspace with the
gating, and the 63 `.fcs` files. All of it is vendored under
`data/raw/collaborator_panels/T47D_cell_cycle/`. The four published SVGs are
byte-identical copies of what that notebook produced, so the match is not a
guess.

The export, notebook, style module and workspace are vendored at
`data/raw/collaborator_panels/T47D_cell_cycle/`; the `.fcs` files (100 MB) are
listed in `data/external/NOT_INCLUDED.csv`. All four panels now redraw
character for character and fill for fill against the published SVGs.

Three things in that package are easy to get wrong and are pinned by tests:

1. **The four gates do not sum to 100%.** G0, G1, S and G2/M are siblings on
   one parent population, so small gaps at the gate boundaries are normal --
   the totals here run 91 to 99%. The stacked panels rescale each condition to
   100; the %S-phase panels use the gate frequency as measured. Doing it the
   other way round changes every bar.
2. **The Bonferroni family is the whole panel, not the treatment.** Welch's
   t-test of each variant against its own non-targeting arm, corrected across
   all four comparisons a panel draws. CDK2 S261G under DMSO has a raw
   p of 0.015: times two it would carry a star, times four it is `ns`, and the
   published panel says `ns`.
3. **The two CDK6 lines are relabelled.** `CDK6_G53S` is printed as
   `CDK6 G53S_M54I` and `CDK6_R46W` as `CDK6 R44C_R46W`, because each line
   carries a second edit the sample name does not record.

The other panels in that set are covered by item 26.

## 29. Figure 6, and the cross-kinase table behind its labels

Figure 6 rendered but never reached `figures/manuscript/`, because
`structures.py` deliberately skips it -- it is a whole figure rather than a
lettered panel, so it was in no registry and `render_figures.py` never copied
it. It is now registered in `klifs_msa.py` as main-figure key `6`, and
writes `main_fig6.svg`, `main_fig6.png` and the background raster the SVG
draws its labels over.

Its mapping row carries a note: *"Additional text added to indicate the
matching aligned kinase positions in other kinases (besides CDKs); these
should be sourced as well; think there's an existing table"*. There is, and it
is already in this repository:

| table | what it holds |
|---|---|
| `resistance_patterns/core7/consensus_sites_homologous_resistance.csv` | one row per consensus site (33), with `homologous_clinical` -- the clinical variants reported at the same KLIFS position in other kinases, with their source |
| `resistance_patterns/core7/matching_sites_cdk_vs_clinical.csv` | the same join by KLIFS position (11 rows): which CDKs, which other kinases, which clinical variants, OncoKB or CIViC |
| `resistance_patterns/core7/consensus_sites.csv` | the 124 alignment columns the figure selects from, with tier, segment and variants |

`figure_6` returns the first of those as its source-data table, so the text on
the figure and the table behind it come from one place. The gatekeeper
question in item 22 is still open.

## 30. Figure 2c: settled, and rebuilt

**Closed.** The mapping recorded 2c only as an older ChimeraX visualisation,
and two candidates in the source trees fitted that description:
`CDK9_cyclin_KLIFS_ribbon.png`, which is reproducible but was drawn in Python
rather than ChimeraX, and three ChimeraX session files, which are binary and
cannot be read without the application. The panel is neither: it is the CDK9
render **with only the variants figure 2b plots**.

That is now built rather than hunted for. `scripts/chimerax/build_cdk9_fig2b_hits.py`
calls figure 2b's own selection -- `screen_figures._cdk9_top_hits` -- and splits
each guide's HGVSp call into residues, which gives five positions from its ten
guides: D104, C106, E107, L156 and I157. `generate_cdk9_fig2b.py` then writes
the same CDK9/cyclin T1 panel every other structure figure uses (8K5R, chain A
and B, KB-0742 as sticks, the shared camera) marking those five and no others,
against the sixteen that figure 1f's `CDK9.png` marks. Four images come out:
plain and labelled, each full-frame and zoomed on the cluster.

Because the residue list is computed from 2b rather than written down, the two
panels cannot drift apart: change 2b's cutoff and 2c follows on the next run.

**What this is not.** It is not a recovery of the original image. The old
session is gone, so the framing, the colours and the label placement are this
repository's conventions, not the published panel's. If the published 2c
differs in any way that matters, say how and it is a one-line change.

## 31. Figure 4g and supplementary 13a: read back out of the published figure

This repository has both panels as published SVGs --
`cdk_variant_T47D_plot_t0_tf.svg` and `cdk_variant_MCF7_plot_t0_tf.svg` -- and
not as tables.

The values were read back out of them. Both are matplotlib line
plots against a labelled linear y axis, so the axis inverts exactly and each
line's two vertices read straight back out; `cdkbe.svg_values` does the
inversion and `scripts/recover_selection_tables.py` writes the tables.

**The inversion is validated, not assumed.** Before writing anything, the
script runs the same inversion over the two T47D cell-cycle panels, whose
published SVGs *and* underlying FlowJo export are both in this repository, and
stops unless the recovered bar heights equal the measured gate frequencies.
They do, to 2e-7 percentage points -- which is the six-decimal precision of the
SVG coordinates, not a limit of the method.

What comes back is fifteen whole-number percentages out of sixteen, the same
form figure 4f's hand-quantified table takes, plus one real fraction: MCF7
CDK2 S261G is 0% before selection and 4.264871% after.

**What this is and is not.** It is an exact reading of what was plotted, which
is not the same thing as the primary quantification of the Sanger traces: a
quantification done from the traces themselves would replace these tables. The
published figures stay in the repository either way.

## 32. Figure 4m cannot be read back the same way, and why

4m is `marker_gene_panel_DMSO_palbo_stacked.svg`: six T47D lines by seventeen
marker genes by two treatments, coloured by z-score of log2 TPM **within arm**.
Like 4g, this repository has it as a published figure rather than a table.

Unlike 4g it is not read back, and the difference is not effort:

* the values are in the cell **fills**, not in geometry. A fill is quantised to
  eight bits per channel and pushed through a colormap, so what comes back is
  an approximation, and `cdkbe.svg_values` deliberately reads no colour.
* even a perfect reading would give z-scores computed within each arm, not the
  TPM matrix they came from. The panel would redraw; nothing else could be done
  with the numbers, and no other analysis could be checked against them.

The experiment behind it is six lines (ABE NT, CBE NT, CDK2 S261G, CDK4 I51V,
CDK6 G53S_M54I, CDK6 R44C_R46W) in DMSO and palbociclib. The T47D RNA-seq this
repository does ship is a different experiment: four CDK4 genotypes, vehicle
only, which is figure 5.

**What would settle it:** the TPM or DESeq2 table behind it. The same table
would close S13e, S13f and S13g, which are the PCA, volcano and GSEA of that
experiment.

## 33. The Sanger traces, and one label that disagrees with its primer

All 288 `.ab1` files this repository ships are now indexed in one table,
`data/raw/sanger/INDEX.csv`, written by `scripts/make_sanger_index.py`: the
sample label and sequencing primer the vendor recorded, the gene, variant,
treatment, selection state and read direction parsed out of that label, and
the run date, instrument, well, order, dye, read length and median Phred read
out of each ABIF header.

The filename is the only record of what a sample is, so the parse is checked
where it can be. Each guide id read out of a filename has to exist in the
screening library: 55 of 56 do.

The one that does not is **`CDK13_2232_selected`** (forward and reverse), in
`SG_2025_8_20_CDK9_12_13_AJ_2nd_selection_30-1229691625_ab1`. There is no
`gRNA_CDK13_targ_2232`. There is a `gRNA_CDK12_targ_2232`, the same run holds
`CDK12_2232_sel_1`, and both files were read with the **CDK12** primer
(`F1_CDK12`, `R1_CDK12`), so the gene in the label is the part that does not
fit.

Nothing is renamed on that inference. The index carries the parse as written
and flags it in `guide_in_library`, and `tests/test_sanger_index.py` pins it so
it cannot quietly disappear.

Two smaller facts recorded in the same place: the cloning QC called a
protospacer for 30 of the 45 colony traces, so fifteen carry no guide id; and
the four Plasmidsaurus deliveries (`*_results`) have no run date in their ABIF
headers, which the capillary vendors' files all do.

## 34. Keeping the source-data workbooks in step with the registry

`source_data_tables/*.xlsx` is what a reader gets instead of this repository,
so it has to follow the code rather than lag it. Three ways it had drifted are
recorded here because each failed silently, and the checks that now catch them
are the useful part.

**The builder stopped part-way through.** `scripts/build_source_data.py` raised
`AttributeError` on `beatseq_qc.PANEL_COMPOSITE`, whose values had become
`(composite, letter)` tuples after the script was written, so the workbooks on
disk were whatever the last complete run produced. `Figure_4.xlsx` had seven
sheets where thirteen panels have source data -- no 4f, no 4c-4e, no 4k/4l --
and still used the pre-relabelling letters. There was no `Figure_6.xlsx`, and
no workbook at all for supplementary 6, 11 or 13.

**Every matrix sheet had lost its row labels.** Sheets were written with
`to_excel(index=False)`, which is right for a tidy frame whose index is just
position and wrong for a matrix: the LFC correlation clustermaps (S4a-S4c), the
Boltz-2 matrices (S6i, S6j) and every other frame carrying labels in the index
shipped as blocks of numbers with nothing naming the rows. The index is now
written whenever it is not an anonymous integer index, and named `label` if it
had no name.

**Three sheets were two tables side by side.** Each cell-count sheet behind
S1b-S1d carries a stale pivoted block to the right of the real columns, under
blank headers, which came through as `Unnamed: 16` onwards. `cpd._sheet` drops
it; the workbook itself is untouched.

**And two sheets did not hold what their panel plots.** Completeness is not
just one sheet per panel. `scripts/audit_source_data.py` draws each panel,
takes every word it prints and reports the ones that appear nowhere in its
sheets; most of what it reports is wording rather than data -- axis titles,
legend headings, counts computed from the sheet -- but two were real:

* **figure 4k and 4l** shipped only their significance table. The bars are
  gate frequencies, which were not in the workbook at all, so the sheet
  reproduced the brackets and nothing else. `t47d_cell_cycle.source` now
  returns the per-replicate frequencies, their means and the significance
  table, and S13c and S13d additionally return the rescaled values the stacked
  panels actually draw. Each frame carries the published cell-line label
  (`CDK6 G53S_M54I`) beside the export's own spelling (`CDK6_G53S`).
* **figure 1f** prints each gene's canonical length beside its track and scales
  the tracks by it, and that length was in neither the sheet nor anything the
  sheet could be joined to. `library.figure_1f` now carries `protein_length`.

The set now holds 148 panel sheets covering every panel in the registry, each
with its row labels and an INDEX line saying what it is -- taken from the
drawer's own docstring where the hand-written list had no entry, so it cannot
go stale against the code. `tests/test_source_data_tables.py` checks the files
on disk rather than the builder, so a workbook falling behind fails the suite
instead of shipping. The label audit stays a script rather than a test: its
match is crude on purpose, and most of what it flags is correct.
