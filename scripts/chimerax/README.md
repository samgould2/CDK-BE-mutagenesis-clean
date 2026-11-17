# CDK-family resistance variants in ChimeraX

Structural figures for the top resistance variants in
`screening_data/03_LFC_FDR_tables/{ABE,CBE}_master_LFC_table_filtered_top_only.csv`,
built to match the style of the earlier sessions in `notebooks/10_PyMOL`
(`viz1.cxs`, `viz2_cdk2.cxs`, `viz3.cxs`, `align_catalytic_v2.py`).

## Run it

```bash
python3 build_hits_table.py      # CSV tables  -> hit_residues.csv
python3 fetch_structures.py      # download structures -> structures/
python3 generate_cxc.py          # -> cxc/*.cxc
```

Then open any panel in ChimeraX:

```
open cxc/00_all_cdks_aligned.cxc
```

or render everything at once from the ChimeraX GUI (`save` needs a GPU context,
so this one does not work under `--nogui`):

```
open cxc/render_all.cxc
```

Optional check that PDB numbering still agrees with the screen's codon numbering:

```bash
/Applications/ChimeraX-1.10.app/Contents/bin/ChimeraX --nogui --exit \
    --cmd "runscript $PWD/validate_structures.py"
```

## Structures

One structure per gene, preferring an experimental structure with the cyclin
partner bound. Chain assignments and the catalytic-domain ranges live in
`cdk_config.py`, and the same table is written out next to the figures as
`figures/chimerax/structure_table.csv` and `.md` every time `generate_cxc.py`
runs, so the provenance travels with the images.

| gene  | structure | CDK chain | cyclin (surface) | ligand | notes |
|-------|-----------|-----------|------------------|--------|-------|
| CDK2  | 1QMZ | A | B — cyclin A | ATP + Mg | alignment reference; substrate peptide (chain E) hidden |
| CDK7  | 6XD3 | J | I — cyclin H | THZ1 | MAT1 (chain H) hidden |
| CDK8  | 5CEI | A | B — cyclin C | compound 22 | |
| CDK19 | AlphaFold Q9BWU1 | A | cyclin C, **borrowed from 5CEI** | none | no experimental CDK19 structure exists; see below |
| CDK9  | 8K5R | A | B — cyclin T1 | KB-0742 | |
| CDK12 | 8BUG | B | C — cyclin K | CR8 | DDB1 (chain A) hidden; copy 1 is the one with the ligand |
| CDK13 | 7NXJ | A | B — cyclin K | THZ531 | |
| CDK4  | 2W96 | B | A — cyclin D1 | none | 2.3 Å, apo active site |
| CDK6  | 2F2C | B | A — KSHV V-cyclin | aminopurvalanol | **no human cyclin D–CDK6 structure exists in the PDB**, so this cyclin is viral |

## Conventions

Carried over from the earlier sessions:

* CDK catalytic domain: cartoon in the per-gene color from `align_catalytic_v2.py`
  (CDK4 and CDK6 are new and were given unused colors).
* Everything outside the catalytic domain: same color at 80% transparency,
  matching `cartoon_transparency 0.8` in the PyMOL script.
* Cyclin partner: molecular surface, wheat (`#f5deb3`).
* Inhibitor / nucleotide: sticks, element-colored on gray carbons.
* White background, `lighting soft`, silhouettes width 2.

CDK19 has no experimental structure and therefore no bound cyclin. Because it is
the CDK8 paralog and shares cyclin C, CDK8 (5CEI) is superposed onto the CDK19
model and 5CEI chain B is shown in its place, at **60% transparency** so it is
visually obvious that the cyclin is inferred rather than observed.

New for these figures:

* Resistance hits are **red CPK spheres**. Every variant in the `*_top_only`
  tables is an enriched hit, so a single uniform color is used rather than the
  red/blue pair from `viz3.cxs`, where those two colors meant something else.
* All nine structures are superposed on CDK2 (1QMZ chain A) with `matchmaker`;
  pruned RMSDs are 0.81–1.12 Å. The reference never moves, so the cyclins and
  ligands are carried into a common frame too.
* Every panel then applies the *same* camera matrix (taken from `viz2_cdk2.cxs`)
  and the same pull-back, so panels are directly comparable side by side.

## Two rendering details worth knowing

* **Glycine hits.** ChimeraX suppresses backbone atom display under a cartoon,
  so a glycine hit (no side chain past CA) would show nothing but a recolored
  ribbon. Each hit selection gets `cartoon <spec> suppressBackboneDisplay false`
  so glycines draw as spheres like every other hit.
* **Hydrogens.** 8BUG (CDK12) carries riding hydrogens while the crystal
  structures do not, which made the CDK12 spheres look bloated relative to the
  rest. Every script ends with `hide @@element=H target ab`.

## Catalytic-domain ranges

CDK7/8/9/12/13/19 use the KLIFS ranges from `align_catalytic_v2.py`. CDK2, CDK4
and CDK6 were not in that script; `map_kd_ranges.py` derived theirs from the
superposition. All six known ranges map onto CDK2 8–149, which was then mapped
back out to CDK4 (10–160) and CDK6 (18–167).

## Coverage

`validation_report.txt` records, per gene, how many hit positions carry the
wild-type residue the HGVSp string predicts. **87 of 91 positions match exactly
and none mismatch**, so the screen's codon numbering and the PDB author
numbering agree throughout.

Four variants cannot be drawn because their residue is not in the model:

| variant | reason |
|---------|--------|
| CDK9 C95R  | disordered in 8K5R |
| CDK12 K1068R | outside the ordered region of 8BUG chain B |
| CDK6 D4N, G5N | N-terminus disordered in 2F2C (chain starts at residue 9) |

One more is deliberately hidden: **CDK19 A412V** falls in the AlphaFold model's
disordered C-terminal tail (pLDDT 41). The cartoon is restricted to residues
1–363 (pLDDT ≥ 50), the same trimming `align_catalytic.cxc` applied.

## Files

| file | what it is |
|------|-----------|
| `cdk_config.py` | structures, chains, colors, camera — the single source of truth |
| `build_hits_table.py` | splits multi-edit HGVSp strings into per-residue rows |
| `hit_residues.csv` | tidy table: gene, editor, position, variant, max LFC, FDR |
| `fetch_structures.py` | downloads `structures/` |
| `generate_cxc.py` | writes every `.cxc` |
| `validate_structures.py` | numbering check; writes `validation_report.txt` |
| `map_kd_ranges.py` | derives catalytic ranges for CDK2/4/6 (run once) |
| `cxc/00_all_cdks_aligned.cxc` | all nine superposed, cyclin surfaces on |
| `cxc/01_all_cdks_no_surface.cxc` | same, surfaces off (they stack up when superposed) |
| `cxc/CDK*.cxc` | one panel per gene |
| `cxc/render_all.cxc` | batch render to `figures/chimerax/` |
| `generate_cdk7_sy5609.py` | writes the CDK7 SY-5609 panel (see below) |
| `cdk7_sy5609_hit_residues.csv` | its hit residues, written by `CDK7_SY5609_hits.ipynb` |
| `cxc/CDK7_SY5609.cxc` | CDK7 panel, SY-5609 hits only, MAT1 shown |
| `cxc/CDK7_SY5609_labeled.cxc` | same panel with the hit residues labeled |
| `cxc/render_cdk7_sy5609.cxc` | renders both panels (wide + zoomed) |
| `../../figures/chimerax/structure_table.{csv,md}` | structure/chain provenance, regenerated by `generate_cxc.py` |

## The CDK7 SY-5609 panel

`generate_cxc.py` draws hits from `hit_residues.csv`, which pools every
condition. `generate_cdk7_sy5609.py` builds a second CDK7 panel restricted to
the **SY-5609** screens, for the resistance figure that accompanies
`figures/CDK7_SY5609/CDK7_SY5609_top_hits.pdf`:

```bash
# hit residues come from the notebook, so run it first
#   notebooks/09_analysis_and_viz/CDK7_SY5609_hits.ipynb
python3 generate_cdk7_sy5609.py
```

Then, from the ChimeraX GUI (`save` needs a GPU context):

```
open cxc/render_cdk7_sy5609.cxc
```

It writes two panels: `CDK7_SY5609.cxc` and `CDK7_SY5609_labeled.cxc`, which
adds a `L18`-style label (wild-type residue + number) at each hit. Labels are
bare black text (`bgColor none`) at 2.2 Å — scene units, so they scale with the
structure rather than the window — and are emitted *after* the `label delete` in
the footer, which would otherwise remove them.

With no background panel, a label that lands on top of its own red sphere is
low-contrast (black on dark red); L65, H71 and I55 are the affected ones in the
canonical orientation. Set `LABEL_COLOR = "white"`, or give those labels an
`offset`, if that matters for a given figure.

Two differences from `cxc/CDK7.cxc`:

* **Hits** are the 7 SY-5609 positions (L18F, I55T, L65P, H71R, H83R, V90A,
  D97G) rather than the 9 all-condition ones. All 7 carry the wild-type residue
  the HGVSp string predicts in 6XD3 chain J, and none are disordered.
* **MAT1 is shown**, as a light-gray molecular surface, where the other panels
  hide it. Only its C-terminal 244–308 is ordered in 6XD3, so it appears as a
  small lobe packed against the kinase rather than a full subunit. The canonical
  camera is kept so the panel stays comparable with the rest, which leaves MAT1
  partly behind the CDK7 cartoon.

Everything else — structure, chains, colors, catalytic-domain fade, THZ1 sticks,
lighting, camera, pull-back — is imported from `cdk_config.py`.

## Gotcha

ChimeraX only treats `#` as a comment at the **start** of a line. An inline
`# comment` after a command is parsed as a model spec and aborts the script.
`generate_cxc.py` keeps all comments on their own lines.

## The CDK12 / CDK13 DDB1 + HQ461 panels

Two separate figures placing this screen's CDK12 and CDK13 resistance variants
against DDB1 and the molecular glue:

```bash
python3 generate_cdk12_13_ddb1.py
```

then, from the ChimeraX GUI:

```
open cxc/render_cdk12_13_ddb1.cxc
```

| panel | structure | what is drawn |
|-------|-----------|---------------|
| `CDK12_DDB1` | 8BUG as deposited | DDB1 (A) + CDK12 (B) + cyclin K (C) + HQ461 |
| `CDK13_DDB1` | 7NXJ superposed onto 8BUG | CDK13 + its cyclin K, against 8BUG's DDB1 and HQ461 |

**8BUG is the HQ461 complex.** Its title is *"Structure of DDB1 bound to
HQ461-engaged CDK12-cyclin K"* (X-ray, 3.53 Å) and its ligand `RPW`
(C15H15N5OS2) **is** HQ461 — so the glue in these panels is the same compound as
the screen's HQ461 arm, not a stand-in. Note that `cdk_config.py`'s note for
8BUG calls the ligand CR8; **that note is wrong**. `generate_cdk12_13_ddb1.py`
does not read it.

**The CDK13 panel is a superposition, not an observed complex.** No CDK13–DDB1
structure exists, so CDK13 from 7NXJ is matched onto CDK12 in 8BUG
(`matchmaker`, RMSD 0.776 Å over 284 pruned atom pairs) and 8BUG's DDB1 and
HQ461 are left in place. Everything about where CDK13 sits relative to DDB1 is
inferred from the CDK12 complex. 8BUG's own CDK12 and cyclin K are hidden there
so only one kinase and one cyclin are drawn. 7NXJ's own ligand (`5I1`, THZ531)
is an ATP-site covalent inhibitor rather than a glue, and is left hidden.

### Hit coloring

| color | meaning | CDK12 | CDK13 |
|-------|---------|-------|-------|
| red | resistance hit | 18 positions | 20 positions |
| gold `#ffd700` | HQ461-specific hit | K743, G822 | E803 |

Hits come from `screening_data/03_LFC_FDR_tables/CDK12_13_hits/CDK12_13_resistance_variants.csv`
(written by `notebooks/09_analysis_and_viz/logo_plots_CDK12_13_canonical.ipynb`),
split per residue into `cdk12_13_ddb1_hit_residues.csv`.

A position is gold only when *every* variant seen at it comes from an
HQ461-specific guide. **CDK12 Y742 is the one place this matters**: it carries
Y742C (HQ461-specific, from the `Y742C_K743E` guide) *and* Y742H
(BSJ-4-116-specific), so it is drawn red while its immediate neighbour K743 is
gold. Adjacent residues in different colors is correct here, not a bug.

### Coverage

Every hit position that is present carries the wild-type residue the HGVSp
string predicts — **0 mismatches in either structure**, so the screen's codon
numbering agrees with both PDB entries. Positions outside the crystallised
construct cannot be drawn:

| gene | chain | modelled range | hits not drawn |
|------|-------|----------------|----------------|
| CDK12 | 8BUG B | 715–1041 | L70, S394, K1068, S1394 |
| CDK13 | 7NXJ A | 695–1023 | D659, L1099 |

All three gold positions are modelled.

### Rendering details

* **Camera.** `compute_ddb1_camera.py` derives the orientation from the complex
  itself: camera x = DDB1 centroid → cyclin K centroid (DDB1 left, cyclin K
  right), camera z = normal of the plane through the three centroids, so the
  subunits sit side by side. Looking *along* the DDB1→CDK axis instead puts DDB1
  in front of the kinase and hides it.
* Only the **rotation** is hard-coded. The panels then call `view <spec>`, which
  moves the camera along its own axis without touching the orientation — far more
  robust than computing a distance, which has to agree with ChimeraX's
  field of view and window aspect ratio.
* Framing is on 8BUG chains **A and C** (DDB1 + cyclin K). Both exist in both
  panels, so the two figures come out in an identical camera; spanning the two
  ends of the assembly also centers the composition.
* `view` measures displayed **atoms**, and DDB1's are hidden under its surface,
  so the scripts show them for that one command and hide them again after.
* **DDB1 is grey (`#b8b8b8`) at 70% transparency.** It was white at 90%, which
  made it invisible: **a white surface on a white background reads as a bare
  silhouette outline no matter what its alpha is.** Renders at 90%, 60%, 30% and
  10% transparency come out *byte-identical* — only full opacity differs, because
  that is the point at which the surface's shading survives. So the lever for
  making DDB1 visible is its **color**, not its transparency; tinting it faintly
  (`#dfe5ee`, `#eef1f5`) is not enough either. Grey at 25% gives it a readable
  body while keeping it subordinate to the kinase, and still lets the interface
  hits that sit inside its envelope (CDK13 E803, CDK12 K743/G822) show through.
  Note `color` resets alpha, so the color command must come **before** the
  transparency command — the generator emits them in that order.
* **HQ461** is element-colored sticks on gray carbons at the default radius —
  the same convention the per-kinase panels use for their inhibitors. It is
  small in the wide panels by design; the `_glue_zoom` panels are what show it.
* **The zoom panels frame the glue (`view HQ461`), not the hits.** Framing the
  hits does not work across both genes: CDK12's hits are a tight cluster
  (730–822) but CDK13's span 711–1004, so `view CDK13_hits` takes in most of the
  kinase and ends up *wider* than the panel it is supposed to zoom into. HQ461
  is the same object in both panels, so this crop is comparable between them.
* **`clip off` must come after the last `view`/`zoom`, not before.** `view <spec>`
  sets the near/far clipping planes to that spec's *own* bounds, so zooming onto
  the small HQ461 ligand leaves a slab only ~20 Å deep and slices any residue
  nearer the camera than the ligand. CDK13 **T715** was being cut in half this
  way: the slab started at 128.80 Å while T715's atoms span 126.65–129.41 Å, so
  five of its seven atoms fell in front of the near plane. `clip off` resets the
  slab to the whole scene (29–243 Å), but only for view commands issued before
  it — having it in the preamble did nothing, because the render script re-aimed
  the camera afterwards. Both the panels and the render script now end with it.
  **If you re-aim the camera by hand, re-issue `clip off`.**
* `validate_ddb1_clipping.py` checks every hit residue against the clipping slab
  in both cameras for both genes and writes `ddb1_clipping_report.txt`
  (currently **0 clipped**). Run it after changing any camera or zoom value.
* Each panel comes in an **`_opaque`** variant with the whole kinase cartoon at
  full opacity, i.e. without the catalytic-domain fade
  (`NONCATALYTIC_TRANSPARENCY`) the other panels in this repo apply. Use those
  when the hits outside the catalytic core matter — most of CDK13's do.
* There are **no subunit name labels** on the panels — DDB1, the kinase and
  cyclin K are identified by their colors and by the caption. The `_labeled`
  variants add a `G731`-style label per hit position, and nothing else.

### Files

| file | what it is |
|------|-----------|
| `generate_cdk12_13_ddb1.py` | writes all four panels + the render script |
| `compute_ddb1_camera.py` | derives the shared camera orientation; writes `ddb1_camera.txt` |
| `validate_ddb1_clipping.py` | checks no hit residue is cut by the clipping planes; writes `ddb1_clipping_report.txt` |
| `cdk12_13_ddb1_hit_residues.csv` | per-residue hit table with the red/gold class |
| `cxc/CDK12_DDB1.cxc`, `cxc/CDK13_DDB1.cxc` | the two panels |
| `cxc/CDK1{2,3}_DDB1_labeled.cxc` | same, with residue labels |
| `cxc/CDK1{2,3}_DDB1_opaque[_labeled].cxc` | same, whole kinase opaque (no catalytic-domain fade) |
| `cxc/render_cdk12_13_ddb1.cxc` | batch render to `figures/chimerax/` |

## The CDK2/4/6 active-site panels

Straight-on views down the ATP cleft for the top resistance variants in each of
CDK2, CDK4 and CDK6, from this screen's own analysis rather than the pooled
`hit_residues.csv`:

```bash
python3 build_cdk2_4_6_hits.py                 # top 10 guides/gene -> cdk2_4_6_hit_residues.csv
ChimeraX --nogui --exit --cmd "runscript $PWD/analyze_cdk2_4_6_active_site.py"
ChimeraX --nogui --exit --cmd "runscript $PWD/compute_active_site_camera.py"
python3 generate_cdk2_4_6_active_site.py       # -> cxc/{GENE}_active_site[_labeled].cxc
```

then, from the ChimeraX GUI:

```
open cxc/render_cdk2_4_6_active_site.cxc
```

Six panels at 3000x3000 (square), two per gene: plain and `_labeled`.

Hits are the top 10 guides per gene by `active_median_LFC` from
`screening_data/03_LFC_FDR_tables/CDK2_4_6_hits/CDK2_4_6_resistance_variants.csv`
(written by `notebooks/09_analysis_and_viz/logo_plots_CDK2_4_6_canonical.ipynb`),
split per residue into `cdk2_4_6_hit_residues.csv`. A compound guide contributes
one row per substituted residue, so the position counts exceed 10: CDK2 11,
CDK4 9, CDK6 14.

### Variants

The chain is parameterised, not duplicated. `CDK246_TAG` suffixes every file it reads
and writes, and `CDK246_FDR` re-applies the hit rule at a tighter FDR:

```bash
# default: the notebook's hit set (FDR < 0.1 & LFC >= 1 in >=1 of 6 conditions)
python3 build_cdk2_4_6_hits.py

# stringent: same rule at FDR < 0.01 -> 36 of 93 guides
CDK246_TAG=_strictFDR CDK246_FDR=0.01 python3 build_cdk2_4_6_hits.py
```

then run `analyze_`, `deconflict_` (twice -- see below), `generate_` and the render
script with the same `CDK246_TAG` set.

| variant | CDK2 | CDK4 | CDK6 | positions drawn |
|---|---|---|---|---|
| default | 20 | 31 | 39 | 20 / 29 / 34 |
| `_strictFDR` | 8 | 12 | 16 | 8 / 12 / 13 |

**`FDR < 0.01` is not comparable across editors, and it shows here.** The surviving
CDK4 set is 11 CBE positions and 1 ABE (`I51`); the large-effect ABE variants `V148A`
(LFC 2.78), `N134G` (2.70) and `H95R` (2.48) are all gone while smaller CBE hits
remain. The strict panels are closer to *which residues the CBE arm called
confidently* than to *the strongest resistance residues*, and are best read beside the
default ones rather than instead of them. Same BH non-comparability the CDK12/13
notebook describes for conditions, biting across editors.

### The ligand-centred zoom panels

`compute_active_site_camera.py` also emits a per-gene `ZOOM_CAMERA_<GENE>`: the
**canonical rotation** with the origin moved in front of that gene's own ligand, at a
20 A half-width (40 A field, camera 74.6 A out). Same orientation as the wide panels,
so they read as a zoom of the same view rather than a different figure. CDK2 and CDK4
share a centre -- CDK4 is apo, so it uses the superposed 1QMZ ATP, the same thing its
wide panel draws -- and CDK6's AP9 sits ~3.5 A off. In frame: 12/20, 25/29, 26/34 of
the drawn hits.

`clip off` is re-issued after every re-aim, per the trap documented for the DDB1
panels: a view onto a small ligand leaves a thin slab that slices anything nearer the
camera.

**Labels need a second deconfliction pass for the zoom.** Label height is in scene
units, so at 3.6x the magnification the wide panel's 3.2 A labels render enormous and
clipped. The zoom panels use `LABEL_HEIGHT_ZOOM = 0.9` (3.2 / 3.6) and their own
offsets, solved against the zoom camera and stored in the `zoom_*` columns; only hits
in frame at that camera are labelled. Run it as a second pass:

```bash
ChimeraX --nogui --exit --cmd "runscript deconflict_labels.py"              # wide
CDK246_ZOOM=1 ChimeraX --nogui --exit --cmd "runscript deconflict_labels.py"  # zoom
python3 generate_cdk2_4_6_active_site.py
```

All six labelled zooms come out with **0 overlapping label pairs**.

Known cosmetic issue, the same one the CDK7 panel has: a label landing on its own red
sphere is black-on-dark-red and low-contrast. The placement cost function penalises
covering a sphere, but at CDK6's density a few still land on one. `LABEL_COLOR =
"white"` is the lever.

### Figures per variant

Three per gene (`_zoom` comes from the *unlabelled* panel by re-aiming; `_zoom_labeled`
is its own `.cxc`), so 12 PNGs per variant, 24 in total:

| | wide | wide + labels | ligand zoom | ligand zoom + labels |
|---|---|---|---|---|
| file | `{GENE}_active_site{TAG}.png` | `..._labeled.png` | `..._zoom.png` | `..._zoom_labeled.png` |

### Coverage

**0 wild-type mismatches** in all three structures, so the screen's codon
numbering agrees with the PDB author numbering. Two positions cannot be drawn:
**CDK6 E21 and Y24**, which are N-terminal to the ordered region of 2F2C
(chain B starts at residue 9 but 21/24 are absent from the model).

### What the panels show, and the three things that had to be measured

* **The view axis was chosen by measuring occlusion, not assumed.** The obvious
  heuristic -- look along kinase centroid -> ATP centroid -- leaves **57 protein
  heavy atoms between the camera and the ATP**, i.e. the pocket is buried behind a
  lobe. `compute_active_site_camera.py` samples 800 directions and counts the
  atoms inside a 5 A cylinder along each ray; the chosen axis has **4**. The roll
  is then set by the kinase -> cyclin axis so the pair reads side-on, which is
  what makes the view look straight-on rather than oblique.
* **Camera: the repo's canonical view.**  These panels use
  `cdk_config.CANONICAL_VIEW` plus the standard `zoom 0.72` pull-back -- byte-identical
  to `cxc/{GENE}.cxc` and the CDK7 SY-5609 panels -- so they line up with the rest of
  the figure set.  Because every structure is matchmaker-superposed onto 1QMZ chain A,
  one matrix orients all three identically.  Verified: the rendered content of
  `CDK6_active_site.png` measures 1506x1297 px against 1427x1257 px for the existing
  `CDK6.png`, i.e. same orientation and same scale.
* **Two bespoke cameras were built first and are not used.**
  `compute_active_site_camera.py` can derive (a) an occlusion-optimised axis -- it
  samples 800 directions and counts protein atoms in a 5 A cylinder along each ray,
  finding axes with 4 blocking atoms against 57 for the naive kinase-centroid -> ATP
  heuristic -- and (b) an axis constrained to put the hinge above the ligand and the
  cyclin on the right.  Both give a defensible close-up and **both were the wrong
  call**: a panel in its own private orientation cannot be compared against the rest of
  the figure set.  They are kept in `active_site_camera.txt` for reference.
* **Aspect.**  The other panels render 3000x1705, because the GUI window is not
  actually square despite `windowsize 1400 1400` and `save width 3000` follows the
  window.  These panels pass `height 3000` as well, giving a square 3000x3000 frame
  with the same content scale and therefore more vertical margin.  Drop the `height`
  argument in the render script to match the others' aspect exactly.

### Conventions specific to these panels

| | |
|---|---|
| hits | red CPK spheres, every modelled position |
| labels | `V96`-style, black, 3.2 A scene units, offsets from `deconflict_labels.py` |
| cyclin | wheat surface, **opaque** (as in the whole-kinase panels) |
| ligand | sticks at radius 0.4, element-coloured on grey carbons |

**CDK4 (2W96) is apo** -- it has no ligand to mark its pocket. The superposed
1QMZ ATP is drawn there instead at 55% transparency, following the convention
`cdk_config` already uses for CDK19's borrowed cyclin: an inferred component is
drawn semi-transparent so it reads as inferred. CDK2's ATP (+Mg) and CDK6's
aminopurvalanol are native and opaque.

### Four traps, for whoever edits this next

* **The framing step hides the ligand.** `view <spec>` measures *displayed* atoms,
  so the panels show the framing spec's atoms, run `view`, then hide them again.
  With `FRAME_ON="complex"` that spec is a whole chain -- which contains the
  ligand's HET residue -- so hiding it removes the ligand from the figure. The
  generator re-issues the hit *and* ligand display commands after the hide. The
  first version only re-issued the hits, and the ligand silently vanished.
* **Label every modelled hit, not a distance-selected subset.** An earlier version
  labelled only hits within 16 A of the ATP, which left hits just past the cutoff
  (CDK6 L33 at 16.5 A) drawn in frame as unlabelled spheres. A label on a hit
  outside the frame simply does not render, so labelling all of them cannot
  under-label.
* **Cyclin transparency.** Early drafts faded it to 25-50% to stop it covering the
  cleft; 50% washed out to nothing against white, the same effect the DDB1 notes
  describe. It is now opaque: with the cleft facing the camera the cyclin sits beside
  the kinase rather than over it. Checked -- **0 of the 20/29/34 drawn hits are fully
  hidden behind the opaque surface** in any of the three genes.
* **Occlusion tests must select the ligand by chain.** 1QMZ holds two copies of
  ATP; averaging both puts the test centroid between them and reports a ligand in
  the open as badly occluded (44 blocking atoms vs the real 1).

### Label overlap

`deconflict_labels.py` places the labels, and the two approaches that failed are worth
knowing about because both *looked* reasonable:

1. **Radial offsets in 3D** (push each label away from the pocket centroid). Overlap
   happens in the *projection*, and a 3D radial direction can project to almost
   nothing, so residues stacked along the view axis were pushed apart in a direction
   the camera cannot see. Going from 2.5 to 4 A only detached labels from their spheres.
2. **Force-directed relaxation in the camera plane, measured in Angstroms.** This
   assumes an orthographic camera. ChimeraX renders with perspective over a complex
   ~80 A deep, so a label near the camera draws larger than one at the back. The
   Angstrom-based check reported 1 overlapping pair for CDK6 where the render plainly
   showed six; with enough spring to hold labels on their spheres, the spring also
   cancelled the separating push (mean shift 0.7 A, and CDK2 went 11 -> 13 pairs).

What works is greedy placement in **normalized screen coordinates** read from the
panel's real camera, so perspective is included by construction: strongest hit first,
try the zero offset then rings of candidates, keep the first that hits no placed label
and covers no hit sphere. Current result:

| gene | labels | overlapping pairs before -> after | mean shift |
|---|---|---|---|
| CDK2 | 20 | 25 -> 1 | 7.1 A |
| CDK4 | 29 | 31 -> 0 | 7.8 A |
| CDK6 | 34 | 44 -> 2 | 7.1 A |

The residual two or three are unavoidable at this density: ChimeraX residue labels have
no leader lines, so past ~13 A of displacement a label stops being attributable to its
own sphere, and a few overlaps are preferable to detached labels. **Run it after
generating the panels and then regenerate** -- it reads the camera from the panel, and
the camera does not depend on the offsets, so one extra pass settles it.

`validate_cdk2_4_6_active_site.py` opens each panel headless and checks that it
runs, that **every labelled hit is inside the camera frustum**, that nothing is cut
by the clipping planes, and that the pocket marker is displayed. Currently all six
pass. Rendering itself cannot be validated headless -- `save` needs OpenGL, which
`--nogui` does not provide.

### Files

| file | what it is |
|------|-----------|
| `build_cdk2_4_6_hits.py` | top-N hits per gene -> `cdk2_4_6_hit_residues.csv` |
| `analyze_cdk2_4_6_active_site.py` | superposes, measures hit-to-ATP distances, derives the pocket and the label offsets; writes `cdk2_4_6_active_site_report.txt`, `cdk2_4_6_hit_distances.csv`, `cdk2_4_6_pocket_residues.csv` |
| `compute_active_site_camera.py` | occlusion-based camera; writes `active_site_camera.txt` |
| `fit_active_site_zoom.py` | measures frame fill; source of the now-unused `PER_GENE_ZOOM` |
| `deconflict_labels.py` | screen-space greedy label placement; rewrites the offset columns |
| `generate_cdk2_4_6_active_site.py` | writes the six panels + the render script |
| `validate_cdk2_4_6_active_site.py` | headless checks; writes `cdk2_4_6_active_site_validation.txt` |
| `cxc/{GENE}_active_site[_labeled].cxc` | the panels |
| `cxc/render_cdk2_4_6_active_site.cxc` | batch render to `figures/chimerax/` |
