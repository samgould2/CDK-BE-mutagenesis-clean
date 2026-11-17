#!/usr/bin/env python3
"""Generate notebooks/04_lfc_fdr/lfc_fdr_tables.ipynb."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _lib.nbwrite import md, code, write, HEADER

ROOT = Path(__file__).resolve().parents[1]

cells = [
md("""# LFC, empirical p-values and FDR

Counts and RPM to log2 fold change per guide per condition, with an empirical
p-value against a non-targeting null and a Benjamini-Hochberg FDR.

**In** `data/processed/count_tables/` and `data/processed/rpm_tables/`.
**Out** `data/processed/lfc_fdr_tables/<flavour>_counts/<baseline>_base/`.

Every screen is built against three baselines, because they answer different
questions:

| baseline | the comparison | used for |
|---|---|---|
| `DMSO` | drug arm vs vehicle, same timepoint | resistance |
| `T0` | end of selection vs start | fitness |
| `Plasmid` | screen vs the cloned library | representation, dropout |

Ported from `LFC_FDR_pvals.ipynb`, which ran
the same 63-line block nine times. Note that as the original notebook was last
saved, its `samples` list had been narrowed to the two CDK2/4/6 screens and the
others commented out — running it top to bottom would have rebuilt two of the
nine tables, not all of them. The loop below covers all of them."""),

md("""## A reproducibility problem worth knowing about

The null distribution is drawn with `np.random.choice`, and **the original set
no seed anywhere**. Two runs of the same code therefore give different p-values
and FDRs.

Measured on `ABE_subpool1` / barcode / DMSO, rebuilding and comparing against
the shipped table:

- the 22 deterministic LFC columns match **exactly**
- FDR columns correlate 0.992–0.999 with the shipped values
- but the largest single difference reaches 0.22, and **98.7–100% of hit calls
  agree** at FDR < 0.1 — so up to about 1 guide in 80 changes hit status

`cdkbe.lfc_fdr` takes an explicit `seed` (default 0), so reruns are now
reproducible with each other. They will not reproduce the shipped tables to the
last digit, because those came from an unseeded run that cannot be recovered.

**Decision: the shipped tables stay the tables of record.** They are what the
manuscript reports, so they are not regenerated here. This notebook is the
reimplementation that documents how they were made and lets you rebuild
equivalent tables; `WRITE = False` below keeps it from overwriting them.

To rebuild anyway, set `WRITE = True` and `OUT_SUFFIX` to something other than
`""` so the rebuild lands beside the originals instead of on top of them."""),

code(HEADER + """
from cdkbe import counts, io, lfc_fdr, paths
"""),

md("""## 1. The guide library

Non-targeting guides supply the null distribution, so they are labelled here."""),

code("""library = counts.library()
library.loc[library["classification"] == "non-targeting control", "Gene"] = "NT"
library["classification"].value_counts()"""),

md("""## 2. What gets built

Nine screens x three count flavours x three baselines. Not every combination
exists — a screen with no plasmid arm has no `Plasmid` baseline — so the loop
skips what it cannot build rather than failing."""),

code("""SCREENS = [s["name"] for s in counts.SCREENS]
print(f"{len(SCREENS)} screens x {len(lfc_fdr.FLAVOUR_SUFFIX)} flavours "
      f"x {len(lfc_fdr.BASELINES)} baselines")
SCREENS"""),

md("""## 3. Build

`seed=0` makes this reproducible. Set `SEED = None` to reproduce the original's
unseeded behaviour.

This is the slow step: each table bootstraps 10,000 null draws per condition."""),

code('''SEED = 0
MIN_BASE_COUNTS = 100     # a guide needs this many baseline reads to inform the null

# The shipped tables are the tables of record (see above). Leave WRITE False to
# rebuild in memory and compare; set it True *and* give OUT_SUFFIX a value to
# write a rebuild beside them.
WRITE = False
OUT_SUFFIX = "_reseeded"

built, skipped = [], []
for flavour in lfc_fdr.FLAVOUR_SUFFIX:
    for baseline in lfc_fdr.BASELINES:
        out_dir = paths.LFC / f"{flavour}_counts" / f"{baseline}_base"
        for screen in SCREENS:
            try:
                table = lfc_fdr.build(screen, flavour, baseline, library,
                                      min_base_counts=MIN_BASE_COUNTS, seed=SEED)
            except (FileNotFoundError, KeyError) as e:
                skipped.append((screen, flavour, baseline, type(e).__name__))
                continue
            if WRITE:
                assert OUT_SUFFIX, "refusing to overwrite the shipped tables"
                paths.ensure(out_dir)
                io.write(table, out_dir / f"{screen}_{flavour}_{baseline}_base{OUT_SUFFIX}.csv")
            built.append((screen, flavour, baseline, len(table)))

print(f"{'wrote' if WRITE else 'rebuilt (not written)'} {len(built)} tables, "
      f"skipped {len(skipped)}")
for s in skipped[:10]:
    print("  skipped:", s)'''),

md("""## 4. Check against the shipped tables

The LFC columns are deterministic and must match exactly. The FDR columns will
not, for the reason above — what matters is that the hit calls barely move."""),

code('''shipped = io.read(paths.LFC / "barcode_counts" / "DMSO_base" /
                  "ABE_subpool1_barcode_DMSO_base.csv")
rebuilt = lfc_fdr.build("ABE_subpool1", "barcode", "DMSO", library, seed=SEED)

det = [c for c in shipped.columns
       if c.startswith(("LFC_avg_", "LFC_median_", "base_")) and c in rebuilt.columns]
bad = [c for c in det
       if not np.allclose(shipped[c], rebuilt[c], rtol=1e-9, atol=1e-12, equal_nan=True)]
print(f"deterministic columns: {len(det)} checked, {len(bad)} differ")

for c in [c for c in shipped.columns if c.startswith("FDR_")][:5]:
    a, b = shipped[c].to_numpy(), rebuilt[c].to_numpy()
    print(f"  {c:<34} corr={np.corrcoef(a, b)[0, 1]:.4f}  "
          f"hit-call agreement={100 * ((a < 0.1) == (b < 0.1)).mean():.1f}%")'''),

md("""## 5. Master tables

One long-format table per editor, stacking every screen and baseline, with each
baseline's columns suffixed so the three can sit side by side. This is what
`fig1_analysis` and the conservation modules read.

Ported from `master_table_gen.ipynb`, where loading nine screens x three
baselines took 169 lines of near-identical `read_csv` calls and 80 more
assigning `Screen_ID` and `Editor` one variable at a time."""),

code("""from cdkbe import master_table as mt

for editor in ("ABE", "CBE"):
    table = mt.build(editor)
    shipped = io.read(mt.path(editor))
    same = (len(table) == len(shipped)
            and list(table.columns) == list(shipped.columns))
    print(f"{editor}: {len(table):,} rows x {len(table.columns)} columns   "
          f"matches shipped: {same}")

mt.build("ABE").groupby(["Screen_ID", "Base_Value"]).size().unstack(fill_value=0)"""),

md("""Next: `notebooks/03_editing/` produces the per-gRNA variant calls that get
joined to these tables, and `notebooks/08_figures/` consumes both."""),
]

p = write(ROOT / "notebooks" / "04_lfc_fdr" / "lfc_fdr_tables.ipynb", cells)
print(f"wrote {p.relative_to(ROOT)} ({len(cells)} cells)")
