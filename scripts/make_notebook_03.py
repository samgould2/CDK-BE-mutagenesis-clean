#!/usr/bin/env python3
"""Generate notebooks/03_editing/editing_tables.ipynb."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _lib.nbwrite import md, code, write, HEADER

ROOT = Path(__file__).resolve().parents[1]

cells = [
md("""# Sensor editing tables

Per-sample, per-read sensor calls to the per-gRNA editing tables that every
residue-level figure depends on.

**In** `data/raw/screen_editing/<screen>/` — one file per sample per sequencing
lane, each row a distinct edit seen on that guide's sensor.
**Out** `data/processed/editing/<screen>_{full,full_canonical,HGVSp,HGVSp_canonical}`.

| table | one row per | used for |
|---|---|---|
| `full` | distinct edit, with its DNA change | anything needing the DNA-level edit |
| `full_canonical` | the same, canonical edits only | the per-gRNA variant call |
| `HGVSp` | (gRNA, protein consequence) | variant-level summaries |
| `HGVSp_canonical` | the same, canonical only | coverage and logo plots |

**Canonical** means `Canonical_edit` is true — the edit is the one this editor
makes (A>G for ABE, C>T for CBE). `Canonical_window` records whether it also
fell in the expected window, and is deliberately *not* part of the filter.

Ported from `sensor_HGVSp_pipeline.ipynb`,
which is a working log rather than a pipeline: it also holds the library design,
the inputs for the cluster CRISPResso jobs and several exploratory sections,
defines `extract_numbers` four times and `HGVSp_combiner` twice, and has the
lines that save these tables commented out. Only the reduction is ported."""),

md("""## The shipped tables came from an earlier CRISPResso run

The per-sample files in `data/raw/screen_editing/` are a **redo** run — their
names say so. The tables shipped in `data/processed/editing/` were built from
the earlier run, whose per-sample outputs are not in the repository. So
rebuilding here does not reproduce them byte for byte:

- shipped `CBE_subpool1_full`: 351,877,689 reads over 1,489,067 rows
- rebuilt from the files present: 345,445,385 reads over 1,431,274 rows

about 2.4% fewer reads. No subset of the available files reproduces the shipped
total, so the difference is a different upstream run, not a selection mistake.

**It does not change any conclusion.** Reducing both to one variant call per
guide:

- **2,439 of 2,439 top variant calls are identical (100%)**
- WT-fraction correlates 0.9999

Read depth differs; what each guide edits does not. The shipped tables remain
the tables of record, and this notebook rebuilds equivalent ones."""),

code(HEADER + """
from cdkbe import counts, editing, editing_pipeline as ep, io, paths
"""),

md("""## 1. Which samples make up each screen's table"""),

code("""key = counts.sample_key()
library = counts.library()

for screen in ep.SCREEN_POOL:
    try:
        files = ep.sample_files(screen, key)
    except FileNotFoundError:
        print(f"{screen:<22} no raw editing files shipped")
        continue
    location, pool = ep.SCREEN_POOL[screen]
    print(f"{screen:<22} {len(files):>3} files from {location:<16} pool {pool}")"""),

md("""## 2. Build

Each screen sums read counts across its samples, keyed on the distinct edit,
then takes percentages within each targeting guide.

This is the slow step — it reads a few hundred MB of gzipped per-read calls."""),

code('''paths.ensure(paths.EDITING)
REBUILD = False     # True overwrites the shipped tables with the rebuild

summary = []
for screen in ep.SCREEN_POOL:
    try:
        tables = ep.build(screen, key, library)
    except FileNotFoundError as e:
        print(f"{screen:<22} skipped: {e}")
        continue
    for kind, df in tables.items():
        if REBUILD:
            io.write(df, paths.EDITING / f"{screen}_{kind}.csv")
        summary.append(dict(screen=screen, table=kind, rows=len(df),
                            reads=int(df["#Reads"].sum())))
    print(f"{screen:<22} " + "  ".join(f"{k}={len(v):,}" for k, v in tables.items()))

pd.DataFrame(summary)'''),

md("""## 3. Check: the variant calls agree

The thing downstream code actually uses is one variant call per guide. Compare
the shipped tables with the rebuild at that level."""),

code('''screen = "CBE_subpool1"

shipped_call = editing.top_hgvsp(editing.hgvsp_from_full_canonical(screen))

rebuilt = ep.build(screen, key, library)["full_canonical"]
ref, alt = editing.EDITS[screen[:3]]
keep = rebuilt["DNA Change"].str.findall(editing._DNA_CHANGE).apply(
    lambda e: all((r == ref) and (a == alt) for _, r, a in e))
h = (rebuilt[keep].groupby(["gRNA_id", "HGVSp"], as_index=False)["#Reads"].sum()
     .sort_values(["gRNA_id", "#Reads"], ascending=[True, False]).reset_index(drop=True))
h["%Reads"] = h["#Reads"] / h.groupby("gRNA_id")["#Reads"].transform("sum") * 100
rebuilt_call = editing.top_hgvsp(h)

m = shipped_call[["gRNA_id", "HGVSp"]].merge(
    rebuilt_call[["gRNA_id", "HGVSp"]], on="gRNA_id", how="outer",
    suffixes=("_shipped", "_rebuilt"))
agree = (m["HGVSp_shipped"] == m["HGVSp_rebuilt"]).sum()
print(f"{agree:,} / {len(m):,} guides get the same top variant "
      f"({100 * agree / len(m):.1f}%)")'''),

md("""Next: `notebooks/04_lfc_fdr/` joins these calls to the LFC tables."""),
]

p = write(ROOT / "notebooks" / "03_editing" / "editing_tables.ipynb", cells)
print(f"wrote {p.relative_to(ROOT)} ({len(cells)} cells)")
