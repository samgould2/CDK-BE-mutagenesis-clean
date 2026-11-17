#!/usr/bin/env python3
"""Generate notebooks/02_counts_and_qc/count_tables.ipynb."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _lib.nbwrite import md, code, write, HEADER

ROOT = Path(__file__).resolve().parents[1]

cells = [
md("""# Count and RPM tables

Per-sample sequencing counts to per-screen count tables, and from those to RPM.

**In** `data/raw/screen_counts/<screen>/` — one file per sample per sequencing
lane, as produced by the demultiplexing step, plus `data/raw/sample_key.xlsx`.

**Out** `data/processed/count_tables/{barcode,proto,matched}_counts/` and the
matching `data/processed/rpm_tables/`.

Three count flavours come out of each screen and they differ in how strictly a
read had to match: `barcode` assigns on the sensor barcode, `proto` on the
protospacer, `matched` requires both to agree. Most analyses use `barcode`; the
compound-mutant screen uses `proto`, because its guides share barcodes.

Ported from `count_tables.ipynb`,
where the same four lines were repeated for each of the nine screen/editor
combinations. They are a table now (`cdkbe.counts.SCREENS`), so adding a screen
is a row rather than a copy-paste.

Verified: all 27 tables this produces are identical to the shipped ones
(`tests/test_counts.py`)."""),

code(HEADER + """
from cdkbe import counts, io, paths
"""),

md("""## 1. The sample key

Every sample, which screen it was sequenced with, and which subpool it targets."""),

code("""key = counts.sample_key()
print(f"{len(key)} samples")
key.groupby(["Screen_data_location", "Editor"]).size().unstack(fill_value=0)"""),

md("""## 2. Load the per-sample counts

A sample split across sequencing lanes has one file per lane; their counts add."""),

code("""raw = counts.load_raw_counts(key)
print(f"loaded {len(raw)} samples")

# every sample must list guides in the same order -- the tables below are
# assembled column-wise, so a reordered sample would misalign every guide
first = next(iter(raw.values()))["Guide_ID"].tolist()
assert all(d["Guide_ID"].tolist() == first for d in raw.values())
print("guide order is consistent across all samples")"""),

md("""## 3. Which samples go into which table

Two screens borrow a sample from another run: subpool 1's plasmid library was
sequenced with the CBE arm, and one CBE replicate was sequenced with the ABE
arm. `cdkbe.counts.SCREENS` records that as an `extra` entry."""),

code("""library = counts.library()

for screen in counts.SCREENS:
    samples = counts.samples_for(screen, key)
    print(f"{screen['name']:<18} {len(samples):>3} samples   pool {screen['pool']}")"""),

md("""## 4. Build and save the count tables"""),

code("""paths.ensure(*(paths.COUNTS / f"{f}_counts" for f in counts.FLAVOURS))

for screen in counts.SCREENS:
    tables = counts.count_tables(counts.samples_for(screen, key), raw,
                                 library, pool=screen["pool"])
    for flavour, df in tables.items():
        suffix = {"barcode": "bc_counts", "proto": "proto_counts",
                  "matched": "matched_counts"}[flavour]
        out = paths.COUNTS / f"{flavour}_counts" / f"{screen['name']}_{suffix}.csv"
        io.write(df, out)
    print(f"{screen['name']:<18} {len(tables['barcode']):>6,} guides x "
          f"{len(tables['barcode'].columns) - 1} samples")"""),

md("""## 5. RPM

No pseudocount here. Adding one at this stage would make it impossible to
filter on genuinely-zero guides downstream, so it is added at the LFC step
instead."""),

code("""paths.ensure(*(paths.RPM / f"{f}_counts" for f in counts.FLAVOURS))

for flavour in counts.FLAVOURS:
    src = paths.COUNTS / f"{flavour}_counts"
    dst = paths.RPM / f"{flavour}_counts"
    n = 0
    for f in sorted(src.glob("*.csv*")):
        name = f.name.split(".")[0]
        io.write(counts.rpm(io.read(f)), dst / f"{name}_RPM.csv")
        n += 1
    print(f"{flavour:<9} {n} tables")"""),

md("""Next: `notebooks/04_lfc_fdr/` turns these into LFC and FDR per guide per
condition."""),
]

p = write(ROOT / "notebooks" / "02_counts_and_qc" / "count_tables.ipynb", cells)
print(f"wrote {p.relative_to(ROOT)} ({len(cells)} cells)")
