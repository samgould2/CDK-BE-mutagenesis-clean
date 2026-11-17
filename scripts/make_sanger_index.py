#!/usr/bin/env python3
"""Index every Sanger trace this repository ships, and build the supplement.

    python3 scripts/make_sanger_index.py            # index + README + workbook
    python3 scripts/make_sanger_index.py --check    # compare, write nothing

The traces arrived as sequencing-vendor deliveries: one directory per order,
named after the order, with the sample labels only in the filenames. That is
enough to work from at the bench and not enough to publish, so this script
reads all of them, parses the labels into fields, pulls the run metadata out of
each ABIF file, and writes:

    data/raw/sanger/INDEX.csv        one row per trace, every field below
    data/raw/sanger/README.md        what each directory is
    data/raw/sanger/basecalls.fasta.gz   the vendor basecalls, so the set is
                                     usable without trace software
    source_data_tables/Supplementary_Data_Sanger_sequencing.xlsx
                                     the index as a two-sheet workbook

**The filename is the only record of what a sample is**, so the parse is
spelled out rather than guessed at. Azenta and Quintara both deliver
`<sample label>-<sequencing primer>.ab1`, and within the sample label this lab
writes the tube number, the gene, the variant, the treatment and the read
direction, separated by underscores. Each field below says what it matches;
anything that does not match is left empty rather than filled in, and
`INDEX.csv` keeps the raw label alongside, so no parse is load-bearing.

Three groups of traces are indexed together:

  endogenous_editing  the edited locus amplified out of treated cells, which
                      is what the validation and post-selection genotyping in
                      figures 2, 3 and 4 rest on
  library_colonies    colonies from the library cloning, read at the U6
                      promoter: the protospacer-barcode pairing QC
  guide_plasmids      the individual validation guide vectors, also read at
                      U6. These ship under data/reference/validation_guides
                      because that is where their designs are; they are
                      indexed here so the three sets have one table.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SANGER = ROOT / "data" / "raw" / "sanger"
GUIDES = ROOT / "data" / "reference" / "validation_guides"
INDEX = SANGER / "INDEX.csv"
README = SANGER / "README.md"
FASTA = SANGER / "basecalls.fasta.gz"
WORKBOOK = (ROOT / "source_data_tables"
            / "Supplementary_Data_Sanger_sequencing.xlsx")

#: where to look, and what to call what is found there
SOURCES = [("endogenous_editing", SANGER / "endogenous_editing"),
           ("library_colonies", SANGER / "library_colonies"),
           ("guide_plasmids", GUIDES)]

# `_` is a word character, so `\b` is no use here: every one of these tokens
# sits between underscores. The boundaries are spelled out instead.
EDGE_L, EDGE_R = r"(?<![A-Za-z0-9])", r"(?![A-Za-z0-9])"
#: a gene token anywhere in the sample label
GENE = re.compile(EDGE_L + r"(CDK\d+)" + EDGE_R, re.IGNORECASE)
#: an amino-acid substitution, e.g. L156F, G53S, E81Q
VARIANT = re.compile(EDGE_L + r"([A-Z]\d{1,4}[A-Z])" + EDGE_R)
#: a validation-guide label, e.g. `10A_CDK4_7925` -> gRNA_CDK4_targ_7925.
#: Three digits at least, so a colony plate label (`CDK1_5`) cannot match.
GUIDE = re.compile(EDGE_L + r"CDK(\d+)[_-](\d{3,5})" + EDGE_R, re.IGNORECASE)

#: treatment tokens this lab writes into sample labels, and what they mean.
#: Matched case-insensitively against the underscore-separated label.
TREATMENTS = {
    "palbo": "palbociclib", "abema": "abemaciclib", "tagto": "tagtociclib",
    "ribo": "ribociclib", "inx": "INX-315",
    "sy": "SY-5609", "kb": "KB-0742", "protac": "KI-CDK9d-32",
    "hq461": "HQ461", "drug": "drug (unspecified)",
}
#: selection state tokens, kept apart from the drug
STATES = {"sel": "selected", "selected": "selected", "select": "selected",
          "unsel": "unselected", "unselected": "unselected",
          "post": "selected", "wt": "wild type", "epo": "post-EPO"}
#: read direction, from the label's own suffix or from the primer name
FORWARD = {"f", "for", "fw", "forward"}
REVERSE = {"r", "rev", "rv", "reverse"}

FIELDS = ["category", "run", "file", "path", "sample_label", "primer",
          "read_part", "gene", "variants", "treatment", "selection", "direction",
          "guide_id", "guide_in_library", "run_date", "instrument", "well", "order",
          "dye", "read_length", "median_phred"]


#: Plasmidsaurus splits a long read into `<name>.1-of-2`, `<name>.2-of-2`.
#: That suffix is not a primer, and its hyphens would otherwise be mistaken
#: for the primer separator.
PART = re.compile(r"\.(\d+)-of-(\d+)$")


def split_primer(stem: str):
    """-> (sample label, sequencing primer, read part).

    Azenta and Quintara deliver `<sample label>-<primer>.ab1`, and the primer
    name can itself contain a hyphen -- `-F1_CDK12_8-22` carries the run date
    -- so the split is on the FIRST hyphen, not the last.
    """
    part = ""
    m = PART.search(stem)
    if m:
        stem, part = stem[:m.start()], f"{m.group(1)} of {m.group(2)}"
    label, sep, primer = stem.partition("-")
    return (label, primer, part) if sep else (stem, "", part)


def tokens(label: str):
    return [t for t in re.split(r"[_\s]+", label) if t]


def direction_of(label: str, primer: str) -> str:
    for source in (tokens(label)[::-1], tokens(primer)[::-1]):
        for t in source:
            low = t.lower().strip("0123456789-")
            if low in FORWARD:
                return "forward"
            if low in REVERSE:
                return "reverse"
    # `...-4F.ab1` / `...-MR.ab1`: the primer name's own last character
    if primer[-1:].upper() == "F":
        return "forward"
    if primer[-1:].upper() == "R":
        return "reverse"
    return ""


def parse_name(stem: str, category: str, colonies: dict) -> dict:
    """Split a trace's filename into the fields the lab wrote into it.

    `category` matters because the same spelling means different things in
    different sets. In `library_colonies` a label like `CDK1_5` is plate 1,
    colony 5 -- not the gene CDK1 -- so nothing is read out of it; the colony's
    identity comes from the protospacer the cloning QC called for it, which
    `colony_seq_identification.csv` holds.
    """
    label, primer, part = split_primer(stem)
    low = [t.lower() for t in tokens(label)]
    row = dict(sample_label=label, primer=primer, read_part=part, gene="",
               variants="", treatment="", selection="",
               direction=direction_of(label, primer), guide_id="")

    if category == "library_colonies":
        called = colonies.get(label, "")
        gene = GENE.search(called)
        row["guide_id"] = called
        row["gene"] = gene.group(1).upper() if gene else ""
        return row

    # the gene is sometimes only in the primer name (`WT_FW-Cdk7_FW`)
    gene = GENE.search(label) or GENE.search(primer)
    row["gene"] = gene.group(1).upper() if gene else ""
    row["variants"] = "_".join(VARIANT.findall(label))
    row["treatment"] = next((TREATMENTS[t] for t in low if t in TREATMENTS), "")
    row["selection"] = next((STATES[t] for t in low if t in STATES), "")
    guide = GUIDE.search(label)
    if guide:
        row["guide_id"] = f"gRNA_CDK{guide.group(1)}_targ_{int(guide.group(2))}"
    return row


def library_guides() -> set:
    """Every gRNA id in the screening library, to check a parsed id against.

    A guide id read out of a filename is a claim about what was sequenced. If
    it is not in the library the claim is wrong somewhere -- a typo in the
    label, or a rule here that is too eager -- and the index says so in
    `guide_in_library` rather than quietly passing it on.
    """
    try:
        from cdkbe import counts
        return set(counts.library()["gRNA_id"])
    except Exception as exc:                       # the library is not shipped
        print(f"  (guide ids not checked: {exc})", file=sys.stderr)
        return set()


def colony_calls() -> dict:
    """{colony sample id: the protospacer the cloning QC called for it}."""
    path = SANGER / "library_colonies" / "colony_seq_identification.csv"
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {r["Sample ID"]: r.get("proto_id", "")
                for r in csv.DictReader(fh) if r.get("Sample ID")}


def read_trace(path: Path) -> dict:
    """Run metadata and read quality, from the ABIF file's own header."""
    from Bio import SeqIO

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        record = SeqIO.read(str(path), "abi")
    raw = record.annotations.get("abif_raw", {})

    def text(key):
        v = raw.get(key, "")
        return v.decode("latin-1").strip() if isinstance(v, bytes) else str(v).strip()

    quality = record.letter_annotations.get("phred_quality", [])
    median = (sorted(quality)[len(quality) // 2] if quality else "")
    return dict(
        run_date=text("RUND1"), instrument=text("MCHN1"), well=text("TUBE1"),
        # the vendor's order number. The order *comment* also carries the
        # ordering scientist's name, which is not indexed here.
        order=text("CTNM1"), dye=text("DySN1"),
        read_length=len(record.seq), median_phred=median), str(record.seq)


def collect() -> tuple:
    """-> (index rows, {trace id: basecalls})"""
    rows, sequences = [], {}
    colonies, known = colony_calls(), library_guides()
    for category, root in SOURCES:
        if not root.is_dir():
            print(f"  {category}: {root} is not there", file=sys.stderr)
            continue
        for path in sorted(root.rglob("*.ab1")):
            run = path.relative_to(root).parts[0] if len(
                path.relative_to(root).parts) > 1 else root.name
            meta, sequence = read_trace(path)
            row = dict(category=category, run=run, file=path.name,
                       path=str(path.relative_to(ROOT)),
                       **parse_name(path.stem, category, colonies), **meta)
            if row["guide_id"] and known:
                row["guide_in_library"] = "yes" if row["guide_id"] in known else "no"
            rows.append({k: row.get(k, "") for k in FIELDS})
            sequences[row["path"]] = sequence
    rows.sort(key=lambda r: (r["category"], r["run"], r["file"]))
    return rows, sequences


def write_index(rows) -> None:
    INDEX.parent.mkdir(parents=True, exist_ok=True)
    with open(INDEX, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {len(rows)} traces -> {INDEX.relative_to(ROOT)}")


def write_fasta(rows, sequences) -> None:
    with gzip.open(FASTA, "wt") as fh:
        for row in rows:
            fh.write(f">{row['path']} {row['sample_label']} primer={row['primer']}\n")
            seq = sequences[row["path"]]
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
    print(f"  basecalls -> {FASTA.relative_to(ROOT)}")


def runs_table(rows) -> list:
    """One row per sequencing run: what it is, and what it covers."""
    out = {}
    for r in rows:
        key = (r["category"], r["run"])
        entry = out.setdefault(key, dict(category=r["category"], run=r["run"],
                                         traces=0, genes=set(), variants=set(),
                                         treatments=set(), dates=set()))
        entry["traces"] += 1
        for field, bucket in (("gene", "genes"), ("variants", "variants"),
                              ("treatment", "treatments"), ("run_date", "dates")):
            if r[field]:
                bucket_set = entry[bucket]
                bucket_set.update(v for v in str(r[field]).split("_") if v) \
                    if field == "variants" else bucket_set.add(r[field])
    table = []
    for entry in out.values():
        table.append(dict(
            category=entry["category"], run=entry["run"], traces=entry["traces"],
            first_run_date=min(entry["dates"]) if entry["dates"] else "",
            last_run_date=max(entry["dates"]) if entry["dates"] else "",
            genes=", ".join(sorted(entry["genes"])),
            variants=", ".join(sorted(entry["variants"])),
            treatments=", ".join(sorted(entry["treatments"]))))
    table.sort(key=lambda r: (r["category"], r["first_run_date"], r["run"]))
    return table


README_TEXT = """\
# Sanger sequencing

Every Sanger trace this repository ships, with one table over all of them.

`INDEX.csv` has one row per `.ab1`: where it lives, the sample label and
sequencing primer the vendor recorded, the gene, variant, treatment, selection
state and read direction parsed out of that label, and the run date,
instrument, well, order number, dye set, read length and median Phred score
read out of the ABIF file itself. `scripts/make_sanger_index.py` writes it and
documents every parsing rule; nothing here is typed by hand.

`basecalls.fasta.gz` carries the vendor basecall for each trace under the same
path, so the set can be searched or aligned without trace-viewing software.
The chromatograms are the primary record -- a base-edited locus is a mixed
peak, which a basecall flattens -- so quantification should go back to the
`.ab1`.

## What is here

| directory | what it is |
|---|---|
{table}

## Reading them

    pip install biopython
    from Bio import SeqIO
    record = SeqIO.read("data/raw/sanger/<path from INDEX.csv>", "abi")

`record.annotations["abif_raw"]` holds the full ABIF header, including the
four dye traces, for quantifying a mixed peak.

## Provenance

The traces are vendor deliveries, copied unchanged -- `data/MANIFEST.csv`
records each one's source path and SHA-256, and `scripts/verify_data.py`
checks them. The sequencing primers are in
`endogenous_editing/cdk_seq_primers.xlsx` and
`endogenous_editing/cdk246_seq_primers.xlsx`, and the SnapGene maps beside
them show where each primer sits in the transcript.
"""

DIRECTORY_NOTES = {
    "endogenous_editing": ("the edited locus amplified out of treated cells: "
                           "the validation and post-selection genotyping behind "
                           "figures 2, 3 and 4. Also holds the sequencing "
                           "primer tables and the SnapGene transcript maps"),
    "library_colonies": ("colonies from the library cloning, read at the U6 "
                         "promoter -- the protospacer/barcode pairing QC. "
                         "`colony_seq_identification.csv` is the call per colony"),
    "guide_plasmids": ("the individual validation guide vectors, read at the "
                       "U6 promoter. These ship under "
                       "`data/reference/validation_guides/`, beside the designs "
                       "they confirm, and are indexed here so that all three "
                       "sets appear in one table"),
}


def write_readme(rows) -> None:
    counts = {}
    for r in rows:
        counts[r["category"]] = counts.get(r["category"], 0) + 1
    lines = []
    for category, _ in SOURCES:
        note = DIRECTORY_NOTES[category]
        lines.append(f"| `{category}/` ({counts.get(category, 0)} traces) | {note} |")
    README.write_text(README_TEXT.format(table="\n".join(lines)))
    print(f"  readme -> {README.relative_to(ROOT)}")


def write_workbook(rows, runs) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    head = Font(name="Arial", bold=True)
    body = Font(name="Arial")

    wb = Workbook()
    ws = wb.active
    ws.title = "Runs"
    run_fields = ["category", "run", "traces", "first_run_date", "last_run_date",
                  "genes", "variants", "treatments"]
    ws.append([f.replace("_", " ") for f in run_fields])
    for i, r in enumerate(runs, start=2):
        ws.append([r[f] for f in run_fields])
        # the count is a formula over the Traces sheet, so the two sheets
        # cannot drift apart if a row is added or removed there
        ws.cell(row=i, column=3).value = (
            f'=COUNTIFS(Traces!$A$2:$A$100000,$A{i},Traces!$B$2:$B$100000,$B{i})')
    widths = {"A": 20, "B": 46, "C": 9, "D": 14, "E": 14, "F": 16, "G": 40, "H": 26}
    for col, width in widths.items():
        ws.column_dimensions[col].width = width

    ts = wb.create_sheet("Traces")
    ts.append([f.replace("_", " ") for f in FIELDS])
    for r in rows:
        ts.append([r[f] for f in FIELDS])
    for i, f in enumerate(FIELDS, start=1):
        ts.column_dimensions[get_column_letter(i)].width = \
            {"path": 60, "file": 38, "sample_label": 30}.get(f, 15)

    for sheet in (ws, ts):
        sheet.freeze_panes = "A2"
        for cell in sheet[1]:
            cell.font = head
            cell.alignment = Alignment(vertical="top", wrap_text=True)
        for row in sheet.iter_rows(min_row=2):
            for cell in row:
                cell.font = body
    WORKBOOK.parent.mkdir(parents=True, exist_ok=True)
    wb.save(WORKBOOK)
    print(f"  workbook -> {WORKBOOK.relative_to(ROOT)}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="compare with INDEX.csv on disk and write nothing")
    args = ap.parse_args()

    rows, sequences = collect()
    if not rows:
        raise SystemExit("no .ab1 files found")
    if args.check:
        if not INDEX.exists():
            print("INDEX.csv is not there yet")
            return 1
        with open(INDEX, newline="") as fh:
            have = [{k: str(v) for k, v in r.items()} for r in csv.DictReader(fh)]
        want = [{k: str(v) for k, v in r.items()} for r in rows]
        same = have == want
        print(f"  INDEX.csv {'matches' if same else 'DIFFERS'} "
              f"({len(have)} rows on disk, {len(want)} read)")
        return 0 if same else 1

    write_index(rows)
    write_fasta(rows, sequences)
    write_readme(rows)
    write_workbook(rows, runs_table(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
