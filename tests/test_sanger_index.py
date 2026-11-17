"""The Sanger trace index.

`data/raw/sanger/INDEX.csv` is the only description of what each of the 288
`.ab1` files is -- the vendor directories carry an order number and the
filenames carry everything else. `scripts/make_sanger_index.py` parses those
filenames, so what is worth testing is that the parse still says the same
thing, that every row points at a file that is really there, and that the
claims it makes which can be checked against something else do check out.

The one that can be checked is `guide_id`: a guide read out of a filename has
to exist in the screening library. One does not, and that is recorded here as
a fact about the data rather than left to be rediscovered -- see
`docs/DISCREPANCIES.md` item 33.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from cdkbe import counts, paths  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INDEX = paths.RAW / "sanger" / "INDEX.csv"
README = paths.RAW / "sanger" / "README.md"
FASTA = paths.RAW / "sanger" / "basecalls.fasta.gz"
WORKBOOK = (ROOT / "source_data_tables"
            / "Supplementary_Data_Sanger_sequencing.xlsx")


@pytest.fixture(scope="module")
def rows():
    with open(INDEX, newline="") as fh:
        return list(csv.DictReader(fh))


def test_the_index_and_its_companions_are_there(rows):
    assert len(rows) == 288
    assert README.is_file() and FASTA.is_file() and WORKBOOK.is_file()


def test_every_row_points_at_a_file_that_exists(rows):
    for row in rows:
        assert (ROOT / row["path"]).is_file(), row["path"]


def test_every_trace_on_disk_has_a_row(rows):
    indexed = {row["path"] for row in rows}
    on_disk = {str(p.relative_to(ROOT))
               for root in (paths.RAW / "sanger",
                            paths.REFERENCE / "validation_guides")
               for p in root.rglob("*.ab1")}
    assert on_disk == indexed


def test_the_three_groups_are_all_present(rows):
    counted = {}
    for row in rows:
        counted[row["category"]] = counted.get(row["category"], 0) + 1
    assert counted == {"endogenous_editing": 187,
                       "guide_plasmids": 56,
                       "library_colonies": 45}


def test_every_guide_id_but_one_is_in_the_screening_library(rows):
    """The one that is not is a mislabelled file, not a parsing slip."""
    known = set(counts.library()["gRNA_id"])
    missing = sorted({row["guide_id"] for row in rows
                      if row["guide_id"] and row["guide_id"] not in known})
    assert missing == ["gRNA_CDK13_targ_2232"]
    # the same run holds CDK12 guide 2232, and this file was read with the
    # CDK12 primer, so the label's gene is what is wrong
    flagged = [r for r in rows if r["guide_in_library"] == "no"]
    assert len(flagged) == 2
    assert all(r["primer"].startswith(("F1_CDK12", "R1_CDK12")) for r in flagged)
    assert "gRNA_CDK12_targ_2232" in known


def test_the_colony_traces_are_identified_by_their_called_protospacer(rows):
    """A colony label is a plate and well, so its guide comes from the QC table."""
    colonies = [r for r in rows if r["category"] == "library_colonies"]
    assert colonies and all(r["sample_label"].startswith("CDK") for r in colonies)
    # the cloning QC called 30 of the 45; the other 15 are not in
    # colony_seq_identification.csv, so their rows carry no guide
    called = [r for r in colonies if r["guide_id"]]
    assert len(called) == 30
    assert all(r["guide_id"].startswith("gRNA_CDK") for r in called)
    # and nothing read a plate number as a gene: where there is a gene it is
    # the one the called protospacer targets, never the plate's own `CDK1`
    assert all(r["gene"] == r["guide_id"].split("_")[1] for r in called)
    assert all(r["gene"] == "" for r in colonies if not r["guide_id"])


def test_the_endogenous_traces_are_paired_forward_and_reverse(rows):
    """Most amplicons were read both ways; the index should show both."""
    endo = [r for r in rows if r["category"] == "endogenous_editing"]
    forward = sum(1 for r in endo if r["direction"] == "forward")
    reverse = sum(1 for r in endo if r["direction"] == "reverse")
    assert abs(forward - reverse) <= 6
    assert forward + reverse > 0.8 * len(endo)


def test_the_guide_plasmid_reads_are_single_reads_off_the_u6_promoter(rows):
    plasmids = [r for r in rows if r["category"] == "guide_plasmids"]
    u6 = [r for r in plasmids if r["primer"].upper() == "U6"]
    assert len(u6) >= 40
    assert all(r["direction"] == "" for r in u6)


def test_the_run_metadata_comes_from_the_files_themselves(rows):
    for row in rows:
        assert int(row["read_length"]) > 0
        assert 0 <= int(row["median_phred"]) <= 93
        if row["run_date"]:
            assert row["run_date"][:2] == "20", row["file"]


def test_only_the_plasmidsaurus_deliveries_carry_no_run_date(rows):
    """Their ABIF files have no RUND field; the capillary vendors' all do."""
    undated = [r for r in rows if not r["run_date"]]
    assert len(undated) == 16
    assert all(r["run"].endswith("_results") for r in undated)
    assert all(r["read_part"] for r in undated)


def test_the_drugs_named_in_the_labels_are_the_ones_in_the_screens(rows):
    named = {r["treatment"] for r in rows if r["treatment"]}
    assert named <= {"palbociclib", "abemaciclib", "tagtociclib", "ribociclib",
                     "INX-315", "SY-5609", "KB-0742", "KI-CDK9d-32", "HQ461",
                     "drug (unspecified)"}
    assert {"palbociclib", "SY-5609", "KB-0742"} <= named


def test_the_script_regenerates_exactly_what_is_on_disk():
    out = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "make_sanger_index.py"), "--check"],
        capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "DIFFERS" not in out.stdout
