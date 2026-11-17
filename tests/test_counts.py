"""The count pipeline must rebuild the shipped tables from the raw per-sample files.

`cdkbe.counts.SCREENS` replaces nine copy-pasted blocks. Two of those blocks
borrowed a sample from another screen's sequencing run -- subpool 1's plasmid
library and one CBE replicate -- and getting that wrong drops a column silently
rather than raising. These tests check the column sets as well as the numbers.
"""
import numpy as np
import pytest

from cdkbe import counts, io, paths

FLAVOUR_SUFFIX = {"barcode": "bc_counts", "proto": "proto_counts",
                  "matched": "matched_counts"}


@pytest.fixture(scope="module")
def built():
    key = counts.sample_key()
    raw = counts.load_raw_counts(key)
    library = counts.library()
    return {s["name"]: counts.count_tables(counts.samples_for(s, key), raw,
                                           library, pool=s["pool"])
            for s in counts.SCREENS}


@pytest.mark.parametrize("screen", [s["name"] for s in counts.SCREENS])
@pytest.mark.parametrize("flavour", sorted(FLAVOUR_SUFFIX))
def test_count_table_matches_shipped(built, screen, flavour):
    got = built[screen][flavour]
    shipped = io.read(paths.COUNTS / f"{flavour}_counts" /
                      f"{screen}_{FLAVOUR_SUFFIX[flavour]}.csv")
    assert list(got.columns) == list(shipped.columns), "sample set differs"
    assert got["gRNA_id"].tolist() == shipped["gRNA_id"].tolist()
    assert np.array_equal(got.iloc[:, 1:].to_numpy(), shipped.iloc[:, 1:].to_numpy())


def test_rpm_sums_to_a_million():
    t = io.read(paths.COUNTS / "barcode_counts" / "ABE_subpool1_bc_counts.csv")
    r = counts.rpm(t)
    assert np.allclose(r.iloc[:, 1:].sum(axis=0), 1e6)


def test_rpm_adds_no_pseudocount():
    """A guide with zero reads must stay zero, so it can be filtered later."""
    t = io.read(paths.COUNTS / "barcode_counts" / "ABE_subpool1_bc_counts.csv")
    r = counts.rpm(t)
    zeros = (t.iloc[:, 1:] == 0).to_numpy()
    if zeros.any():
        assert (r.iloc[:, 1:].to_numpy()[zeros] == 0).all()
