"""LFC must be deterministic; FDR cannot be, and the test says why.

The original pipeline drew its null distribution with an unseeded
`np.random.choice`, so the shipped FDR values came from a run that cannot be
reproduced. What *can* be checked is that the LFC columns -- which involve no
randomness -- come out identical, and that the seeded FDR is stable across runs.
"""
import numpy as np
import pytest

from cdkbe import counts, io, lfc_fdr, paths

SCREEN, FLAVOUR, BASELINE = "ABE_subpool1", "barcode", "DMSO"


@pytest.fixture(scope="module")
def library():
    lib = counts.library()
    lib.loc[lib["classification"] == "non-targeting control", "Gene"] = "NT"
    return lib


@pytest.fixture(scope="module")
def rebuilt(library):
    return lfc_fdr.build(SCREEN, FLAVOUR, BASELINE, library, seed=0)


@pytest.fixture(scope="module")
def shipped():
    return io.read(paths.LFC / f"{FLAVOUR}_counts" / f"{BASELINE}_base" /
                   f"{SCREEN}_barcode_{BASELINE}_base.csv")


def test_lfc_columns_match_exactly(rebuilt, shipped):
    cols = [c for c in shipped.columns
            if c.startswith(("LFC_avg_", "LFC_median_", "base_")) and c in rebuilt.columns]
    assert cols, "no LFC columns found to compare"
    for c in cols:
        assert np.allclose(shipped[c], rebuilt[c], rtol=1e-9, atol=1e-12, equal_nan=True), c


def test_seeded_fdr_is_reproducible(library):
    a = lfc_fdr.build(SCREEN, FLAVOUR, BASELINE, library, seed=7)
    b = lfc_fdr.build(SCREEN, FLAVOUR, BASELINE, library, seed=7)
    col = next(c for c in a.columns if c.startswith("FDR_"))
    assert np.array_equal(a[col].to_numpy(), b[col].to_numpy())


def test_hit_calls_are_robust_to_the_seed(library):
    """Different seeds must not meaningfully change which guides are hits."""
    a = lfc_fdr.build(SCREEN, FLAVOUR, BASELINE, library, seed=0)
    b = lfc_fdr.build(SCREEN, FLAVOUR, BASELINE, library, seed=1)
    for col in [c for c in a.columns if c.startswith("FDR_")][:4]:
        agree = ((a[col] < 0.1) == (b[col] < 0.1)).mean()
        assert agree > 0.97, f"{col}: only {agree:.1%} of hit calls agree between seeds"
