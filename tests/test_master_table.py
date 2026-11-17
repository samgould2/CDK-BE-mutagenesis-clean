"""The master tables must rebuild exactly.

`master_table.SOURCES` replaces 169 lines of `read_csv` plus 80 lines of
per-variable `Screen_ID`/`Editor` assignment. Two details in there are easy to
get wrong and neither raises: subpool 1's columns already lack the `ABE_`/`CBE_`
prefix that the other screens carry, and the compound-mutant screen uses proto
counts rather than barcode counts.
"""
import numpy as np
import pandas as pd
import pytest

from cdkbe import io, master_table as mt


@pytest.mark.parametrize("editor", ["ABE", "CBE"])
def test_master_table_matches_shipped(editor):
    built = mt.build(editor).reset_index(drop=True)
    shipped = io.read(mt.path(editor)).reset_index(drop=True)
    assert len(built) == len(shipped)
    assert list(built.columns) == list(shipped.columns)
    built = built[shipped.columns]
    for k in ("gRNA_id", "Screen_ID", "Editor", "Base_Value"):
        assert built[k].astype(str).tolist() == shipped[k].astype(str).tolist(), k
    for c in shipped.columns:
        if pd.api.types.is_numeric_dtype(shipped[c]):
            assert np.allclose(built[c].to_numpy(float), shipped[c].to_numpy(float),
                               rtol=1e-9, atol=1e-12, equal_nan=True), c


def test_compound_mutant_uses_proto_counts():
    """Its guides share barcodes, so barcode counts would be wrong."""
    entry = next(e for e in mt.SOURCES if e[0] == "Iterative")
    assert entry[3] == "proto"


def test_only_subpool1_keeps_its_column_prefix():
    strip = {e[2]: e[4] for e in mt.SOURCES}
    assert strip["ABE_subpool1"] is False and strip["CBE_subpool1"] is False
    assert all(v for k, v in strip.items() if k.startswith(("CDK", "SY_")))
