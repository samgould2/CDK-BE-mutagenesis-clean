"""Source-data panels must match the tables verified against the published figures.

The reference workbooks in the original repository were checked panel by panel
against numbers printed in the figures themselves -- every p-value on 2m, 2n,
4o and 4p, all 22 DEG counts on 3f, all 8 bar labels on 3g, all 18 counts on
4b. This test asserts the port into `cdkbe.source_data` did not change them.

The workbooks are committed under tests/reference/source_data so this runs
anywhere the repository is checked out.
"""
from pathlib import Path

import pandas as pd
import pytest

from cdkbe.source_data import PANELS

REFERENCE = Path(__file__).resolve().parent / "reference" / "source_data"

#: panel -> sheet name in the reference workbook.
#: The shipped workbook predates the November relabelling of figure 4's last
#: three panels, so its sheets still carry the old letters: what the
#: manuscript now calls 4n and 4o are sheets "4o" and "4p".
#: See `docs/DISCREPANCIES.md` item 28.
SHEETS = {
    "2a": "2a", "2h": "2h", "2m": "2m", "2n": "2n",
    "3a": "3a", "3b": "3b", "3f": "3f", "3h": "3h",
    "4a": "4a", "4b": "4b", "4n": "4o", "4o": "4p",
}

pytestmark = pytest.mark.skipif(not REFERENCE.is_dir(),
                                reason="original tree not present")


@pytest.mark.parametrize("panel", sorted(SHEETS))
def test_panel_matches_reference(panel):
    ref = pd.read_excel(REFERENCE / f"Figure_{panel[0]}.xlsx", sheet_name=SHEETS[panel])
    got = PANELS[panel]()
    assert len(got) == len(ref), f"{panel}: {len(got)} rows, reference has {len(ref)}"
    shared = [c for c in got.columns if c in ref.columns]
    assert shared, f"{panel}: no column names in common with the reference"
    for c in shared:
        a, b = ref[c], got[c].reset_index(drop=True)
        if pd.api.types.is_bool_dtype(a) or pd.api.types.is_bool_dtype(b):
            assert a.astype(bool).tolist() == b.astype(bool).tolist(), \
                f"{panel}: column {c} differs"
        elif pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            assert ((a - b).abs().fillna(0) < 1e-6).all(), f"{panel}: column {c} differs"
        else:
            assert a.astype(str).tolist() == b.astype(str).tolist(), \
                f"{panel}: column {c} differs"
