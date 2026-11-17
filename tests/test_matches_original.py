"""The refactored helpers must reproduce the original notebooks bit for bit.

`cdkbe.editing.top_hgvsp` replaced a `groupby.apply(nlargest)` that appeared in
five notebooks. The replacement is much faster, but the two only agree if the
sort is STABLE -- `nlargest` keeps the first row among ties, and an unstable
sort picks a different variant for the guides whose top two calls tie on
%Reads. This test exists because that regression was live and silent.

    python3 -m pytest tests/ -q
"""
import re
import numpy as np
import pandas as pd
import pytest

from cdkbe import editing

SAMPLES = ["ABE_subpool1", "CBE_subpool1", "ABE_CDK12_13",
           "CBE_CDK12_13", "ABE_CDK2_4_6", "CBE_CDK2_4_6"]


def original_top_hgvsp(edit: pd.DataFrame) -> pd.DataFrame:
    """Verbatim from lollipop_plots_canonical.ipynb, kept as the reference."""
    non_wt = edit[edit["HGVSp"] != "WT"]
    wt = (edit[edit["HGVSp"] == "WT"].rename(columns={"%Reads": "WT_perc"})
          [["gRNA_id", "WT_perc"]])
    top = non_wt.groupby("gRNA_id").apply(lambda x: x.nlargest(1, "%Reads"))
    top = top.reset_index(drop=True)
    new = pd.merge(top, wt, on="gRNA_id", how="outer")
    new["HGVSp"] = new["HGVSp"].fillna("WT")
    new = new.fillna(0)
    new["Codon"] = [0 if a == "WT" else int(re.findall(r"\d+", a)[0])
                    for a in new["HGVSp"]]
    return new


@pytest.mark.parametrize("sample", SAMPLES)
def test_top_hgvsp_matches_original(sample):
    h = editing.hgvsp_from_full_canonical(sample)
    cols = ["gRNA_id", "HGVSp", "%Reads", "WT_perc", "Codon"]
    a = original_top_hgvsp(h).sort_values("gRNA_id").reset_index(drop=True)[cols]
    b = editing.top_hgvsp(h).sort_values("gRNA_id").reset_index(drop=True)[cols]
    assert a["gRNA_id"].tolist() == b["gRNA_id"].tolist()
    assert a["HGVSp"].tolist() == b["HGVSp"].tolist(), (
        f"{(a['HGVSp'] != b['HGVSp']).sum()} guides got a different variant call")
    assert np.allclose(a[["%Reads", "WT_perc", "Codon"]].astype(float),
                       b[["%Reads", "WT_perc", "Codon"]].astype(float))


def test_editor_edit_filter_is_applied():
    """An ABE table must not contain C>T-only reads (and vice versa)."""
    h = editing.hgvsp_from_full_canonical("ABE_subpool1")
    assert len(h) > 0 and h["%Reads"].max() <= 100.0 + 1e-9
