"""Xenograft panels.

invivo_xenograft_CDK9_CDK7_combined.ipynb
  cell 3  builds `fc` and `ki`;  cell 5 -> figure 2m;  cell 7 -> figure 2n
invivo_xenograft_CDK6_palbociclib.ipynb
  cell 3  builds the CDK6 table;  cell 20 -> figure 4n;  cell 22 -> figure 4o

Verified: the p-values the notebooks compute from these frames reproduce all
four annotations on 2m, all eight on 2n, and all four on each of 4n and 4o.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

from .. import io, paths

RAW = paths.RAW / "invivo_xenograft"

COHORTS = {
    "CDK9 I157V": dict(endpoint="Day5", drug="Atuveciclib"),
    "CDK7 L18F":  dict(endpoint="Day7", drug="SY-5609"),
}
#: sacrificed for illness after Day 3, tissue collected off-schedule
KI67_EXCLUDE = {"CDK9-ABE-16"}

_FC_COLS = ["MouseID", "Tag", "Cohort", "Genotype", "Genotype_label", "Drug",
            "Condition", "Treatment", "Day", "FoldChange"]
_KI_COLS = ["MouseID", "Tag", "Cohort", "Genotype", "Genotype_label", "Drug",
            "Condition", "Treatment", "Num Positive", "Num Detections", "Ki67 %"]


def _fold_changes() -> pd.DataFrame:
    fc = io.read(RAW / "fold_changes.csv")
    fc["Cohort"] = fc["Cohort"].str.replace("Cdk", "CDK")
    fc["Treatment"] = np.where(fc["Condition"] == "Control", "Control", fc["Drug"])
    fc["Genotype_label"] = np.where(fc["Genotype"] == "WT", "WT", fc["Cohort"])
    return fc


def _ki67() -> pd.DataFrame:
    fc = _fold_changes()
    ki = io.read_excel(RAW / "Ki67_unblinded_final.xlsx", "Ki67 Unblinded")
    ki["Tag"] = ki["Mouse Tag"].astype(int)
    ki["Genotype"] = ki["Genotype (Flank)"].str.split(" ").str[0]
    ki["Ki67 %"] = 100 * ki["Num Positive"] / ki["Num Detections"]
    info = fc.drop_duplicates(["Tag", "Genotype"])[
        ["Tag", "Genotype", "Cohort", "MouseID", "Drug", "Condition",
         "Treatment", "Genotype_label"]]
    ki = ki.merge(info, on=["Tag", "Genotype"], how="left", validate="one_to_one")
    assert ki["Cohort"].notna().all(), "a Ki67 tumour has no cohort"
    return ki[~ki["MouseID"].isin(KI67_EXCLUDE)].copy()


def figure_2m() -> pd.DataFrame:
    """Tumour volume fold change at each cohort's endpoint."""
    fc = _fold_changes()
    parts = [fc[(fc["Cohort"] == c) & (fc["Day"] == i["endpoint"])]
             for c, i in COHORTS.items()]
    return pd.concat(parts)[_FC_COLS]


def figure_2n() -> pd.DataFrame:
    """Ki67 positivity, both CDK9 and CDK7 cohorts."""
    ki = _ki67()
    return ki[ki["Cohort"].isin(COHORTS)][_KI_COLS]


def _cdk6() -> pd.DataFrame:
    return io.read(paths.PROCESSED / "invivo" / "C018_CDK6_per_tumor_annotated.csv")


def figure_4n() -> pd.DataFrame:
    """CDK6 R46W tumour volume fold change, palbociclib, Day 7."""
    d = _cdk6()
    d = d[(d["Cohort"] == "CDK6 R46W") & (d["Day"] == "Day7")]
    return d[["MouseID", "Tag", "Cohort", "Genotype", "Genotype_label", "Drug",
              "Condition", "Treatment", "Day", "V_day1", "V_day7",
              "FoldChange", "log2FC"]]


def figure_4o() -> pd.DataFrame:
    """CDK6 R46W Ki67, palbociclib."""
    d = _cdk6()
    d = d[d["Cohort"] == "CDK6 R46W"]
    return d[["MouseID", "Tag", "Cohort", "Genotype", "Genotype_label", "Drug",
              "Condition", "Treatment", "Viable cells", "Ki67+ cells", "Ki67 %"]]
