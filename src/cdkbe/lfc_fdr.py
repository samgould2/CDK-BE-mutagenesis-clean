"""Counts and RPM -> log2 fold change, empirical p-value and BH FDR per guide.

From LFC_FDR_pvals.ipynb, which ran the same
63-line block nine times (3 count flavours x 3 baselines).

Three things about the original are worth knowing, because they affect results:

1. **The bootstrap was unseeded.** `np.random.choice` draws the null
   distribution, and no seed was ever set, so two runs give different p-values
   and FDRs. `seed` is now an explicit argument with a default, which makes a
   rerun reproducible but means recomputed FDRs will not match the shipped
   tables to the last digit. See `notebooks/04_lfc_fdr/` for how far apart they
   land.

2. **`one_sided_FDR_v2` is two-sided.** It computes a high-tail and a low-tail
   p-value and keeps `min` of the two. The name is kept here as
   `empirical_fdr` to stop the confusion spreading.

3. **The function was defined twice.** A second, incompatible version appears
   later in the notebook for an exploratory section; rerunning an earlier cell
   after it raises a TypeError. Only the first version is ported.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests

from . import io, paths

PLASMID_NAMES = ("PLASMID", "Plasmid", "PLASMID_LIBRARY")
DEFAULT_SEED = 0
MIN_BASE_COUNTS = 100
NUM_BOOT = 10_000

FLAVOUR_SUFFIX = {"barcode": "bc", "proto": "proto", "matched": "matched"}
BASELINES = ("DMSO", "Plasmid", "T0")


def pick_columns(rpm: pd.DataFrame, base: str):
    """-> (base-mean columns, {condition: replicate columns}).

    Replicates are named `<prefix><condition>_REP{1,2,3}`, where the prefix is
    empty, `CBE_` or `ABE_` depending on the screen. The plasmid arm has a
    single column and no replicates.
    """
    conditions = list(rpm.columns[1:])
    if base not in PLASMID_NAMES:
        for prefix in ("", "CBE_", "ABE_"):
            if f"{prefix}{base}_REP1" in conditions:
                break
        else:
            raise KeyError(f"no replicate columns for baseline {base!r} in {conditions[:6]}...")
        base_cols = [f"{prefix}{base}_REP{i}" for i in (1, 2, 3)]
    else:
        base_cols = [base]

    rest = [c for c in conditions if c not in base_cols and c not in PLASMID_NAMES]
    uniq = sorted({c[:-5] for c in rest})            # strip _REPn
    samp = {c: [f"{c}_REP{i}" for i in (1, 2, 3)] for c in uniq}
    if base not in PLASMID_NAMES:
        for c in conditions:
            if c in PLASMID_NAMES:
                samp["Plasmid"] = [c]
    return base_cols, samp


def log_fold_change(rpm: pd.DataFrame, raw: pd.DataFrame, base_cols, samp: dict,
                    pseudocount: int = 1, base_method: str = "median") -> pd.DataFrame:
    """LFC per replicate, plus the per-condition average and median.

    `base_RPM`/`base_RAW` are carried through so downstream filters can drop
    guides whose baseline had almost no reads.
    """
    if base_method not in ("median", "average"):
        raise ValueError("base_method must be 'median' or 'average'")
    agg = np.median if base_method == "median" else np.average
    if len(base_cols) > 1:
        base_rpm = agg(rpm[base_cols], axis=1)
        base_raw = agg(raw[base_cols], axis=1)
    else:                                            # plasmid: a single arm
        base_rpm = rpm[base_cols[0]].to_numpy()
        base_raw = raw[base_cols[0]].to_numpy()

    out = []
    for condition, cols in samp.items():
        t = pd.DataFrame({c: np.log2((rpm[c] + pseudocount) / (base_rpm + pseudocount))
                          for c in cols})
        t[f"LFC_avg_{condition}"] = np.average(t[cols], axis=1)
        t[f"LFC_median_{condition}"] = np.median(t[cols], axis=1)
        t["base_RPM"] = base_rpm
        t["base_RAW"] = base_raw
        out.append(pd.concat((rpm[["gRNA_id"]], t), axis=1))
    combined = pd.concat(out, axis=1)
    return combined.loc[:, ~combined.columns.duplicated()].copy()


def empirical_fdr(null_ids, samp: dict, condition: str, lfc_df: pd.DataFrame,
                  num_boot: int = NUM_BOOT, seed: int | None = DEFAULT_SEED) -> pd.DataFrame:
    """Two-sided empirical p-value against a non-targeting null, then BH FDR.

    The null is `num_boot` draws with replacement from the pooled per-replicate
    LFCs of the non-targeting guides that cleared the base-count filter. Each
    replicate of each guide is scored against that null, the three are combined
    with Fisher's method, and BH is applied across guides. The stored p and FDR
    are the smaller of the high and low tails -- that is, two-sided.
    """
    cols = samp[condition]
    rng = np.random.default_rng(seed)
    nt = lfc_df[lfc_df["gRNA_id"].isin(null_ids)]
    data = nt[cols].to_numpy().ravel()
    if data.size == 0:
        raise ValueError(f"{condition}: no non-targeting guides passed the base-count filter")
    null = rng.choice(data, size=num_boot, replace=True)
    null_sorted = np.sort(null)
    n = null.size

    obs = lfc_df[cols].to_numpy()
    # p_high = P(null >= observed), p_low = P(null <= observed), +1 smoothing
    idx_ge = n - np.searchsorted(null_sorted, obs, side="left")
    idx_le = np.searchsorted(null_sorted, obs, side="right")
    p_high = (idx_ge + 1) / (n + 1)
    p_low = (idx_le + 1) / (n + 1)

    comb_high = np.array([combine_pvalues(r, method="fisher")[1] for r in p_high])
    comb_low = np.array([combine_pvalues(r, method="fisher")[1] for r in p_low])
    fdr_high = multipletests(comb_high, method="fdr_bh")[1]
    fdr_low = multipletests(comb_low, method="fdr_bh")[1]

    out = lfc_df.copy()
    out[f"p_{condition}"] = np.minimum(comb_high, comb_low)
    out[f"FDR_{condition}"] = np.minimum(fdr_high, fdr_low)
    return out


def resolve_plasmid_column(rpm: pd.DataFrame) -> str:
    for name in PLASMID_NAMES:
        if name in rpm.columns:
            return name
    raise KeyError(f"no plasmid column among {list(rpm.columns[:8])}...")


def build(screen: str, flavour: str, baseline: str, library: pd.DataFrame,
          min_base_counts: int = MIN_BASE_COUNTS, num_boot: int = NUM_BOOT,
          seed: int | None = DEFAULT_SEED) -> pd.DataFrame:
    """One screen x count flavour x baseline -> the LFC/FDR table."""
    suffix = FLAVOUR_SUFFIX[flavour]
    rpm = io.read(paths.RPM / f"{flavour}_counts" / f"{screen}_{suffix}_counts_RPM.csv")
    raw = io.read(paths.COUNTS / f"{flavour}_counts" / f"{screen}_{suffix}_counts.csv")
    base = resolve_plasmid_column(rpm) if baseline == "Plasmid" else baseline
    base_cols, samp = pick_columns(rpm, base)
    lfc = log_fold_change(rpm, raw, base_cols, samp)
    lfc = pd.merge(library[["gRNA_id", "Gene", "classification"]], lfc, on="gRNA_id")
    null_ids = lfc.loc[(lfc["classification"] == "non-targeting control")
                       & (lfc["base_RAW"] >= min_base_counts), "gRNA_id"]
    for condition in samp:
        lfc = empirical_fdr(null_ids, samp, condition, lfc, num_boot=num_boot, seed=seed)
    return lfc
