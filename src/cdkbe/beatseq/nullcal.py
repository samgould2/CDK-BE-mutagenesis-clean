"""Size-matched calibration of the pseudo-guide null. Shared by severity.py and norm_cc.py.

WHY THIS MODULE EXISTS
----------------------
Both `severity.py` and `norm_cc.py` standardised the cross-fitted statistic D against a
null POOLED across cell counts:

    stats = null.groupby("arm")["D"].agg(["mean", "std"])
    severity = (D - stats["mean"]) / stats["std"]

The stated justification, in `severity.build_null`'s docstring, was "D is n-unbiased, so one
pooled null per arm suffices". The verification actually run was

    spearmanr(null.n_cells, null.D)

which tests whether the null's MEAN depends on n. It does not -- that is what cross-fitting
buys. But the SD was never checked, and D is a split-half inner product whose sampling
variance falls steeply with cell number: on the shipped null it varies about ELEVEN-FOLD
between the smallest and largest size drawn. Dividing every guide by one pooled SD is
therefore conservative for well-covered guides and anti-conservative for sparse ones, and
the empirical p inherits the same error more severely, because the pooled tail is set by the
small-n draws.

A second, independent defect: the sizes were hardcoded ([40, 80, 150, 300, 500] for A549,
[25, 40, 70, 110, 160] for T47D) and do not span the real guide range, so 54%, 69%, 24% and
0% of guides in the four arms fell outside any size the null was drawn at.

WHAT THIS MODULE DOES
---------------------
`sizes_for_guides`  draw sizes from the observed guide cell counts, so the null spans them
`fit_size_model`    per-size mean and SD, plus a size-standardised pooled draw set
`standardise`       score each guide against ITS OWN size, and compute p and q from the
                    size-standardised draws

    severity = (D - mu(n_g)) / sigma(n_g)          mu, sigma interpolated in log n
    p        = fraction of size-standardised null draws >= severity
    q        = Benjamini-Hochberg within arm

`calibration_report` runs the four checks that should have been run originally:

  1. does the null SD depend on n?            (the check that was missed)
  2. var(severity) across NT guides           (must be ~1; 95% interval 0.46-3.67 at n=9)
  3. rho(severity, cells per guide)           (must be ~0 after correction)
  4. NT p-value uniformity                    (fraction below 0.05 must be ~5%)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kstest, binomtest, mannwhitneyu

# 95% interval for a variance estimate from 9 observations (chi-square, df=8)
CHI2_LO, CHI2_HI = 0.46, 3.67


def sizes_for_guides(n_cells, k: int = 9) -> list[int]:
    """Draw sizes spanning the observed guide cell counts, including both extremes."""
    n = np.asarray(n_cells, dtype=int)
    q = np.percentile(n, np.linspace(0, 100, k)).astype(int)
    return sorted({int(x) for x in q if x >= 20})


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    o = np.argsort(p)
    n = len(p)
    adj = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        running = min(running, p[o[rank]] * n / (rank + 1))
        adj[o[rank]] = running
    return adj


def fit_size_model(null: pd.DataFrame, d_col: str = "D"):
    """Per-arm: per-size moments, and the size-standardised pooled draws.

    Standardising each draw by its own size band makes the pooled set size-free, so it can
    legitimately be used as a common reference for guides of any n. The assumption is that
    the null SHAPE does not change with n; `calibration_report` checks it.
    """
    models, std_draws = {}, {}
    for arm, d in null.groupby("arm"):
        g = d.groupby("n_cells")[d_col].agg(["mean", "std", "count"])
        g = g[g["count"] >= 30]
        models[arm] = g
        std_draws[arm] = np.concatenate([
            (d.loc[d.n_cells == k, d_col].to_numpy() - g.loc[k, "mean"]) / g.loc[k, "std"]
            for k in g.index])
    return models, std_draws


def standardise(m: pd.DataFrame, models: dict, std_draws: dict, d_col: str = "D",
                prefix: str = "") -> pd.DataFrame:
    """Add <prefix>severity, <prefix>p, <prefix>q, scoring each guide at its own size."""
    m = m.copy()
    mu, sd = [], []
    for _, r in m.iterrows():
        g = models[r.arm]
        ln = np.log(g.index.to_numpy(float))
        x = float(np.log(np.clip(r.n_cells, g.index.min(), g.index.max())))
        mu.append(float(np.interp(x, ln, g["mean"].to_numpy())))
        sd.append(float(np.exp(np.interp(x, ln, np.log(g["std"].to_numpy())))))
    m[prefix + "null_mean"] = mu
    m[prefix + "null_sd"] = sd
    m[prefix + "severity"] = (m[d_col] - m[prefix + "null_mean"]) / m[prefix + "null_sd"]

    m[prefix + "p"] = [
        float((1 + (std_draws[r.arm] >= r[prefix + "severity"]).sum())
              / (1 + len(std_draws[r.arm]))) for _, r in m.iterrows()]
    m[prefix + "q"] = np.nan
    for arm, idx in m.groupby("arm").groups.items():
        m.loc[idx, prefix + "q"] = bh(m.loc[idx, prefix + "p"].to_numpy())
    return m


def calibration_report(m: pd.DataFrame, null: pd.DataFrame, models: dict,
                       d_col: str = "D", prefix: str = "",
                       old_sev: str | None = None, old_p: str | None = None) -> pd.DataFrame:
    """The four checks. Returns a per-arm table and prints a verdict."""
    rows = []
    for arm, g in models.items():
        d = null[null.arm == arm]
        sd_ratio = float(g["std"].max() / g["std"].min())
        nt = m[(m.arm == arm) & m.is_control]
        tg = m[(m.arm == arm) & (~m.is_control)]
        var_nt = float(nt[prefix + "severity"].var())
        rho_n = float(spearmanr(tg.n_cells, tg[prefix + "severity"]).statistic)
        nt_fp = float((nt[prefix + "p"] < 0.05).mean())
        ks = float(kstest(nt[prefix + "p"], "uniform").pvalue)
        row = dict(arm=arm, null_sd_range=sd_ratio,
                   pooled_sd=float(d[d_col].std()),
                   var_nt_severity=var_nt,
                   calibrated=bool(CHI2_LO <= var_nt <= CHI2_HI),
                   rho_severity_vs_ncells=rho_n,
                   nt_frac_p_below_05=nt_fp,
                   nt_p_uniform_ks_p=ks)
        if old_sev is not None:
            row["var_nt_severity_old"] = float(nt[old_sev].var())
            row["rho_severity_vs_ncells_old"] = float(
                spearmanr(tg.n_cells, tg[old_sev]).statistic)
        if old_p is not None:
            row["nt_frac_p_below_05_old"] = float((nt[old_p] < 0.05).mean())
        rows.append(row)
    rep = pd.DataFrame(rows)

    print("\n  CALIBRATION CHECKS")
    print("  1. null SD range across sizes (the check that was originally missed):")
    for _, r in rep.iterrows():
        print(f"       {r.arm:9s} SD varies {r.null_sd_range:5.1f}x across sizes  "
              f"(pooled SD {r.pooled_sd:.2e})")
    print("  2. var(severity) across the 9 NT guides — target 1, 95% interval "
          f"{CHI2_LO}-{CHI2_HI}:")
    for _, r in rep.iterrows():
        old = f"  (was {r.var_nt_severity_old:.3f})" if old_sev else ""
        print(f"       {r.arm:9s} {r.var_nt_severity:6.3f}"
              f"  {'OK' if r.calibrated else 'STILL OFF'}{old}")
    print("  3. rho(severity, cells per guide) — must be ~0:")
    for _, r in rep.iterrows():
        old = f"  (was {r.rho_severity_vs_ncells_old:+.3f})" if old_sev else ""
        print(f"       {r.arm:9s} {r.rho_severity_vs_ncells:+.3f}{old}")
    print("  4. NT p-values below 0.05 — target 5%:")
    for _, r in rep.iterrows():
        old = f"  (was {100*r.nt_frac_p_below_05_old:.0f}%)" if old_p else ""
        print(f"       {r.arm:9s} {100*r.nt_frac_p_below_05:5.1f}%   "
              f"KS uniformity p={r.nt_p_uniform_ks_p:.3f}{old}")
    n_ok = int(rep.calibrated.sum())
    print(f"\n  VERDICT: {n_ok} of {len(rep)} arms calibrated on check 2.")
    return rep


def nt_referenced(m: pd.DataFrame, sev_col: str = "severity",
                  n_perm: int = 50_000, seed: int = 0) -> pd.DataFrame:
    """Hit count referenced to the NON-TARGETING guides, using no null model at all.

    WHY THIS IS NEEDED EVEN AFTER SIZE-MATCHING. Pseudo-guides are random draws from the
    control cell pool, so they capture sampling noise but NOT guide-specific technical
    variation -- integration site, clonal expansion, representation. Real non-targeting
    guides capture both. Two of the nine NT guides in A549/ABE reach p<0.01 against the
    size-matched pseudo-guide null, which means that null is still anti-conservative for
    real guides.

    The assumption-free alternative: count targeting guides whose severity exceeds the
    MAXIMUM of that arm's nine NT guides. Under exchangeability of all n+k guides the
    expected count is n/(k+1), i.e. n/10 for k=9.

    THE TEST IS NOT BINOMIAL -- this was a bug, fixed 2026-08-19.
    ------------------------------------------------------------
    The original implementation used ``binomtest(k, n, 1/10)``, which assumes the n
    targeting guides are independent Bernoulli trials against a FIXED threshold. The
    threshold is not fixed: it is the maximum of nine noisy NT draws, and every targeting
    guide is compared to the SAME draw. That shared threshold induces strong positive
    correlation -- if the NT maximum lands low, all targeting guides clear it together --
    so the true variance is far larger than binomial.

    Exactly: let R = #{targeting > max(NT)} under exchangeability of the n+k values. R >= r
    iff the top r of the pooled ranking are all targeting, so

        P(R >= r) = prod_{i=0}^{r-1} (n - i) / (n + k - i)

    a negative-hypergeometric tail. Its mean is n/(k+1) -- the binomial got the EXPECTATION
    right, which is why the error went unnoticed -- but for n=82, k=9 its SD is 7.9 against
    the binomial's 2.7. The published p-values were therefore ~100x too small:

        arm         k/n      binomial p   exact p
        A549/ABE     5/82      0.92         0.59
        A549/CBE     7/62      0.43         0.37
        T47D/ABE    18/82      0.0011       0.124
        T47D/CBE    20/62      0.00002      0.041

    `binom_p` is retained for traceability ONLY. Use `exact_p`. `perm_p` is a label
    permutation that makes no distributional assumption at all and agrees with `exact_p`
    to three decimals; `mwu_p` is the corresponding location test over the whole
    distribution rather than the tail.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for arm, a in m.groupby("arm"):
        nt, t = a[a.is_control], a[~a.is_control]
        if len(nt) < 5 or len(t) < 10:
            continue
        n, k_nt = len(t), len(nt)
        mx = float(nt[sev_col].max())
        k = int((t[sev_col] > mx).sum())

        # exact negative-hypergeometric tail: P(top k of the pooled ranking all targeting)
        exact = 1.0
        for i in range(k):
            exact *= (n - i) / (n + k_nt - i)

        # label permutation, assumption-free
        v = a[sev_col].to_numpy(float)
        hits = 0
        for _ in range(n_perm):
            idx = rng.permutation(n + k_nt)
            if (v[idx[:n]] > v[idx[n:]].max()).sum() >= k:
                hits += 1
        perm = (hits + 1) / (n_perm + 1)

        rows.append(dict(
            arm=arm, nt_max=mx, n_targeting=n, n_nt=k_nt, above_nt_max=k,
            expected=n / (k_nt + 1),
            exact_p=float(exact),
            perm_p=float(perm),
            mwu_p=float(mannwhitneyu(t[sev_col], nt[sev_col],
                                     alternative="greater").pvalue),
            binom_p=float(binomtest(k, n, 1 / (k_nt + 1),
                                    alternative="greater").pvalue),  # WRONG - traceability only
        ))
    return pd.DataFrame(rows)
