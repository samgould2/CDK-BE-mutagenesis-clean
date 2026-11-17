"""Foundation: held-out control normalization, a common gene universe, and an
n-unbiased perturbation statistic.

Three problems with the delivered pseudobulk, all fixed here.

1. LEAVE-ONE-OUT BIAS. The shipped `control_norm` layer Z-scores against the whole
   clean-control pool, and each NT guide's own cells are inside that pool. Its
   pseudobulk is shrunk toward zero by (1 - n_g/N_ctl) while targeting guides are
   not, which manufactures a ~14% "targeting > control" difference out of nothing.
   Here the control pool for guide g always excludes g's own cells.

2. INCOMPATIBLE GENE UNIVERSES. The per-arm sets shipped with the data span
   5,563-10,331 genes, a 1.9x spread on one platform. Cross-arm statistics need one
   universe; we take genes detected in >=5% of cells in every arm.

3. n-DEPENDENT MAGNITUDE. ||mean||^2 carries a sampling floor that scales as 1/n,
   and cells-per-guide varies 27-836 here, so naive magnitude ranks guides by how
   few cells they have. The cross-fitted statistic

       D(g) = <m_A, m_B> / G

   splits a guide's cells in half and takes the inner product of the two
   half-means. The noise in the two halves is independent, so the cross term has
   expectation zero and D is unbiased for ||mu||^2/G with no floor to subtract.
   For A549 the split is by gem_group, which also removes replicate-specific
   technical variance; T47D has one gem group per arm, so splits are random and
   averaged over repeats (unbiased, but it cannot remove clonal variance -- state
   that in any methods text).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import paths, singlet_io as io  # noqa: E402

MIN_DETECT_FRAC = 0.05
N_RANDOM_SPLITS = 50


def _cp10k_log(counts: sp.csr_matrix) -> sp.csr_matrix:
    tot = np.asarray(counts.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    out = counts.multiply(1e4 / tot[:, None]).tocsr()
    out.data = np.log1p(out.data)
    return out


def common_gene_universe(min_detect: float = MIN_DETECT_FRAC):
    """Genes detected in >= min_detect of cells in EVERY arm."""
    keep = None
    names = None
    for cl in ["A549", "T47D"]:
        obs = io.read_obs(cl)
        X = io.read_matrix(cl, layer="counts")
        names = io.gene_names(cl) if names is None else names
        binx = X.copy()
        binx.data = np.ones_like(binx.data)
        for arm, idx in obs.groupby("arm", observed=True).indices.items():
            frac = np.asarray(binx[idx].sum(axis=0)).ravel() / len(idx)
            k = frac >= min_detect
            keep = k if keep is None else (keep & k)
        del X, binx
    genes = names[keep]
    # Ribosomal / mitochondrial genes are short, highly expressed and depth-driven;
    # they dominate any length regression and add nothing to a perturbation signal.
    mask = ~pd.Series(genes).str.startswith(paths.EXCLUDE_PREFIXES).to_numpy()
    return genes[mask]


def arm_matrix(cell_line: str, arm: str, genes):
    """Dense log-CP10K matrix for one arm restricted to `genes`, plus its obs."""
    obs = io.read_obs(cell_line)
    sel = (obs["arm"].astype(str) == arm).to_numpy()
    X = io.read_matrix(cell_line, layer="counts")[sel]
    X = _cp10k_log(X)
    names = io.gene_names(cell_line)
    idx = pd.Series(np.arange(len(names)), index=names).reindex(genes).dropna().astype(int)
    M = np.asarray(X[:, idx.to_numpy()].todense(), dtype=np.float32)
    del X
    return obs[sel].copy(), M, np.asarray(idx.index)


def control_moments(M, gem, ctrl_mask, exclude_mask=None):
    """Per-gem-group control mean/SD, optionally holding out a set of cells.

    Only NT guides are ever inside the control pool, so for the 144 targeting
    guides the held-out moments are identical to the full-pool ones -- we exploit
    that below rather than recomputing per guide.
    """
    out = {}
    for g in np.unique(gem):
        pool = (gem == g) & ctrl_mask
        if exclude_mask is not None:
            pool = pool & ~exclude_mask
            if pool.sum() < 20:                       # guard: never leave too few
                pool = (gem == g) & ctrl_mask
        mu = M[pool].mean(axis=0)
        sd = M[pool].std(axis=0)
        sd[sd < 1e-8] = np.inf                        # zero-variance genes -> z = 0
        out[g] = (mu, sd)
    return out


def _z_mean(M, rows, gem, moments):
    """Mean z-profile of a set of cells, standardised within each gem group.

    Never materialises a full Z matrix: standardising the group mean is identical
    to averaging the standardised cells, and costs O(genes) instead of O(cells x genes).
    """
    idx = np.where(rows)[0]
    if len(idx) == 0:
        return None
    tot, acc = 0, None
    for g in np.unique(gem[idx]):
        sub = idx[gem[idx] == g]
        mu, sd = moments[g]
        v = (M[sub].mean(axis=0) - mu) / sd
        acc = v * len(sub) if acc is None else acc + v * len(sub)
        tot += len(sub)
    return acc / tot


def cross_fitted_stats(cell_line: str, arm: str, genes, rng_seed: int = 0):
    """Per-guide cross-fitted D and the pseudobulk profile, with held-out controls."""
    obs, M, genes_used = arm_matrix(cell_line, arm, genes)
    ctrl = obs["clean_control"].to_numpy(bool)
    gid = obs["guide_identity"].astype(str).to_numpy()
    gem = obs["gem_group"].astype(int).to_numpy()
    gems = np.unique(gem)
    rng = np.random.default_rng(rng_seed)
    G = M.shape[1]

    full = control_moments(M, gem, ctrl)          # reused by every targeting guide
    rows, profiles = [], {}
    for g in np.unique(gid):
        gm = gid == g
        n = int(gm.sum())
        if n < 10:
            continue
        # Holding a guide out only changes the pool if that guide IS a control.
        is_ctrl_guide = bool(ctrl[gm].any())
        mom = control_moments(M, gem, ctrl, gm) if is_ctrl_guide else full

        prof = _z_mean(M, gm, gem, mom)
        profiles[g] = prof

        if len(gems) > 1:
            a, b = gm & (gem == gems[0]), gm & (gem == gems[1])
            if a.sum() > 2 and b.sum() > 2:
                d = float(_z_mean(M, a, gem, mom) @ _z_mean(M, b, gem, mom) / G)
            else:
                d = np.nan
        else:
            idx = np.where(gm)[0]
            ds = []
            for _ in range(N_RANDOM_SPLITS):
                perm = rng.permutation(idx)
                h = len(perm) // 2
                if h < 3:
                    break
                ma = np.zeros(len(gm), bool); ma[perm[:h]] = True
                mb = np.zeros(len(gm), bool); mb[perm[h:2 * h]] = True
                ds.append(_z_mean(M, ma, gem, mom) @ _z_mean(M, mb, gem, mom) / G)
            d = float(np.mean(ds)) if ds else np.nan
        rows.append(dict(arm=arm, guide_identity=g, n_cells=n, D=d,
                         naive=float(prof @ prof / G), is_control=is_ctrl_guide))
    return pd.DataFrame(rows), pd.DataFrame(profiles, index=genes_used).T


def pseudo_guide_null(cell_line: str, arm: str, genes, sizes, n_draws: int = 200,
                      rng_seed: int = 1):
    """Null distribution of D from NT cells partitioned into synthetic guides.

    Sampling is matched to each real guide's cell count and gem-group composition,
    and runs through the identical held-out-normalisation and split-half pipeline.
    """
    obs, M, _ = arm_matrix(cell_line, arm, genes)
    ctrl = obs["clean_control"].to_numpy(bool)
    gem = obs["gem_group"].astype(int).to_numpy()
    gems = np.unique(gem)
    rng = np.random.default_rng(rng_seed)
    ctrl_idx = np.where(ctrl)[0]
    G = M.shape[1]

    out = []
    for n in sizes:
        n = int(min(n, len(ctrl_idx) // 2))
        for _ in range(n_draws):
            pick = rng.choice(ctrl_idx, n, replace=False)
            fake = np.zeros(len(obs), bool); fake[pick] = True
            mom = control_moments(M, gem, ctrl, fake)
            if len(gems) > 1:
                a, b = fake & (gem == gems[0]), fake & (gem == gems[1])
                if a.sum() < 3 or b.sum() < 3:
                    continue
                d = _z_mean(M, a, gem, mom) @ _z_mean(M, b, gem, mom) / G
            else:
                perm = rng.permutation(pick); h = len(perm) // 2
                if h < 3:
                    continue
                ma = np.zeros(len(obs), bool); ma[perm[:h]] = True
                mb = np.zeros(len(obs), bool); mb[perm[h:2*h]] = True
                d = _z_mean(M, ma, gem, mom) @ _z_mean(M, mb, gem, mom) / G
            out.append(dict(arm=arm, n_cells=n, D=float(d)))
    return pd.DataFrame(out)


if __name__ == "__main__":
    paths.ensure_dirs()
    genes = common_gene_universe()
    np.save(paths.RESULTS / "gene_universe.npy", genes)
    print(f"common gene universe: {len(genes)} genes detected in >={MIN_DETECT_FRAC:.0%} of cells in all 4 arms")

    stats, profs = [], {}
    for cl, arm in [("A549", "A549/ABE"), ("A549", "A549/CBE"),
                    ("T47D", "T47D/ABE"), ("T47D", "T47D/CBE")]:
        s, p = cross_fitted_stats(cl, arm, genes)
        stats.append(s)
        profs[arm] = p
        rho = s[["D", "n_cells"]].corr(method="spearman").iloc[0, 1]
        rho_naive = s[["naive", "n_cells"]].corr(method="spearman").iloc[0, 1]
        print(f"  {arm:11s} {len(s):3d} guides   rho(D, n)={rho:+.3f}   "
              f"[naive statistic: {rho_naive:+.3f}]")
    S = pd.concat(stats, ignore_index=True)
    S.to_csv(paths.RESULTS / "cross_fitted_D.csv", index=False)
    for arm, p in profs.items():
        # pyarrow/fastparquet are absent in scanpy_env; npz keeps this dependency-free
        np.savez_compressed(paths.RESULTS / f"profile_{arm.replace('/', '_')}.npz",
                            values=p.to_numpy().astype(np.float32),
                            guides=p.index.to_numpy().astype(str),
                            genes=p.columns.to_numpy().astype(str))
    print("wrote cross_fitted_D.csv and 4 profile matrices")


def load_profile(arm: str):
    """Guide x gene control-normalised profile matrix for one arm."""
    z = np.load(paths.RESULTS / f"profile_{arm.replace('/', '_')}.npz", allow_pickle=False)
    return pd.DataFrame(z["values"], index=z["guides"], columns=z["genes"])
