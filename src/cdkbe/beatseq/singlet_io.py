"""Fast loaders for the CDK BEAT-seq singlet objects.

Named `singlet_io` rather than `io` as in the original tree: inside this package
`from . import io` would shadow `cdkbe.io`, the compression-agnostic table
reader every other module here uses, and the two do completely different jobs.

Everything in this module opens one of the ~34 GB singlet `.h5ad` objects, which
are NOT in this repository (`data/external/NOT_INCLUDED.csv`). Nothing that
draws a figure or recomputes a shipped statistic reaches it: the derived tables
those need are all copied in. It is here so the upstream derivation is readable,
and so re-deriving from the raw objects is a matter of fetching them.

Design note: we never materialise the 15 GB norm_pop objects. Reading the raw
counts out of the singlet h5ad with h5py and recomputing normalization on the fly
takes a few seconds, and it is the only way to get a *held-out* control pool
(see nulls.held_out_control_norm) -- the shipped `control_norm` layer includes each
NT guide's own cells in the pool used to centre it, which shrinks NT pseudobulks
toward zero and fakes a targeting-vs-control difference.
"""

from __future__ import annotations

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import paths


def _decode(arr: np.ndarray) -> np.ndarray:
    """h5py returns variable-length strings as object arrays of bytes."""
    arr = np.asarray(arr)
    if arr.dtype.kind in "SO" and arr.size and isinstance(arr.flat[0], bytes):
        return np.array([v.decode() for v in arr.ravel()]).reshape(arr.shape)
    if arr.dtype.kind == "S":
        return arr.astype(str)
    return arr


def _read_dataframe(grp: h5py.Group) -> pd.DataFrame:
    """Decode an AnnData-encoded obs/var group (handles categoricals + nullable)."""
    index_key = grp.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    cols = {}
    for key in grp.keys():
        if key == index_key:
            continue
        node = grp[key]
        if isinstance(node, h5py.Group):
            enc = node.attrs.get("encoding-type", "")
            if enc == "categorical":
                cats = _decode(node["categories"][:])
                codes = node["codes"][:]
                cols[key] = pd.Categorical.from_codes(codes, categories=cats)
            elif enc in ("nullable-integer", "nullable-boolean"):
                vals = node["values"][:]
                mask = node["mask"][:] if "mask" in node else None
                arr = pd.array(vals, dtype="Int64" if enc == "nullable-integer" else "boolean")
                if mask is not None:
                    arr[mask] = pd.NA
                cols[key] = arr
        else:
            cols[key] = _decode(node[:])
    idx = _decode(grp[index_key][:])
    return pd.DataFrame(cols, index=pd.Index(idx, name="cell_barcode"))


def read_obs(cell_line: str) -> pd.DataFrame:
    """Cell metadata only -- no matrix touched."""
    with h5py.File(paths.SINGLETS[cell_line], "r") as f:
        return _read_dataframe(f["obs"])


def read_obsm(cell_line: str, key: str = "X_umap") -> np.ndarray:
    """A cell embedding from obsm (e.g. 'X_umap', 'X_pca') without touching the matrix.

    The delivered objects carry X_umap and X_pca computed by the original pipeline. Rows
    are in the same order as `read_obs`, so boolean masks built from obs apply directly.
    """
    with h5py.File(paths.SINGLETS[cell_line], "r") as f:
        return np.asarray(f["obsm"][key][:])


def read_var(cell_line: str) -> pd.DataFrame:
    with h5py.File(paths.SINGLETS[cell_line], "r") as f:
        return _read_dataframe(f["var"])


def gene_names(cell_line: str) -> np.ndarray:
    with h5py.File(paths.SINGLETS[cell_line], "r") as f:
        idx = f["var"][f["var"].attrs.get("_index", "_index")][:]
    return _decode(idx)


def read_matrix(cell_line: str, layer: str | None = "counts") -> sp.csr_matrix:
    """Load the full CSR matrix.

    layer='counts' gives raw CellRanger counts; layer=None gives X, which is
    log1p(CP10000). A549 is 84,466 x 38,606 with 292M nonzeros (~3.5 GB in memory);
    T47D is 17,339 x 38,606 with 84M nonzeros. Both fit comfortably in 64 GB.
    """
    with h5py.File(paths.SINGLETS[cell_line], "r") as f:
        grp = f["X"] if layer is None else f["layers"][layer]
        shape = tuple(grp.attrs["shape"])
        mat = sp.csr_matrix(
            (grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape
        )
    return mat


def read_genes(cell_line: str, genes, layer: str | None = "counts") -> pd.DataFrame:
    """Read a handful of named genes without densifying the whole matrix.

    Used by the NMD readout, which needs only the 9 CDK transcripts. Reads the CSR
    once and slices columns; ~10 s for A549, seconds for T47D.
    """
    names = gene_names(cell_line)
    lookup = pd.Series(np.arange(len(names)), index=names)
    genes = [g for g in genes if g in lookup.index]
    if not genes:
        return pd.DataFrame(index=read_obs(cell_line).index)
    cols = lookup.loc[genes].to_numpy()
    mat = read_matrix(cell_line, layer=layer)
    sub = np.asarray(mat[:, cols].todense())
    del mat
    obs_index = read_obs(cell_line).index
    return pd.DataFrame(sub, columns=genes, index=obs_index)


def guide_annotation() -> pd.DataFrame:
    """Per-(guide, arm) annotation joined from the QC tables.

    Keyed on (id, editor) -- never on sgRNA or id alone, because four targeting
    protospacers are shared across editor arms while installing different alleles.
    """
    ue = pd.read_csv(paths.T_UNTREATED_EXPECT)
    ge = pd.read_csv(paths.T_GUIDE_EDITING)
    ann = ue.merge(
        ge[["id", "arm", "low_edit"]], on=["id", "arm"], how="left", validate="1:1"
    )
    ann = ann.rename(columns={"arm": "editor"})
    ann["variant_class"] = np.where(
        ann["Mutation"].astype(str).str.contains(r"\*", regex=True), "nonsense",
        np.where(ann["Mutation"].astype(str).str.match(r"^M1[A-Z]"), "start_loss", "missense"),
    )
    return ann


def load_arm(cell_line: str, layer: str | None = "counts"):
    """Return (obs, matrix, gene_names) for one cell line."""
    return read_obs(cell_line), read_matrix(cell_line, layer=layer), gene_names(cell_line)
