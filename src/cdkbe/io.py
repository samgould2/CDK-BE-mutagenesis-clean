"""Reading that does not care how a table was stored.

Large CSVs in this repository are shipped gzipped (roughly 7x smaller). `read`
accepts a logical path without an extension suffix and finds whichever of
`.csv`, `.csv.gz` or `.zip` is actually present, so the notebooks never mention
compression and re-compressing a file does not require touching any code.
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable
import pandas as pd

_SUFFIXES = ("", ".gz", ".zip")

def resolve(path) -> Path:
    """Find the real file behind a logical path."""
    p = Path(path)
    if p.exists():
        return p
    for s in _SUFFIXES:
        q = Path(str(p) + s)
        if q.exists():
            return q
    # caller gave `foo.csv.gz` but only `foo.csv` exists (or vice versa)
    if p.suffix in (".gz", ".zip"):
        q = p.with_suffix("")
        if q.exists():
            return q
    # `foo.csv` asked for, but the file on disk is `foo.zip` / `foo.gz`
    for s in (".zip", ".gz"):
        q = p.with_suffix(s)
        if q.exists():
            return q
    raise FileNotFoundError(f"no file for {path} (tried {', '.join(repr(str(p)+s) for s in _SUFFIXES)})")

def read(path, **kw) -> pd.DataFrame:
    """read_csv that transparently handles .csv / .csv.gz / .zip."""
    return pd.read_csv(resolve(path), **kw)

def read_excel(path, sheet=0, **kw) -> pd.DataFrame:
    return pd.read_excel(resolve(path), sheet_name=sheet, **kw)

#: a written table larger than this is gzipped (matches how data/ was built)
GZIP_MIN_BYTES = 1_000_000


def write(df: pd.DataFrame, path, compress: bool | None = None, **kw) -> Path:
    """Write a table, gzipping it when the result is big enough to be worth it.

    The decision is made on the size of the written file, not on row count: the
    RPM tables are only a few thousand rows but are wide and full of floats, so
    a row-count rule leaves megabyte CSVs uncompressed.

    Writing `foo.csv` where `foo.csv.gz` already exists (or the reverse) would
    leave two files that `read` could pick between, so the stale one is
    removed. That is the only file this function ever deletes.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    index = kw.pop("index", False)
    plain = Path(str(p)[:-3]) if str(p).endswith(".gz") else p

    df.to_csv(plain, index=index, **kw)
    if compress is None:
        compress = plain.stat().st_size >= GZIP_MIN_BYTES
    if compress:
        import gzip, shutil
        gz = Path(str(plain) + ".gz")
        with open(plain, "rb") as a, gzip.open(gz, "wb", compresslevel=6) as b:
            shutil.copyfileobj(a, b, 1 << 20)
        plain.unlink()
        return gz
    stale = Path(str(plain) + ".gz")
    if stale.exists():
        stale.unlink()
    return plain


def exists(path) -> bool:
    try:
        resolve(path); return True
    except FileNotFoundError:
        return False

def concat_read(paths: Iterable, **kw) -> pd.DataFrame:
    return pd.concat([read(p, **kw) for p in paths], ignore_index=True)
