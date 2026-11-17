"""Canonical locations inside the repository.

Everything is resolved from the repo root, so notebooks work regardless of the
directory they are launched from -- the old tree relied on `../../` prefixes,
which broke whenever a notebook moved.
"""
from __future__ import annotations
import os
from pathlib import Path

def repo_root() -> Path:
    env = os.environ.get("CDKBE_ROOT")
    if env:
        return Path(env).resolve()
    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "data").is_dir() and (p / "src" / "cdkbe").is_dir():
            return p
    return here.parents[2]

ROOT       = repo_root()
DATA       = ROOT / "data"
RAW        = DATA / "raw"            # as close to the sequencer as this repo goes
PROCESSED  = DATA / "processed"      # everything derived, regenerable by the notebooks
REFERENCE  = DATA / "reference"      # external annotation: KLIFS, library design, gene sets
EXTERNAL   = DATA / "external"       # manifest + fetch script for data too large to ship

COUNTS       = PROCESSED / "count_tables"
RPM          = PROCESSED / "rpm_tables"
LFC          = PROCESSED / "lfc_fdr_tables"
EDITING      = PROCESSED / "editing"
RAW_COUNTS   = RAW / "screen_counts"
RAW_EDITING  = RAW / "screen_editing"

FIGURES      = ROOT / "figures"
FIG_MAIN     = FIGURES / "main"
FIG_SUPP     = FIGURES / "supplementary"
SOURCE_DATA  = ROOT / "source_data_tables"

def ensure(*paths: Path) -> None:
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)


def rel(p) -> str:
    """`p` written relative to the repository root, for printing.

    A message naming an absolute path is baked into the output of any notebook
    that prints it, and that path is only right on the machine it ran on.
    """
    try:
        return str(Path(p).relative_to(ROOT))
    except ValueError:
        return str(p)
