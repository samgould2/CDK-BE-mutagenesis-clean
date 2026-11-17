"""Where the conservation modules read and write.

The originals each rebuilt a `ROOT` from `__file__/../..` and then joined
`source_data/22_conservation`, `screening_data/...` and `figures/` onto it.
Those locations are resolved from the repository root now, so the modules work
from any working directory and a moved file does not break them.
"""
from __future__ import annotations
import os

from .. import paths

CONS = str(paths.REFERENCE / "conservation")     # was source_data/22_conservation
BLAST = os.path.join(CONS, "blast")
ORTH = os.path.join(CONS, "orthologs")
ALN = os.path.join(CONS, "alignments")
OUT = CONS                                       # the originals wrote back here
KLIFS = str(paths.REFERENCE / "KLIFS")           # was source_data/16_KLIFS
LFC = str(paths.LFC)                             # was screening_data/03_LFC_FDR_tables
EDITING = str(paths.EDITING)                     # was screening_data/04_editing
FIG = str(paths.FIGURES)
ROOT = str(paths.ROOT)

for _d in (ORTH, ALN):
    os.makedirs(_d, exist_ok=True)
