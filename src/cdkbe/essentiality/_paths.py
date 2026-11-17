"""Where the essentiality and KLIFS-detail figure modules read and write.

The originals sat in `notebooks/09_analysis_and_viz/` and each rebuilt a
`ROOT` from `__file__/../..`, then joined `source_data/...`,
`screening_data/...` and `figures/` onto it. Those locations are resolved from
the repository root now, so the modules work from any working directory.

`DBNSFP` is the one location that is not in this repository: the dbNSFP dumps
are 259 MB and only `dbnsfp_recheck.py` reads them. See
`data/external/NOT_INCLUDED.csv`; that module raises a normal
FileNotFoundError naming the path.
"""
from __future__ import annotations
import os

from .. import paths

ROOT = str(paths.ROOT)
LIBRARY = str(paths.REFERENCE / "library")
LIBRARY_GENERATION = str(paths.REFERENCE / "library_generation")   # was source_data/01_library_generation_files         # was source_data/02_library
KLIFS_DIR = str(paths.REFERENCE / "KLIFS")             # was source_data/16_KLIFS
DISORDER = str(paths.REFERENCE / "disorder")       # was source_data/20_metapredict_disorder
ALPHAMISSENSE = str(paths.REFERENCE / "alphamissense")   # was source_data/12_alphamissense_scores
STRUCTURES = str(paths.ROOT / "scripts" / "chimerax" / "structures")
LFC_DIR = str(paths.LFC)                               # was screening_data/03_LFC_FDR_tables
EDITING = str(paths.EDITING)                       # was screening_data/04_editing
ESSENTIALITY = os.path.join(LFC_DIR, "T0_plasmid_essentiality")
FIG_DIR = str(paths.FIGURES)
RENDERS = str(paths.REFERENCE / "structure_renders")

#: not vendored -- see the module docstring
DBNSFP = str(paths.EXTERNAL / "dbNSFP")
