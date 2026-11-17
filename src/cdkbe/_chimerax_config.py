"""Re-export the ChimeraX structure configuration as an importable module.

`scripts/chimerax/cdk_config.py` is the single place the PDB/AlphaFold entry,
chain assignment and colour for each CDK are declared, and both the ChimeraX
generators and the essentiality figures read it. The generators run as plain
scripts from that directory, so the file stays there; this shim makes the same
declarations importable as `cdkbe._chimerax_config` without a second copy.
"""
from __future__ import annotations

import importlib.util
import sys

from . import paths

_spec = importlib.util.spec_from_file_location(
    "cdkbe._cdk_config_impl",
    paths.ROOT / "scripts" / "chimerax" / "cdk_config.py")
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)

globals().update({k: v for k, v in vars(_mod).items() if not k.startswith("_")})
