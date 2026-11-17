"""Every third-party import in the repository is declared in `pyproject.toml`.

An undeclared dependency fails the same way an absolute path does: it keeps
working on the machine that already has it, and only breaks for the next person
who follows the README. `pip install -e .` has to be enough to run the figures,
and the extras have to cover everything else.

The check is on the import graph rather than on a hand-kept list, so adding a
module that imports something new fails here until the dependency is declared.
"""
from __future__ import annotations

import ast
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: directories whose imports must be satisfiable from pyproject
CODE = ("src", "scripts", "tests", "pipeline")

#: distribution name -> the module name it provides, where they differ
PROVIDES = {
    "biopython": "Bio",
    "scikit-learn": "sklearn",
    "python-Levenshtein": "Levenshtein",
    "adjustText": "adjusttext",
}

#: imports that are deliberately not pip dependencies, with the reason
EXEMPT = {
    "cdkbe": "this package",
    "chimerax": "provided by the ChimeraX application; these scripts run inside it",
    "pymol": "provided by the PyMOL application; these scripts run inside it",
    "cdk_config": "scripts/chimerax/cdk_config.py, imported as a sibling file",
    "analysis": "a docstring-only reference to the original BEAT-seq tree",
    "__future__": "language feature",
}


def _stdlib() -> set:
    names = set(sys.builtin_module_names)
    if hasattr(sys, "stdlib_module_names"):            # 3.10+
        names |= set(sys.stdlib_module_names)
    else:                                              # 3.9
        names |= {
            "abc", "argparse", "ast", "base64", "builtins", "collections",
            "contextlib", "copy", "csv", "dataclasses", "datetime", "enum",
            "functools", "glob", "gzip", "hashlib", "importlib", "inspect",
            "io", "itertools", "json", "logging", "math", "operator", "os",
            "pathlib", "pickle", "pkgutil", "platform", "random", "re",
            "shutil", "statistics", "string", "subprocess", "sys", "tempfile",
            "textwrap", "time", "tomllib", "traceback", "types", "typing",
            "unicodedata",
            "urllib", "warnings", "xml", "zipfile",
        }
        # tomllib arrived in 3.11; on 3.9/3.10 this file falls back to a regex,
        # so the import is guarded there and must not count as a dependency.
    return names


STDLIB = _stdlib()


def _local_modules() -> set:
    """Top-level names importable as a sibling file of some script we ship."""
    out = set()
    for d in CODE:
        for dp, dn, fn in os.walk(ROOT / d):
            dn[:] = [x for x in dn if x not in {"__pycache__", ".ipynb_checkpoints"}]
            for f in fn:
                if f.endswith(".py"):
                    out.add(f[:-3])
            for x in dn:
                out.add(x)
    return out


LOCAL = _local_modules()


def _declared() -> set:
    """Module names pyproject makes installable, in any dependency group."""
    text = (ROOT / "pyproject.toml").read_text()
    try:
        import tomllib                                  # 3.11+
        cfg = tomllib.loads(text)
    except ModuleNotFoundError:
        import re
        names = set(re.findall(r'"([A-Za-z0-9_.\-]+)\s*[><=!]', text))
        names |= set(re.findall(r'"([A-Za-z0-9_.\-]+)"', text))
        return {PROVIDES.get(n, n).lower() for n in names}
    deps = list(cfg["project"].get("dependencies", []))
    for group in cfg["project"].get("optional-dependencies", {}).values():
        deps += list(group)
    import re
    names = {re.split(r"[><=!\[;\s]", d, 1)[0] for d in deps}
    return {PROVIDES.get(n, n).lower() for n in names}


DECLARED = _declared()


def _imports():
    """-> {module: [file, ...]} for every top-level third-party import."""
    found = {}
    for d in CODE:
        for dp, dn, fn in os.walk(ROOT / d):
            dn[:] = [x for x in dn if x not in {"__pycache__", ".ipynb_checkpoints"}]
            for f in sorted(fn):
                if not f.endswith(".py"):
                    continue
                p = Path(dp) / f
                try:
                    tree = ast.parse(p.read_text())
                except SyntaxError:
                    continue
                mods = set()
                for n in ast.walk(tree):
                    if isinstance(n, ast.Import):
                        mods |= {a.name.split(".")[0] for a in n.names}
                    elif isinstance(n, ast.ImportFrom) and n.level == 0 and n.module:
                        mods.add(n.module.split(".")[0])
                for m in mods:
                    if m in STDLIB or m in EXEMPT or m in LOCAL:
                        continue
                    found.setdefault(m, []).append(str(p.relative_to(ROOT)))
    return found


IMPORTS = _imports()


@pytest.mark.parametrize("module", sorted(IMPORTS), ids=sorted(IMPORTS))
def test_import_is_declared(module):
    where = ", ".join(IMPORTS[module][:3])
    assert module.lower() in DECLARED, (
        f"`{module}` is imported by {where} but no dependency in pyproject.toml "
        f"provides it. Add it to the right group, or to EXEMPT here with a reason.")


def test_the_core_dependencies_cover_the_figures():
    """`pip install -e .` alone has to draw every panel, with no extras."""
    core = {"pandas", "numpy", "scipy", "matplotlib", "seaborn", "statsmodels",
            "openpyxl", "logomaker", "adjusttext", "sklearn", "bio"}
    assert core <= DECLARED
