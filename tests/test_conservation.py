"""The conservation modules must run from the packaged locations.

They were standalone scripts that rebuilt a ROOT from `__file__/../..`; moving
them into a package would silently break every path if `_paths` were wrong, and
the failure mode is a FileNotFoundError deep inside a long run. These tests are
cheap and catch that immediately.
"""
import os

import pytest

from cdkbe.conservation import _paths


def test_every_module_imports():
    import importlib
    import cdkbe.conservation as pkg
    here = os.path.dirname(pkg.__file__)
    mods = [f[:-3] for f in sorted(os.listdir(here))
            if f.endswith(".py") and not f.startswith("_")]
    assert len(mods) >= 10
    for m in mods:
        importlib.import_module(f"cdkbe.conservation.{m}")


@pytest.mark.parametrize("name", ["CONS", "BLAST", "ORTH", "ALN", "KLIFS"])
def test_data_directories_exist(name):
    assert os.path.isdir(getattr(_paths, name)), f"{name} -> {getattr(_paths, name)}"


def test_kinome_conservation_runs():
    from cdkbe.conservation import kinome_conservation as kc
    master, per_group, cdk_tab, mat, info = kc.main()
    assert mat.shape[1] == 85, "the KLIFS pocket is 85 positions"
    assert mat.shape[0] > 400, "expected the whole human kinome"
    assert len(master) == 85


def test_top_hgvsp_delegates_to_the_shared_helper():
    """The conservation copy of this reduction was the sixth; it must not drift."""
    from cdkbe import editing, io, paths
    from cdkbe.conservation import screen_conservation_integration as sci
    e = io.read(paths.EDITING / "ABE_subpool1_HGVSp_canonical.csv")
    a, b = sci.top_HGVSp_compiler(e), editing.top_hgvsp(e)
    assert a["HGVSp"].tolist() == b["HGVSp"].tolist()
