"""Evolutionary conservation of the kinase pocket.

Ported from the standalone scripts that sat beside the original 07_KLIFS notebook.
They were already modules; what changed is that paths now resolve from the
repository root (see `_paths.py`) and the cross-module imports are relative,
so they work from any working directory.

    from cdkbe.conservation import kinome_conservation as kc
    master, per_group, cdk_tab, mat, info = kc.main()
"""
