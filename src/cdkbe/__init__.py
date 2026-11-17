"""Shared code for the CDK base-editing screens.

    from cdkbe import paths, io, editing, lfc, plotting

`paths` resolves everything from the repository root, so a notebook works from
wherever it is run.
"""
from . import paths, io, editing, lfc, plotting   # noqa: F401

__all__ = ["paths", "io", "editing", "lfc", "plotting"]
__version__ = "0.1.0"
