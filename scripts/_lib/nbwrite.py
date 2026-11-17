"""Build a .ipynb from a list of (kind, source) pairs."""
from __future__ import annotations
import json
from pathlib import Path

HEADER = '''%matplotlib inline
import sys, warnings
from pathlib import Path

sys.path.insert(0, str(next(p for p in Path.cwd().resolve().parents
                            if (p / "src" / "cdkbe").is_dir()) / "src"))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
'''


def md(src): return ("markdown", src)
def code(src): return ("code", src)


def write(path, cells, kernel="Python 3"):
    out = []
    for i, (kind, src) in enumerate(cells):
        src = src.rstrip("\n")
        # nbformat >=4.5 wants a stable id per cell; without one, execution warns
        c = dict(cell_type=kind, id=f"c{i:03d}", metadata={},
                 source=src.splitlines(keepends=True))
        if kind == "code":
            c.update(execution_count=None, outputs=[])
        out.append(c)
    nb = dict(cells=out,
              metadata=dict(kernelspec=dict(display_name=kernel, language="python", name="python3"),
                            language_info=dict(name="python", version="3.9")),
              nbformat=4, nbformat_minor=5)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    json.dump(nb, open(p, "w"), indent=1)
    return p
