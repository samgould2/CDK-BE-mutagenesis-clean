#!/usr/bin/env python3
"""Execute notebooks in dependency order and report what each produced.

    python3 scripts/run_notebooks.py                  # the figure notebooks
    python3 scripts/run_notebooks.py --all            # every notebook
    python3 scripts/run_notebooks.py --inplace        # keep the rendered output

Without `--inplace` the executed copies go to a temporary directory, so running
this does not fill the repository with rendered output. With it, the notebooks
keep their figures and tables, which is what you want before sharing them.
"""
from __future__ import annotations
import argparse, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: dependency order; figures last, since they only read processed tables
ORDER = ["02_counts_and_qc", "03_editing", "04_lfc_fdr", "05_rnaseq",
         "06_structure_and_conservation", "07_other_analyses", "08_figures",
         "09_supplementary"]

#: the figure notebooks, run by default -- the processing ones above them are
#: slow and only need rerunning when the tables they build change
FIGURE_DIRS = ["08_figures", "09_supplementary"]


def notebooks(which: str):
    base = ROOT / "notebooks"
    dirs = ORDER if which == "all" else FIGURE_DIRS
    out = []
    for d in dirs:
        out += sorted((base / d).glob("*.ipynb"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--inplace", action="store_true")
    ap.add_argument("--timeout", type=int, default=1800)
    a = ap.parse_args()

    import nbformat
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    tmp = Path(tempfile.mkdtemp(prefix="nbrun_"))
    failures = 0
    for nb_path in notebooks("all" if a.all else "figures"):
        nb = nbformat.read(nb_path, as_version=4)
        client = NotebookClient(nb, timeout=a.timeout,
                                resources={"metadata": {"path": str(nb_path.parent)}})
        t = time.time()
        try:
            client.execute()
            status = "ok"
        except CellExecutionError as e:
            status = f"FAILED ({str(e).splitlines()[-1][:60]})"
            failures += 1
        images = sum(1 for c in nb.cells for o in c.get("outputs", [])
                     if "image/png" in o.get("data", {}))
        out = nb_path if a.inplace else tmp / nb_path.name
        nbformat.write(nb, out)
        print(f"  {nb_path.parent.name}/{nb_path.name:<26} {status:<10} "
              f"{images:>2} figure(s)  {time.time() - t:5.0f}s")
    if not a.inplace:
        print(f"\nexecuted copies in {tmp} (pass --inplace to keep them in the repo)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
