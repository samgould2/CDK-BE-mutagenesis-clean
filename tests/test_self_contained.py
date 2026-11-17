"""Nothing in the repository may reach outside it for data.

The point of this repository is that someone who clones it can rebuild every
figure without having any of the original trees on their machine.  That breaks
silently: an absolute path to `~/Documents/GitHub/CDK-BE-mutagenesis-v2` keeps
working on the machine it was written on, and only fails for the next person.
So the check is mechanical rather than a convention.

One legitimate exception, explicit: `cdkbe.beatseq.paths._ABSENT` documents
the ~34 GB of .h5ad objects that are deliberately not vendored.  Nothing a
figure touches resolves through it.
"""
import ast
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: directories whose code must resolve only inside the repository
CODE = ("src", "tests", "scripts", "notebooks", "pipeline")

#: files allowed to name an outside tree, with the reason
ALLOWED = {
    "src/cdkbe/beatseq/paths.py": "_ABSENT documents the un-vendored .h5ad objects",
    "tests/test_self_contained.py": "this file, whose patterns are the outside paths",
}

#: a path that leaves the repository: a home directory, a tree we migrated from,
#: or an absolute path that is not /tmp
OUTSIDE = re.compile(
    r"""(?x)
    /Users/ | ~/Documents | \$HOME/
    | CDK-BE-mutagenesis-v2
    | CDK-BEAT-seq/CDK_Beatseq       # the delivered bundle, not our copy of it
    | (?<![\w./])screening_data/     # the old repository's data root
    """
)


def _files():
    for d in CODE:
        for p in sorted((ROOT / d).rglob("*")):
            if p.suffix not in {".py", ".ipynb", ".cxc", ".sh"}:
                continue
            parts = p.relative_to(ROOT).parts
            if "/.ipynb_checkpoints/" in str(p) or "__pycache__" in str(p):
                continue
            # Generated ChimeraX command files: they are gitignored, and they
            # must carry an absolute `cd` to this checkout, because ChimeraX
            # silently aborts a render on a relative path.
            if len(parts) > 2 and parts[1] == "chimerax" and parts[2].startswith("cxc"):
                continue
            yield p


def _strings(path: Path):
    """Every string literal in the file, with its line number.

    Comments and docstrings are skipped on purpose: a comment saying a table
    "was screening_data/03_LFC_FDR_tables" is documentation, not a path the
    code resolves.  Only literals that could actually be opened are checked.
    """
    text = path.read_text(errors="replace")
    if path.suffix == ".ipynb":
        import json

        try:
            nb = json.loads(text)
        except json.JSONDecodeError:
            return
        for i, cell in enumerate(nb.get("cells", [])):
            if cell.get("cell_type") != "code":
                continue
            src = "".join(cell.get("source", []))
            for ln, line in enumerate(src.splitlines(), 1):
                if line.lstrip().startswith("#"):
                    continue
                yield f"cell {i}, line {ln}", line
        return
    if path.suffix in {".cxc", ".sh"}:
        for ln, line in enumerate(text.splitlines(), 1):
            if line.lstrip().startswith("#"):
                continue
            yield f"line {ln}", line
        return
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return
    # Collect the docstring *Constant* nodes, not the Expr wrapping them:
    # ast.walk yields the child independently, so skipping the Expr alone
    # still lets the docstring text through.
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if ast.get_docstring(node, clean=False) is not None and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr):
                    docstrings.add(id(first.value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if id(node) in docstrings:
                continue
            yield f"line {node.lineno}", node.value


@pytest.mark.parametrize("path", list(_files()), ids=lambda p: str(p.relative_to(ROOT)))
def test_no_path_outside_the_repository(path):
    rel = str(path.relative_to(ROOT))
    if rel in ALLOWED:
        pytest.skip(ALLOWED[rel])
    hits = [(where, s.strip()[:120]) for where, s in _strings(path) if OUTSIDE.search(s)]
    assert not hits, "\n".join(f"  {rel} {w}: {s}" for w, s in hits)


def test_every_declared_data_path_exists():
    """`cdkbe.paths` and `cdkbe.beatseq.paths` must resolve inside the tree."""
    import cdkbe.paths as P
    from cdkbe.beatseq import paths as B

    missing = []
    for mod, names in ((P, [n for n in dir(P) if n.isupper()]),
                       (B, [n for n in dir(B) if n.startswith("T_")])):
        for n in names:
            v = getattr(mod, n)
            if isinstance(v, Path) and not v.exists():
                missing.append(f"{mod.__name__}.{n} -> {v}")
    assert not missing, "declared paths that do not exist:\n  " + "\n  ".join(missing)


def test_absent_objects_are_declared():
    """The un-vendored .h5ad objects stay listed, so their absence is on record."""
    note = ROOT / "data" / "external" / "NOT_INCLUDED.csv"
    assert note.is_file(), "data/external/NOT_INCLUDED.csv is missing"
    body = note.read_text()
    assert "h5ad" in body, "NOT_INCLUDED.csv no longer records the single-cell objects"


def _executed_notebooks():
    return sorted((ROOT / "notebooks").rglob("*.ipynb"))


def _output_text(cell):
    """Every string a cell's stored output would show a reader."""
    for out in cell.get("outputs", []):
        kind = out.get("output_type")
        if kind == "stream":
            yield "".join(out.get("text", []))
        elif kind in ("execute_result", "display_data"):
            data = out.get("data", {})
            for mime in ("text/plain", "text/html", "image/svg+xml"):
                if mime in data:
                    yield "".join(data[mime])
        elif kind == "error":
            yield "\n".join(out.get("traceback", []))


@pytest.mark.parametrize("path", _executed_notebooks(),
                         ids=lambda p: str(p.relative_to(ROOT)))
def test_notebook_output_names_no_outside_path(path):
    """A notebook ships with its figures, so it ships with its printed text too.

    An absolute path printed by a drawer is then baked into the file, correct
    only on the machine it ran on.  `cdkbe.paths.rel` is how a module prints a
    location; a hit here means something printed a raw path instead.
    """
    import json

    nb = json.loads(path.read_text())
    hits = []
    for i, cell in enumerate(nb.get("cells", [])):
        for text in _output_text(cell):
            for line in text.splitlines():
                if OUTSIDE.search(line):
                    hits.append(f"cell {i}: {line.strip()[:120]}")
    rel = path.relative_to(ROOT)
    assert not hits, f"{rel} stores an outside path in its output:\n  " + "\n  ".join(hits[:10])


def test_notebooks_store_no_errors():
    """An executed notebook that ships must have run clean."""
    import json

    failed = []
    for path in _executed_notebooks():
        nb = json.loads(path.read_text())
        errs = [o for c in nb.get("cells", []) for o in c.get("outputs", [])
                if o.get("output_type") == "error"]
        if errs:
            failed.append(f"{path.relative_to(ROOT)}: {len(errs)} "
                          f"({errs[0].get('ename')}: {str(errs[0].get('evalue'))[:60]})")
    assert not failed, "notebooks stored with an error output:\n  " + "\n  ".join(failed)
