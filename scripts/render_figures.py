#!/usr/bin/env python3
"""Render every panel this repository can draw, named after the panel it is.

    python scripts/render_figures.py              # everything, to figures/manuscript/
    python scripts/render_figures.py 5 S12        # only figure 5 and supplementary 12
    python scripts/render_figures.py 5j --png     # one panel, with a PNG beside the PDF

Files are named for their place in the manuscript, not for the function that drew
them:

    figures/manuscript/main_fig3a.pdf
    figures/manuscript/main_fig5b_umap_target.pdf     a panel built from two parts
    figures/manuscript/main_fig4a_CDK2.pdf            a panel drawn once per gene
    figures/manuscript/supp_fig12c.pdf

A panel that is assembled in Illustrator from more than one drawn piece gets one
file per piece, suffixed with what that piece is. `figures/manuscript/INDEX.csv`
records every file written, which panel it belongs to, and the module and
function that drew it, so a file on disk can always be traced back to its code.

Panels listed in the registries but not drawn here are reported at the end rather
than skipped silently; `docs/COVERAGE.md` says what each is waiting on.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from cdkbe import paths, source_data as S  # noqa: E402

OUT = paths.FIGURES / "manuscript"

#: `pdf.fonttype` is read when the file is written, not when the figure is drawn,
#: so forcing it here is what keeps text live in Illustrator no matter which
#: module drew the panel. See the note in `source_data/rnaseq_t47d_figures.py`.
LIVE_TEXT = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}


def stem(key: str) -> str:
    """Registry key -> file stem. "3a" -> "main_fig3a"; "S12c" -> "supp_fig12c".

    The panel letter is lower-cased. Two supplementary-5 panels are keyed "S5A"
    and "S5B" in `supp_figures.SUPP_SPEC` while every other key in both
    registries -- and every panel letter in the manuscript's own figures -- is
    lower case, so without this the same figure would produce `supp_fig5A.pdf`
    beside `supp_fig5e.pdf`. The registry keys are left alone; only the file
    name is normalised. Flagged in docs/DISCREPANCIES.md for confirmation.
    """
    body = key[1:] if key.startswith("S") else key
    body = re.sub(r"\(.*\)$", "", body)      # "6(whole figure)" -> "6"
    num = body.rstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
    letter = body[len(num):].lower()
    return ("supp_fig" if key.startswith("S") else "main_fig") + num + letter


def figures_of(result, names=None):
    """Normalise what a drawing function returned to {suffix: figure}.

    The registries predate this script and return whatever suited each panel: a
    bare figure, a figure with its source rows beside it, or a dict keyed by gene
    or by sub-panel. All four shapes end up here as {suffix: figure}, with an
    empty suffix for a panel that is a single picture.
    """
    def one(v):
        if isinstance(v, Figure):
            return v
        if isinstance(v, (tuple, list)):
            for x in v:
                if isinstance(x, Figure):
                    return x
        return None

    f = one(result)
    if f is not None:
        return {"": f}
    if isinstance(result, dict):
        out = {}
        for k, v in result.items():
            g = one(v)
            if g is not None:
                out[str(k)] = g
        if out:
            # a panel that happens to have one part is still one picture, so it
            # keeps the bare `main_fig2k.pdf` rather than gaining a suffix
            return out if len(out) > 1 else {"": next(iter(out.values()))}
    if isinstance(result, (list, tuple)):
        # a positional list of panels, one per cohort/gene. `names` comes from the
        # drawing module's own PARTS table, never from the list's order alone.
        out = {}
        for i, v in enumerate(result):
            g = one(v)
            if g is None:
                continue
            out[names[i] if names and i < len(names) else str(i + 1)] = g
        if out:
            return out if len(out) > 1 else {"": next(iter(out.values()))}
    raise TypeError(f"no figure in a {type(result).__name__}")


def part_names(key: str, fn):
    """The names a drawing module gives the entries of a list it returns."""
    import importlib
    mod = importlib.import_module(fn.__module__) if fn.__module__ else None
    return (getattr(mod, "PARTS", {}) or {}).get(key)


def copies_of(result):
    """Normalise a result that is already-rendered FILES to {suffix: path}.

    The ChimeraX panels are images on disk rather than matplotlib figures --
    rendering one needs the application and an OpenGL context -- so they are
    copied under their manuscript name instead of being re-drawn and re-encoded.
    """
    paths_ = result[0] if isinstance(result, tuple) and result else result
    if isinstance(paths_, (str, Path)):
        paths_ = [paths_]
    if not isinstance(paths_, (list, tuple)) or not paths_:
        return None
    if not all(isinstance(p, (str, Path)) for p in paths_):
        return None
    out = {Path(p).stem: Path(p) for p in paths_}
    if len(out) < len(paths_):
        # one panel delivered in several formats -- the composited PNG and the
        # layered SVG of a leader panel share a stem, so key them by extension
        out = {Path(p).suffix.lstrip("."): Path(p) for p in paths_}
    return out if len(out) > 1 else {"": next(iter(out.values()))}


def copy_referenced(src: Path, outdir: Path) -> list:
    """Copy the files an SVG links to, so the copied SVG still renders.

    The leader panels are a raster structure with vector labels over it, and
    the SVG points at that raster by a relative name. Copy the SVG on its own
    and it comes out as labels floating on nothing.
    """
    if src.suffix.lower() != ".svg":
        return []
    import re as _re
    import shutil as _shutil

    out = []
    for href in _re.findall(r'xlink:href="([^"]+)"', src.read_text()):
        if "://" in href or href.startswith("data:"):
            continue
        sib = src.parent / href
        if sib.exists():
            dst = outdir / sib.name
            _shutil.copy2(sib, dst)
            out.append(dst)
    return out


def extra_renders(key: str, index) -> list:
    """Copy any ChimeraX images this panel uses, beside whatever it draws.

    Figure 1f is a plot and nine structure renders. The plot has a drawer; the
    renders do not, and without this they would be the only published artwork
    in the repository with no manuscript name.
    """
    import shutil

    from cdkbe.source_data import structures as _st

    if key in _st.FIGURES:            # the panel is renders only; handled there
        return []
    out = []
    for src in _st.images(key):
        dst = OUT / f"{stem(key)}_{src.stem}{src.suffix}"
        shutil.copy2(src, dst)
        out.append(dst)
        index.append(dict(panel=key, part=src.stem, file=dst.name,
                          module="scripts.chimerax", drawn_by="render.py"))
    return out


def render(key: str, fn, formats, index):
    result = fn()

    made = copies_of(result)
    if made is not None:
        import shutil

        written = []
        for suffix, src in made.items():
            # an extension key means the parts are formats of one panel, not
            # separate pictures, so they all keep the panel's own name
            name = stem(key) + ("" if suffix == src.suffix.lstrip(".")
                                else f"_{suffix}" if suffix else "")
            dst = OUT / f"{name}{src.suffix}"
            shutil.copy2(src, dst)
            written.append(dst)
            written += copy_referenced(src, OUT)
            index.append(dict(panel=key, part=suffix, file=dst.name,
                              module=fn.__module__ or "",
                              drawn_by=getattr(fn, "__qualname__", repr(fn))))
        return written

    figs = figures_of(result, part_names(key, fn))
    written = extra_renders(key, index)
    with matplotlib.rc_context(LIVE_TEXT):
        for suffix, fig in figs.items():
            name = stem(key) + (f"_{suffix}" if suffix else "")
            for ext in formats:
                p = OUT / f"{name}.{ext}"
                fig.savefig(p, dpi=400, bbox_inches="tight",
                            transparent=(ext != "png"))
                written.append(p)
            index.append(dict(panel=key, part=suffix, file=f"{name}.pdf",
                              module=fn.__module__ or "", drawn_by=getattr(
                                  fn, "__qualname__", repr(fn))))
            plt.close(fig)
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("only", nargs="*", help="figure numbers or panel keys, "
                    "e.g. 5, S12, 5j, S12c (default: everything)")
    ap.add_argument("--png", action="store_true", help="also write a PNG")
    ap.add_argument("--svg", action="store_true", help="also write an SVG")
    ap.add_argument("--outdir", type=Path, default=OUT)
    a = ap.parse_args()

    formats = ["pdf"] + (["svg"] if a.svg else []) + (["png"] if a.png else [])
    globals()["OUT"] = a.outdir
    OUT.mkdir(parents=True, exist_ok=True)

    registry = {**S.DRAWN, **S.SUPP_DRAWN}
    all_panels = set(S.PANELS) | set(S.SUPP_PANELS)

    def wanted(k: str) -> bool:
        if not a.only:
            return True
        return any(k == sel or k.rstrip("abcdefghijklmnopqrstuvwxyz") == sel
                   for sel in a.only)

    keys = sorted(k for k in registry if wanted(k))
    if not keys:
        print(f"nothing matches {a.only}", file=sys.stderr)
        return 2

    clash = {}
    for k in keys:
        clash.setdefault(stem(k), []).append(k)
    dupes = {s: v for s, v in clash.items() if len(v) > 1}
    if dupes:
        print("two registry keys would write the same file: "
              + "; ".join(f"{s} <- {v}" for s, v in dupes.items()), file=sys.stderr)
        return 2

    index, failed, n = [], [], 0
    for k in keys:
        try:
            files = render(k, registry[k], formats, index)
            n += len(files)
            shown = [p.name for p in files if p.suffix in (".pdf", ".png", ".svg")]
            print(f"  {k:6s} -> {', '.join(shown)}")
        except Exception as exc:
            failed.append((k, exc))
            print(f"  {k:6s} -> FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
            traceback.print_exc(limit=3, file=sys.stderr)

    if index:
        ip = OUT / "INDEX.csv"
        # A partial run (`render_figures.py 4k 4l`) must not shrink the index
        # to the panels it happened to draw -- the other files are still on
        # disk, and the index is what traces each of them back to its code.
        # Rows this run produced replace the matching old rows; the rest stay.
        if ip.exists():
            fresh = {(r["panel"], r["part"]) for r in index}
            with open(ip, newline="") as fh:
                index += [r for r in csv.DictReader(fh)
                          if (r["panel"], r["part"]) not in fresh
                          and (OUT / r["file"]).exists()]
        with open(ip, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["panel", "part", "file",
                                               "module", "drawn_by"])
            w.writeheader()
            w.writerows(sorted(index, key=lambda r: (r["panel"].startswith("S"),
                                                     r["panel"], r["part"])))
        print(f"\n{n} file(s) written; {len(index)} listed in {ip.name}")

    undrawn = sorted(k for k in all_panels - set(registry) if wanted(k))
    if undrawn:
        print(f"\nnot drawn here ({len(undrawn)}): {', '.join(undrawn)}"
              "\n  these have source data but no plotting code; "
              "see docs/DISCREPANCIES.md item 13")
    if failed:
        print(f"\n{len(failed)} FAILED: {', '.join(k for k, _ in failed)}",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
