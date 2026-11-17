"""Read the plotted values back out of a matplotlib-written SVG.

Some panels are available only as the published figure. Where that figure is
a vector file written by matplotlib, the numbers are still in it: the axes are drawn in figure points at full
precision, each tick is a mark at a known position with its label spelled out
next to it, and every mark, bar and line vertex is an absolute coordinate in
that same space. Two ticks fix the affine map from points to data units, and
every other number on the page follows from it.

That makes this an exact inversion, not an estimate, for anything drawn
against a labelled linear axis -- bar heights, line vertices, marker centres.
It is **not** a way to read a heatmap: there the value is encoded in the fill
colour, quantised to eight bits per channel and pushed through a colormap, and
what comes back is an approximation of a derived quantity. Nothing here reads
colour.

`tests/test_svg_values.py` checks the inversion against the four T47D
cell-cycle panels, whose published SVGs *and* underlying FlowJo export are both
in this repository: the recovered bar heights equal the measured gate
frequencies exactly. That is the warrant for using it on the panels this
repository holds only as a figure.

Two spellings of text exist in matplotlib's SVG output and both are handled:
real `<text>` elements (`svg.fonttype: "none"`) and references into a `<defs>
of glyph outlines (the default), whose ids end in the character's hex
codepoint.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

SVG = "{http://www.w3.org/2000/svg}"
XLINK = "{http://www.w3.org/1999/xlink}href"

#: a glyph id such as `Helvetica-4b` or `DejaVuSans-2212`
GLYPH = re.compile(r"-([0-9a-f]{2,4})$")
NUMBER = re.compile(r"-?\d+\.?\d*(?:[eE]-?\d+)?")


def load(source) -> ET.Element:
    """Parse a path, or anything else ElementTree can read (a file object)."""
    if isinstance(source, (str, Path)):
        return ET.parse(str(source)).getroot()
    return ET.parse(source).getroot()


def text_of(group: ET.Element) -> str:
    """The characters a matplotlib text group spells out."""
    real = "".join("".join(t.itertext()) for t in group.iter(SVG + "text"))
    if real.strip():
        return real
    out = []
    for use in group.iter(SVG + "use"):
        m = GLYPH.search(use.get(XLINK, ""))
        if m:
            out.append(chr(int(m.group(1), 16)))
    return "".join(out)


def _mark_xy(group: ET.Element):
    """Where a tick's mark sits: the first `<use>` under it, with transforms."""
    for node in group.iter():
        tr = node.get("transform")
        m = tr and re.match(r"translate\(([-\d.]+)[ ,]+([-\d.]+)\)", tr)
        if m:
            use = node.find(".//" + SVG + "use")
            dx, dy = float(m.group(1)), float(m.group(2))
            if use is None:
                return dx, dy
            return dx + float(use.get("x", 0)), dy + float(use.get("y", 0))
    use = group.find(".//" + SVG + "use")
    if use is None:
        return None
    return float(use.get("x", 0)), float(use.get("y", 0))


def ticks(root: ET.Element, axis: str) -> list:
    """[(position in points, label value)] for every numeric tick on an axis."""
    want = "xtick_" if axis == "x" else "ytick_"
    out = []
    for group in root.iter(SVG + "g"):
        if not (group.get("id") or "").startswith(want):
            continue
        mark = label = None
        for sub in group:
            sid = sub.get("id") or ""
            if sid.startswith("line2d_"):
                mark = _mark_xy(sub)
            elif sid.startswith("text_"):
                label = text_of(sub)
        if mark is None or not label:
            continue
        try:                       # U+2212 MINUS SIGN is what matplotlib writes
            value = float(label.replace("−", "-").strip())
        except ValueError:
            continue
        out.append((mark[0] if axis == "x" else mark[1], value))
    return sorted(out)


def scale(root: ET.Element, axis: str):
    """-> f(points) = data units, from the axis's own tick labels.

    Raises if the axis carries fewer than two numeric ticks, or if any
    intermediate tick does not land on the line the outer two define -- a log
    or symlog axis, where nothing this module returns would mean anything.
    """
    found = ticks(root, axis)
    if len(found) < 2:
        raise ValueError(f"{axis} axis: need two numeric ticks, found {len(found)}")
    (p0, v0), (p1, v1) = found[0], found[-1]
    if p1 == p0:
        raise ValueError(f"{axis} axis: ticks are on top of each other")
    slope = (v1 - v0) / (p1 - p0)
    for p, v in found:
        if abs((v0 + slope * (p - p0)) - v) > 1e-6 * max(1.0, abs(v)):
            raise ValueError(f"{axis} axis is not linear: tick {v} is off the line")
    return lambda p: v0 + slope * (p - p0)


def vertices(d: str) -> list:
    """Every coordinate pair in a path's `d`, in the order it was written."""
    nums = [float(x) for x in NUMBER.findall(d)]
    return list(zip(nums[0::2], nums[1::2]))


def paths(root: ET.Element, prefix: str = "") -> dict:
    """{group id: (style, [(x, y), ...])} for each group holding one path."""
    out = {}
    for group in root.iter(SVG + "g"):
        gid = group.get("id") or ""
        if prefix and not gid.startswith(prefix):
            continue
        path = group.find(SVG + "path")
        if path is None or not path.get("d"):
            continue
        out[gid] = (path.get("style") or "", vertices(path.get("d")))
    return out


def stroke(style: str):
    """The stroke colour a style string sets, or None."""
    m = re.search(r"stroke:\s*(#[0-9a-fA-F]{6})", style)
    return m.group(1).lower() if m else None


def bars(root: ET.Element, axes: str = "axes_1") -> list:
    """[(left edge, height in data units)] for each axis-aligned bar, left to right.

    A bar is a *stroked rectangle*: four or more vertices over exactly two
    distinct x and two distinct y values, with a stroke. That excludes the
    axes background (a rectangle, but unstroked) and the spines (stroked, but
    two vertices over one x). What survives is then kept only if it stands on
    the baseline every other candidate shares.
    """
    to_data = scale(root, "y")
    found = []
    for gid, (style, pts) in paths(root, "patch_").items():
        xs = {round(x, 6) for x, _ in pts}
        ys = {round(y, 6) for _, y in pts}
        if len(pts) < 4 or len(xs) != 2 or len(ys) != 2 or stroke(style) is None:
            continue
        found.append((min(xs), max(ys), min(ys)))
    if not found:
        return []
    baseline = max(b for _, b, _ in found)          # bottom of the axes, in points
    out = [(left, to_data(top)) for left, base, top in found
           if abs(base - baseline) < 1e-6]
    return sorted(out)


def series(root: ET.Element, axes: str = "axes_1") -> list:
    """[(stroke colour, [(x, y) in data units, ...])] for each line in the axes.

    Only lines drawn directly in the axes group, so a legend's sample lines --
    which live under `legend_1` -- are left out.

    A categorical x axis carries no numeric ticks to invert against -- the two
    columns of a pre/post panel are labelled, not numbered -- so x is left in
    points when it cannot be converted and y, which is what such a panel is
    about, is converted as usual.
    """
    node = root.find(f".//{SVG}g[@id='{axes}']")
    if node is None:
        raise ValueError(f"no {axes} in this file")
    fy = scale(root, "y")
    try:
        fx = scale(root, "x")
    except ValueError:
        fx = lambda p: p
    out = []
    for sub in node:
        if not (sub.get("id") or "").startswith("line2d_"):
            continue
        path = sub.find(SVG + "path")
        if path is None or not path.get("d"):
            continue
        colour = stroke(path.get("style") or "")
        if colour is None:
            continue
        out.append((colour, [(fx(x), fy(y)) for x, y in vertices(path.get("d"))]))
    return out


def legend(root: ET.Element, group: str = "legend_1") -> list:
    """[(stroke colour or None, label)] in the order the legend draws them."""
    node = root.find(f".//{SVG}g[@id='{group}']")
    if node is None:
        return []
    out, pending = [], None
    for sub in node:
        sid = sub.get("id") or ""
        if sid.startswith("line2d_"):
            path = sub.find(SVG + "path")
            pending = stroke(path.get("style") or "") if path is not None else None
        elif sid.startswith("text_"):
            out.append((pending, text_of(sub).strip()))
            pending = None
    return out


def labelled_series(path) -> list:
    """[(label, [(x, y), ...])] for a line plot with a colour legend.

    The convenience case: read a published line panel straight into named
    series, by matching each line's stroke colour to its legend entry.
    """
    root = load(path)
    names = {colour: label for colour, label in legend(root) if colour}
    return [(names.get(colour, colour), pts) for colour, pts in series(root)]
