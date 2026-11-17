#!/usr/bin/env python3
"""Check data/ against data/MANIFEST.csv, or refresh the manifest's checksums.

    python3 scripts/verify_data.py            # report drift
    python3 scripts/verify_data.py --refresh  # record the current state
    python3 scripts/verify_data.py --prune    # also drop rows whose file is gone
    python3 scripts/verify_data.py --adopt    # also add rows for untracked files

Rerunning a notebook rewrites files it produces, so drift under
`data/processed/` is expected after a rebuild and is not by itself a problem.
Drift under `data/raw/` is -- nothing should ever write there.
A missing file is an error by default and `--refresh` leaves its row alone,
because forgetting a file that went missing by accident is exactly what this
script exists to prevent.  `--prune` is the explicit way to record a deliberate
removal.
"""
from __future__ import annotations
import argparse, csv, hashlib, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "MANIFEST.csv"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--prune", action="store_true",
                    help="with --refresh, drop the rows whose file is gone")
    ap.add_argument("--adopt", action="store_true",
                    help="with --refresh, add a row for every file under data/ "
                         "that has none, so the manifest covers the whole tree")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(MANIFEST)))
    missing, changed, ok = [], [], 0
    for r in rows:
        p = ROOT / r["dest"]
        if not p.exists():
            missing.append(r["dest"]); continue
        digest = sha256(p)
        if r["sha256"] and digest != r["sha256"]:
            changed.append(r["dest"])
            if a.refresh:
                r["sha256"] = digest
                r["bytes_out"] = p.stat().st_size
        else:
            ok += 1

    raw_changed = [d for d in changed if d.startswith("data/raw/")]
    print(f"{len(rows)} files in manifest: {ok} unchanged, {len(changed)} changed, "
          f"{len(missing)} missing")
    if raw_changed:
        print(f"\n  WARNING: {len(raw_changed)} file(s) under data/raw/ changed; "
              "nothing should write there")
        for d in raw_changed[:10]:
            print("   ", d)
    if changed and not raw_changed:
        print("  (all drift is under data/processed/, i.e. notebook output)")
    for d in missing[:10]:
        print("   missing:", d)

    adopted = []
    if a.refresh and a.adopt:
        tracked = {r["dest"] for r in rows}
        data = ROOT / "data"
        for f in sorted(data.rglob("*")):
            if not f.is_file() or f.name.startswith(".") or f == MANIFEST:
                continue
            dest = str(f.relative_to(ROOT))
            if dest in tracked:
                continue
            size = f.stat().st_size
            adopted.append({
                "source": "", "dest": dest, "bytes_in": size, "bytes_out": size,
                "gzipped": str(f.suffix == ".gz"),
                "note": "added after the ingest; checksummed from the file on disk",
                "sha256": sha256(f)})
        rows = rows + adopted
        print(f"\n  adopted {len(adopted)} untracked file(s)")

    if a.refresh:
        keep = rows
        if a.prune and missing:
            gone = set(missing)
            keep = [r for r in rows if r["dest"] not in gone]
            print(f"\n  pruned {len(rows) - len(keep)} row(s) whose file is gone")
        with open(MANIFEST, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(keep)
        print(f"\nrefreshed {MANIFEST.relative_to(ROOT)}")
        if a.prune:
            missing = []
    return 1 if (missing or raw_changed) else 0


if __name__ == "__main__":
    sys.exit(main())
