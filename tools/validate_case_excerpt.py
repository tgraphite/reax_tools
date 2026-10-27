#!/usr/bin/env python3
"""Re-run a bounded XYZ excerpt without changing the source case.

Accepts four-column element/xyz or five-column id/element/xyz rows, retaining
the original cell comment, coordinates and row ordering. No MD is performed.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import time

from prepublication_validate import inspect, sha, xyz_frames


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--core", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--frames", type=int, default=200)
    p.add_argument("--scale", type=float, default=1.0)
    a = p.parse_args()
    if a.frames < 1:
        p.error("--frames must be positive")
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    selected = []
    for block in xyz_frames(a.input):
        selected.append(block)
        if len(selected) == a.frames:
            break
    if not selected:
        raise ValueError("empty source")
    atoms = int(selected[0][0])
    inventory = Counter()
    for line in selected[0][2:]:
        fields = line.split()
        if len(fields) not in (4, 5):
            raise ValueError("only four/five column XYZ supported")
        inventory[fields[0] if len(fields) == 4 else fields[1]] += 1
    excerpt = root / "excerpt.xyz"
    excerpt.write_text("".join("".join(block) for block in selected))
    out = root / "analysis"
    cmd = [str(a.core.resolve()), "-f", str(excerpt), "-o", str(out), "-nt", "1", "-r", str(a.scale)]
    start = time.perf_counter()
    run = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    elapsed = time.perf_counter()-start
    (root / "run.log").write_text(run.stdout+run.stderr)
    if run.returncode:
        raise RuntimeError(f"core exit {run.returncode}; inspect {root}/run.log")
    record = inspect(out, atoms, inventory)
    if record["frames"] != len(selected):
        record["bookkeeping_errors"].append("frame count mismatch")
    record.update(input=str(a.input.resolve()), source_sha256=sha(a.input),
                  excerpt_sha256=sha(excerpt), selected_frames=len(selected),
                  selection="first consecutive frames; no coordinate/cell/ID modification",
                  command=cmd, wall_seconds=elapsed, core_sha256=sha(a.core),
                  source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                  scope="Table bookkeeping and bounded rerun only; not independent chemical accuracy.")
    (root / "results.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))
    return 1 if record["bookkeeping_errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
