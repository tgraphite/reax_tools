#!/usr/bin/env python3
"""Reproducible software pilots; synthetic coordinate cases are NOT MD simulations.

Writes only to a new --output directory. Records limitations without turning them
into passing chemistry tests. Uses independent formula/overlap bookkeeping.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    if not path.exists():
        return []
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def composition(formula):
    result = Counter()
    for element, count in re.findall(r"([A-Z][a-z]*)(\d*)", formula):
        result[element] += int(count or 1)
    return result


def event_signature(out):
    return Counter((r["frame"], tuple(sorted(r["reactant_formulas"].split("+"))),
                    tuple(sorted(r["product_formulas"].split("+"))))
                   for r in rows(out / "reaction_events.csv"))


def inspect(out, expected_atoms, expected_composition):
    events = rows(out / "reaction_events.csv")
    pairs = rows(out / "reaction_event_pairs.csv")
    flow = rows(out / "transfer_flow.csv")
    errors = []
    event_by_id = {r["event_id"]: r for r in events}
    if len(event_by_id) != len(events):
        errors.append("duplicate event IDs")
    pair_count, pair_atoms = Counter(), Counter()
    per_event_atoms = Counter()
    for r in pairs:
        event = event_by_id.get(r["event_id"])
        if not event or event["frame"] != r["frame"]:
            errors.append("orphan or wrong-frame pair")
        elif (r["source_id"] not in event["reactant_hashes"].split("+") or
              r["target_id"] not in event["product_hashes"].split("+")):
            errors.append("pair species not in event")
        key = (r["source_id"], r["target_id"])
        pair_count[key] += 1
        pair_atoms[key] += int(r["atom_overlap"])
        per_event_atoms[r["event_id"]] += int(r["atom_overlap"])
    flow_count, flow_atoms = Counter(), Counter()
    for r in flow:
        key = (r["source_id"], r["target_id"])
        flow_count[key] += int(r["count"])
        flow_atoms[key] += int(r["atom_transfer"])
    if pair_count != flow_count or pair_atoms != flow_atoms:
        errors.append("pair/flow mismatch")
    for r in events:
        left, right = Counter(), Counter()
        for formula in r["reactant_formulas"].split("+"):
            left.update(composition(formula))
        for formula in r["product_formulas"].split("+"):
            right.update(composition(formula))
        if left != right:
            errors.append("event formula imbalance")
        if per_event_atoms[r["event_id"]] != sum(left.values()):
            errors.append("event overlap total differs from reactant atoms")
    frame_compositions = defaultdict(Counter)
    for r in rows(out / "species_count.csv"):
        for element, number in composition(r["formula"]).items():
            frame_compositions[r["frame"]][element] += number * int(r["count"])
    if not frame_compositions or any(c != expected_composition for c in frame_compositions.values()):
        errors.append("species element inventory mismatch")
    manifest = json.loads((out / "reax_tools_manifest.json").read_text())
    return {"events": len(events), "pairs": len(pairs), "edges": len(flow),
            "atom_transfer": sum(flow_atoms.values()), "frames": len(frame_compositions),
            "input_atoms": expected_atoms,
            "bookkeeping_errors": sorted(set(errors)),
            "manifest_missing_files": [f for f in manifest["files"] if not (out / f).exists()],
            "output_sha256": {f.name: sha(f) for f in sorted(out.glob("*.csv"))}}


def xyz_frames(path):
    with path.open() as stream:
        while True:
            first = stream.readline()
            if not first:
                return
            n = int(first)
            block = [first, stream.readline()]
            block.extend(stream.readline() for _ in range(n))
            if any(not s for s in block):
                raise ValueError("truncated XYZ")
            yield block


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    core, target = args.core.resolve(), args.output.resolve()
    target.mkdir(parents=True, exist_ok=False)
    inputs = target / "inputs"
    inputs.mkdir()
    results = []

    def run(name, source, options):
        out = target / name
        cmd = [str(core), "-f", str(source), "-o", str(out), *options]
        # A fresh supervisor isolates ru_maxrss from previous child processes.
        supervisor = ("import pathlib,resource,subprocess,sys,time; t=time.perf_counter(); "
                      "p=subprocess.run(sys.argv[2:]); "
                      "pathlib.Path(sys.argv[1]).write_text(str(time.perf_counter()-t)+' '+"
                      "str(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)); sys.exit(p.returncode)")
        timed_cmd = [sys.executable, "-c", supervisor, str(target / (name + ".time")), *cmd]
        start = time.perf_counter()
        completed = subprocess.run(timed_cmd, capture_output=True, text=True, timeout=180)
        (target / (name + ".log")).write_text(completed.stdout + completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"{name} failed ({completed.returncode}); inspect log")
        first = next(xyz_frames(source))
        elements = Counter(line.split()[0] for line in first[2:])
        record = inspect(out, int(first[0]), elements)
        timing = (target / (name + ".time")).read_text().strip().split()
        record.update(name=name, command=cmd, input_sha256=sha(source),
                      wall_seconds=time.perf_counter() - start,
                      time_wall_seconds=float(timing[0]), peak_rss_kib=int(timing[1]))
        results.append(record)
        print(json.dumps({k: record[k] for k in ("name", "events", "pairs", "frames", "bookkeeping_errors")}), flush=True)
        return out

    # All cells are orthorhombic and comfortably larger than the cutoff.
    # H-H distance 0.7 is bonded, distance >= 3 is unbonded at default scale.
    geometries = {
        "split": [[(1,1,1),(1.7,1,1)], [(1,1,1),(5,1,1)]],
        "jitter": [[(1,1,1),(1.7,1,1)], [(1,1,1),(5,1,1)], [(1,1,1),(1.7,1,1)]],
        "exchange": [[(1,1,1),(1.7,1,1),(6,1,1),(6.7,1,1)],
                     [(1,1,1),(6,1,1),(1.7,1,1),(6.7,1,1)]],
        "null": [[(1,1,1),(1.7,1,1)]] * 3,
        "pbc": [[(.2,1,1),(19.5,1,1)], [(.2,1,1),(.9,1,1)]],
    }
    for name, frames in geometries.items():
        source = inputs / (name + ".xyz")
        source.write_text("".join(str(len(frame)) + '\nLattice="20 0 0 0 20 0 0 0 20"\n' +
                                  "".join(f"H {x} {y} {z}\n" for x,y,z in frame) for frame in frames))
        run(name, source, ["-nt", "1", "--no-rings"])
    run("jitter_stable10", inputs / "jitter.xyz", ["-nt", "1", "--no-rings", "--stable-time", "10"])
    source = ROOT / "test/energetic.xyz"
    for threads in (1, 4):
        for repeat in range(3):
            run(f"energetic_nt{threads}_rep{repeat}", source, ["-nt", str(threads)])
    for scale in (1.0, 1.1):
        run(f"energetic_r{scale}", source, ["-nt", "1", "-r", str(scale)])
    for stride in (2, 4):
        sample = inputs / f"energetic_stride{stride}.xyz"
        sample.write_text("".join("".join(block) for i, block in enumerate(xyz_frames(source)) if i % stride == 0))
        run(f"energetic_stride{stride}", sample, ["-nt", "1"])
    ref = target / "energetic_nt1_rep0"
    comparisons = {}
    for threads in (1, 4):
        for repeat in range(3):
            other = target / f"energetic_nt{threads}_rep{repeat}"
            comparisons[other.name] = {
                "event_multiset_equal": event_signature(ref) == event_signature(other),
                "species_byte_equal": sha(ref / "species_count.csv") == sha(other / "species_count.csv"),
                "flow_byte_equal": sha(ref / "transfer_flow.csv") == sha(other / "transfer_flow.csv"),
            }
    observations = {
        "split_expected_one_event": len(event_signature(target / "split")) == 1,
        "same_formula_exchange_detected_events": sum(event_signature(target / "exchange").values()),
        "jitter_event_count": sum(event_signature(target / "jitter").values()),
        "stable_time_default3_vs_10_same_events": event_signature(target / "jitter") == event_signature(target / "jitter_stable10"),
        "null_event_count": sum(event_signature(target / "null").values()),
        "pbc_event_count": sum(event_signature(target / "pbc").values()),
    }
    report = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "core_sha256": sha(core), "platform": platform.platform(),
              "cpu": next((s.split(":",1)[1].strip() for s in Path("/proc/cpuinfo").read_text().splitlines() if s.startswith("model name")), "unknown"),
              "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "runs": results, "thread_comparisons": comparisons, "observations": observations,
              "scope": "Software pilots only; no physical reaction ground truth or competitor timing."}
    (target / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(observations, indent=2))
    return 1 if any(r["bookkeeping_errors"] for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
