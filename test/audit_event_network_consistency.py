#!/usr/bin/env python3
"""Audit that transfer_flow.csv is the aggregation of reaction_event_pairs.csv."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fp:
        return list(csv.DictReader(fp))


def aggregate_pairs(rows: list[dict[str, str]]) -> tuple[Counter[tuple[str, str]], Counter[tuple[str, str]]]:
    counts: Counter[tuple[str, str]] = Counter()
    atom_transfer: Counter[tuple[str, str]] = Counter()
    for row in rows:
        key = (row["source_id"], row["target_id"])
        counts[key] += 1
        atom_transfer[key] += int(row["atom_overlap"])
    return counts, atom_transfer


def aggregate_flow(rows: list[dict[str, str]]) -> tuple[Counter[tuple[str, str]], Counter[tuple[str, str]]]:
    counts: Counter[tuple[str, str]] = Counter()
    atom_transfer: Counter[tuple[str, str]] = Counter()
    for row in rows:
        key = (row["source_id"], row["target_id"])
        counts[key] += int(row["count"])
        atom_transfer[key] += int(row["atom_transfer"])
    return counts, atom_transfer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--expect-event-pairs", type=int)
    args = parser.parse_args()

    out = args.output_dir
    pairs = read_rows(out / "reaction_event_pairs.csv")
    flow = read_rows(out / "transfer_flow.csv")

    if args.expect_event_pairs is not None and len(pairs) != args.expect_event_pairs:
        raise SystemExit(f"event-pair count mismatch: got {len(pairs)}, expected {args.expect_event_pairs}")

    pair_atom_transfer = sum(int(row["atom_overlap"]) for row in pairs)
    flow_atom_transfer = sum(int(row["atom_transfer"]) for row in flow)
    if pair_atom_transfer != flow_atom_transfer:
        raise SystemExit(
            "reaction_event_pairs atom_overlap total does not match transfer_flow atom_transfer total: "
            f"{pair_atom_transfer} != {flow_atom_transfer}"
        )

    pair_counts, pair_atoms = aggregate_pairs(pairs)
    flow_counts, flow_atoms = aggregate_flow(flow)
    if pair_counts != flow_counts or pair_atoms != flow_atoms:
        mismatches = []
        for key in sorted(set(pair_counts) | set(flow_counts)):
            if pair_counts[key] != flow_counts[key] or pair_atoms[key] != flow_atoms[key]:
                mismatches.append((key, pair_counts[key], flow_counts[key], pair_atoms[key], flow_atoms[key]))
        preview = "; ".join(
            f"{src}->{tgt}: pair_count={pc}, flow_count={fc}, pair_atoms={pa}, flow_atoms={fa}"
            for (src, tgt), pc, fc, pa, fa in mismatches[:10]
        )
        raise SystemExit(f"transfer_flow.csv is not an exact aggregation of reaction_event_pairs.csv: {preview}")

    self_loop_count = sum(int(row["count"]) for row in flow if row.get("self_loop") == "1")
    self_loop_atoms = sum(int(row["atom_transfer"]) for row in flow if row.get("self_loop") == "1")
    print("Event-network consistency audit passed")
    print(f"  reaction event pairs: {len(pairs)}")
    print(f"  transfer edges: {len(flow)}")
    print(f"  atom transfer total: {flow_atom_transfer}")
    print(f"  self-loop count: {self_loop_count}")
    print(f"  self-loop atom transfer: {self_loop_atoms}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
