#!/usr/bin/env python3
"""Compare the new structural-identity output (v5) with the old formula output (v4).

The question: if we manually merge v5 structural species back to formulas,
do the counts/events/flows agree with the old formula-level run?
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V4 = ROOT / "test" / "energetic_v4"
V5 = ROOT / "test" / "energetic_v5"
OUT = Path(__file__).resolve().parent / "v4_v5_compare"
OUT.mkdir(exist_ok=True)


def rows(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def split(value, sep):
    return tuple(sorted(part for part in str(value).split(sep) if part))


def cancel(reactants, products):
    left, right = Counter(reactants), Counter(products)
    for key in set(left) & set(right):
        n = min(left[key], right[key])
        left[key] -= n
        right[key] -= n
    return (
        tuple(sorted(k for k, n in left.items() for _ in range(n))),
        tuple(sorted(k for k, n in right.items() for _ in range(n))),
    )


def compare_text(name):
    a = (V4 / name).read_text()
    b = (V5 / name).read_text()
    return a == b


def species_tables():
    v4 = defaultdict(lambda: defaultdict(int))
    for r in rows(V4 / "species_count.csv"):
        v4[r["formula"]][int(r["frame"])] += int(r["count"])

    v5_formula = {r["id"]: r["formula"] for r in rows(V5 / "molecules.csv")}
    header = next(csv.reader((V5 / "species_count.csv").open()))
    ids = header[1:]
    v5 = defaultdict(lambda: defaultdict(int))
    with (V5 / "species_count.csv").open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            frame = int(row[0])
            for molecule_id, value in zip(ids, row[1:]):
                count = int(value)
                if count:
                    v5[v5_formula[molecule_id]][frame] += count
    return v4, v5


def compare_species(v4, v5):
    formulas = sorted(set(v4) | set(v5))
    mismatches = []
    for formula in formulas:
        frames = sorted(set(v4[formula]) | set(v5[formula]))
        for frame in frames:
            a = v4[formula].get(frame, 0)
            b = v5[formula].get(frame, 0)
            if a != b:
                mismatches.append((formula, frame, a, b))
    return mismatches


def events(path, sep, structural):
    table = Counter()
    for r in rows(path):
        table[(int(r["frame"]),
               split(r["reactant_formulas"], sep),
               split(r["product_formulas"], sep))] += 1
    return table


def event_cancel(table):
    out = Counter()
    for (frame, react, prod), n in table.items():
        r, p = cancel(react, prod)
        if not r and not p:
            # fully cancelling same-formula event: keep a marker for reporting
            out[(frame, ("<cancelled>",), ())] += n
        else:
            out[(frame, r, p)] += n
    return out


def compare_events():
    v4_raw = events(V4 / "reaction_events.csv", "+", False)
    v5_raw = events(V5 / "reaction_events.csv", ";", True)
    v4_c = event_cancel(v4_raw)
    v5_c = event_cancel(v5_raw)
    return {
        "v4_raw": len(v4_raw), "v5_raw": len(v5_raw),
        "v4_cancelled": len(v4_c), "v5_cancelled": len(v5_c),
        "raw_extra": v5_raw - v4_raw, "raw_missing": v4_raw - v5_raw,
        "cancel_extra": v5_c - v4_c, "cancel_missing": v4_c - v5_c,
        "same_formula_v5": sum(n for (f, r, p), n in v5_raw.items() if r == p),
    }


def flows(path, source_col, target_col):
    counts, atoms = Counter(), Counter()
    for r in rows(path):
        key = (r[source_col], r[target_col])
        counts[key] += int(r["count"])
        atoms[key] += int(r["atom_transfer"])
    return counts, atoms


def main():
    global V4, V5, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v4", type=Path, default=V4)
    parser.add_argument("--v5", type=Path, default=V5)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    V4, V5, OUT = args.v4, args.v5, args.out
    OUT.mkdir(parents=True, exist_ok=True)

    summary = {}
    summary["count_tables_identical"] = {
        name: compare_text(name)
        for name in ["bond_count.csv", "atom_bonded_num_count.csv", "ring_count.csv"]
    }

    v4_species, v5_species = species_tables()
    mismatches = compare_species(v4_species, v5_species)
    summary["species"] = {
        "v4_formulas": len(v4_species),
        "v5_formulas": len(v5_species),
        "formula_sets_equal": set(v4_species) == set(v5_species),
        "mismatched_formula_frame_cells": len(mismatches),
        "first_mismatches": mismatches[:20],
    }

    v4_molecules = json.loads((V4 / "molecules.json").read_text())["molecules"]
    v4_mol_formulas = {m["formula"] for m in v4_molecules}
    v5_mol_rows = rows(V5 / "molecules.csv")
    v5_mol_formulas = {m["formula"] for m in v5_mol_rows}
    summary["molecules"] = {
        "v4_rows": len(v4_molecules),
        "v5_rows": len(v5_mol_rows),
        "v4_formulas": len(v4_mol_formulas),
        "v5_formulas": len(v5_mol_formulas),
        "formula_sets_equal": v4_mol_formulas == v5_mol_formulas,
    }

    ev = compare_events()
    summary["events"] = {
        k: (len(v) if isinstance(v, Counter) else v)
        for k, v in ev.items()
    }
    # Save the actual extra/missing signatures for inspection.
    with (OUT / "event_extra_raw.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["frame", "reactants", "products", "count"])
        for (frame, r, p), n in sorted(ev["raw_extra"].items()):
            writer.writerow([frame, ";".join(r), ";".join(p), n])
    with (OUT / "event_missing_raw.csv").open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["frame", "reactants", "products", "count"])
        for (frame, r, p), n in sorted(ev["raw_missing"].items()):
            writer.writerow([frame, ";".join(r), ";".join(p), n])

    v4_counts, v4_atoms = flows(V4 / "transfer_flow.csv", "source_label", "target_label")
    v5_counts, v5_atoms = flows(V5 / "transfer_flow.csv", "source_formula", "target_formula")
    keys = set(v4_counts) | set(v5_counts)
    summary["transfer_flow"] = {
        "v4_edges": len(v4_counts),
        "v5_collapsed_edges": len(v5_counts),
        "edge_sets_equal": set(v4_counts) == set(v5_counts),
        "count_mismatches": sum(1 for k in keys if v4_counts.get(k, 0) != v5_counts.get(k, 0)),
        "atom_mismatches": sum(1 for k in keys if v4_atoms.get(k, 0) != v5_atoms.get(k, 0)),
        "top_count_diffs": sorted(
            [(k, v4_counts.get(k, 0), v5_counts.get(k, 0)) for k in keys
             if v4_counts.get(k, 0) != v5_counts.get(k, 0)],
            key=lambda item: abs(item[1] - item[2]), reverse=True)[:15],
    }

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    print("Count tables identical:", summary["count_tables_identical"])
    print("Species:", summary["species"])
    print("Molecules:", summary["molecules"])
    print("Events:", summary["events"])
    print("Transfer flow:", {k: v for k, v in summary["transfer_flow"].items() if k != "top_count_diffs"})
    print("Top flow count diffs:")
    for item in summary["transfer_flow"]["top_count_diffs"]:
        print("  ", item)
    print("Wrote", OUT / "summary.json")


if __name__ == "__main__":
    main()
