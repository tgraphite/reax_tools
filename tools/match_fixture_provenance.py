#!/usr/bin/env python3
"""Match an element/x/y/z XYZ fixture to orthorhombic LAMMPS frames.

Read-only on inputs. Rounds coordinates to 5 decimals for candidate matching,
then records source IDs and maximum coordinate differences for every match.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path


def checksum(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atom_key(symbol, coords):
    return (symbol, *(round(x, 5) for x in coords))


def frame_key(atoms):
    return tuple(sorted(atom_key(s, c) for s, c in atoms))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xyz", type=Path, required=True)
    p.add_argument("--dump", type=Path, required=True)
    p.add_argument("--types", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    symbols = a.types.split(",")
    fixtures = []
    with a.xyz.open() as f:
        while first := f.readline():
            n = int(first)
            comment = f.readline().strip()
            atoms = []
            for _ in range(n):
                row = f.readline().split()
                atoms.append((row[0], tuple(map(float, row[1:4]))))
            fixtures.append((atoms, comment))
    wanted = defaultdict(list)
    for i, (atoms, _) in enumerate(fixtures):
        wanted[frame_key(atoms)].append(i)
    matches = [[] for _ in fixtures]
    raw_frames = 0
    with a.dump.open() as f:
        while first := f.readline():
            if first.strip() != "ITEM: TIMESTEP":
                raise ValueError("Expected TIMESTEP")
            step = int(f.readline())
            f.readline()
            n = int(f.readline())
            box_header = f.readline().strip()
            if any(v in box_header.split() for v in ("xy", "xz", "yz")):
                raise ValueError("This provenance matcher accepts orthorhombic dumps only")
            bounds = [list(map(float, f.readline().split())) for _ in range(3)]
            fields = f.readline().split()[2:]
            ix = {v: i for i, v in enumerate(fields)}
            atoms, by_key = [], defaultdict(list)
            for _ in range(n):
                r = f.readline().split()
                symbol = symbols[int(r[ix["type"]]) - 1]
                coords = tuple(float(r[ix[k]]) for k in ("x", "y", "z"))
                atoms.append((symbol, coords))
                by_key[atom_key(symbol, coords)].append((int(r[ix["id"]]), coords))
            raw_frames += 1
            for i in wanted.get(frame_key(atoms), []):
                ids, diffs = [], []
                for s, c in fixtures[i][0]:
                    candidates = by_key[atom_key(s, c)]
                    if len(candidates) != 1:
                        raise ValueError("Ambiguous rounded atom match")
                    source_id, source_coords = candidates[0]
                    ids.append(source_id)
                    diffs.extend(abs(x-y) for x,y in zip(c, source_coords))
                matches[i].append({"source_step": step, "box_bounds": bounds,
                                   "row_source_ids": ids, "max_coordinate_delta_angstrom": max(diffs)})
    unique = all(len(m) == 1 for m in matches)
    mapping = matches[0][0]["row_source_ids"] if unique else None
    stable = unique and all(m[0]["row_source_ids"] == mapping for m in matches)
    report = {"fixture": str(a.xyz.resolve()), "source": str(a.dump.resolve()),
              "fixture_sha256": checksum(a.xyz), "source_sha256": checksum(a.dump),
              "types": symbols, "raw_frames": raw_frames, "fixture_frames": len(fixtures),
              "all_frames_uniquely_matched": unique, "row_to_source_id_stable": stable,
              "row_ids_equal_source_ids": stable and mapping == list(range(1, len(mapping)+1)),
              "matching_precision_decimal_places": 5,
              "frame_matches": [[{k:v for k,v in m.items() if k != "row_source_ids"} for m in mm] for mm in matches],
              "source_id_mapping_sha256": hashlib.sha256(repr(mapping).encode()).hexdigest() if stable else None,
              "scope": "Coordinates, elements and ID ordering only; not force-field or chemical accuracy validation."}
    with a.output.open("x") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    print(json.dumps({k:v for k,v in report.items() if k != "frame_matches"}, indent=2))
    return 0 if unique and stable else 1


if __name__ == "__main__":
    raise SystemExit(main())
