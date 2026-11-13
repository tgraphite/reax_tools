#!/usr/bin/env python3
"""Triclinic trajectory, --dump and --dumpxyz regression tests."""

from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORE = ROOT / "bin" / "reax_tools_core"
COMBUSTION = ROOT / "test" / "combustion.xyz"
SMALL = ROOT / "test" / "fixtures" / "small_molecules.xyz"
MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999}


def read_extxyz(path: Path):
    lines = path.read_text().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].strip())
        header = lines[i + 1]
        lattice_match = re.search(r'Lattice="([^"]+)"', header)
        if not lattice_match:
            raise ValueError(f"{path}: missing Lattice at frame {len(frames)}")
        lattice = [float(x) for x in lattice_match.group(1).split()]
        properties_match = re.search(r"Properties=([^\s]+)", header)
        if not properties_match:
            raise ValueError(f"{path}: missing Properties")
        prop_parts = properties_match.group(1).split(":")
        columns = []
        offset = 0
        for j in range(0, len(prop_parts) - 2, 3):
            name = prop_parts[j]
            count = int(prop_parts[j + 2])
            columns.append((name, offset, count))
            offset += count
        atoms = []
        for line in lines[i + 2 : i + 2 + n]:
            tokens = line.split()
            values = {}
            for name, start, count in columns:
                if count == 1:
                    values[name] = tokens[start]
                elif name in {"pos", "positions"}:
                    values["pos"] = tuple(float(tokens[start + k]) for k in range(3))
                else:
                    values[name] = tuple(tokens[start + k] for k in range(count))
            atoms.append(values)
        frames.append({"lattice": lattice, "header": header, "atoms": atoms})
        i += 2 + n
    return frames


def write_extxyz(path: Path, frames):
    out = []
    for frame in frames:
        out.append(f"{len(frame['atoms'])}\n")
        out.append(
            'Lattice="' + " ".join(f"{v:.10f}" for v in frame["lattice"]) +
            '" Properties=species:S:1:pos:R:3 pbc="T T T"\n'
        )
        for atom in frame["atoms"]:
            species = atom["species"]
            x, y, z = atom["pos"]
            out.append(f"{species} {x:.10f} {y:.10f} {z:.10f}\n")
    path.write_text("".join(out))


def write_extxyz_2d(path: Path, frames):
    out = []
    for frame in frames:
        a, b, c = columns_from_lattice(frame["lattice"])
        matrix = [
            [a[0], b[0], c[0]],
            [a[1], b[1], c[1]],
            [a[2], b[2], c[2]],
        ]
        lattice = "[[" + "], [".join(", ".join(f"{value:.10f}" for value in row) for row in matrix) + "]]"
        out.append(f"{len(frame['atoms'])}\n")
        out.append(f'Lattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for atom in frame["atoms"]:
            x, y, z = atom["pos"]
            out.append(f"{atom['species']} {x:.10f} {y:.10f} {z:.10f}\n")
    path.write_text("".join(out))


def columns_from_lattice(lattice):
    a = lattice[0:3]
    b = lattice[3:6]
    c = lattice[6:9]
    return a, b, c


def cart(columns, frac):
    a, b, c = columns
    return tuple(a[i] * frac[0] + b[i] * frac[1] + c[i] * frac[2] for i in range(3))


def cell_for_shear(lattice, xy=4.0, xz=3.0, yz=2.0):
    a, b, c = columns_from_lattice(lattice)
    # Restricted cell columns: a=(lx,0,0), b=(xy,ly,0), c=(xz,yz,lz).
    return ((a[0], 0.0, 0.0), (xy, b[1], 0.0), (xz, yz, c[2]))


def write_lammpstrj(path: Path, frames, cell, scaled=False):
    lx = cell[0][0]
    ly = cell[1][1]
    lz = cell[2][2]
    xy = cell[1][0]
    xz = cell[2][0]
    yz = cell[2][1]
    xlo = ylo = zlo = 0.0
    xhi, yhi, zhi = lx, ly, lz
    xlo_b = xlo + min(0.0, xy, xz, xy + xz)
    xhi_b = xhi + max(0.0, xy, xz, xy + xz)
    ylo_b = ylo + min(0.0, yz)
    yhi_b = yhi + max(0.0, yz)

    out = []
    for frame_index, frame in enumerate(frames):
        out.append(f"ITEM: TIMESTEP\n{frame_index}\n")
        out.append(f"ITEM: NUMBER OF ATOMS\n{len(frame['atoms'])}\n")
        out.append("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        out.append(f"{xlo_b} {xhi_b} {xy}\n{ylo_b} {yhi_b} {xz}\n{zlo} {zhi} {yz}\n")
        out.append("ITEM: ATOMS id type xs ys zs\n" if scaled else "ITEM: ATOMS id type x y z\n")
        for atom in frame["atoms"]:
            species = atom["species"]
            type_id = {"C": 1, "H": 2, "O": 3}[species]
            if scaled:
                values = atom["frac"]
            else:
                values = atom["pos"]
            out.append(
                f"{atom['id']} {type_id} {values[0]:.10f} {values[1]:.10f} {values[2]:.10f}\n"
            )
    path.write_text("".join(out))


def frames_from_combustion(limit: int):
    lines = COMBUSTION.read_text().splitlines(keepends=True)
    frames = []
    i = 0
    while len(frames) < limit:
        n = int(lines[i].strip())
        header = lines[i + 1]
        match = re.search(r'Lattice="([^"]+)"', header)
        lattice = [float(x) for x in match.group(1).split()]
        atoms = []
        for atom_id, line in enumerate(lines[i + 2 : i + 2 + n], 1):
            tokens = line.split()
            atoms.append(
                {
                    "id": atom_id,
                    "species": tokens[0],
                    "pos": (float(tokens[1]), float(tokens[2]), float(tokens[3])),
                }
            )
        frames.append({"lattice": lattice, "atoms": atoms})
        i += 2 + n
    return frames


def parse_data(path: Path):
    lines = path.read_text().splitlines()
    result = {"atoms": [], "bonds": [], "angles": [], "dihedrals": [], "box": []}
    section = None
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        match = re.match(r"(\d+) (atoms|bonds|angles|dihedrals|atom types|bond types|angle types|dihedral types)$", stripped)
        if match:
            result["n_" + match.group(2).replace(" ", "_")] = int(match.group(1))
            continue
        if stripped.startswith("Atoms #"):
            section = "atoms"
            continue
        if stripped in {"Bonds", "Angles", "Dihedrals"}:
            section = stripped.lower()
            continue
        if stripped.startswith("Masses"):
            section = "masses"
            continue
        if "xlo" in stripped or "ylo" in stripped or "zlo" in stripped:
            result["box"].append(stripped)
            continue
        if section in {"atoms", "bonds", "angles", "dihedrals"}:
            values = stripped.split()
            if section == "atoms":
                if len(values) >= 7 and values[0].isdigit():
                    result["atoms"].append(
                        {
                            "id": int(values[0]),
                            "molid": int(values[1]),
                            "type": int(values[2]),
                            "pos": (float(values[4]), float(values[5]), float(values[6])),
                        }
                    )
            elif section == "bonds":
                if len(values) >= 4 and values[0].isdigit():
                    result["bonds"].append(tuple(int(x) for x in values[:4]))
            elif section == "angles":
                if len(values) >= 5 and values[0].isdigit():
                    result["angles"].append(tuple(int(x) for x in values[:5]))
            else:
                if len(values) >= 6 and values[0].isdigit():
                    result["dihedrals"].append(tuple(int(x) for x in values[:6]))
    return result


def expected_topology(bonds):
    adjacency = defaultdict(set)
    for _, _, i, j in bonds:
        adjacency[i].add(j)
        adjacency[j].add(i)
    angles = set()
    for center, neighbors in adjacency.items():
        ordered = sorted(neighbors)
        for i in range(len(ordered)):
            for j in range(i + 1, len(ordered)):
                angles.add((ordered[i], center, ordered[j]))
    dihedrals = set()
    for _, _, i, j in bonds:
        left = sorted(adjacency[i] - {j})
        right = sorted(adjacency[j] - {i})
        for a in left:
            for b in right:
                forward = (a, i, j, b)
                reverse = (b, j, i, a)
                dihedrals.add(min(forward, reverse))
    return angles, dihedrals


def parse_formula(formula: str):
    counts = Counter()
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[element] += int(number) if number else 1
    return counts


class ReaxToolsTriclinicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument("--core", default=str(DEFAULT_CORE))
        args, _ = parser.parse_known_args()
        cls.core = Path(args.core)
        if not cls.core.exists():
            raise unittest.SkipTest(f"core binary not found: {cls.core}")

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="reax_triclinic_")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_core(self, traj: Path, out: Path, extra=None):
        out.mkdir(parents=True, exist_ok=True)
        command = [
            str(self.core), "-f", str(traj), "-o", str(out),
            "--no-track-reactions", "--no-reactions",
        ] + list(extra or [])
        completed = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if completed.returncode != 0:
            raise AssertionError(f"core failed: {' '.join(command)}\n{completed.stdout}\n{completed.stderr}")
        return completed

    def compare_outputs(self, first: Path, second: Path):
        for name in ["species_count.csv", "molecules.csv", "bond_count.csv", "ring_count.csv"]:
            self.assertEqual(
                (first / name).read_text(), (second / name).read_text(),
                f"{name} differs between {first} and {second}",
            )

    def test_xz_permuted_lattice(self):
        frames = frames_from_combustion(3)
        original = self.root / "original.xyz"
        permuted = self.root / "permuted.xyz"
        write_extxyz(original, frames)
        permuted_frames = []
        for frame in frames:
            a, b, c = columns_from_lattice(frame["lattice"])
            permuted_frames.append(
                {
                    "lattice": [a[2], a[1], a[0], b[2], b[1], b[0], c[2], c[1], c[0]],
                    "atoms": [
                        {"species": atom["species"], "pos": (atom["pos"][2], atom["pos"][1], atom["pos"][0])}
                        for atom in frame["atoms"]
                    ],
                }
            )
        write_extxyz(permuted, permuted_frames)
        out_original = self.root / "out_original"
        out_permuted = self.root / "out_permuted"
        self.run_core(original, out_original, ["--dump", "1"])
        self.run_core(permuted, out_permuted, ["--dump", "1"])
        self.compare_outputs(out_original, out_permuted)

        # A reflected xz flip is left-handed.  The data writer must still
        # produce the same orthorhombic box lengths, not fall back to a bbox.
        data = parse_data(out_permuted / "frame_1.data")
        lengths = [float(line.split()[1]) for line in data["box"]]
        self.assertAlmostEqual(lengths[0], 49.019968, places=5)
        self.assertAlmostEqual(lengths[1], 49.019816, places=5)
        self.assertAlmostEqual(lengths[2], 49.044123, places=5)

    def test_triclinic_xyz_and_lammpstrj(self):
        original_frames = frames_from_combustion(3)
        lattice = original_frames[0]["lattice"]
        cell = cell_for_shear(lattice)
        frames = []
        for frame in original_frames:
            a, b, c = columns_from_lattice(frame["lattice"])
            atoms = []
            for atom in frame["atoms"]:
                frac = (atom["pos"][0] / a[0], atom["pos"][1] / b[1], atom["pos"][2] / c[2])
                atoms.append(
                    {
                        "id": atom["id"],
                        "species": atom["species"],
                        "frac": frac,
                        "pos": cart(cell, frac),
                    }
                )
            frames.append({"lattice": list(cell[0]) + list(cell[1]) + list(cell[2]), "atoms": atoms})

        xyz = self.root / "triclinic.xyz"
        xyz_2d = self.root / "triclinic_2d.xyz"
        dump = self.root / "triclinic.lammpstrj"
        scaled_dump = self.root / "triclinic_scaled.lammpstrj"
        permuted = self.root / "triclinic_permuted.xyz"
        write_extxyz(xyz, frames)
        write_extxyz_2d(xyz_2d, frames)
        write_lammpstrj(dump, frames, cell)
        write_lammpstrj(scaled_dump, frames, cell, scaled=True)

        permuted_frames = []
        for frame in frames:
            a, b, c = columns_from_lattice(frame["lattice"])
            permuted_frames.append(
                {
                    "lattice": [a[2], a[1], a[0], b[2], b[1], b[0], c[2], c[1], c[0]],
                    "atoms": [
                        {"species": atom["species"], "pos": (atom["pos"][2], atom["pos"][1], atom["pos"][0])}
                        for atom in frame["atoms"]
                    ],
                }
            )
        write_extxyz(permuted, permuted_frames)

        out_xyz = self.root / "out_xyz"
        out_2d = self.root / "out_2d"
        out_dump = self.root / "out_dump"
        out_scaled = self.root / "out_scaled"
        out_permuted = self.root / "out_permuted"
        self.run_core(xyz, out_xyz)
        self.run_core(xyz_2d, out_2d)
        self.run_core(dump, out_dump, ["-t", "C,H,O"])
        self.run_core(scaled_dump, out_scaled, ["-t", "C,H,O"])
        self.run_core(permuted, out_permuted)
        self.compare_outputs(out_xyz, out_2d)
        self.compare_outputs(out_xyz, out_dump)
        self.compare_outputs(out_xyz, out_scaled)
        self.compare_outputs(out_xyz, out_permuted)

        out_data = self.root / "out_triclinic_data"
        self.run_core(xyz, out_data, ["--dump", "2"])
        self.assertTrue((out_data / "frame_1.data").exists())
        self.assertTrue((out_data / "frame_3.data").exists())
        self.assertFalse((out_data / "frame_2.data").exists())
        data = parse_data(out_data / "frame_1.data")
        self.assertEqual(len(data["box"]), 3)
        box_values = [[float(value) for value in line.split()[:3]] for line in data["box"]]
        self.assertAlmostEqual(box_values[0][2], 4.0)
        self.assertAlmostEqual(box_values[1][2], 3.0)
        self.assertAlmostEqual(box_values[2][2], 2.0)
        self.assertGreater(data["n_angles"], 0)
        self.assertGreater(data["n_dihedrals"], 0)

    def test_dump_topology_molid_and_ringsize(self):
        out_default = self.root / "out_default"
        out_mark = self.root / "out_mark"
        out_noring = self.root / "out_noring"
        self.run_core(SMALL, out_default, ["--dump", "1"])
        self.run_core(SMALL, out_mark, ["--dump", "1", "--dumpxyz", "1", "--mark-ring-atoms"])
        self.run_core(SMALL, out_noring, ["--dump", "1", "--dumpxyz", "1", "--mark-ring-atoms", "--no-ring"])

        data = parse_data(out_default / "frame_1.data")
        self.assertEqual(len({atom["id"] for atom in data["atoms"]}), len(data["atoms"]))
        self.assertEqual(len(data["atoms"]), data["n_atoms"])
        self.assertEqual(len(data["bonds"]), data["n_bonds"])
        self.assertEqual(len(data["angles"]), data["n_angles"])
        self.assertEqual(len(data["dihedrals"]), data["n_dihedrals"])
        expected_angles, expected_dihedrals = expected_topology(data["bonds"])
        self.assertEqual({tuple(x[2:]) for x in data["angles"]}, expected_angles)
        self.assertEqual({tuple(x[2:]) for x in data["dihedrals"]}, expected_dihedrals)

        sizes = defaultdict(int)
        for atom in data["atoms"]:
            sizes[atom["molid"]] += 1
        self.assertEqual(sorted(sizes.values(), reverse=True), list(sizes.values()))
        self.assertEqual(sizes[1], 12)  # benzene
        self.assertEqual(sizes[2], 8)   # ethane
        self.assertEqual(sizes[3], 3)   # water

        marked = parse_data(out_mark / "frame_1.data")
        self.assertEqual(len(marked["atoms"]), marked["n_atoms"])
        self.assertEqual({atom["molid"] for atom in marked["atoms"]}, {0, 6})
        self.assertEqual(sum(atom["molid"] == 6 for atom in marked["atoms"]), 6)

        marked_xyz = read_extxyz(out_mark / "small_molecules.marked.xyz")[0]
        ring_sizes = Counter(atom["ringsize"] for atom in marked_xyz["atoms"])
        self.assertEqual(ring_sizes, Counter({"6": 6, "0": 17}))

        no_ring = parse_data(out_noring / "frame_1.data")
        self.assertEqual({atom["molid"] for atom in no_ring["atoms"]}, {0})
        no_ring_xyz = read_extxyz(out_noring / "small_molecules.marked.xyz")[0]
        self.assertEqual({atom["ringsize"] for atom in no_ring_xyz["atoms"]}, {"0"})

        try:
            import warnings
            warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
            from ovito.io import import_file
        except Exception as exc:  # pragma: no cover - depends on local OVITO
            self.skipTest(f"OVITO Python not available: {exc}")

        pipeline = import_file(str(out_mark / "small_molecules.marked.xyz"))
        data = pipeline.compute(0)
        self.assertEqual(data.particles.count, 23)
        ring = [int(data.particles["ringsize"][i]) for i in range(23)]
        self.assertEqual(ring[:6], [6] * 6)
        self.assertEqual(set(ring[6:]), {0})
        carbon = [int(data.particles["molC"][i]) for i in range(23)]
        self.assertEqual(carbon[:12], [6] * 12)
        self.assertEqual(carbon[12:20], [2] * 8)
        self.assertEqual(carbon[20:], [0] * 3)

    def test_origin_preserved(self):
        source = self.root / "origin.xyz"
        source.write_text(
            "3\n"
            'Lattice="10 0 0 0 10 0 0 0 10" Origin="5 5 5" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n'
            "O 6 5 5\nH 7 5 5\nH 6 6 5\n"
        )
        out = self.root / "out_origin"
        self.run_core(source, out, ["--dumpxyz", "1"])
        text = (out / "origin.marked.xyz").read_text()
        self.assertIn('Origin="5.00000000 5.00000000 5.00000000"', text.splitlines()[1])
        first_atom = text.splitlines()[2].split()
        self.assertEqual(first_atom[0], "O")
        self.assertEqual(first_atom[1:4], ["6.00000000", "5.00000000", "5.00000000"])

    def test_dumpxyz_fields_and_frames(self):
        frames = frames_from_combustion(5)
        trajectory = self.root / "five_frames.xyz"
        write_extxyz(trajectory, frames)
        out = self.root / "out_dumpxyz"
        self.run_core(trajectory, out, ["--dumpxyz", "2"])
        marked = out / "five_frames.marked.xyz"
        self.assertTrue(marked.exists())
        parsed_frames = read_extxyz(marked)
        self.assertEqual(len(parsed_frames), 3)

        with (out / "species_count.csv").open() as fp:
            species_rows = list(csv.DictReader(fp))
        molecule_formula = {}
        with (out / "molecules.csv").open() as fp:
            for row in csv.DictReader(fp):
                molecule_formula[row["id"]] = row["formula"]

        for parsed_index, source_index in enumerate([0, 2, 4]):
            parsed = parsed_frames[parsed_index]
            self.assertEqual(len(parsed["atoms"]), 1500)
            groups = defaultdict(list)
            for atom in parsed["atoms"]:
                groups[int(atom["molecule"])].append(atom)
            self.assertTrue(groups)
            formula_counts = Counter()
            for atoms in groups.values():
                formulas = {atom["formula"] for atom in atoms}
                self.assertEqual(len(formulas), 1)
                formula = formulas.pop()
                expected = parse_formula(formula)
                for atom in atoms:
                    self.assertEqual(int(atom["molsize"]), len(atoms))
                    for element, count in expected.items():
                        self.assertEqual(int(atom[f"mol{element}"]), count)
                self.assertAlmostEqual(float(atoms[0]["molmass"]),
                    sum(MASSES[element] * count for element, count in expected.items()), places=4)
                formula_counts[formula] += 1

            expected_counts = Counter()
            for molecule_id, count in species_rows[source_index].items():
                if molecule_id == "frame" or int(count) == 0:
                    continue
                expected_counts[molecule_formula[molecule_id]] += int(count)
            self.assertEqual(formula_counts, expected_counts)

            # Numeric fields must remain readable by OVITO even though OVITO
            # 3.15 skips the optional string formula column.
            self.assertEqual({atom["ringsize"] for atom in parsed["atoms"]}, {"0"})

        try:
            import warnings
            warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
            from ovito.io import import_file
        except Exception as exc:  # pragma: no cover - depends on local OVITO
            self.skipTest(f"OVITO Python not available: {exc}")

        pipeline = import_file(str(marked))
        self.assertEqual(pipeline.source.num_frames, 3)
        for index in range(pipeline.source.num_frames):
            data = pipeline.compute(index)
            self.assertEqual(data.particles.count, 1500)
            mol_ids = {int(data.particles["Molecule Identifier"][i]) for i in range(data.particles.count)}
            self.assertEqual(len(mol_ids), len({int(a["molecule"]) for a in parsed_frames[index]["atoms"]}))
            self.assertIn("molC", [prop.name for prop in data.particles.properties])
            self.assertIn("molH", [prop.name for prop in data.particles.properties])
            self.assertIn("ringsize", [prop.name for prop in data.particles.properties])


if __name__ == "__main__":
    unittest.main()
