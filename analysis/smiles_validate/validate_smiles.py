"""Validate the native C++ canonical SMILES writer against RDKit.

Builds explicit-atom graphs, sends them to bin/smiles_cli, then checks:
  1. RDKit can parse the produced SMILES;
  2. the parsed graph is isomorphic to the input graph (element labels);
  3. molecular formula is unchanged;
  4. the SMILES is invariant to input atom order;
  5. molecules above the 200-atom cutoff fall back to formula identity.

Usage:
    python3 analysis/smiles_validate/validate_smiles.py
"""
from __future__ import annotations

import csv
import random
import subprocess
import sys
from pathlib import Path

import networkx as nx
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "hash_bench"))

from datasets import (  # noqa: E402
    alkane_families,
    curated_families,
    regular_classes,
)

RDLogger.DisableLog("rdApp.*")
CLI = ROOT / "bin" / "smiles_cli"
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)


def explicit_graph(mol: Chem.Mol):
    mol = Chem.AddHs(mol)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    return atoms, bonds


def run_cli(atoms, bonds):
    lines = [f"{len(atoms)} {len(bonds)}"]
    lines.extend(atoms)
    lines.extend(f"{i} {j}" for i, j in bonds)
    proc = subprocess.run(
        [str(CLI)],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip())
    out = proc.stdout.strip()
    if out.startswith("SMILES "):
        return out[len("SMILES ") :], False
    if out.startswith("FALLBACK"):
        return "", True
    raise RuntimeError(f"unexpected CLI output: {out!r}")


def to_nx(atoms, bonds):
    graph = nx.Graph()
    for i, element in enumerate(atoms):
        graph.add_node(i, element=element)
    for i, j in bonds:
        graph.add_edge(i, j)
    return graph


def all_single_reference(atoms, bonds):
    rw = Chem.RWMol()
    for element in atoms:
        atom = Chem.Atom(element)
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    for i, j in bonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def validate_case(name, mol: Chem.Mol):
    atoms, bonds = explicit_graph(mol)
    original_formula = rdMolDescriptors.CalcMolFormula(Chem.AddHs(mol))
    original_graph = to_nx(atoms, bonds)

    smiles, fallback = run_cli(atoms, bonds)
    if fallback:
        return {
            "name": name,
            "atoms": len(atoms),
            "status": "fallback",
            "smiles": "",
            "reason": "over cutoff" if len(atoms) > 200 else "unexpected",
        }
    if len(atoms) > 200:
        return {
            "name": name,
            "atoms": len(atoms),
            "status": "fail",
            "smiles": smiles,
            "reason": "expected fallback above 200 atoms",
        }

    parsed = Chem.MolFromSmiles(smiles)
    if parsed is None:
        return {
            "name": name,
            "atoms": len(atoms),
            "status": "fail",
            "smiles": smiles,
            "reason": "RDKit parse failed",
        }

    parsed_h = Chem.AddHs(parsed)
    p_atoms = [atom.GetSymbol() for atom in parsed_h.GetAtoms()]
    p_bonds = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in parsed_h.GetBonds()
    ]
    parsed_graph = to_nx(p_atoms, p_bonds)

    reasons = []
    if not nx.is_isomorphic(
        original_graph, parsed_graph, node_match=lambda a, b: a["element"] == b["element"]
    ):
        reasons.append("graph mismatch")
    parsed_formula = rdMolDescriptors.CalcMolFormula(parsed)
    if parsed_formula != original_formula:
        reasons.append(f"formula {parsed_formula} != {original_formula}")

    reference = all_single_reference(atoms, bonds)
    roundtrip = Chem.MolToSmiles(parsed_h)
    if reference is not None and roundtrip != reference:
        reasons.append("round-trip canonical SMILES mismatch")

    # Atom-order invariance.
    rng = random.Random(20261027)
    for _ in range(3):
        order = list(range(len(atoms)))
        rng.shuffle(order)
        inverse = [0] * len(order)
        for new, old in enumerate(order):
            inverse[old] = new
        s_atoms = [atoms[old] for old in order]
        s_bonds = [(inverse[i], inverse[j]) for i, j in bonds]
        s_smiles, s_fallback = run_cli(s_atoms, s_bonds)
        if s_fallback or s_smiles != smiles:
            reasons.append("permutation mismatch")
            break

    return {
        "name": name,
        "atoms": len(atoms),
        "status": "pass" if not reasons else "fail",
        "smiles": smiles,
        "reason": "; ".join(reasons),
    }


def build_corpus():
    cases = []
    for formula, mols in curated_families().items():
        for i, mol in enumerate(mols):
            cases.append((f"{formula}_{i}", mol))
    for formula, mols in alkane_families([7, 8, 9, 10]).items():
        for i, mol in enumerate(mols):
            cases.append((f"{formula}_{i}", mol))
    for n, d in [(6, 3), (8, 3), (10, 3), (12, 3)]:
        for i, mol in enumerate(regular_classes(n, d, range(0, 200))):
            cases.append((f"regular_d{d}_n{n}_{i}", mol))

    radical_smiles = [
        "[CH3]", "[CH2]", "[CH]", "[C]", "[OH]", "[NH2]", "[NH]", "[N]",
        "[O]", "[H]", "[H][H]", "[O][O]", "[N][N]", "[NO]", "[NO2]",
    ]
    for i, smiles in enumerate(radical_smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is not None:
            cases.append((f"radical_{i}", mol))

    ring_smiles = [
        "c1ccccc1", "C1CCCCC1", "C12C3C4C1C5C2C3C45", "c1ccc2ccccc2c1",
        "C1CCC2CCCCC2C1", "C1CC2CCC1C2", "C1CC1", "C1CCC1",
    ]
    for i, smiles in enumerate(ring_smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            cases.append((f"ring_{i}", mol))

    # Size cutoff cases: C66H134 = 200 atoms, C67H136 = 203 atoms.
    cases.append(("cutoff_200_atoms", Chem.MolFromSmiles("C" * 66)))
    cases.append(("cutoff_203_atoms", Chem.MolFromSmiles("C" * 67)))

    return cases


def main():
    if not CLI.exists():
        raise SystemExit(f"missing {CLI}; build the smiles_cli target first")

    rows = []
    for name, mol in build_corpus():
        try:
            rows.append(validate_case(name, mol))
        except Exception as exc:  # pragma: no cover - diagnostic path
            rows.append(
                {
                    "name": name,
                    "atoms": mol.GetNumAtoms() if mol else 0,
                    "status": "error",
                    "smiles": "",
                    "reason": str(exc),
                }
            )

    with open(OUT / "validation.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "atoms", "status", "smiles", "reason"])
        writer.writeheader()
        writer.writerows(rows)

    counts = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print("cases:", len(rows), counts)
    failures = [row for row in rows if row["status"] != "pass"]
    for row in failures[:40]:
        print(f"  {row['status']:8s} {row['name']:28s} {row['reason']}")
    print(f"results: {OUT / 'validation.csv'}")
    return 0 if not any(row["status"] in ("fail", "error") for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
