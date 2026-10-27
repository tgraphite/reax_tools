"""Readers for ReaxTools raw output directories."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from .chemistry import build_molecule_labels
from .models import MoleculeRecord, RunBundle


REQUIRED_FILES = {
    "reax_tools.log",
    "molecules.csv",
    "species_count.csv",
    "bond_count.csv",
    "atom_bonded_num_count.csv",
    "ring_count.csv",
    "reaction_events.csv",
    "reaction_event_pairs.csv",
    "transfer_flow.csv",
}


def parse_run_log(path: str | Path) -> dict[str, dict[str, str]]:
    """Parse the human-readable run log into sections of key/value pairs."""
    sections: dict[str, dict[str, str]] = {}
    current: str | None = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1]
            sections.setdefault(current, {})
            continue
        if current is not None and "=" in line:
            key, value = line.split("=", 1)
            sections[current][key.strip()] = value.strip()
    return sections


def load_bundle(output_dir: str | Path, validate: bool = True) -> RunBundle:
    root = Path(output_dir)
    log = parse_run_log(root / "reax_tools.log")

    molecules: dict[str, MoleculeRecord] = {}
    with (root / "molecules.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            molecule_id = str(row["id"])
            if molecule_id in molecules:
                raise ValueError(f"molecules.csv contains duplicate id {molecule_id}")
            molecules[molecule_id] = MoleculeRecord(
                id=molecule_id,
                formula=str(row.get("formula", "")),
                smiles=str(row.get("smiles", "")),
            )

    bundle = RunBundle(
        root=root,
        log=log,
        molecules=molecules,
        labels=build_molecule_labels(molecules),
    )
    if validate:
        validate_bundle(bundle)
    return bundle


def validate_bundle(bundle: RunBundle) -> None:
    missing = [name for name in REQUIRED_FILES if not (bundle.root / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing ReaxTools output files: {', '.join(missing)}")

    molecule_ids = set(bundle.molecules)

    species = pd.read_csv(bundle.root / "species_count.csv", dtype=str)
    if "frame" not in species.columns:
        raise ValueError("species_count.csv must start with a frame column")
    unknown = set(species.columns[1:]) - molecule_ids
    if unknown:
        raise ValueError(f"species_count.csv references unknown molecule ids: {sorted(unknown)[:10]}")

    flow = pd.read_csv(
        bundle.root / "transfer_flow.csv", dtype={"source_id": str, "target_id": str}
    )
    unknown = (set(flow["source_id"]) | set(flow["target_id"])) - molecule_ids
    if unknown:
        raise ValueError(f"transfer_flow.csv references unknown molecule ids: {sorted(unknown)[:10]}")

    events = pd.read_csv(bundle.root / "reaction_events.csv", dtype=str)
    for column in [
        "reactant_ids",
        "product_ids",
        "reactant_formulas",
        "product_formulas",
        "reactant_smiles",
        "product_smiles",
    ]:
        if column not in events.columns:
            raise ValueError(f"reaction_events.csv missing {column}")

    unknown_event_ids: set[str] = set()
    for column in ["reactant_ids", "product_ids"]:
        for value in events[column].dropna():
            unknown_event_ids.update(part for part in str(value).split(";") if part and part not in molecule_ids)
    if unknown_event_ids:
        raise ValueError(f"reaction_events.csv references unknown molecule ids: {sorted(unknown_event_ids)[:10]}")


def load_molecule_labels(root: str | Path) -> dict[str, str]:
    """Load only molecules.csv and return ID -> display label mapping."""
    path = Path(root) / "molecules.csv"
    if not path.exists():
        return {}
    molecules: dict[str, MoleculeRecord] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            molecule_id = str(row["id"])
            molecules[molecule_id] = MoleculeRecord(
                id=molecule_id,
                formula=str(row.get("formula", "")),
                smiles=str(row.get("smiles", "")),
            )
    return build_molecule_labels(molecules)


def read_count_table(bundle: RunBundle, name: str) -> pd.DataFrame:
    return pd.read_csv(bundle.root / name)


def read_transfer_flow(bundle: RunBundle | None) -> pd.DataFrame:
    if bundle is None:
        return pd.read_csv("transfer_flow.csv", dtype={"source_id": str, "target_id": str})

    flow = pd.read_csv(
        bundle.root / "transfer_flow.csv",
        dtype={"source_id": str, "target_id": str},
    )
    flow["source_label"] = flow["source_id"].map(bundle.labels).fillna(
        flow.get("source_formula", flow["source_id"])
    )
    flow["target_label"] = flow["target_id"].map(bundle.labels).fillna(
        flow.get("target_formula", flow["target_id"])
    )
    return flow


def read_transfer_flow_file(path: str | Path) -> pd.DataFrame:
    flow = pd.read_csv(path, dtype={"source_id": str, "target_id": str})
    if "source_label" not in flow.columns:
        flow["source_label"] = flow["source_formula"] if "source_formula" in flow.columns else flow["source_id"]
    if "target_label" not in flow.columns:
        flow["target_label"] = flow["target_formula"] if "target_formula" in flow.columns else flow["target_id"]
    return flow


def read_reaction_events(bundle: RunBundle) -> pd.DataFrame:
    return pd.read_csv(bundle.root / "reaction_events.csv", dtype=str)
