"""Data models for ReaxTools Python processing.

These classes describe loaded C++ outputs. They intentionally contain no file
parsing, plotting, or filtering logic.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class MoleculeRecord:
    id: str
    formula: str
    smiles: str


@dataclass(frozen=True)
class TransferEdge:
    source_id: str
    target_id: str
    source_label: str
    target_label: str
    source_formula: str
    target_formula: str
    source_smiles: str
    target_smiles: str
    count: int
    atom_transfer: int
    self_loop: bool


@dataclass(frozen=True)
class ReactionEvent:
    event_id: int
    frame: int
    reactant_ids: tuple[str, ...]
    product_ids: tuple[str, ...]
    reactant_formulas: tuple[str, ...]
    product_formulas: tuple[str, ...]
    reactant_smiles: tuple[str, ...]
    product_smiles: tuple[str, ...]


@dataclass(frozen=True)
class RunBundle:
    root: Path
    log: dict[str, dict[str, str]]
    molecules: dict[str, MoleculeRecord]
    labels: dict[str, str] = field(default_factory=dict)
