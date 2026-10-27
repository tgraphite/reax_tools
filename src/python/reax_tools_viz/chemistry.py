"""Small chemistry-facing formatting helpers."""

from __future__ import annotations

import re
from collections import defaultdict


SUBSCRIPT_TRANS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def formula_to_subscript(formula: object) -> str:
    """Render only digits in formula-like labels as unicode subscripts."""
    return str(formula).translate(SUBSCRIPT_TRANS)


def build_molecule_labels(molecules: dict) -> dict[str, str]:
    """Map structure IDs to display labels.

    A formula with one structure keeps the formula. A formula shared by
    several constitutional isomers gets a letter suffix, e.g. C2H6O (a).
    """
    by_formula: dict[str, list[str]] = defaultdict(list)
    for molecule_id, record in molecules.items():
        by_formula[str(record.formula)].append(str(molecule_id))

    labels: dict[str, str] = {}
    for formula, ids in by_formula.items():
        ids.sort(key=lambda value: int(value) if value.isdigit() else value)
        if len(ids) == 1:
            labels[ids[0]] = formula
            continue
        for index, molecule_id in enumerate(ids):
            suffix = chr(ord("a") + index) if index < 26 else f"i{index + 1}"
            labels[molecule_id] = f"{formula} ({suffix})"
    return labels


def reaction_to_subscript(label: object) -> str:
    return " + ".join(formula_to_subscript(part) for part in str(label).split("+"))


def natural_count_key(label: object) -> tuple:
    text = str(label)
    parts = re.split(r"(\d+)", text)
    return tuple(int(part) if part.isdigit() else part for part in parts)


def parse_atom_degree(label: str) -> tuple[str, int] | None:
    match = re.fullmatch(r"([A-Z][a-z]?)-(\d+)", str(label))
    if not match:
        return None
    return match.group(1), int(match.group(2))


def parse_bond(label: str) -> tuple[str, str] | None:
    match = re.fullmatch(r"([A-Z][a-z]?)-([A-Z][a-z]?)", str(label))
    if not match:
        return None
    return match.group(1), match.group(2)
