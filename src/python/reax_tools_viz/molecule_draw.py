"""Optional 2D molecule drawing from canonical SMILES.

RDKit is an optional dependency of the plotting layer.  It is not required by
the C++ core or by the count/network/event commands.
"""

from __future__ import annotations

from pathlib import Path


def _rdkit():
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "RDKit is required for the 'molecules' command: pip install rdkit"
        ) from exc
    return rdMolDraw2D


def draw_molecule_record(record: dict, output: str | Path) -> None:
    smiles = str(record.get("smiles", "") or "")
    if not smiles:
        raise ValueError(
            f"molecule {record.get('id')} has no SMILES (formula-level fallback)"
        )

    from .drawing import prepare_drawing_mol

    rdMolDraw2D = _rdkit()
    mol = prepare_drawing_mol(smiles)
    if mol is None:
        raise ValueError(f"molecule {record.get('id')} has an unparsable SMILES: {smiles}")

    drawer = rdMolDraw2D.MolDraw2DCairo(900, 900)
    drawer.drawOptions().padding = 0.06
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    Path(output).write_bytes(drawer.GetDrawingText())


def draw_selected_molecules(
    molecules: list[dict],
    output_dir: str | Path,
    species: list[str] | None = None,
) -> list[Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    wanted = {str(item) for item in (species or [])}
    written: list[Path] = []

    for record in molecules:
        molecule_id = str(record.get("id", ""))
        formula = str(record.get("formula", ""))
        smiles = str(record.get("smiles", "") or "")
        if wanted and (molecule_id not in wanted and formula not in wanted and smiles not in wanted):
            continue
        if not smiles:
            print(f"Skipping {molecule_id}: formula-level fallback has no SMILES")
            continue
        filename = f"{formula or molecule_id}_{molecule_id}.png"
        path = out / filename
        draw_molecule_record(record, path)
        written.append(path)
    return written
