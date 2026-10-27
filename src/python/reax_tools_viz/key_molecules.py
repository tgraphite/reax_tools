"""Importance ranking and 2D montage for key molecular structures.

Importance combines the time-averaged species population with the weighted
degree in the reaction-transfer graph. Labels come from the run-global
formula/isomer assignment (formula (a), formula (b), ...).
"""

from __future__ import annotations

import io
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

from .drawing import prepare_drawing_mol
from .io import read_transfer_flow


def _rank_normalize(values: pd.Series) -> pd.Series:
    if values.empty:
        return values
    if len(values) == 1:
        return pd.Series([1.0], index=values.index)
    ranks = values.rank(method="average", ascending=False)
    return 1.0 - (ranks - 1.0) / (len(values) - 1.0)


def select_key_molecules(
    bundle,
    *,
    n: int = 48,
    weights: dict | None = None,
    degree_weight: str = "atom_transfer",
    abundance: str = "mean",
    exclude_self_loops: bool = False,
) -> list[dict]:
    """Return the top-n structures by abundance + reaction-graph degree."""
    species = pd.read_csv(bundle.root / "species_count.csv", dtype=str)
    if species.empty or "frame" not in species.columns:
        return []
    counts = species.drop(columns=["frame"]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if abundance == "sum":
        abundance_values = counts.sum(axis=0)
    elif abundance == "max":
        abundance_values = counts.max(axis=0)
    elif abundance == "presence":
        abundance_values = (counts > 0).mean(axis=0)
    else:
        abundance_values = counts.mean(axis=0)

    degree = pd.Series(0.0, index=abundance_values.index)
    flow = read_transfer_flow(bundle)
    weight_column = degree_weight if degree_weight in flow.columns else "atom_transfer"
    for _, row in flow.iterrows():
        source = str(row["source_id"])
        target = str(row["target_id"])
        if source == target:
            if exclude_self_loops:
                continue
        elif exclude_self_loops and str(row.get("self_loop", "0")) == "1":
            continue
        weight = float(row.get(weight_column, 0.0) or 0.0)
        if source in degree.index:
            degree[source] += weight
        if target in degree.index:
            degree[target] += weight

    candidates = [
        molecule_id
        for molecule_id in abundance_values.index
        if molecule_id in bundle.molecules and bundle.molecules[molecule_id].smiles
    ]
    if not candidates:
        return []

    abundance_values = abundance_values[candidates]
    degree = degree[candidates]
    abundance_rank = _rank_normalize(abundance_values)
    degree_rank = _rank_normalize(degree)
    weights = dict(weights or {})
    weight_abundance = float(weights.get("abundance", 0.6))
    weight_degree = float(weights.get("degree", 0.4))
    score = weight_abundance * abundance_rank + weight_degree * degree_rank

    ordered = sorted(
        candidates,
        key=lambda molecule_id: (
            -float(score[molecule_id]),
            -float(abundance_values[molecule_id]),
            -float(degree[molecule_id]),
            int(molecule_id) if molecule_id.isdigit() else molecule_id,
        ),
    )[: int(n)]

    records = []
    for rank, molecule_id in enumerate(ordered, start=1):
        molecule = bundle.molecules[molecule_id]
        records.append(
            {
                "rank": rank,
                "id": molecule_id,
                "label": bundle.labels.get(molecule_id, molecule.formula),
                "formula": molecule.formula,
                "smiles": molecule.smiles,
                "abundance": float(abundance_values[molecule_id]),
                "degree": float(degree[molecule_id]),
                "score": float(score[molecule_id]),
            }
        )
    return records


def _draw_smiles_png(smiles: str, size: int = 420) -> bytes | None:
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = prepare_drawing_mol(smiles)
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
    drawer.drawOptions().padding = 0.06
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def plot_key_molecules(bundle, output_dir: str | Path, config: dict | None = None) -> list[Path]:
    """Write key_molecules_page1.png ... for the top-N structures.

    RDKit is an optional dependency. If it is missing, print a note and return
    an empty list so the rest of plot still succeeds.
    """
    config = dict(config or {})
    n = int(config.get("n", 48))
    per_page = int(config.get("per_page", 16))
    weights = config.get("weights", {}) or {}
    records = select_key_molecules(
        bundle,
        n=n,
        weights=weights,
        degree_weight=str(config.get("degree_weight", "atom_transfer")),
        abundance=str(config.get("abundance", "mean")),
        exclude_self_loops=bool(config.get("exclude_self_loops", False)),
    )
    if not records:
        print("Key molecules: no structures with SMILES; skipping")
        return []

    try:
        import rdkit  # noqa: F401
    except ImportError:
        print("Key molecules: RDKit is not installed; skipping 2D montage")
        return []

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    pages = math.ceil(len(records) / per_page)
    columns = 4
    rows = math.ceil(per_page / columns)
    written: list[Path] = []

    for page_index in range(pages):
        page_records = records[page_index * per_page : (page_index + 1) * per_page]
        fig, axes = plt.subplots(rows, columns, figsize=(3.2 * columns, 3.4 * rows))
        axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
        for ax in axes:
            ax.axis("off")
        for ax, record in zip(axes, page_records):
            png = _draw_smiles_png(record["smiles"])
            if png is None:
                ax.text(0.5, 0.5, "unparsable SMILES", ha="center", va="center", fontsize=9)
            else:
                ax.imshow(mpimg.imread(io.BytesIO(png), format="png"))
            ax.set_title(
                f"{record['label']}\n#{record['rank']}  score {record['score']:.2f}",
                fontsize=10,
            )
        fig.suptitle(f"Key molecules (page {page_index + 1}/{pages})", fontsize=14)
        fig.tight_layout()
        path = out / f"key_molecules_page{page_index + 1}.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        written.append(path)
    return written

