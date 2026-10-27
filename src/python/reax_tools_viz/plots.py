"""Plotting helpers for ReaxTools count and event tables."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .chemistry import formula_to_subscript, natural_count_key, parse_atom_degree, parse_bond


COMMON_ELEMENT_ORDER = ["C", "O", "N", "H", "S", "P", "F", "Cl", "Br", "I"]


def _element_rank(element: str) -> tuple[int, str]:
    try:
        return COMMON_ELEMENT_ORDER.index(element), element
    except ValueError:
        return len(COMMON_ELEMENT_ORDER), element


def _importance_order(table: pd.DataFrame) -> list[str]:
    if table.empty:
        return []
    mean = table.mean(axis=0)
    std = table.std(axis=0).fillna(0)
    trend = (table.iloc[-1] - table.iloc[0]).abs() if len(table) > 1 else std
    span = (table.max(axis=0) - table.min(axis=0)).abs()
    variation = std + trend + 0.25 * span
    score = np.log1p(mean.clip(lower=0)) * np.log1p(variation.clip(lower=0))
    return score.sort_values(ascending=False).index.to_list()


def _valuable_species_order(table: pd.DataFrame) -> list[str]:
    ordered = _importance_order(table)
    if table.empty:
        return ordered
    mean = table.mean(axis=0)
    std = table.std(axis=0).fillna(0)
    trend = (table.iloc[-1] - table.iloc[0]).abs() if len(table) > 1 else std
    span = (table.max(axis=0) - table.min(axis=0)).abs()
    variation = std + trend + 0.25 * span
    selected = [column for column in ordered if mean[column] > 0 and variation[column] > 0]
    return selected or ordered


def _paginate(items: Sequence, min_per_page: int = 4, max_per_page: int = 8) -> list[list]:
    items = list(items)
    if not items:
        return []
    if len(items) <= max_per_page:
        return [items]
    pages = math.ceil(len(items) / max_per_page)
    while pages > 1 and len(items) / pages < min_per_page:
        pages -= 1
    base = len(items) // pages
    remainder = len(items) % pages
    result = []
    start = 0
    for idx in range(pages):
        size = base + (1 if idx < remainder else 0)
        result.append(items[start : start + size])
        start += size
    return result


def _chunk_pages(items: Sequence, max_per_page: int) -> list[list]:
    items = list(items)
    return [items[idx : idx + max_per_page] for idx in range(0, len(items), max_per_page)]


def _split_targets(targets: Sequence[str] | str | None) -> list[str] | None:
    if targets is None:
        return None
    if isinstance(targets, str):
        raw = targets.split(",")
    else:
        raw = []
        for item in targets:
            raw.extend(str(item).split(","))
    result = [item.strip() for item in raw if item.strip()]
    return result or None


def _page_output(base_output: Path, page_index: int, page_count: int) -> Path:
    if page_count == 1:
        return base_output
    return base_output.with_name(f"{base_output.stem}_page{page_index + 1}{base_output.suffix}")


def _cleanup_stale_page_outputs(base_output: Path, written: Sequence[Path]) -> None:
    expected = {path.resolve() for path in written}
    if len(written) > 1 and base_output.exists() and base_output.resolve() not in expected:
        base_output.unlink(missing_ok=True)
    for candidate in base_output.parent.glob(f"{base_output.stem}_page*{base_output.suffix}"):
        if candidate.resolve() not in expected:
            candidate.unlink(missing_ok=True)


def _bottom_legend(ax, labels: list[str], columns: int) -> None:
    handles, _ = ax.get_legend_handles_labels()
    if not handles:
        return
    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=max(1, min(columns, len(labels))),
        frameon=False,
        fontsize=9,
    )


def _plot_line_page(
    table: pd.DataFrame,
    frame_values: pd.Series,
    columns: list[str],
    labels: list[str],
    output: Path,
    *,
    figure_width: float,
    figure_height: float,
    line_width: float,
    legend_columns: int,
    ylabel: str = "Count",
) -> None:
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))
    for column, label in zip(columns, labels):
        ax.plot(frame_values, table[column], linewidth=line_width, label=label)
    ax.set_xlabel("Frame", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.tick_params(labelsize=10)
    _bottom_legend(ax, labels, legend_columns)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.94, bottom=0.30)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def _prepare_wide_table(df: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    frame_col = "frame" if "frame" in df.columns else df.columns[0]
    value_cols = [col for col in df.columns if col != frame_col]
    numeric = df[value_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    return df[frame_col], numeric


def _default_count_pages(file_name: str, table: pd.DataFrame, targets: list[str] | None, config: dict) -> list[list[tuple[str, str]]]:
    columns = table.columns.to_list()
    max_per_page = int(config.get("max_per_page", 8))
    min_per_page = int(config.get("min_per_page", 4))

    if targets:
        selected = [target for target in targets if target in columns]
        return [[(column, column) for column in page] for page in _paginate(selected, min_per_page, max_per_page)]

    if file_name == "species_count.csv":
        ordered = _valuable_species_order(table)
        return [[(column, column) for column in page] for page in _paginate(ordered, min_per_page, max_per_page)]

    if file_name == "ring_count.csv":
        ordered = sorted(columns, key=natural_count_key)
        return [[(column, column) for column in ordered]]

    if file_name == "atom_bonded_num_count.csv":
        parsed = [(column, parse_atom_degree(column)) for column in columns]
        if len(parsed) <= max_per_page:
            ordered = [column for column, _ in sorted(parsed, key=lambda item: natural_count_key(item[0]))]
            return [[(column, column) for column in ordered]]
        nonzero = [(column, data) for column, data in parsed if data and data[1] > 0]
        zero = [(column, data) for column, data in parsed if data and data[1] == 0]
        pages: list[list[tuple[str, str]]] = []
        elements = sorted({data[0] for _, data in nonzero}, key=_element_rank)
        for element in elements:
            group = [(column, column) for column, data in nonzero if data[0] == element]
            group.sort(key=lambda item: parse_atom_degree(item[0])[1], reverse=True)
            pages.extend(_paginate(group, min_per_page, max_per_page))
        zero_group = [(column, column) for column, _ in sorted(zero, key=lambda item: _element_rank(item[1][0]))]
        if zero_group:
            pages.extend(_paginate(zero_group, min_per_page, max_per_page))
        return pages

    if file_name == "bond_count.csv":
        parsed_bonds = {column: parse_bond(column) for column in columns}
        if len(columns) <= max_per_page:
            ordered = sorted(columns, key=lambda col: tuple(_element_rank(x) for x in parsed_bonds.get(col, ("", ""))))
            return [[(column, column) for column in ordered]]
        elements = sorted({element for pair in parsed_bonds.values() if pair for element in pair}, key=_element_rank)
        pages = []
        for element in elements:
            group = []
            for other in elements:
                candidates = [f"{element}-{other}", f"{other}-{element}"]
                column = next((candidate for candidate in candidates if candidate in parsed_bonds), None)
                if column is not None:
                    group.append((column, f"{element}-{other}"))
            if group:
                pages.extend(_paginate(group, min_per_page, max_per_page))
        return pages

    ordered = _importance_order(table)
    return [[(column, column) for column in page] for page in _paginate(ordered, min_per_page, max_per_page)]


def plot_count_pages(
    df: pd.DataFrame,
    output: str | Path,
    *,
    file_name: str,
    config: dict,
    targets: Sequence[str] | str | None = None,
    ylabel: str = "Count",
) -> list[Path]:
    output = Path(output)
    frames, table = _prepare_wide_table(df)
    if table.empty:
        return []

    target_list = _split_targets(targets)
    pages = _default_count_pages(file_name, table, target_list, config)
    max_pages = config.get("max_pages")
    if max_pages:
        pages = pages[: int(max_pages)]
    written: list[Path] = []
    label_subscripts = bool(config.get("label_subscripts", file_name == "species_count.csv"))
    for idx, page in enumerate(pages):
        columns = [column for column, _ in page]
        labels = [label for _, label in page]
        display_labels = [formula_to_subscript(label) if label_subscripts else label for label in labels]
        page_output = _page_output(output, idx, len(pages))
        _plot_line_page(
            table,
            frames,
            columns,
            display_labels,
            page_output,
            figure_width=float(config.get("figure_width", 6.0)),
            figure_height=float(config.get("figure_height", 4.0)),
            line_width=float(config.get("line_width", 1.8)),
            legend_columns=int(config.get("legend_columns", 4)),
            ylabel=ylabel,
        )
        written.append(page_output)
    _cleanup_stale_page_outputs(output, written)
    return written


def plot_species_count_long(species_count: pd.DataFrame, output: str | Path, top_n: int = 8) -> list[Path]:
    config = {"max_per_page": top_n, "min_per_page": 4, "max_pages": 5, "label_subscripts": True}
    return plot_count_pages(species_count, output, file_name="species_count.csv", config=config)


def plot_wide_count_table(df: pd.DataFrame, output: str | Path, top_n: int = 8, ylabel: str = "Count") -> None:
    config = {"max_per_page": top_n, "min_per_page": 4, "label_subscripts": False}
    plot_count_pages(df, output, file_name=Path(output).name.replace(".png", ".csv"), config=config, ylabel=ylabel)


def plot_reaction_events(events: pd.DataFrame, output: str | Path, top_n: int | None = None, config: dict | None = None) -> list[Path]:
    if events.empty:
        return []
    config = config or {}
    output = Path(output)
    max_per_page = int(config.get("max_per_page", 20))
    max_pages = int(config.get("max_pages", 5))
    total = top_n if top_n is not None else int(config.get("max_reactions", max_per_page * max_pages))
    total = min(int(total), max_per_page * max_pages)
    top = events.sort_values("frequency", ascending=False).head(total).copy()
    rows = list(top.iterrows())
    pages = _chunk_pages(rows, max_per_page)[:max_pages]
    written: list[Path] = []
    for page_index, page in enumerate(pages):
        labels = [
            f"{row.get('reactants', row.get('reactant_labels', row['reactant_ids']))} -> "
            f"{row.get('products', row.get('product_labels', row['product_ids']))}"
            for _, row in page
        ]
        labels = [formula_to_subscript(label) for label in labels]
        values = np.array([int(row["frequency"]) for _, row in page])

        fig_height = max(float(config.get("min_height", 4.0)), float(config.get("row_height", 0.34)) * len(labels) + 1.0)
        fig, ax = plt.subplots(figsize=(float(config.get("figure_width", 7.0)), fig_height))
        y = np.arange(len(labels))
        cmap = plt.get_cmap(str(config.get("color_map", "YlGnBu")))
        norm = values / values.max() if values.max() else values
        colors = [cmap(0.28 + 0.55 * value) for value in norm]
        alpha = float(config.get("alpha", 0.72))
        ax.barh(y, values, color=colors, alpha=alpha)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Frequency", fontsize=11)
        ax.set_ylabel("")
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=9)
        fig.subplots_adjust(left=0.48, right=0.96, top=0.96, bottom=0.12)
        page_output = _page_output(output, page_index, len(pages))
        fig.savefig(page_output, dpi=300)
        plt.close(fig)
        written.append(page_output)
    _cleanup_stale_page_outputs(output, written)
    return written
