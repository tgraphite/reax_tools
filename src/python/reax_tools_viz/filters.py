"""Filtering and cleaning helpers for ReaxTools tables."""

from __future__ import annotations

from collections import Counter

import networkx as nx
import pandas as pd


def filter_transfer_flow(
    flow: pd.DataFrame,
    *,
    include_self_loops: bool = False,
    top_n: int | None = None,
    min_count: int | None = None,
    min_atom_transfer: int | None = None,
    sort_by: str = "count",
    cancel_reverse: bool = True,
    max_molecules: int | None = None,
    max_subgraphs: int | None = None,
) -> pd.DataFrame:
    result = flow.copy()
    if not include_self_loops and "self_loop" in result.columns:
        result = result[result["self_loop"].astype(int) == 0]
    for column in ["count", "atom_transfer"]:
        result[column] = result[column].astype(int)
    if sort_by not in result.columns:
        sort_by = "count"
    if cancel_reverse:
        result = cancel_reverse_transfer_edges(result, direction_by=sort_by)
    if min_count is not None:
        result = result[result["count"].astype(int) >= min_count]
    if min_atom_transfer is not None:
        result = result[result["atom_transfer"].astype(int) >= min_atom_transfer]
    result = result.sort_values(sort_by, ascending=False)
    if top_n is not None:
        result = result.head(top_n)
    if max_molecules is not None:
        result = filter_flow_by_top_molecules(result, max_molecules, weight=sort_by)
    if max_subgraphs is not None:
        result = filter_flow_by_largest_subgraphs(result, max_subgraphs)
    return result.reset_index(drop=True)


def prepare_transfer_flow(
    flow: pd.DataFrame,
    *,
    include_self_loops: bool = False,
    cancel_reverse: bool = True,
    direction_by: str = "count",
) -> pd.DataFrame:
    result = flow.copy()
    if not include_self_loops and "self_loop" in result.columns:
        result = result[result["self_loop"].astype(int) == 0]
    for column in ["count", "atom_transfer"]:
        result[column] = result[column].astype(int)
    if cancel_reverse:
        result = cancel_reverse_transfer_edges(result, direction_by=direction_by)
    return result.reset_index(drop=True)


def cancel_reverse_transfer_edges(flow: pd.DataFrame, *, direction_by: str = "count") -> pd.DataFrame:
    if direction_by not in {"count", "atom_transfer"}:
        direction_by = "count"
    fallback_by = "atom_transfer" if direction_by == "count" else "count"
    records = {}
    for _, row in flow.iterrows():
        source = str(row["source_id"])
        target = str(row["target_id"])
        if source == target:
            continue
        key = tuple(sorted((source, target)))
        sign = 1 if (source, target) == key else -1
        if key not in records:
            flip = source != key[0]
            def pick(column: str, default: str = "") -> str:
                return str(row.get(column, default))
            records[key] = {
                "source_id": key[0],
                "target_id": key[1],
                "source_label": pick("source_label") if not flip else pick("target_label"),
                "target_label": pick("target_label") if not flip else pick("source_label"),
                "source_formula": pick("source_formula") if not flip else pick("target_formula"),
                "target_formula": pick("target_formula") if not flip else pick("source_formula"),
                "source_smiles": pick("source_smiles") if not flip else pick("target_smiles"),
                "target_smiles": pick("target_smiles") if not flip else pick("source_smiles"),
                "count": 0,
                "atom_transfer": 0,
                "self_loop": 0,
            }
        records[key]["count"] += sign * int(row["count"])
        records[key]["atom_transfer"] += sign * int(row["atom_transfer"])
    rows = []
    for rec in records.values():
        direction_value = rec[direction_by] or rec[fallback_by]
        if direction_value == 0:
            continue
        if direction_value < 0:
            rec = {
                "source_id": rec["target_id"],
                "target_id": rec["source_id"],
                "source_label": rec["target_label"],
                "target_label": rec["source_label"],
                "source_formula": rec.get("target_formula", ""),
                "target_formula": rec.get("source_formula", ""),
                "source_smiles": rec.get("target_smiles", ""),
                "target_smiles": rec.get("source_smiles", ""),
                "count": abs(rec["count"]),
                "atom_transfer": abs(rec["atom_transfer"]),
                "self_loop": 0,
            }
        else:
            rec["count"] = abs(rec["count"])
            rec["atom_transfer"] = abs(rec["atom_transfer"])
        rows.append(rec)
    columns = [column for column in flow.columns if column in {
        "source_id",
        "target_id",
        "source_label",
        "target_label",
        "source_formula",
        "target_formula",
        "source_smiles",
        "target_smiles",
        "count",
        "atom_transfer",
        "self_loop",
    }]
    return pd.DataFrame(rows, columns=columns or None)


def _flow_graph(flow: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    for _, row in flow.iterrows():
        graph.add_edge(str(row["source_id"]), str(row["target_id"]), weight=int(row["count"]))
        graph.nodes[str(row["source_id"])]["label"] = row["source_label"]
        graph.nodes[str(row["target_id"])]["label"] = row["target_label"]
    return graph


def _resolve_node(graph: nx.Graph, value: str) -> str | None:
    value = str(value)
    if value in graph:
        return value
    matches = [node for node, data in graph.nodes(data=True) if str(data.get("label")) == value]
    return matches[0] if matches else None


def filter_flow_by_top_molecules(flow: pd.DataFrame, max_molecules: int, *, weight: str = "count") -> pd.DataFrame:
    if flow.empty:
        return flow
    if weight not in flow.columns:
        weight = "count"
    scores = Counter()
    for _, row in flow.iterrows():
        value = int(row[weight])
        scores[str(row["source_id"])] += value
        scores[str(row["target_id"])] += value
    keep = {node for node, _ in scores.most_common(max_molecules)}
    # Induced subgraph semantics: max_molecules is an AND filter on endpoints.
    mask = flow["source_id"].astype(str).isin(keep) & flow["target_id"].astype(str).isin(keep)
    return flow[mask].copy()


def filter_flow_by_largest_subgraphs(flow: pd.DataFrame, max_subgraphs: int) -> pd.DataFrame:
    if flow.empty:
        return flow
    graph = _flow_graph(flow)
    components = sorted(nx.connected_components(graph), key=len, reverse=True)[:max_subgraphs]
    keep = set().union(*components) if components else set()
    mask = flow["source_id"].astype(str).isin(keep) & flow["target_id"].astype(str).isin(keep)
    return flow[mask].copy()


def filter_flow_by_center(flow: pd.DataFrame, center: str, depth: int) -> pd.DataFrame:
    if flow.empty:
        return flow
    graph = _flow_graph(flow)
    center_id = _resolve_node(graph, center)
    if center_id is None:
        return flow.iloc[0:0].copy()
    lengths = nx.single_source_shortest_path_length(graph, center_id, cutoff=depth)
    keep = set(lengths)
    mask = flow["source_id"].astype(str).isin(keep) & flow["target_id"].astype(str).isin(keep)
    return flow[mask].copy()


def top_flow_centers(flow: pd.DataFrame, n: int) -> list[str]:
    if flow.empty:
        return []
    graph = _flow_graph(flow)
    centrality = nx.degree_centrality(graph)
    weighted = Counter()
    for _, row in flow.iterrows():
        weight = int(row["count"])
        weighted[str(row["source_id"])] += weight
        weighted[str(row["target_id"])] += weight
    ranked = sorted(graph.nodes, key=lambda node: (centrality.get(node, 0), weighted[node]), reverse=True)
    return ranked[:n]


def _split_ids(value: str) -> list[str]:
    if not value or value == "nan":
        return []
    # Semicolon is the structural-event separator.  Fall back to "+" for
    # hand-written legacy tables that do not contain charged SMILES.
    text = str(value)
    separator = ";" if ";" in text else "+"
    return [part for part in text.split(separator) if part]


def cancel_common_ids(reactants: list[str], products: list[str]) -> tuple[list[str], list[str]]:
    left = Counter(reactants)
    right = Counter(products)
    for molecule_id in set(left) & set(right):
        n = min(left[molecule_id], right[molecule_id])
        left[molecule_id] -= n
        right[molecule_id] -= n
    new_reactants = sorted(id_ for id_, count in left.items() for _ in range(count))
    new_products = sorted(id_ for id_, count in right.items() for _ in range(count))
    return new_reactants, new_products


def aggregate_reaction_events(
    events: pd.DataFrame,
    molecules: dict[str, object] | None = None,
    *,
    labels: dict[str, str] | None = None,
    cancel_common: bool = True,
) -> pd.DataFrame:
    records: dict[tuple[tuple[str, ...], tuple[str, ...]], dict] = {}
    for _, row in events.iterrows():
        reactant_column = "reactant_ids" if "reactant_ids" in events.columns else "reactant_hashes"
        product_column = "product_ids" if "product_ids" in events.columns else "product_hashes"
        reactants = _split_ids(row[reactant_column])
        products = _split_ids(row[product_column])
        if cancel_common:
            reactants, products = cancel_common_ids(reactants, products)
        if not reactants and not products:
            continue
        key = (tuple(sorted(reactants)), tuple(sorted(products)))
        frame = int(row["frame"])
        if key not in records:
            def label_for(molecule_id: str) -> str:
                if labels and molecule_id in labels:
                    return labels[molecule_id]
                if molecules and molecule_id in molecules:
                    return str(getattr(molecules[molecule_id], "formula", molecule_id))
                return molecule_id

            reactant_labels = " + ".join(label_for(id_) for id_ in key[0])
            product_labels = " + ".join(label_for(id_) for id_ in key[1])
            records[key] = {
                "frequency": 0,
                "first_frame": frame,
                "last_frame": frame,
                "reactant_ids": ";".join(key[0]),
                "product_ids": ";".join(key[1]),
                "reactants": reactant_labels,
                "products": product_labels,
                "reactant_smiles": ";".join(
                    str(getattr(molecules[id_], "smiles", "")) if molecules and id_ in molecules else ""
                    for id_ in key[0]
                ),
                "product_smiles": ";".join(
                    str(getattr(molecules[id_], "smiles", "")) if molecules and id_ in molecules else ""
                    for id_ in key[1]
                ),
            }
        rec = records[key]
        rec["frequency"] += 1
        rec["first_frame"] = min(rec["first_frame"], frame)
        rec["last_frame"] = max(rec["last_frame"], frame)
    result = pd.DataFrame(records.values())
    if result.empty:
        return result
    return result.sort_values("frequency", ascending=False).reset_index(drop=True)
