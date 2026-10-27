"""Command-line entry point for the ReaxTools Python layer."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from .export import export_transfer_flow_dot
from .filters import (
    aggregate_reaction_events,
    filter_transfer_flow,
)
from .io import (
    load_bundle,
    load_molecule_labels,
    read_reaction_events,
    read_transfer_flow,
    read_transfer_flow_file,
)
from .molecule_draw import draw_selected_molecules
from .network import build_transfer_graph, draw_transfer_graph
from .key_molecules import plot_key_molecules
from .plots import plot_count_pages, plot_reaction_events
from .sankey import draw_focus_sankey_set, draw_role_flow
from .templates import load_plot_template


COUNT_FILES = [
    "species_count.csv",
    "bond_count.csv",
    "atom_bonded_num_count.csv",
    "ring_count.csv",
]


def _input_path(value: str) -> Path:
    return Path(value)


def _load_bundle_or_none(path: Path):
    if path.is_dir():
        return load_bundle(path)
    return None


def _count_output_name(path: Path) -> str:
    return path.name.replace(".csv", ".png")


def _count_config_key(file_name: str) -> str:
    return file_name.replace(".csv", "")


def _with_cli_override(config: dict, **values) -> dict:
    result = dict(config)
    for key, value in values.items():
        if value is not None:
            result[key] = value
    return result


def _resolve_config(path: Path, template: str | None) -> dict:
    output_dir = path if path.is_dir() else path.parent
    return load_plot_template(output_dir, template)


def _parse_center(values: list[str] | None, default_depth: int) -> tuple[str | None, int]:
    if not values:
        return None, default_depth
    if len(values) == 1:
        return values[0], default_depth
    return values[0], int(values[1])


def _parse_multiple_centers(values: list[str] | None, default_n: int, default_depth: int) -> tuple[int, int] | None:
    if values is None:
        return None
    if len(values) == 0:
        return default_n, default_depth
    if len(values) == 1:
        return int(values[0]), default_depth
    return int(values[0]), int(values[1])


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def plot_count_file(
    path: Path,
    output_dir: Path,
    config: dict,
    targets=None,
    top: int | None = None,
    bundle=None,
) -> list[Path]:
    if path.name == "species_count.csv":
        df = pd.read_csv(path, dtype=str)
        labels = bundle.labels if bundle is not None else load_molecule_labels(path.parent)
        if labels:
            df = df.rename(columns=labels)
    else:
        df = pd.read_csv(path)
    section = dict(config.get(_count_config_key(path.name), {}))
    if top is not None:
        section["max_per_page"] = top
    output = output_dir / _count_output_name(path)
    return plot_count_pages(df, output, file_name=path.name, config=section, targets=targets)


def _network_config(config: dict, args) -> dict:
    section = dict(config.get("transfer_flow", {}))
    result = _with_cli_override(
        section,
        layout=args.layout,
        filter_max_molecules=args.max_molecules,
        filter_max_subgraphs=args.max_subgraphs,
    )
    if str(result.get("layout", "dot")).lower() not in {"dot", "fdp", "sfdp", "neato"}:
        result["layout"] = "dot"
    return result


def _max_reactions(config: dict, args) -> int | None:
    if args.max_reactions is not None:
        return args.max_reactions
    if args.top is not None:
        return args.top
    value = config.get("transfer_flow", {}).get("filter_max_reactions", False)
    return int(value) if value else None


def _flow_sort_by(network_config: dict) -> str:
    value = str(network_config.get("weight", network_config.get("sort_by", "atom_transfer")))
    return value if value in {"count", "atom_transfer"} else "atom_transfer"


def _filtered_transfer_view(
    raw_flow: pd.DataFrame,
    network_config: dict,
    *,
    include_self_loops: bool,
    max_reactions: int | None = None,
    max_molecules: int | None = None,
    max_subgraphs: int | None = None,
    min_count: int | None = None,
    min_atom_transfer: int | None = None,
) -> pd.DataFrame:
    return filter_transfer_flow(
        raw_flow,
        include_self_loops=include_self_loops,
        top_n=max_reactions,
        max_molecules=max_molecules,
        max_subgraphs=max_subgraphs,
        min_count=min_count,
        min_atom_transfer=min_atom_transfer,
        sort_by=_flow_sort_by(network_config),
    )


def _draw_transfer_view(
    raw_flow: pd.DataFrame,
    output_dir: Path,
    network_config: dict,
    *,
    stem: str,
    title: str,
    include_self_loops: bool = False,
    max_reactions: int | None = None,
    max_molecules: int | None = None,
    max_subgraphs: int | None = None,
    min_count: int | None = None,
    min_atom_transfer: int | None = None,
    use_colors: bool = True,
) -> list[Path]:
    filtered = _filtered_transfer_view(
        raw_flow,
        network_config,
        include_self_loops=include_self_loops,
        max_reactions=max_reactions,
        max_molecules=max_molecules,
        max_subgraphs=max_subgraphs,
        min_count=min_count,
        min_atom_transfer=min_atom_transfer,
    )
    csv_path = output_dir / f"{stem}_filtered.csv"
    png_path = output_dir / f"{stem}.png"
    filtered.to_csv(csv_path, index=False)
    draw_transfer_graph(
        build_transfer_graph(filtered),
        png_path,
        layout=str(network_config.get("layout", "dot")),
        config=network_config,
        use_colors=use_colors,
        title=title,
    )
    return [csv_path, png_path]


def _draw_default_transfer_views(raw_flow: pd.DataFrame, output_dir: Path, network_config: dict, *, use_colors: bool = True) -> list[Path]:
    max_molecules = int(network_config.get("default_max_molecules", 24))
    max_reactions = int(network_config.get("default_max_reactions", 48))
    written: list[Path] = []
    written.extend(
        _draw_transfer_view(
            raw_flow,
            output_dir,
            network_config,
            stem=f"molecule_transfer_top{max_molecules}",
            title=f"Top {max_molecules} molecules by atom transfer",
            max_molecules=max_molecules,
            max_subgraphs=int(network_config.get("filter_max_subgraphs")) if network_config.get("filter_max_subgraphs") else None,
            use_colors=use_colors,
        )
    )
    written.extend(
        _draw_transfer_view(
            raw_flow,
            output_dir,
            network_config,
            stem=f"reaction_transfer_top{max_reactions}",
            title=f"Top {max_reactions} reactions by atom transfer",
            max_reactions=max_reactions,
            max_subgraphs=int(network_config.get("filter_max_subgraphs")) if network_config.get("filter_max_subgraphs") else None,
            use_colors=use_colors,
        )
    )
    return written


def _draw_network_command(args, *, legacy_flow: bool = False) -> int:
    path = _input_path(args.file)
    bundle = _load_bundle_or_none(path)
    config = _resolve_config(path, args.template)
    network_config = _network_config(config, args)
    raw_flow = read_transfer_flow(bundle) if bundle else read_transfer_flow_file(path)
    output_dir = bundle.root if bundle else path.parent
    cli_filtered = any(
        value is not None
        for value in [args.max_reactions, args.max_molecules, args.max_subgraphs, args.top, args.min_count, args.min_atom_transfer]
    ) or args.include_self_loops
    if not cli_filtered and not legacy_flow:
        for output in _draw_default_transfer_views(raw_flow, output_dir, network_config, use_colors=not args.no_colors):
            print(f"Wrote {output}")
        return 0

    max_reactions = _max_reactions(config, args)
    max_molecules = int(network_config.get("filter_max_molecules")) if network_config.get("filter_max_molecules") else None
    stem = f"reaction_transfer_top{max_reactions}" if max_reactions else f"molecule_transfer_top{max_molecules}" if max_molecules else "molecule_transfer"
    title = (
        f"Top {max_reactions} reactions by atom transfer"
        if max_reactions
        else f"Top {max_molecules} molecules by atom transfer"
        if max_molecules
        else "Molecule transfer"
    )
    written = _draw_transfer_view(
        raw_flow,
        output_dir,
        network_config,
        stem=stem,
        title=title,
        include_self_loops=args.include_self_loops,
        max_reactions=max_reactions,
        max_molecules=max_molecules,
        max_subgraphs=int(network_config.get("filter_max_subgraphs")) if network_config.get("filter_max_subgraphs") else None,
        min_count=args.min_count,
        min_atom_transfer=args.min_atom_transfer,
        use_colors=not args.no_colors,
    )
    for output in written:
        print(f"Wrote {output}")
    if legacy_flow:
        filtered = pd.read_csv(written[0])
        filtered.to_csv(output_dir / "transfer_flow_filtered.csv", index=False)
        draw_transfer_graph(
            build_transfer_graph(filtered),
            output_dir / "transfer_flow.png",
            layout=str(network_config.get("layout", "dot")),
            config=network_config,
            use_colors=not args.no_colors,
            title=title,
        )
        print(f"Wrote {output_dir / 'transfer_flow_filtered.csv'}")
        print(f"Wrote {output_dir / 'transfer_flow.png'}")
    if args.dot:
        filtered = pd.read_csv(written[0])
        export_transfer_flow_dot(filtered, output_dir / f"{stem}_filtered.dot")
        print(f"Wrote {output_dir / f'{stem}_filtered.dot'}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot and filter ReaxTools raw outputs")
    sub = parser.add_subparsers(dest="command", required=True)

    summary = sub.add_parser("summary", help="Validate and summarize a ReaxTools output directory")
    summary.add_argument("-f", "--file", required=True, help="ReaxTools output directory")

    plot = sub.add_parser("plot", help="Generate default plots from a ReaxTools output directory")
    plot.add_argument("-f", "--file", required=True, help="ReaxTools output directory")
    plot.add_argument("--template", help="YAML plot template")
    plot.add_argument("--top", type=int, help="Compatibility alias for count max objects per page")

    counts = sub.add_parser("counts", help="Plot count tables")
    counts.add_argument("-f", "--file", required=True, help="Count CSV file or ReaxTools output directory")
    counts.add_argument("-t", "--targets", nargs="*", help="Comma-separated or space-separated objects to plot")
    counts.add_argument("--template", help="YAML plot template")
    counts.add_argument("--top", type=int, help="Compatibility alias for count max objects per page")

    network = sub.add_parser("network", help="Draw molecule/reaction transfer diagrams")
    network.add_argument("-f", "--file", required=True, help="transfer_flow.csv or ReaxTools output directory")
    network.add_argument("--template", help="YAML plot template")
    network.add_argument("--top", type=int, help="Compatibility alias for --max-reactions")
    network.add_argument("--max-reactions", type=int, help="Keep the strongest N reactions")
    network.add_argument("--max-molecules", type=int, help="Keep the strongest N molecules and their direct transfers")
    network.add_argument("--max-subgraphs", type=int, help="Keep the largest N connected molecule groups")
    network.add_argument("--layout", choices=["dot", "fdp", "sfdp", "neato"])
    network.add_argument("--no-colors", action="store_true")
    network.add_argument("--min-count", type=int)
    network.add_argument("--min-atom-transfer", type=int)
    network.add_argument("--include-self-loops", action="store_true")
    network.add_argument("--dot", action="store_true", help="Also export the filtered transfer diagram as DOT source")

    flow = sub.add_parser("flow", help="Draw role-layer atom-transfer flow diagram")
    flow.add_argument("-f", "--file", required=True, help="transfer_flow.csv or ReaxTools output directory")
    flow.add_argument("--template", help="YAML plot template")
    flow.add_argument("--max-species-per-level", type=int, help="Maximum molecule nodes shown in each role layer")
    flow.add_argument("--max-links", type=int, help="Maximum forward inter-layer links shown")
    flow.add_argument(
        "--role-breaks",
        nargs="+",
        type=float,
        help="Sorted role-ratio cut points for layers, e.g. -0.65 -0.25 0.25 0.65",
    )
    flow.add_argument("--include-self-loops", action="store_true")

    focus = sub.add_parser("focus", help="Draw centered molecule inflow/outflow diagrams")
    focus.add_argument("-f", "--file", required=True, help="transfer_flow.csv or ReaxTools output directory")
    focus.add_argument("--template", help="YAML plot template")
    focus.add_argument("--centers", nargs="*", help="Formula or molecule ids to draw; default uses key molecules")
    focus.add_argument("-n", "--num-centers", type=int, help="Number of default key molecules")
    focus.add_argument("--max-in", type=int, help="Maximum incoming sources per focus graph")
    focus.add_argument("--max-out", type=int, help="Maximum outgoing targets per focus graph")
    focus.add_argument("--include-self-loops", action="store_true")

    events = sub.add_parser("events", help="Clean and plot reaction events")
    events.add_argument("-f", "--file", required=True, help="reaction_events.csv or ReaxTools output directory")
    events.add_argument("--template", help="YAML plot template")
    events.add_argument("-max", "--max", dest="max_events", type=int, help="Number of frequent reactions to plot")
    events.add_argument("--top", type=int, help="Compatibility alias for --max")
    events.add_argument("-mol", "--mol", action="store_true", help="Reserve space for molecule embedding")
    events.add_argument("--no-cancel", action="store_true")

    molecules = sub.add_parser("molecules", help="Draw structures from molecules.csv SMILES")
    molecules.add_argument("-f", "--file", required=True, help="molecules.csv or ReaxTools output directory")
    molecules.add_argument("--species", nargs="*", help="Formula, SMILES or molecule id to draw")

    args = parser.parse_args()

    if args.command == "summary":
        bundle = load_bundle(args.file)
        print(f"Output directory: {bundle.root}")
        print(f"Identity model: {bundle.log.get('run', {}).get('identity_model', 'unknown')}")
        print(f"Molecules: {len(bundle.molecules)}")
        return 0

    if args.command in {"plot", "counts"}:
        path = _input_path(args.file)
        bundle = _load_bundle_or_none(path)
        output_dir = path if path.is_dir() else path.parent
        config = _resolve_config(path, args.template)
        files = [bundle.root / name for name in COUNT_FILES] if bundle else [path]
        for count_file in files:
            if count_file.exists():
                written = plot_count_file(
                    count_file,
                    output_dir,
                    config,
                    targets=getattr(args, "targets", None),
                    top=args.top,
                    bundle=bundle,
                )
                for output in written:
                    print(f"Wrote {output}")
        if args.command == "plot":
            flow_args = argparse.Namespace(
                layout=None,
                max_molecules=None,
                max_subgraphs=None,
                max_reactions=None,
                top=None,
            )
            network_config = _network_config(config, flow_args)
            raw_flow = read_transfer_flow(bundle)
            transfer_outputs = _draw_default_transfer_views(raw_flow, bundle.root, network_config)
            raw_events = read_reaction_events(bundle)
            cleaned = aggregate_reaction_events(
                raw_events,
                bundle.molecules,
                labels=bundle.labels,
                cancel_common=True,
            )
            cleaned.to_csv(bundle.root / "reaction_events_cleaned.csv", index=False)
            event_config = config.get("reaction_events", {})
            event_outputs = plot_reaction_events(
                cleaned,
                bundle.root / "reaction_events.png",
                top_n=int(event_config.get("max_reactions", 15)),
                config=event_config,
            )
            for output in transfer_outputs:
                print(f"Wrote {output}")
            print(f"Wrote {bundle.root / 'reaction_events_cleaned.csv'}")
            for output in event_outputs:
                print(f"Wrote {output}")
            key_outputs = plot_key_molecules(
                bundle,
                bundle.root,
                config.get("key_molecules", {}),
            )
            for output in key_outputs:
                print(f"Wrote {output}")
        return 0

    if args.command == "flow":
        path = _input_path(args.file)
        bundle = _load_bundle_or_none(path)
        config = _resolve_config(path, args.template)
        flow_config = dict(config.get("flow", {}))
        raw_flow = read_transfer_flow(bundle) if bundle else read_transfer_flow_file(path)
        output_dir = bundle.root if bundle else path.parent
        for key, value in {
            "max_species_per_level": args.max_species_per_level,
            "max_links": args.max_links,
            "role_breaks": args.role_breaks,
        }.items():
            if value is not None:
                flow_config[key] = value
        if draw_role_flow(raw_flow, output_dir / "transfer_flow.png", config=flow_config, include_self_loops=args.include_self_loops):
            print(f"Wrote {output_dir / 'transfer_flow.png'}")
        return 0

    if args.command == "network":
        return _draw_network_command(args)

    if args.command == "focus":
        path = _input_path(args.file)
        bundle = _load_bundle_or_none(path)
        config = _resolve_config(path, args.template)
        focus_config = dict(config.get("center_flow", {}))
        if args.num_centers is not None:
            focus_config["centers"] = args.num_centers
        if args.max_in is not None:
            focus_config["max_in"] = args.max_in
        if args.max_out is not None:
            focus_config["max_out"] = args.max_out
        raw_flow = read_transfer_flow(bundle) if bundle else read_transfer_flow_file(path)
        output_dir = bundle.root if bundle else path.parent
        written = draw_focus_sankey_set(
            raw_flow,
            output_dir,
            centers=args.centers,
            config=focus_config,
            include_self_loops=args.include_self_loops,
        )
        for output in written:
            print(f"Wrote {output}")
        return 0

    if args.command == "events":
        path = _input_path(args.file)
        bundle = _load_bundle_or_none(path)
        config = _resolve_config(path, args.template)
        event_config = dict(config.get("reaction_events", {}))
        if args.mol:
            event_config["embed_molecules"] = True
        raw_events = read_reaction_events(bundle) if bundle else pd.read_csv(path, dtype=str)
        output_dir = bundle.root if bundle else path.parent
        if {"frequency", "reactants", "products"}.issubset(raw_events.columns):
            cleaned = raw_events.copy()
            cleaned["frequency"] = cleaned["frequency"].astype(int)
        else:
            cleaned = aggregate_reaction_events(
                raw_events,
                bundle.molecules if bundle else None,
                labels=bundle.labels if bundle else None,
                cancel_common=not args.no_cancel,
            )
        cleaned.to_csv(output_dir / "reaction_events_cleaned.csv", index=False)
        top_n = args.max_events or args.top or int(event_config.get("max_reactions", 15))
        event_outputs = plot_reaction_events(cleaned, output_dir / "reaction_events.png", top_n=top_n, config=event_config)
        print(f"Wrote {output_dir / 'reaction_events_cleaned.csv'}")
        for output in event_outputs:
            print(f"Wrote {output}")
        return 0

    if args.command == "molecules":
        path = _input_path(args.file)
        if path.is_dir():
            molecule_path = path / "molecules.csv"
            output_dir = path / "molecule_pictures"
        else:
            molecule_path = path
            output_dir = path.parent / "molecule_pictures"
        records = pd.read_csv(molecule_path, dtype=str).to_dict("records")
        written = draw_selected_molecules(records, output_dir, args.species)
        for item in written:
            print(f"Wrote {item}")
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
