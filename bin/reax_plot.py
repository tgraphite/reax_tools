#!/usr/bin/env python3

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from dataclasses import dataclass
import numpy as np
import pandas as pd
import os
import re
import math


@dataclass
class PlotType:
    pattern: str
    xlabel: str
    ylabel: str
    timestep: float
    timeunit: str
    suptitle: str
    to_filter: bool = True
    filter_trivial: bool = False
    filter_zeros: bool = False
    filter_stable: bool = False
    max_subplots_in_page: int = 6
    max_pages: int = 3
    sharex = False
    sharey = False


def formula_to_subscript(formula: str) -> str:
    """Convert chemical formula to subscript format (e.g., H2O -> H₂O)"""
    if not formula:
        return formula
    # Replace numbers with subscript characters
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    return formula.translate(subscript_map)


def filter_data_series(
    df: pd.DataFrame,
    filter_trivial: bool = False,
    filter_zeros: bool = True,
    filter_stable: bool = True,
) -> tuple[list[str], str]:

    notif_dict = {"Trivial": [], "Zeros": [], "Stable": []}

    df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)
    df_abs = df.abs()

    averages = df.mean()
    trivial_threshold = np.percentile(averages, 10)

    selected_columns = {c: averages[c] for c in df.columns}

    for column in df.columns:
        average = averages[column]

        stdev = np.std(df[column])

        if average != 0:
            rel_stdev = float(stdev) / float(average)
        else:
            rel_stdev = float(stdev)

        if filter_trivial and average < trivial_threshold:
            selected_columns.pop(column)
            notif_dict["Trivial"].append(column)
            continue
        if filter_zeros and np.abs(average) < (0.1 * len(df[column])) and stdev < 1.0:
            selected_columns.pop(column)
            notif_dict["Zeros"].append(column)
            continue
        if filter_stable and rel_stdev < 0.02:
            selected_columns.pop(column)
            notif_dict["Stable"].append(column)
            continue

    notif_string = ""
    for key, filtered_columns in notif_dict.items():
        size = len(notif_dict[key])
        repr_columns = ", ".join(notif_dict[key][: min(size, 5)])

        if size > 0:
            string = f"Hide {size} data series ({repr_columns} ...), reason: {key}\n"
            notif_string += string

    return list(selected_columns.keys()), notif_string


def filter_data_series_for_comparison(
    dfs: list[pd.DataFrame],
    filter_trivial: bool = False,
    filter_zeros: bool = True,
    filter_stable: bool = True,
    min_exist: float = 1.0,
) -> list[str]:

    compare_columns = dict()
    min_count = math.ceil(min_exist * len(dfs))

    for df in dfs:
        curr_selected_columns, _ = filter_data_series(
            df, filter_trivial, filter_zeros, filter_stable
        )
        for column in curr_selected_columns:
            if column not in compare_columns.keys():
                compare_columns[column] = 1
            else:
                compare_columns[column] += 1

    selected_columns = list()

    for column, count in compare_columns.items():
        if count >= min_count:
            selected_columns.append(column)

    return selected_columns


def plot_line_series(
    df: pd.DataFrame,
    selected_columns: list[str],
    output_basename: str,
    plot_type: PlotType,
) -> None:

    if len(selected_columns) == 0:
        print(f"No important data need to plot, skip.")
        return

    required_pages = min(
        math.ceil(len(selected_columns) / plot_type.max_subplots_in_page),
        plot_type.max_pages,
    )
    figsize = (10, 6)
    columns_in_pages = list()

    for page_id in range(required_pages):
        start_column_id = page_id * plot_type.max_subplots_in_page
        end_column_id = min(
            (page_id + 1) * plot_type.max_subplots_in_page, len(selected_columns)
        )
        columns = selected_columns[start_column_id:end_column_id]
        columns_in_pages.append(columns)

    flatten = list()
    for p in columns_in_pages:
        flatten.extend(p)
    size = len(flatten)
    repr_columns = ", ".join(flatten[: min(size, 5)])
    print(
        f"Trying to plot {size} data series ({repr_columns} ...)  out of {len(selected_columns)}"
    )

    # Color palette for multiple series
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    for page_id, columns in enumerate(columns_in_pages):
        fig, ax = plt.subplots(figsize=figsize)
        x = df.index * plot_type.timestep

        for idx, column in enumerate(columns):
            y = df[column]
            color = colors[idx % len(colors)]
            ax.plot(x, y, label=column, color=color, linewidth=1.5)

        ax.set_xlabel(f"{plot_type.xlabel} ({plot_type.timeunit})")
        ax.set_ylabel(plot_type.ylabel)
        ax.xaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_major_locator(MaxNLocator(6))

        # Add legend with subscript formatting
        handles, labels = ax.get_legend_handles_labels()
        subscript_labels = [formula_to_subscript(label) for label in labels]
        ax.legend(handles, subscript_labels, loc='best', framealpha=0.9)

        # Add grid for better readability
        ax.grid(True, alpha=0.3, linestyle='--')

        series_name = os.path.dirname(output_basename)
        series_name = re.sub("_.*$", "", series_name)
        if series_name == ".":
            series_name = "default"
        if plot_type.suptitle and len(columns_in_pages) == 1:
            fig.suptitle(f"{plot_type.suptitle}, workdir: {series_name}", fontsize=12)
        elif plot_type.suptitle and len(columns_in_pages) > 1:
            curr_suptitle = f"{plot_type.suptitle}, workdir: {series_name} (page {page_id+1}/{len(columns_in_pages)})"
            fig.suptitle(curr_suptitle, fontsize=12)

        if len(columns_in_pages) == 1:
            output_file = output_basename + ".png"
        else:
            output_file = output_basename + f".page{page_id+1}.png"

        plt.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Plot and saved file {output_file}")

        plt.close()


def plot_reaction_track(input_file: str, output_file: str = None, top_n: int = 20):
    """
    Plot reaction track events as horizontal bar chart.
    Reaction equation on left, frequency bar on right.
    """
    if not os.path.isfile(input_file):
        print(f"{input_file} not found, skip reaction track plot.")
        return

    df = pd.read_csv(input_file)
    
    # Get top N reactions by frequency
    df_top = df.head(top_n).copy()
    
    if len(df_top) == 0:
        print(f"No reaction data to plot.")
        return

    # Create reaction equation string
    def format_reaction(row):
        reactants = row['reactants']
        products = row['products']
        # Apply subscript formatting
        reactants = formula_to_subscript(str(reactants))
        products = formula_to_subscript(str(products))
        return f"{reactants} → {products}"

    df_top['reaction'] = df_top.apply(format_reaction, axis=1)
    
    # Reverse order for top-to-bottom display
    df_top = df_top.iloc[::-1]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, top_n * 0.4)))
    
    # Colors based on frequency
    frequencies = df_top['frequency'].values
    colors = plt.cm.viridis(frequencies / frequencies.max())
    
    # Plot horizontal bars
    y_pos = range(len(df_top))
    bars = ax.barh(y_pos, frequencies, color=colors, edgecolor='black', linewidth=0.5)
    
    # Set y-axis labels with reaction equations
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_top['reaction'].values, fontsize=9)
    
    # Set labels and title
    ax.set_xlabel('Frequency (count)', fontsize=11)
    ax.set_ylabel('Reaction', fontsize=11)
    ax.set_title(f'Top {top_n} Most Frequent Reactions', fontsize=13, fontweight='bold')
    
    # Add frequency labels at the end of bars
    for i, (bar, freq) in enumerate(zip(bars, frequencies)):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f'{int(freq)}', va='center', ha='left', fontsize=9)
    
    # Adjust layout to fit long labels
    plt.tight_layout()
    
    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + ".png"
    
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Reaction track plot saved to {output_file}")
    plt.close()


def get_plot_type(input_file: str, plot_types: list[PlotType]) -> PlotType:
    for pt in plot_types:
        if pt.pattern in input_file:
            return pt

    all_pt_patterns = [pt.pattern for pt in plot_types]
    all_pt_patterns = ", ".join(all_pt_patterns)
    print(
        f"Error: can not determine plot type for file: {input_file}\n\
        file name should have one of the patterns in below:\
\
        {all_pt_patterns}"
    )


def plot_single_file(input_file: str, plot_types: list[PlotType]):
    print("-" * 40)
    print(f"ReaxTools post-process plot kit, working on file: {input_file}\n")

    if not os.path.isfile(input_file):
        print(f"{input_file} not found, skip.")
        return

    # Handle reaction track events file separately
    if "reaction_track_events" in input_file:
        plot_reaction_track(input_file)
        return

    output_basename = os.path.splitext(input_file)[0]
    df = pd.read_csv(input_file)
    plot_type = get_plot_type(input_file, plot_types)

    if plot_type.to_filter:
        selected_columns, notif_string = filter_data_series(df, False, True, True)
        print(notif_string)
    else:
        selected_columns = list(df.columns)

    plot_line_series(df, selected_columns, output_basename, plot_type)


def plot_comparison(input_files: list[str], plot_types: list[PlotType]):
    print("-" * 40)
    files_string = ", ".join(input_files)
    print(
        f"ReaxTools post-process plot kit, comparison mode, working on files: {files_string}\n"
    )

    dfs = list()
    series_names = list()
    basedir = os.path.dirname(input_files[0])

    plot_type = get_plot_type(input_files[0], plot_types)
    raw_suptitle = plot_type.suptitle
    plot_type.max_subplots_in_page = 6  # Single plot with up to 6 series
    plot_type.to_filter = False
    plot_type.sharey = True

    for input_file in input_files:
        df = pd.read_csv(input_file)
        dfs.append(df)
        # Use parent directory name + file basename to create unique series name
        abs_path = os.path.abspath(input_file)
        dir_name = os.path.basename(os.path.dirname(abs_path))
        base_name = os.path.basename(input_file)
        base_name = re.sub("_.*$", "", base_name)
        series_name = f"{dir_name}_{base_name}"
        series_names.append(series_name)

    selected_columns = filter_data_series_for_comparison(dfs, False, False, False)
    sel_columns_string = ", ".join(selected_columns)
    print(f"Found {len(selected_columns)} properties to compare: {sel_columns_string}")

    for column in selected_columns:
        curr_compare_data = dict()
        max_len = 0

        for series_name, df in zip(series_names, dfs):
            if column in df.columns:
                curr_compare_data[series_name] = df[column].to_numpy()
                max_len = max(max_len, len(df[column]))

        # Handle different length arrays by padding with NaN
        for key in curr_compare_data:
            arr = curr_compare_data[key]
            if len(arr) < max_len:
                padded = np.full(max_len, np.nan)
                padded[:len(arr)] = arr
                curr_compare_data[key] = padded

        curr_compare_df = pd.DataFrame(curr_compare_data)

        plot_type.suptitle = f"{raw_suptitle}, comparison for {formula_to_subscript(column)}."

        output_basename = os.path.join(basedir, f"compare_{plot_type.pattern}_{column}")
        plot_line_series(
            curr_compare_df, curr_compare_df.columns, output_basename, plot_type
        )

        compare_file_path = f"{output_basename}.csv"
        curr_compare_df.to_csv(compare_file_path, float_format="%.4f", index=False)


def main():
    import argparse
    import glob

    usage = "ReaxTools post-process plot kit for bond count, ring count, and species count.\npython reax_plot.py -d <workdir> -m <max_species>"

    parser = argparse.ArgumentParser(usage=usage)

    parser.add_argument("--dir", "-d", type=str, help="Work directory")
    parser.add_argument("--file", "-f", type=str, help="Plot a file")
    parser.add_argument(
        "--compare",
        "-c",
        type=str,
        nargs="+",
        help="Compare a series of files, should be same data type",
    )
    parser.add_argument(
        "--timestep", "-t", type=float, help="Timestep between frames", default=1
    )
    parser.add_argument(
        "--timeunit", "-u", type=str, help="Time unit of timestep", default="ps"
    )
    parser.add_argument(
        "--reaction-track", "-r", type=str, help="Plot reaction track events CSV file"
    )
    parser.add_argument(
        "--top-n", "-n", type=int, help="Number of top reactions to plot (default=20)", default=20
    )

    dir = parser.parse_args().dir
    single_file = parser.parse_args().file
    compare_files = parser.parse_args().compare
    timestep = parser.parse_args().timestep
    timeunit = parser.parse_args().timeunit
    reaction_track_file = parser.parse_args().reaction_track
    top_n = parser.parse_args().top_n

    # Handle reaction track plot directly
    if reaction_track_file:
        plot_reaction_track(reaction_track_file, top_n=top_n)
        return

    plot_type_species_count = PlotType(
        pattern="species_count",
        xlabel="Time",
        ylabel="Count / Weight",
        suptitle="species count plot (y axis = weight in group-merge mode)",
        timestep=timestep,
        timeunit=timeunit,
        to_filter=True,
        filter_trivial=True,
        filter_zeros=True,
        filter_stable=False,
    )

    plot_type_bond_count = PlotType(
        pattern="bond_count",
        xlabel="Time",
        ylabel="Count",
        suptitle="ReaxTools autoplot for bond count",
        timestep=timestep,
        timeunit=timeunit,
        to_filter=False,
        filter_trivial=False,
        filter_zeros=True,
        filter_stable=False,
    )

    plot_type_atom_bonded_num_count = PlotType(
        pattern="atom_bonded_num_count",
        xlabel="Time",
        ylabel="Count",
        suptitle="ReaxTools autoplot for atom type count",
        timestep=timestep,
        timeunit=timeunit,
        to_filter=False,
        filter_trivial=True,
        filter_zeros=True,
        filter_stable=False,
    )

    plot_type_ring_count = PlotType(
        pattern="ring_count",
        xlabel="Time",
        ylabel="Count",
        suptitle="ReaxTools autoplot for ring size count",
        timestep=timestep,
        timeunit=timeunit,
        to_filter=False,
        filter_trivial=True,
        filter_zeros=True,
        filter_stable=False,
    )

    plot_types = [
        plot_type_species_count,
        plot_type_bond_count,
        plot_type_atom_bonded_num_count,
        plot_type_ring_count,
    ]

    patterns = [pt.pattern for pt in plot_types]
    patterns.append("reaction_track_events")  # Add reaction track pattern

    if dir:
        files = list()
        for pattern in patterns:
            matched_files = glob.glob(os.path.join(dir, f"*{pattern}.csv"))
            files.extend(matched_files)

        for input_file in files:
            plot_single_file(input_file, plot_types)

    elif single_file:
        plot_single_file(single_file, plot_types)

    elif compare_files:
        plot_comparison(compare_files, plot_types)


if __name__ == "__main__":
    main()
