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
    fig_nrow = math.ceil(plot_type.max_subplots_in_page / 2)
    fig_ncol = 2
    figsize = (fig_ncol * 5, fig_nrow * 2)
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

    for page_id, columns in enumerate(columns_in_pages):
        fig, axs = plt.subplots(
            fig_nrow,
            fig_ncol,
            figsize=figsize,
            sharex=plot_type.sharex,
            sharey=plot_type.sharey,
        )
        x = df.index * plot_type.timestep

        for column_id, column in enumerate(columns):
            y = df[column]
            ax = axs.flat[column_id]
            ax.plot(x, y)
            ax.set_xlabel(f"{plot_type.xlabel} ({plot_type.timeunit})")
            ax.set_ylabel(plot_type.ylabel)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))

            ax.text(0, 1.02, column, ha="left", va="bottom", transform=ax.transAxes)

            p_5 = np.percentile(y, 5)
            p_95 = np.percentile(y, 95)
            avg = np.mean(y)
            brief_report = (
                f"Higher 5%: {p_95:.1f} | Avg: {avg:.1f} | Lower 5%: {p_5:.1f} "
            )
            ax.text(
                1, 1.02, brief_report, ha="right", va="bottom", transform=ax.transAxes
            )

        if len(columns) < plot_type.max_subplots_in_page:
            for i in range(len(columns), plot_type.max_subplots_in_page):
                if i > 1:
                    fig.delaxes(axs.flat[i])
                    # The canvas will be cutted to half the width if only one subplot.
                    # In this case, temporary keep a placeholder in the last page.

        if len(columns) == 1:
            axs.flat[1].text(
                0.5,
                0.5,
                "Placeholder, no data to plot.",
                transform=axs.flat[1].transAxes,
            )

        series_name = os.path.basename(output_basename)
        series_name = re.sub("_.*$", "", series_name)

        if plot_type.suptitle and len(columns_in_pages) == 1:
            fig.suptitle(f"{series_name} {plot_type.suptitle}")
        elif plot_type.suptitle and len(columns_in_pages) > 1:
            curr_suptitle = f"{series_name} {plot_type.suptitle} (page {page_id+1}/{len(columns_in_pages)})"
            fig.suptitle(curr_suptitle)

        if len(columns_in_pages) == 1:
            output_file = output_basename + ".png"
        else:
            output_file = output_basename + f".page{page_id+1}.png"

        plt.tight_layout()
        fig.savefig(output_file, dpi=600, bbox_inches="tight")
        print(f"Plot and saved file {output_file}")

        plt.close()


def get_plot_type(input_file: str, plot_types: list[PlotType]) -> PlotType:
    for pt in plot_types:
        if pt.pattern in input_file:
            return pt

    all_pt_patterns = [pt.pattern for pt in plot_types]
    all_pt_patterns = ", ".join(all_pt_patterns)
    print(
        f"Error: can not determine plot type for file: {input_file}\n\
        file name should have one of the patterns in below:\n\
        {all_pt_patterns}"
    )


def plot_single_file(input_file: str, plot_types: list[PlotType]):
    print("-" * 40)
    print(f"ReaxTools post-process plot kit, working on file: {input_file}\n")

    if not os.path.isfile(input_file):
        print(f"{input_file} not found, skip.")
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
    plot_type.max_subplots_in_page = 8
    plot_type.to_filter = False
    plot_type.sharey = True

    for input_file in input_files:
        df = pd.read_csv(input_file)
        dfs.append(df)
        series_name = os.path.basename(input_file)
        series_name = re.sub("_.*$", "", series_name)
        series_names.append(series_name)

    selected_columns = filter_data_series_for_comparison(dfs, False, False, False)
    sel_columns_string = ", ".join(selected_columns)
    print(f"Found {len(selected_columns)} properties to compare: {sel_columns_string}")

    for column in selected_columns:
        curr_compare_data = dict()

        for series_name, df in zip(series_names, dfs):
            if column in df.columns:
                curr_compare_data[series_name] = df[column].to_numpy()

        curr_compare_df = pd.DataFrame(curr_compare_data, index=None)

        plot_type.suptitle = f"{raw_suptitle}, comparison for {column}."

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

    dir = parser.parse_args().dir
    single_file = parser.parse_args().file
    compare_files = parser.parse_args().compare
    timestep = parser.parse_args().timestep
    timeunit = parser.parse_args().timeunit

    plot_type_species_count = PlotType(
        pattern="species_count",
        xlabel="Time",
        ylabel="Count / Weight",
        suptitle="species count plot",
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
        suptitle="bond count plot",
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
        ylabel="Count / Weight",
        suptitle="atom type count plot",
        timestep=timestep,
        timeunit=timeunit,
        to_filter=True,
        filter_trivial=True,
        filter_zeros=True,
        filter_stable=False,
    )

    plot_type_ring_count = PlotType(
        pattern="ring_count",
        xlabel="Time",
        ylabel="Count / Weight",
        suptitle="ring size count plot",
        timestep=timestep,
        timeunit=timeunit,
        to_filter=True,
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
