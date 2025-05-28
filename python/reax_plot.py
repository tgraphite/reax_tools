#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os


def plot_bond_count(dir):
    file = os.path.join(dir, "bond_count.csv")

    if not os.path.exists(file):
        print(f"Bond count file bond_count.csv not found in directory: {dir}.")
        return

    output = os.path.join(dir, "bond_count.png")
    df = pd.read_csv(file)

    fig_columns = 3
    fig_rows = (df.columns.size + fig_columns - 1) // fig_columns

    fig, axs = plt.subplots(
        fig_rows, fig_columns, figsize=(4 * fig_columns, 3 * fig_rows)
    )

    for i, col in enumerate(df.columns):
        axs.flat[i].plot(df[col])
        axs.flat[i].set_title(col)
        axs.flat[i].set_xlabel("Frame")
        axs.flat[i].set_ylabel("Bond count")

    for i in range(len(df.columns), fig_rows * fig_columns):
        fig.delaxes(axs.flat[i])

    fig.suptitle(
        f"Reax_tools auto-plot for bond count, workdir: {os.path.basename(dir)}"
    )

    plt.tight_layout()
    fig.savefig(output, dpi=600)
    print(f"Bond count plot saved to {output}.")


def plot_ring_count(dir):
    file = os.path.join(dir, "ring_count.csv")

    if not os.path.exists(file):
        print(f"Ring count file ring_count.csv not found in directory: {dir}.")
        return

    output = os.path.join(dir, "ring_count.png")
    df = pd.read_csv(file)

    fig_columns = 3
    fig_rows = (df.columns.size + fig_columns - 1) // fig_columns

    fig, axs = plt.subplots(
        fig_rows, fig_columns, figsize=(4 * fig_columns, 3 * fig_rows)
    )

    for i, col in enumerate(df.columns):
        axs.flat[i].plot(df[col])
        axs.flat[i].set_title(col)
        axs.flat[i].set_xlabel("Frame")
        axs.flat[i].set_ylabel("Independent ring count")

    # Remove empty subplots
    for i in range(len(df.columns), fig_rows * fig_columns):
        fig.delaxes(axs.flat[i])

    fig.suptitle(f"Reax_tools auto-plot for ring count")

    plt.tight_layout()
    fig.savefig(output, dpi=600)
    print(f"Ring count plot saved to {output}.")


def plot_atom_bonded_num_count(dir):
    file = os.path.join(dir, "atom_bonded_num_count.csv")

    if not os.path.exists(file):
        print(
            f"Atom bonded num count file atom_bonded_num_count not found in directory: {dir}."
        )
        return

    output = os.path.join(dir, "atom_bonded_num_count.png")
    df = pd.read_csv(file)

    fig_columns = 3
    fig_rows = (df.columns.size + fig_columns - 1) // fig_columns

    fig, axs = plt.subplots(
        fig_rows, fig_columns, figsize=(4 * fig_columns, 3 * fig_rows)
    )

    for i, col in enumerate(df.columns):
        axs.flat[i].plot(df[col])
        axs.flat[i].set_title(col)
        axs.flat[i].set_xlabel("Frame")
        axs.flat[i].set_ylabel("Count")

    # Remove empty subplots
    for i in range(len(df.columns), fig_rows * fig_columns):
        fig.delaxes(axs.flat[i])

    fig.suptitle(
        f"Reax_tools auto-plot for atom bonded num count\nE.g.: C(3) means carbon with 3 bonded atoms"
    )

    plt.tight_layout()
    fig.savefig(output, dpi=600)
    print(f"Atom bonded num count plot saved to {output}.")


def plot_species_count(dir, max_species=9):
    file = os.path.join(dir, "species_count.csv")

    if not os.path.exists(file):
        print(f"Species count file species_count.csv not found in directory: {dir}.")
        return

    output = os.path.join(dir, "species_count.png")
    df = pd.read_csv(file)

    fig_columns = 3
    fig_rows = max_species // fig_columns

    fig, axs = plt.subplots(
        fig_rows, fig_columns, figsize=(4 * fig_columns, 3 * fig_rows)
    )

    # sort df by mean of each column
    df = df.reindex(df.mean().sort_values(ascending=False).index, axis=1)
    df = df.iloc[:, :max_species]

    for i, col in enumerate(df.columns):
        axs.flat[i].plot(df[col])
        axs.flat[i].set_title(col)
        axs.flat[i].set_xlabel("Frame")
        axs.flat[i].set_ylabel("Species count or weight")

    for i in range(len(df.columns), fig_rows * fig_columns):
        fig.delaxes(axs.flat[i])

    fig.suptitle(
        f"Reax_tools auto-plot for species count\nNote: When using -rc option, the species count is rescaled in weight."
    )

    plt.tight_layout()
    fig.savefig(output, dpi=600)
    print(f"Species count plot saved to {output}.")


if __name__ == "__main__":
    import argparse

    usage = "Reax_tools auto-plot for bond count, ring count, and species count.\npython reax_plot.py -d <workdir> -m <max_species>"

    parser = argparse.ArgumentParser(usage=usage)

    parser.add_argument("--dir", "-d", type=str, required=True, help="Work directory")
    parser.add_argument(
        "--max_species", "-m", type=int, default=9, help="Max species to plot"
    )

    dir = parser.parse_args().dir
    max_species = parser.parse_args().max_species

    plot_bond_count(dir)
    plot_ring_count(dir)
    plot_species_count(dir, max_species)
    plot_atom_bonded_num_count(dir)
