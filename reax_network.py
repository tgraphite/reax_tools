#!/usr/bin/env python3

import networkx as nx
import pygraphviz
from dataclasses import dataclass
import re
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import time
import signal
from cairosvg import svg2png
import argparse
import glob


@dataclass
class Node:
    name: str
    smiles = ""
    img_path: str = ""
    img_generated: bool = False  # Track if image was successfully generated


class ReaxNetwork:
    def __init__(self, file, smiles_file=None):
        """
        Initialize ReaxNetwork by parsing a dot file

        Args:
            file (str): Path to the dot file
            smiles_file (str): Path to the SMILES CSV file
        """
        self.nx = nx.DiGraph()  # Create directed graph
        self.nodes = {}  # Store node objects
        self.smiles_dict = {}  # Store SMILES mapping
        self.dot_file_path = file  # Store the dot file path

        # Load SMILES data if provided
        if smiles_file:
            self._load_smiles(smiles_file)

        self._parse_dot_file(file)
        self._generate_molecular_images()

    def _load_smiles(self, smiles_file):
        """
        Load SMILES data from CSV file

        Args:
            smiles_file (str): Path to the SMILES CSV file
        """
        try:
            df = pd.read_csv(smiles_file, header=None, names=["formula", "smiles"])
            self.smiles_dict = dict(zip(df["formula"], df["smiles"]))
        except Exception as e:
            print(f"Warning: Could not load SMILES file: {e}")
            self.smiles_dict = {}

    def _convert_svg_to_png(self, svg_path, png_path):
        """
        Convert SVG file to PNG for use with PyGraphviz

        Args:
            svg_path (str): Path to SVG file
            png_path (str): Path to output PNG file
        """
        try:
            svg2png(
                url=svg_path, write_to=png_path, output_width=100, output_height=100
            )
            return True
        except Exception as e:
            print(f"Warning: Could not convert SVG to PNG for {svg_path}: {e}")
            return False

    def _generate_molecular_images(self):
        """
        Generate molecular images for nodes that have SMILES data
        First convert all SVG files to PNG, then process nodes
        """
        # Get the parent directory of the dot file
        dot_file_dir = (
            os.path.dirname(self.dot_file_path)
            if hasattr(self, "dot_file_path")
            else "."
        )
        molecule_pictures_dir = os.path.join(dot_file_dir, "molecule_pictures")

        # Create molecule_pictures directory if it doesn't exist
        os.makedirs(molecule_pictures_dir, exist_ok=True)

        # First, convert all SVG files to PNG files
        self._convert_all_svg_to_png(molecule_pictures_dir)

        def timeout_handler(signum, frame):
            raise TimeoutError("Image generation timeout")

        for node_id, node in self.nodes.items():
            formula = node.name
            if formula in self.smiles_dict:
                smiles = self.smiles_dict[formula]
                node.smiles = smiles

                # Look for PNG file (either converted from SVG or generated)
                png_path = os.path.join(
                    molecule_pictures_dir, f"molecule_{formula}.png"
                )

                # Check if PNG exists
                if os.path.exists(png_path):
                    node.img_path = png_path
                    node.img_generated = True
                    print(f"Using existing PNG for {formula}: {png_path}")
                    continue

                # Generate new PNG if it doesn't exist
                try:
                    # Set timeout for image generation (500ms)
                    signal.signal(signal.SIGALRM, timeout_handler)
                    signal.alarm(1)  # 1 second timeout (1000ms)

                    start_time = time.time()

                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        # Generate PNG image
                        img = Draw.MolToImage(mol, size=(100, 100))
                        img.save(png_path)

                        # Check if generation took more than 500ms
                        generation_time = (time.time() - start_time) * 1000
                        if generation_time > 500:
                            print(
                                f"Warning: Image generation for {formula} took {generation_time:.1f}ms (>500ms), using formula instead"
                            )
                            node.img_generated = False
                            # Remove the generated file
                            if os.path.exists(png_path):
                                os.remove(png_path)
                        else:
                            node.img_path = png_path
                            node.img_generated = True
                            print(
                                f"Successfully generated PNG for {formula} in {generation_time:.1f}ms: {png_path}"
                            )
                    else:
                        print(
                            f"Warning: Could not parse SMILES for {formula}: {smiles}"
                        )
                        node.img_generated = False

                    # Cancel the alarm
                    signal.alarm(0)

                except TimeoutError:
                    print(
                        f"Warning: Image generation timeout for {formula}, using formula instead"
                    )
                    node.img_generated = False
                    # Remove partial file if exists
                    if os.path.exists(png_path):
                        os.remove(png_path)
                    signal.alarm(0)

                except Exception as e:
                    print(f"Warning: Could not generate image for {formula}: {e}")
                    node.img_generated = False
                    # Remove partial file if exists
                    if os.path.exists(png_path):
                        os.remove(png_path)
                    signal.alarm(0)

    def _convert_all_svg_to_png(self, molecule_pictures_dir):
        """
        Convert all SVG files in the molecule_pictures directory to PNG files

        Args:
            molecule_pictures_dir (str): Path to the molecule_pictures directory
        """
        if not os.path.exists(molecule_pictures_dir):
            return

        svg_files = [f for f in os.listdir(molecule_pictures_dir) if f.endswith(".svg")]

        if svg_files:
            print(f"Converting {len(svg_files)} SVG files to PNG...")

            for svg_file in svg_files:
                svg_path = os.path.join(molecule_pictures_dir, svg_file)
                png_file = svg_file.replace(".svg", ".png")
                png_path = os.path.join(molecule_pictures_dir, png_file)

                # Only convert if PNG doesn't already exist
                if not os.path.exists(png_path):
                    try:
                        if self._convert_svg_to_png(svg_path, png_path):
                            print(f"  Converted: {svg_file} -> {png_file}")
                        else:
                            print(f"  Failed to convert: {svg_file}")
                    except Exception as e:
                        print(f"  Error converting {svg_file}: {e}")
                else:
                    print(f"  PNG already exists: {png_file}")
        else:
            print("No SVG files found in molecule_pictures directory")

    def _parse_dot_file(self, file_path):
        """
        Parse the dot file and populate the network

        Args:
            file_path (str): Path to the dot file
        """
        with open(file_path, "r") as f:
            content = f.read()

        # Parse nodes
        node_pattern = r'node(\d+)\s*\[label="([^"]+)"\];'
        node_matches = re.findall(node_pattern, content)

        for node_id, label in node_matches:
            # Create Node object and add to graph
            node = Node(name=label)
            self.nodes[node_id] = node
            self.nx.add_node(node_id, label=label, node_obj=node)

        # Parse edges
        edge_pattern = r'node(\d+)\s*->\s*node(\d+)\s*\[label="([^"]+)"(?:,\s*penwidth=([^,]+))?(?:,\s*color=([^\]]+))?\];'
        edge_matches = re.findall(edge_pattern, content)

        for source, target, label, penwidth, color in edge_matches:
            # Parse reaction count from label (format: "count : count")
            reaction_count = int(label.split(":")[0].strip())

            # Create edge attributes
            edge_attrs = {"label": label, "reaction_count": reaction_count}

            if penwidth:
                edge_attrs["penwidth"] = float(penwidth)
            if color:
                edge_attrs["color"] = color

            # Add edge to graph
            self.nx.add_edge(source, target, **edge_attrs)

    def draw_network(
        self,
        output_file="reaction_network.png",
        layout="dot",
        rankdir="LR",
        dpi=300,
        node_size=1.0,
        font_size=12,
    ):
        """
        Draw the reaction network using PyGraphviz for professional visualization

        Args:
            output_file (str): Output file path for the network image
            layout (str): Layout engine ('dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo')
            rankdir (str): Direction of layout ('TB', 'LR', 'BT', 'RL')
            dpi (int): DPI for the output image
            node_size (float): Size multiplier for nodes
            font_size (int): Font size for labels
        """
        # Create PyGraphviz graph
        G = pygraphviz.AGraph(directed=True, strict=False)
        G.graph_attr.update(
            rankdir=rankdir,
            layout=layout,
            dpi=str(dpi),
            fontsize=str(font_size),
            fontname="Arial",
            splines="ortho",  # Use orthogonal edges for cleaner look
            nodesep="0.5",
            ranksep="1.0",
        )

        # Add nodes with appropriate styling
        for node_id, node in self.nodes.items():
            if node.img_generated and node.img_path and os.path.exists(node.img_path):
                # Node with molecular image
                G.add_node(
                    node_id,
                    label=node.name,
                    shape="box",
                    style="filled",
                    fillcolor="lightblue",
                    fontsize=str(font_size),
                    width=str(node_size * 2.0),
                    height=str(node_size * 2.0),
                    fixedsize="true",
                    image=node.img_path,
                    imagescale="true",
                    labelloc="b",  # Label below image
                )
            else:
                # Node with chemical formula only
                G.add_node(
                    node_id,
                    label=node.name,
                    shape="box",
                    style="filled",
                    fillcolor="lightgray",
                    fontsize=str(font_size),
                    width=str(node_size * 1.5),
                    height=str(node_size * 0.8),
                    fixedsize="true",
                )

        # Add edges with reaction count information
        for u, v, data in self.nx.edges(data=True):
            original_label = data.get("label", "")  # Use original label from dot file
            original_penwidth = data.get(
                "penwidth", "1.0"
            )  # Use original penwidth from dot file
            original_color = data.get(
                "color", "gray"
            )  # Use original color from dot file

            G.add_edge(
                u,
                v,
                label=original_label,  # Use original label from dot file
                color=original_color,  # Use original color from dot file
                penwidth=original_penwidth,  # Use original penwidth from dot file
                fontsize="8",
            )

        # Layout and render the graph
        try:
            G.layout(prog=layout)
            G.draw(output_file, format="png")
            print(f"Network visualization saved to {output_file}")

        except Exception as e:
            print(f"Error generating network visualization: {e}")
            # Fallback to basic layout
            try:
                G.layout(prog="dot")
                G.draw(output_file, format="png")
                print(f"Network visualization (fallback) saved to {output_file}")
            except Exception as e2:
                print(f"Failed to generate network visualization: {e2}")

    def draw_network_simple(
        self,
        output_file="reaction_network_simple.png",
        layout="dot",
        rankdir="LR",
    ):
        """
        Draw a simplified version of the network without molecular images

        Args:
            output_file (str): Output file path
            layout (str): Layout engine
            rankdir (str): Direction of layout
        """
        G = pygraphviz.AGraph(directed=True, strict=False)
        G.graph_attr.update(
            rankdir=rankdir,
            layout=layout,
            fontsize="12",
            fontname="Arial",
            splines="ortho",
            nodesep="0.8",
            ranksep="1.2",
        )

        # Add all nodes as simple boxes
        for node_id, node in self.nodes.items():
            G.add_node(
                node_id,
                label=node.name,
                shape="box",
                style="filled",
                fillcolor="lightblue" if node.smiles else "lightgray",
                fontsize="10",
                width="1.2",
                height="0.6",
                fixedsize="true",
            )

        # Add edges
        for u, v, data in self.nx.edges(data=True):
            original_label = data.get("label", "")  # Use original label from dot file
            original_penwidth = data.get(
                "penwidth", "1.0"
            )  # Use original penwidth from dot file
            original_color = data.get(
                "color", "gray"
            )  # Use original color from dot file

            G.add_edge(
                u,
                v,
                label=original_label,  # Use original label from dot file
                color=original_color,  # Use original color from dot file
                penwidth=original_penwidth,  # Use original penwidth from dot file
                fontsize="8",
            )

        try:
            G.layout(prog=layout)
            G.draw(output_file, format="png")
            print(f"Simple network visualization saved to {output_file}")
        except Exception as e:
            print(f"Error generating simple network: {e}")

    def get_network_stats(self):
        """
        Get basic statistics about the network

        Returns:
            dict: Network statistics
        """
        stats = {
            "total_nodes": self.nx.number_of_nodes(),
            "total_edges": self.nx.number_of_edges(),
            "nodes_with_smiles": sum(1 for node in self.nodes.values() if node.smiles),
            "nodes_without_smiles": sum(
                1 for node in self.nodes.values() if not node.smiles
            ),
            "nodes_with_images": sum(
                1 for node in self.nodes.values() if node.img_generated
            ),
            "nodes_with_formulas_only": sum(
                1 for node in self.nodes.values() if not node.img_generated
            ),
            "total_reactions": sum(
                data.get("reaction_count", 0) for _, _, data in self.nx.edges(data=True)
            ),
        }
        return stats


def process_directory(
    directory_path,
    include_full=False,
    layout="dot",
    rankdir="LR",
    dpi=300,
    node_size=1.0,
    font_size=12,
):
    """
    Process all dot files in a directory and generate network visualizations

    Args:
        directory_path (str): Path to directory containing dot files and molecule_pictures
        include_full (bool): Whether to include reactions_full.dot file
        layout (str): Layout engine for network visualization
        rankdir (str): Direction of layout
        dpi (int): DPI for output images
        node_size (float): Size multiplier for nodes
        font_size (int): Font size for labels
    """
    print(f"Processing directory: {directory_path}")

    # Check if directory exists
    if not os.path.exists(directory_path):
        print(f"Error: Directory {directory_path} does not exist")
        return

    # Find all dot files in the directory
    dot_files = glob.glob(os.path.join(directory_path, "*.dot"))

    if not dot_files:
        print(f"No dot files found in {directory_path}")
        return

    print(f"Found {len(dot_files)} dot files:")
    for dot_file in dot_files:
        print(f"  - {os.path.basename(dot_file)}")

    # Filter out reactions_full.dot unless include_full is True
    if not include_full:
        dot_files = [
            f for f in dot_files if os.path.basename(f) != "reactions_full.dot"
        ]
        print(f"Excluding reactions_full.dot (use --full to include)")

    print(f"Processing {len(dot_files)} dot files...")

    # Find SMILES file
    smiles_file = os.path.join(directory_path, "molecules_smiles.csv")
    if not os.path.exists(smiles_file):
        print(f"Warning: SMILES file not found at {smiles_file}")
        smiles_file = None

    # Process each dot file
    for dot_file in dot_files:
        print(f"\nProcessing: {os.path.basename(dot_file)}")

        try:
            # Initialize network
            network = ReaxNetwork(file=dot_file, smiles_file=smiles_file)

            # Get network statistics
            stats = network.get_network_stats()
            print(f"  Nodes: {stats['total_nodes']}, Edges: {stats['total_edges']}")
            print(f"  Nodes with images: {stats['nodes_with_images']}")
            print(f"  Total reactions: {stats['total_reactions']}")

            # Generate output filename
            base_name = os.path.splitext(os.path.basename(dot_file))[0]
            output_file = os.path.join(directory_path, f"{base_name}_network.png")

            # Generate network visualization
            print(f"  Generating visualization: {os.path.basename(output_file)}")
            network.draw_network(
                output_file=output_file,
                layout=layout,
                rankdir=rankdir,
                dpi=dpi,
                node_size=node_size,
                font_size=font_size,
            )

            # Also generate simple version
            simple_output = os.path.join(
                directory_path, f"{base_name}_network_simple.png"
            )
            print(
                f"  Generating simple visualization: {os.path.basename(simple_output)}"
            )
            network.draw_network_simple(
                output_file=simple_output, layout=layout, rankdir=rankdir
            )

        except Exception as e:
            print(f"  Error processing {os.path.basename(dot_file)}: {e}")
            continue

    print(
        f"\nProcessing completed! Check the output directory for generated network images."
    )


def main():
    """
    Command line interface for ReaxNetwork
    """
    parser = argparse.ArgumentParser(
        description="Generate reaction network visualizations from dot files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reax_network.py /path/to/directory
  python reax_network.py /path/to/directory --full
  python reax_network.py /path/to/directory --layout neato --rankdir TB
  python reax_network.py /path/to/directory --dpi 600 --node-size 1.5
        """,
    )

    parser.add_argument(
        "directory", help="Directory containing dot files and molecule_pictures folder"
    )

    parser.add_argument(
        "--full",
        action="store_true",
        help="Include reactions_full.dot file in processing",
    )

    parser.add_argument(
        "--layout",
        default="dot",
        choices=["dot", "neato", "fdp", "sfdp", "twopi", "circo"],
        help="Layout engine for network visualization (default: dot)",
    )

    parser.add_argument(
        "--rankdir",
        default="LR",
        choices=["TB", "LR", "BT", "RL"],
        help="Direction of layout (default: LR)",
    )

    parser.add_argument(
        "--dpi", type=int, default=300, help="DPI for output images (default: 300)"
    )

    parser.add_argument(
        "--node-size",
        type=float,
        default=1.0,
        help="Size multiplier for nodes (default: 1.0)",
    )

    parser.add_argument(
        "--font-size", type=int, default=12, help="Font size for labels (default: 12)"
    )

    args = parser.parse_args()

    # Process the directory
    process_directory(
        directory_path=args.directory,
        include_full=args.full,
        layout=args.layout,
        rankdir=args.rankdir,
        dpi=args.dpi,
        node_size=args.node_size,
        font_size=args.font_size,
    )


if __name__ == "__main__":
    main()
