#!/usr/bin/env python3

import networkx as nx
import sys
import os
from matplotlib import pyplot as plt
from dataclasses import dataclass
# import pygraphviz

@dataclass 
class DotFileNode:
    hash: str
    formula: str
    degree: int
    in_degree: int
    out_degree: int

@dataclass
class DotFileEdge:
    source: DotFileNode
    target: DotFileNode
    degree: int

class ReaxNetworkSolver:
    def __init__(self, graph_file):
        self.graph_file = graph_file
        self.graph = None

    def load_gml(self):
        if not os.path.exists(self.graph_file):
            raise FileNotFoundError(f"Graph file {self.graph_file} does not exist.")
        self.graph = nx.read_gml(self.graph_file)

    def solve_network(self):
        if self.graph is None:
            raise ValueError("Graph not loaded. Call load_graph() first.")
        # Example: Find shortest path between two nodes
        start_node = list(self.graph.nodes)[0]
        end_node = list(self.graph.nodes)[-1]
        shortest_path = nx.shortest_path(self.graph, source=start_node, target=end_node)
        return shortest_path

    def visualize_graph(self, output_image):
        if self.graph is None:
            raise ValueError("Graph not loaded. Call load_graph() first.")
        pos = nx.spring_layout(self.graph)
        nx.draw(self.graph, pos, with_labels=True)
        plt.savefig(output_image)
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python reax_network_solver.py <graph_file.gml> <output_image.png>")
        sys.exit(1)

    graph_file = sys.argv[1]
    output_image = sys.argv[2]

    solver = ReaxNetworkSolver(graph_file)
    solver.load_gml()
    shortest_path = solver.solve_network()
    print(f"Shortest path: {shortest_path}")
    solver.visualize_graph(output_image)
    print(f"Graph visualization saved to {output_image}")