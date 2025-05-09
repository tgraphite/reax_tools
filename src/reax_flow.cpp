#include "reax_flow.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>

#include "draw_molecule.h"
#include "fmt/format.h"
#include "universe.h"

ReaxFlow::Node *ReaxFlow::get_node_from_molecule(Molecule *mol) {
    for (auto &node : nodes) {
        if (node->molecule == *mol) {
            return node;
        }
    }
    return nullptr;
}

ReaxFlow::Edge *ReaxFlow::get_edge_from_molecules(Molecule *source, Molecule *target) {
    for (auto &edge : edges) {
        if (edge->source_node->molecule == *source && edge->target_node->molecule == *target) {
            return edge;
        }
    }
    return nullptr;
}

void ReaxFlow::add_reaction(int frame, Molecule *source_mol, Molecule *target_mol) {
    // already locked in upper calling.

    // Get or create source node
    Node *current_source_node = nullptr;
    Node *current_target_node = nullptr;

    Node *tmp_source_node = get_node_from_molecule(source_mol);
    if (tmp_source_node) {
        current_source_node = tmp_source_node;
    } else {
        // Create new node
        Node *source_node = new Node(*source_mol);  // Use copy constructor
        source_node->formula = source_mol->formula;

        nodes.push_back(source_node);
        current_source_node = source_node;
    }

    Node *tmp_target_node = get_node_from_molecule(target_mol);
    if (tmp_target_node) {
        current_target_node = tmp_target_node;
    } else {
        // Create new node
        Node *target_node = new Node(*target_mol);  // Use copy constructor
        target_node->formula = target_mol->formula;

        nodes.push_back(target_node);
        current_target_node = target_node;
    }

    // Check if the edge already exists
    Edge *tmp_edge = get_edge_from_molecules(source_mol, target_mol);
    if (tmp_edge) {
        tmp_edge->reaction_count++;
    } else {
        // If the edge does not exist, create a new edge
        Edge *new_edge = new Edge();
        new_edge->source_node = current_source_node;
        new_edge->target_node = current_target_node;
        new_edge->reaction_count = 1;
        edges.push_back(new_edge);
    }
}

void ReaxFlow::reduce_graph(int max_reactions) {
    // Sort edges by reaction count in descending order
    std::vector<Edge *> sorted_edges = edges;
    std::sort(sorted_edges.begin(), sorted_edges.end(),
              [](const Edge *a, const Edge *b) { return a->reaction_count > b->reaction_count; });

    int reaction_count_threshold;

    if (sorted_edges.size() < max_reactions) {
        return;  // No need to reduce
    } else {
        reaction_count_threshold = sorted_edges[max_reactions - 1]->reaction_count;
    }

    // Create sets to track which nodes and edges to keep
    std::unordered_set<Node *> nodes_to_keep;
    std::unordered_set<Edge *> edges_to_keep;

    // First, identify edges to keep
    for (auto &edge : edges) {
        if (edge->reaction_count >= reaction_count_threshold) {
            edges_to_keep.insert(edge);
            nodes_to_keep.insert(edge->source_node);
            nodes_to_keep.insert(edge->target_node);
        }
    }

    // Delete edges that are not in edges_to_keep
    for (auto it = edges.begin(); it != edges.end();) {
        if (edges_to_keep.find(*it) == edges_to_keep.end()) {
            delete *it;
            it = edges.erase(it);
        } else {
            ++it;
        }
    }

    // Delete nodes that are not in nodes_to_keep
    for (auto it = nodes.begin(); it != nodes.end();) {
        if (nodes_to_keep.find(*it) == nodes_to_keep.end()) {
            delete *it;
            it = nodes.erase(it);
        } else {
            ++it;
        }
    }
}

// Generate reaction flow report
void ReaxFlow::brief_report() {
    std::cout << "=== Reaction Flow Report ===" << std::endl;
    std::cout << "Total nodes (species): " << nodes.size() << std::endl;
    std::cout << "Total edges (reactions): " << edges.size() << std::endl;

    // Sort edges by reaction count
    std::vector<std::pair<int, int>> sorted_edges;
    for (size_t i = 0; i < edges.size(); i++) {
        sorted_edges.push_back({i, edges[i]->reaction_count});
    }

    std::sort(sorted_edges.begin(), sorted_edges.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    // Display the top 10 most frequent reactions
    std::cout << "\nTop reactions:" << std::endl;
    int count = 0;
    for (const auto &pair : sorted_edges) {
        if (count >= 10) break;

        const Edge *edge = edges[pair.first];
        std::cout << fmt::format("{}: {} -> {} (count: {})", count + 1, edge->source_node->formula,
                                 edge->target_node->formula, pair.second)
                  << std::endl;
        count++;
    }
}

// Save reaction flow graph as DOT format for Graphviz.
void ReaxFlow::save_graph(const std::string &raw_file_path, int &max_reactions) {
    std::string save_basedir = raw_file_path.substr(0, raw_file_path.find_last_of(".")) + "_graph";

    // Create directory if not exists
    if (!std::filesystem::exists(save_basedir)) {
        std::filesystem::create_directory(save_basedir);
    } else {
        // Delete all files in directory if exists
        for (const auto &entry : std::filesystem::directory_iterator(save_basedir)) {
            std::filesystem::remove(entry.path());
        }
    }

    if (max_reactions > 0) {
        reduce_graph(max_reactions);
    }

    // Create node and edge id maps
    std::map<Node *, int> node_id_map;
    std::map<Edge *, int> edge_id_map;

    for (size_t i = 0; i < nodes.size(); i++) {
        node_id_map[nodes[i]] = i;
    }

    for (size_t i = 0; i < edges.size(); i++) {
        edge_id_map[edges[i]] = i;
    }

    std::string save_graph_path = save_basedir + "/reaction_flow.dot";
    std::ofstream file(save_graph_path);

    if (!file) {
        std::cerr << "Error: Could not open dot file {} for writing.";
        return;
    }

    // DOT file header
    file << "digraph ReactionFlow {\n";
    file << "  rankdir=LR;\n";
    file << "  node [shape=box, style=filled, fillcolor=lightblue];\n\n";

    // Write nodes
    for (size_t i = 0; i < nodes.size(); i++) {
        file << "  node" << i << " [label=\"" << nodes[i]->formula << "\"];\n";
    }

    file << "\n";

    // Write edges
    float penwidth = 1.0f;
    for (const auto &edge : edges) {
        penwidth = std::min(5.0, 1.0 + log(edge->reaction_count));
        file << " node" << node_id_map[edge->source_node] << " -> node" << node_id_map[edge->target_node]
             << " [label=\"" << edge->reaction_count << "\", penwidth=" << penwidth << "];\n";
    }

    file << "}\n";
    file.close();

    fmt::print("Reaction graph saved to {}\n", save_graph_path);

    // Draw molecules
    for (const auto &node : nodes) {
        std::string save_path = save_basedir + "/node_" + std::to_string(node_id_map[node]) + ".svg";
        draw_molecule(node->molecule, save_path);
    }
}
