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
        if (node->formula == mol->formula) {
            return node;
        }
    }
    return nullptr;
}

ReaxFlow::Edge *ReaxFlow::get_edge_from_molecules(Molecule *source, Molecule *target) {
    for (auto &edge : edges) {
        if (edge->source_node->formula == source->formula && edge->target_node->formula == target->formula) {
            return edge;
        }
    }
    return nullptr;
}

void ReaxFlow::add_reaction(int frame, Molecule *source_mol, Molecule *target_mol) {
    // already locked in upper calling.
    if (source_mol == target_mol) {
        return;
    }

    // If the reaction (edge) already exists, increment the reaction count and return.
    Edge *tmp_reaction = get_edge_from_molecules(source_mol, target_mol);

    if (tmp_reaction) {
        tmp_reaction->reaction_count++;
        tmp_reaction->source_node->reaction_count++;
        tmp_reaction->target_node->reaction_count++;
        return;
    }

    // ------------------------------------------------------------
    // Create new reaction.
    Node *current_source_node = nullptr;
    Node *current_target_node = nullptr;

    // Insertion is unsafe without lock, use class lock.
    std::lock_guard<std::mutex> lock(reaxflow_mutex);

    Node *tmp_source_node = get_node_from_molecule(source_mol);
    if (tmp_source_node) {
        current_source_node = tmp_source_node;
    } else {
        // Create new node, the default reaction count is 0.
        Node *source_node = new Node(*source_mol);  // Use copy constructor
        source_node->formula = source_mol->formula;

        nodes.push_back(source_node);
        current_source_node = source_node;
    }

    Node *tmp_target_node = get_node_from_molecule(target_mol);
    if (tmp_target_node) {
        current_target_node = tmp_target_node;
    } else {
        // Create new node, the default reaction count is 0.
        Node *target_node = new Node(*target_mol);  // Use copy constructor
        target_node->formula = target_mol->formula;

        nodes.push_back(target_node);
        current_target_node = target_node;
    }

    // Create new edge, the default reaction count is 0.
    Edge *new_edge = new Edge();
    new_edge->source_node = current_source_node;
    new_edge->target_node = current_target_node;

    // Initial increment.
    current_source_node->reaction_count++;
    current_target_node->reaction_count++;
    new_edge->reaction_count++;

    edges.push_back(new_edge);
}

void ReaxFlow::reduce_graph(int max_molecules) {
    // Only keep the top max_molecules nodes.
    std::unordered_set<Node *> nodes_to_remove;
    std::unordered_set<Edge *> edges_to_remove;

    // Sort edges by reaction count in descending order

    std::sort(nodes.begin(), nodes.end(),
              [](const Node *a, const Node *b) { return a->reaction_count > b->reaction_count; });

    // drop nodes with reaction count less than max_molecules
    if (max_molecules < nodes.size()) {
        for (int i = max_molecules; i < nodes.size(); i++) {
            nodes_to_remove.insert(nodes[i]);
        }
    } else {
        // keep all nodes
        return;
    }

    // drop edges connected to nodes to remove
    for (auto &edge : edges) {
        if (nodes_to_remove.count(edge->source_node) || nodes_to_remove.count(edge->target_node)) {
            edges_to_remove.insert(edge);
        }
    }

    // save delete edges and nodes
    for (auto it = edges.begin(); it != edges.end();) {
        if (edges_to_remove.count(*it)) {
            (*it)->source_node->reaction_count -= (*it)->reaction_count;
            (*it)->target_node->reaction_count -= (*it)->reaction_count;
            delete *it;
            it = edges.erase(it);
        } else {
            ++it;
        }
    }

    for (auto it = nodes.begin(); it != nodes.end();) {
        if (nodes_to_remove.count(*it)) {
            delete *it;
            it = nodes.erase(it);
        } else {
            ++it;
        }
    }

    edges_to_remove.clear();
    nodes_to_remove.clear();

    // Remove isolated nodes
    // In the process above, X-A-B-C-Y, when A and C are trivial and removed, B is isolated.
    // We need to remove B.
    for (auto &node : nodes) {
        if (node->reaction_count == 0) {
            nodes_to_remove.insert(node);
        }
    }

    for (auto it = nodes.begin(); it != nodes.end();) {
        if (nodes_to_remove.count(*it)) {
            delete *it;
            it = nodes.erase(it);
        } else {
            ++it;
        }
    }

    nodes_to_remove.clear();
}

// Generate reaction flow report
void ReaxFlow::brief_report() {
    fmt::print("=== Reaction Flow Report ===\n");
    fmt::print("Total key molecules: {}\n", nodes.size());
    fmt::print("Total reactions: {}\n", edges.size());

    // Sort edges by reaction count
    std::vector<std::pair<int, int>> sorted_edges;
    for (size_t i = 0; i < edges.size(); i++) {
        sorted_edges.push_back({i, edges[i]->reaction_count});
    }

    std::sort(sorted_edges.begin(), sorted_edges.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    // Display the top 10 most frequent reactions
    fmt::print("Top reactions:\n");
    int count = 0;
    for (const auto &pair : sorted_edges) {
        if (count >= 10) break;

        const Edge *edge = edges[pair.first];
        fmt::print("{}: {} -> {} (count: {})\n", count + 1, edge->source_node->formula, edge->target_node->formula,
                   pair.second);
        count++;
    }
}

// Save reaction flow graph as DOT format for Graphviz.
void ReaxFlow::save_graph(const std::string &output_dir, int &max_molecules, bool draw_molecules) {
    if (max_molecules > 0) {
        reduce_graph(max_molecules);
    }

    int key_molecules_count = 5;
    int key_reactions_count = 10;

    // Sort nodes and edges by reaction count in descending order
    // For highlighting the most important reactions and molecules.
    std::sort(nodes.begin(), nodes.end(),
              [](const Node *a, const Node *b) { return a->reaction_count > b->reaction_count; });
    std::sort(edges.begin(), edges.end(),
              [](const Edge *a, const Edge *b) { return a->reaction_count > b->reaction_count; });

    // Create node and edge id maps
    std::map<Node *, int> node_id_map;
    std::map<Edge *, int> edge_id_map;

    for (size_t i = 0; i < nodes.size(); i++) node_id_map[nodes[i]] = i;
    for (size_t i = 0; i < edges.size(); i++) edge_id_map[edges[i]] = i;

    // Draw molecules
    if (draw_molecules) {
        std::string formulas_to_draw;

        for (const auto &node : nodes) {
            formulas_to_draw += fmt::format("{} ", node->formula);
        }

        fmt::print("Molecules to draw: {}\n", formulas_to_draw);

        for (const auto &node : nodes) {
            std::string save_path = fmt::format("{}/molecule_{}_{}.svg", output_dir, node_id_map[node], node->formula);
            draw_molecule(node->molecule, save_path);
        }
    }

    // -------------------------------------------------------------------
    // The main graph file.
    std::string save_graph_path = output_dir + "reaction_flow.dot";

    FILE *fp = fopen(save_graph_path.c_str(), "w");

    if (!fp) {
        fmt::print(stderr, "Error: Could not open dot file {} for writing.", save_graph_path);
        return;
    }

    // DOT file header
    fmt::print(fp, "digraph ReactionFlow {{\n");
    fmt::print(fp, "  rankdir=LR;\n");
    fmt::print(fp, "  layout=circo;\n");
    fmt::print(fp, "  node [shape=box, style=filled, fillcolor=azure2, height=0.5, width=1.5];\n");
    fmt::print(fp, "  edge [color=dimgray];\n");
    fmt::print(fp, "\n");

    // Write nodes
    for (size_t i = 0; i < nodes.size(); i++) {
        // Highlight the key molecules, but if there are not many molecules, just draw simple.
        if (i < key_molecules_count && key_molecules_count < nodes.size()) {
            fmt::print(fp, "  node{} [label=\"[{}] {}\", fillcolor=cornsilk1];\n", i, i, nodes[i]->formula);
        } else {
            fmt::print(fp, "  node{} [label=\"[{}] {}\"];\n", i, i, nodes[i]->formula);
        }
    }

    fmt::print(fp, "\n");

    // Write edges
    float penwidth = 1.0f;
    for (size_t i = 0; i < edges.size(); i++) {
        penwidth = std::min(5.0, 2.0 + log(edges[i]->reaction_count));

        // Highlight the key reactions, but if there are not many reactions, just draw simple.
        if (i < key_reactions_count && key_reactions_count < edges.size()) {
            fmt::print(fp, " node{} -> node{} [label=\"{}\", penwidth={}, color=goldenrod];\n",
                       node_id_map[edges[i]->source_node], node_id_map[edges[i]->target_node], edges[i]->reaction_count,
                       penwidth);
        } else {
            fmt::print(fp, " node{} -> node{} [label=\"{}\", penwidth={}];\n", node_id_map[edges[i]->source_node],
                       node_id_map[edges[i]->target_node], edges[i]->reaction_count, penwidth);
        }
    }

    fmt::print(fp, "}}\n");
    fclose(fp);

    fmt::print("\nReaction graph saved to {}\n", save_graph_path);
}
