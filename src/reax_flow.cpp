#include "reax_flow.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <unordered_map>

#ifdef ENABLE_DRAW_MOLECULE
#include "draw_molecule.h"
#endif
#include "fmt/format.h"
#include "universe.h"

int ReaxFlow::add_molecule(Molecule *mol) {
    // If already exists
    // Use hash map for O(1) lookup
    auto it = molecule_map.find(mol->formula);
    if (it != molecule_map.end()) {
        return it->second;
    }

    // create new by copy
    // the current molecule in system will be destructed with system later.
    Molecule *new_mol = new Molecule(*mol);
    int new_id = nodes.size();
    nodes.emplace_back(new_mol);
    molecule_map[mol->formula] = new_id;
    return new_id;
}

std::pair<int, int> ReaxFlow::add_reaction(const int &frame, const int &atom_transfer_count, Molecule *source,
                                           Molecule *target) {
    // lock it to avoid race condition when parallel execution
    std::lock_guard<std::mutex> lock(reaxflow_mutex);

    int source_id = add_molecule(source);
    int target_id = add_molecule(target);
    std::pair<int, int> edge = {source_id, target_id};

    // If already exists, optimize by unordered_set later.
    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        if (edge == edges[edge_id]) {
            edge_reaction_counts[edge_id]++;
            edge_atom_transfer_counts[edge_id] += atom_transfer_count;
            return edge;
        }
    }

    // create new
    int new_edge_id = edges.size();
    edges.push_back(edge);
    edge_reaction_counts.push_back(1);
    edge_atom_transfer_counts.push_back(atom_transfer_count);
    return edge;
}

void ReaxFlow::reduce_graph() {
    // Reduce reversive reactions.
    std::vector<int> edges_visited;
    int reaction_count_diff = 0;
    int atom_transfer_count_diff = 0;

    size_t reverse_edge_id = 0;
    std::pair<int, int> reverse_edge = {0, 0};

    // Verify data consistency
    if (edges.size() != edge_reaction_counts.size()) {
        fmt::print(stderr, "Error: Inconsistent data in reduce_graph: edges size {} != edge_reaction_counts size {}\n",
                   edges.size(), edge_reaction_counts.size());
        return;
    }

    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        if (std::find(edges_visited.begin(), edges_visited.end(), edge_id) != edges_visited.end()) continue;

        edges_visited.push_back(edge_id);
        reverse_edge = {edges[edge_id].second, edges[edge_id].first};

        reverse_edge_id = std::find(edges.begin(), edges.end(), reverse_edge) - edges.begin();

        // no reverse reaction
        if (reverse_edge_id == edges.size()) continue;

        // found reverse reaction
        reaction_count_diff = edge_reaction_counts[edge_id] - edge_reaction_counts[reverse_edge_id];
        atom_transfer_count_diff = edge_atom_transfer_counts[edge_id] - edge_atom_transfer_counts[reverse_edge_id];

        if (reaction_count_diff > 0) {
            // Keep the net forward reaction count
            edge_reaction_counts[edge_id] = reaction_count_diff;
            edge_reaction_counts[reverse_edge_id] = 0;
        } else if (reaction_count_diff < 0) {
            // Keep the net reverse reaction count
            edge_reaction_counts[edge_id] = 0;
            edge_reaction_counts[reverse_edge_id] = -reaction_count_diff;
        } else {
            // cancel each other
            edge_reaction_counts[edge_id] = 0;
            edge_reaction_counts[reverse_edge_id] = 0;
        }

        if (atom_transfer_count_diff > 0) {
            // Keep the net forward reaction count
            edge_atom_transfer_counts[edge_id] = atom_transfer_count_diff;
            edge_atom_transfer_counts[reverse_edge_id] = 0;
        } else if (atom_transfer_count_diff < 0) {
            // Keep the net reverse reaction count
            edge_atom_transfer_counts[edge_id] = 0;
            edge_atom_transfer_counts[reverse_edge_id] = -atom_transfer_count_diff;
        } else {
            // cancel each other
            edge_atom_transfer_counts[edge_id] = 0;
            edge_atom_transfer_counts[reverse_edge_id] = 0;
        }
    }

    // Create new vectors for non-zero reactions
    std::vector<std::pair<int, int>> new_edges;
    std::vector<int> new_edge_reaction_counts;
    std::vector<int> new_edge_atom_transfer_counts;

    new_edges.reserve(edges.size());
    new_edge_reaction_counts.reserve(edge_reaction_counts.size());
    new_edge_atom_transfer_counts.reserve(edge_atom_transfer_counts.size());

    for (size_t tmp_id = 0; tmp_id < edge_reaction_counts.size(); tmp_id++) {
        if (edge_reaction_counts[tmp_id] == 0 || edge_atom_transfer_counts[tmp_id] == 0) continue;
        new_edges.push_back(edges[tmp_id]);
        new_edge_reaction_counts.push_back(edge_reaction_counts[tmp_id]);
        new_edge_atom_transfer_counts.push_back(edge_atom_transfer_counts[tmp_id]);
    }

    edges.clear();
    edge_reaction_counts.clear();
    edge_atom_transfer_counts.clear();

    edges = std::move(new_edges);
    edge_reaction_counts = std::move(new_edge_reaction_counts);
    edge_atom_transfer_counts = std::move(new_edge_atom_transfer_counts);
}

void ReaxFlow::brief_report() { return; }

void ReaxFlow::write_dot_file(const std::string &output_file, const std::vector<int> &edge_indices,
                              bool write_atom_transfer_count, std::string layout) {
    FILE *fp = fopen(output_file.c_str(), "w");
    std::set<int> node_indices;

    for (const auto &edge_id : edge_indices) {
        node_indices.insert(edges[edge_id].first);
        node_indices.insert(edges[edge_id].second);
    }

    if (!fp) {
        fmt::print(stderr, "Error: Could not open dot file {} for writing.", output_file);
        return;
    }

    // Dot file header.
    fmt::print(fp, "digraph ReactionFlow {{\n");
    fmt::print(fp, "  rankdir=LR;\n");
    fmt::print(fp, "  layout={};\n", layout);
    fmt::print(fp, "  node [shape=box, style=filled, fillcolor=azure2, height=0.5, width=1.5];\n");
    fmt::print(fp, "  edge [color=dimgray];\n");
    fmt::print(fp, "\n");

    // Write nodes
    for (const auto &node_id : node_indices) {
        fmt::print(fp, "  node{} [label=\"[{}] {}\"];\n", node_id, node_id, nodes[node_id]->formula);
    }
    fmt::print(fp, "\n");

    // Write edges, some of them are high-lighted because of bigger counts
    float penwidth = 1.0f;
    int max_highlights = 10;
    int curr_highlights = 0;
    int reaction_weight = 0;

    // Write edges based on sorted reaction counts
    for (const auto &edge_id : edge_indices) {
        if (write_atom_transfer_count) {
            reaction_weight = edge_atom_transfer_counts[edge_id];
            penwidth = std::min(5.0, 1.0 + log(reaction_weight));
        } else {
            reaction_weight = edge_reaction_counts[edge_id];
            penwidth = std::min(5.0, 2.0 + log(reaction_weight));
        }

        if (curr_highlights < max_highlights) {
            fmt::print(fp, " node{} -> node{} [label=\"{}\", penwidth={}, color=goldenrod];\n", edges[edge_id].first,
                       edges[edge_id].second, reaction_weight, penwidth);
        } else {
            fmt::print(fp, " node{} -> node{} [label=\"{}\", penwidth={}];\n", edges[edge_id].first,
                       edges[edge_id].second, reaction_weight, penwidth);
        }

        curr_highlights++;
    }

    fmt::print(fp, "}}\n");
    fclose(fp);

    fmt::print("\nReaction graph saved to {}\n", output_file);
}

// Save reaction flow graph as DOT format for Graphviz.
void ReaxFlow::save_graph(const std::string &output_dir, int &max_reactions, bool draw_molecules,
                          bool reduce_reactions) {
    std::string save_path_reactions = output_dir + "reactions.dot";
    std::string save_path_reactions_full = output_dir + "reactions_full.dot";

    std::string save_path_atom_transfer = output_dir + "atom_transfer.dot";
    std::string save_path_atom_transfer_full = output_dir + "atom_transfer_full.dot";

    if (reduce_reactions) {
        reduce_graph();
    }

    // Sort edges by reaction count
    std::vector<std::pair<int, int>> sorted_edge_id_count;
    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        if (edge_reaction_counts[edge_id] > 0) {  // only consider reactions with positive counts
            sorted_edge_id_count.push_back({edge_id, edge_reaction_counts[edge_id]});
        }
    }
    std::sort(sorted_edge_id_count.begin(), sorted_edge_id_count.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    // Resize the graph by max reactions
    if (max_reactions > 0 && max_reactions < sorted_edge_id_count.size()) {
        sorted_edge_id_count.resize(max_reactions);
    }

    std::vector<int> selected_edge_indices;
    for (const auto &pair : sorted_edge_id_count) {
        selected_edge_indices.push_back(pair.first);
    }

    // Display the topmost frequent reactions
    fmt::print("Note: To avoid too much screen output, only the top 20 reactions are printed.");
    fmt::print("=== Reaction Flow Report ===\n");
    fmt::print("Total key molecules: {}\n", nodes.size());
    fmt::print("Total reactions: {}\n", edges.size());

    int max_display = std::min(20, int(selected_edge_indices.size()));
    int tmp_edge_id = 0;
    int tmp_reaction_count = 0;
    int tmp_atom_transfer_count = 0;
    std::string tmp_source_formula = "";
    std::string tmp_target_formula = "";

    std::string display_header = fmt::format("{:<12s}{:<12s}{:<8s}{:<15s}", "from", "to", "count", "atom transfered");
    fmt::print("{}\n", display_header);

    for (size_t tmp_display_id = 0; tmp_display_id < max_display; tmp_display_id++) {
        tmp_edge_id = selected_edge_indices[tmp_display_id];
        tmp_source_formula = nodes[edges[tmp_edge_id].first]->formula;
        tmp_target_formula = nodes[edges[tmp_edge_id].second]->formula;
        tmp_reaction_count = edge_reaction_counts[tmp_edge_id];
        tmp_atom_transfer_count = edge_atom_transfer_counts[tmp_edge_id];

        fmt::print("{:<12s}{:<12s}{:<8d}{:<15d}\n", tmp_source_formula, tmp_target_formula, tmp_reaction_count,
                   tmp_atom_transfer_count);
    }

    // Write main dot file.
    write_dot_file(save_path_reactions, selected_edge_indices);
    write_dot_file(save_path_atom_transfer, selected_edge_indices, true);

    // If
    if (selected_edge_indices.size() < edges.size()) {
        std::vector<int> all_edge_indices;
        for (int tmp_id = 0; tmp_id < edges.size(); tmp_id++) {
            all_edge_indices.push_back(tmp_id);
        }

        fmt::print("Note: Reaction graph too complex, write full graph (*full.dot) and default subgraph seperately.\n");

        write_dot_file(save_path_reactions_full, all_edge_indices);
        write_dot_file(save_path_atom_transfer_full, selected_edge_indices, true);
    }

#ifdef ENABLE_DRAW_MOLECULE
    // Draw molecules
    if (draw_molecules) {
        for (const auto &node_id : nodes_to_output) {
            std::string save_path = output_dir + fmt::format("molecule_{}_{}.svg", node_id, nodes[node_id]->formula);
            draw_molecule(*nodes[node_id], save_path);
        }
        fmt::print("\n");
    }
#endif
}
