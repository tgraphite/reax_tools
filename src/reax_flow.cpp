#include "reax_flow.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>

#ifdef ENABLE_DRAW_MOLECULE
#include "draw_molecule.h"
#endif
#include "fmt/format.h"
#include "universe.h"

int ReaxFlow::add_molecule(Molecule *mol) {
    // already exists
    for (size_t i = 0; i < nodes.size(); i++) {
        if (nodes[i]->formula == mol->formula) {
            return i;
        }
    }
    // create new
    Molecule *new_mol = new Molecule(*mol);
    nodes.emplace_back(new_mol);
    return nodes.size() - 1;
}

std::pair<int, int> ReaxFlow::add_reaction(int frame, Molecule *source, Molecule *target) {
    // lock it to avoid race condition when parallel execution
    std::lock_guard<std::mutex> lock(reaxflow_mutex);

    int source_id = add_molecule(source);
    int target_id = add_molecule(target);
    std::pair<int, int> edge = {source_id, target_id};

    // already exists
    for (size_t i = 0; i < edges.size(); i++) {
        if (edges[i] == edge) {
            edge_reaction_counts[i]++;
            return edge;
        }
    }
    // create new
    edges.push_back(edge);
    edge_reaction_counts.push_back(1);
    return edge;
}

void ReaxFlow::reduce_graph() {
    // Reduce reversive reactions.
    std::vector<int> edges_visited;
    int reaction_count_diff = 0;
    int reaction_count_to_reduce = 0;
    size_t reverse_edge_id = 0;
    std::pair<int, int> reverse_edge = {0, 0};

    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        if (std::find(edges_visited.begin(), edges_visited.end(), edge_id) != edges_visited.end()) continue;

        edges_visited.push_back(edge_id);
        reverse_edge = {edges[edge_id].second, edges[edge_id].first};

        reverse_edge_id = std::find(edges.begin(), edges.end(), reverse_edge) - edges.begin();

        // no reverse reaction
        if (reverse_edge_id == edges.size()) continue;

        // found reverse reaction
        reaction_count_diff = edge_reaction_counts[edge_id] - edge_reaction_counts[reverse_edge_id];
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
    }
}

void ReaxFlow::brief_report() { return; }
// Generate reaction flow report
// void ReaxFlow::brief_report() {
//     fmt::print("=== Reaction Flow Report ===\n");
//     fmt::print("Total key molecules: {}\n", nodes.size());
//     fmt::print("Total reactions: {}\n", edges.size());

//     // Sort edges by reaction count
//     std::vector<std::pair<int, int>> sorted_edge_id_count;
//     for (size_t i = 0; i < edges.size(); i++) {
//         sorted_edge_id_count.push_back({i, edge_reaction_counts[i]});
//     }
//     std::sort(sorted_edge_id_count.begin(), sorted_edge_id_count.end(),
//               [](const auto &a, const auto &b) { return a.second > b.second; });

//     // Display the top 10 most frequent reactions
//     fmt::print("Top reactions:\n");
//     int count = 0;
//     for (const auto &pair : sorted_edge_id_count) {
//         if (count >= 10) break;

//         const std::pair<int, int> edge = edges[pair.first];
//         fmt::print("{}: {} -> {} (count: {})\n", count + 1, nodes[edge.first]->formula, nodes[edge.second]->formula,
//                    pair.second);
//         count++;
//     }
// }

// Save reaction flow graph as DOT format for Graphviz.
void ReaxFlow::save_graph(const std::string &output_dir, int &max_reactions, bool draw_molecules,
                          bool reduce_reactions) {
    std::string edge_operator = "->";

    if (reduce_reactions) {
        reduce_graph();
        // edge_operator = "--";
        // Undirected graph
        // Try still use diagraph, 2025.05.29
    }
    int edge_id = 0;
    int source_node_id = 0;
    int target_node_id = 0;

    // get reactions priority by reaction count
    std::vector<std::pair<int, int>> sorted_edge_id_count;
    // Fill the vector with edge indices and their counts
    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        if (edge_reaction_counts[edge_id] > 0) {  // only consider reactions with positive counts
            sorted_edge_id_count.push_back({edge_id, edge_reaction_counts[edge_id]});
        }
    }
    // Sort by reaction count in descending order
    std::sort(sorted_edge_id_count.begin(), sorted_edge_id_count.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });

    // if max_reactions is specified, limit the number of reactions
    if (max_reactions > 0 && max_reactions < sorted_edge_id_count.size()) {
        sorted_edge_id_count.resize(max_reactions);
    }

    std::set<int> nodes_to_output;

    // Add selected reactions and their nodes
    for (const auto &pair : sorted_edge_id_count) {
        edge_id = pair.first;

        source_node_id = edges[edge_id].first;
        target_node_id = edges[edge_id].second;

        nodes_to_output.insert(source_node_id);
        nodes_to_output.insert(target_node_id);
    }

    // Display the top most frequent reactions
    fmt::print("=== Reaction Flow Report ===\n");
    fmt::print("Total key molecules: {}\n", nodes.size());
    fmt::print("Total reactions: {}\n", edges.size());

    for (const auto &pair : sorted_edge_id_count) {
        int edge_id = pair.first;
        int reaction_count = pair.second;
        fmt::print("{}: {} {} {} (count: {})\n", edge_id, nodes[edges[edge_id].first]->formula, edge_operator,
                   nodes[edges[edge_id].second]->formula, reaction_count);
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

    // -------------------------------------------------------------------
    // The main graph file.
    std::string save_graph_path = output_dir + "reaction_flow.dot";

    FILE *fp = fopen(save_graph_path.c_str(), "w");

    if (!fp) {
        fmt::print(stderr, "Error: Could not open dot file {} for writing.", save_graph_path);
        return;
    }

    // DOT file header
    // if (reduce_reactions) {
    //     fmt::print(fp, "strict graph ReactionFlow {{\n");
    // } else {
    //
    // }
    // Try still use diagraph, 2025.05.29
    fmt::print(fp, "digraph ReactionFlow {{\n");

    fmt::print(fp, "  rankdir=LR;\n");
    fmt::print(fp, "  layout=circo;\n");
    fmt::print(fp, "  node [shape=box, style=filled, fillcolor=azure2, height=0.5, width=1.5];\n");
    fmt::print(fp, "  edge [color=dimgray];\n");
    fmt::print(fp, "\n");

    // Write nodes
    for (const auto &node_id : nodes_to_output) {
        fmt::print(fp, "  node{} [label=\"[{}] {}\"];\n", node_id, node_id, nodes[node_id]->formula);
    }

    fmt::print(fp, "\n");

    // Write edges
    float penwidth = 1.0f;
    int reaction_count = 0;
    int max_key_reactions_count = 10;
    int key_reactions_count = 0;
    // Write edges based on sorted reaction counts
    for (const auto &pair : sorted_edge_id_count) {
        int edge_id = pair.first;
        reaction_count = pair.second;
        penwidth = std::min(5.0, 2.0 + log(reaction_count));

        if (key_reactions_count < max_key_reactions_count) {
            fmt::print(fp, " node{} {} node{} [label=\"{}\", penwidth={}, color=goldenrod];\n", edges[edge_id].first,
                       edge_operator, edges[edge_id].second, reaction_count, penwidth);
        } else {
            fmt::print(fp, " node{} {} node{} [label=\"{}\", penwidth={}];\n", edges[edge_id].first, edge_operator,
                       edges[edge_id].second, reaction_count, penwidth);
        }

        key_reactions_count++;
    }

    fmt::print(fp, "}}\n");
    fclose(fp);

    fmt::print("\nReaction graph saved to {}\n", save_graph_path);
}
