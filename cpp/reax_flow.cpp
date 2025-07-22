#include "reax_flow.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <unordered_map>

#include "argparser.h"
#include "fmt/format.h"
#include "string_tools.h"
#include "universe.h"

#ifndef WASM_MODE
#include "rdkit_utils.h"
#endif

ReaxFlow::~ReaxFlow() {
    // for (const auto& instance : nodes) {
    //     if (instance != nullptr)
    //         delete instance;
    // }

    for (const auto& pair : nodes) {
        if (pair.second != nullptr) {
            delete pair.second;
        }
    }

    nodes.clear();
}

/**
 * @brief Adds a molecule to the reaction flow system
 * @param mol Pointer to the molecule to be added
 * @return int The ID of the molecule in the system
 * @note If the molecule already exists, returns its existing ID
 */
inline int ReaxFlow::add_molecule(Molecule* mol) {
    // If already exists
    // Use hash map for O(1) lookup
    auto it = molecule_map.find(mol->hash);
    if (it != molecule_map.end()) {
        return it->second;
    }

    // create new by copy
    // the current molecule in system will be destructed with system later.
    Molecule* new_mol = new Molecule(*mol);
    int new_id = nodes.size();
    // nodes.emplace_back(new_mol);
    nodes[new_id] = new_mol;

    molecule_map[mol->hash] = new_id;
    return new_id;
}

/**
 * @brief Adds a reaction between two molecules to the system
 * @param frame Frame number of the reaction
 * @param atom_transfer_count Number of atoms transferred in the reaction
 * @param source Source molecule of the reaction
 * @param target Target molecule of the reaction
 * @return std::pair<int, int> Pair of molecule IDs representing the reaction edge
 * @note Thread-safe implementation with mutex lock
 */
std::pair<int, int> ReaxFlow::add_reaction(const int& frame_id, const int& atom_transfer_count, Molecule* source,
                                           Molecule* target) {
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
    // edges.push_back(edge);
    // edge_reaction_counts.push_back(1);
    // edge_atom_transfer_counts.push_back(atom_transfer_count);

    edges[new_edge_id] = edge;
    edge_reaction_counts[new_edge_id] = 1;
    edge_atom_transfer_counts[new_edge_id] = atom_transfer_count;

    return edge;
}

/**
 * @brief Reduces the reaction graph by handling reversible reactions
 * @details Combines forward and reverse reactions, keeping only the net reaction counts
 * and atom transfer counts. Removes zero-count reactions from the graph.
 */
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
        if (std::find(edges_visited.begin(), edges_visited.end(), edge_id) != edges_visited.end())
            continue;

        edges_visited.push_back(edge_id);
        reverse_edge = {edges[edge_id].second, edges[edge_id].first};

        // reverse_edge_id = std::find(edges.begin(), edges.end(), reverse_edge) - edges.begin();

        int reverse_edge_id = -1;
        for (const auto& pair : edges) {
            if (pair.second == reverse_edge)
                reverse_edge_id = pair.first;
        }

        // no reverse reaction
        if (reverse_edge_id == -1)
            continue;

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

    // std::vector<std::pair<int, int>> new_edges;
    // std::vector<int> new_edge_reaction_counts;
    // std::vector<int> new_edge_atom_transfer_counts;

    std::unordered_map<int, std::pair<int, int>> new_edges;
    std::unordered_map<int, int> new_edge_reaction_counts;
    std::unordered_map<int, int> new_edge_atom_transfer_counts;

    new_edges.reserve(edges.size());
    new_edge_reaction_counts.reserve(edge_reaction_counts.size());
    new_edge_atom_transfer_counts.reserve(edge_atom_transfer_counts.size());

    for (size_t tmp_id = 0; tmp_id < edge_reaction_counts.size(); tmp_id++) {
        if (edge_reaction_counts[tmp_id] == 0 || edge_atom_transfer_counts[tmp_id] == 0)
            continue;
        new_edges[tmp_id] = edges[tmp_id];
        new_edge_reaction_counts[tmp_id] = edge_reaction_counts[tmp_id];
        new_edge_atom_transfer_counts[tmp_id] = edge_atom_transfer_counts[tmp_id];
    }

    edges.clear();
    edge_reaction_counts.clear();
    edge_atom_transfer_counts.clear();

    edges = std::move(new_edges);
    edge_reaction_counts = std::move(new_edge_reaction_counts);
    edge_atom_transfer_counts = std::move(new_edge_atom_transfer_counts);
}

/**
 * @brief Generates a brief report of the reaction flow
 */
void ReaxFlow::brief_report() { return; }

/**
 * @brief Writes the reaction flow graph to a DOT file
 * @param output_file Path to the output DOT file
 * @param edge_indices Indices of edges to include in the graph
 * @param write_atom_transfer Whether to use atom transfer counts as edge weights
 * @param layout Graph layout algorithm to use
 */
void ReaxFlow::write_dot_file(std::string basename, const std::vector<int>& edge_indices, bool write_atom_transfer,
                              std::string layout) {
    FILE* fp = create_file(basename);
    std::set<int> node_indices;

    for (const auto& edge_id : edge_indices) {
        node_indices.insert(edges[edge_id].first);
        node_indices.insert(edges[edge_id].second);
    }

    // Dot file header.
    fmt::print(fp, "digraph ReactionFlow {{\n");
    fmt::print(fp, "  rankdir=LR;\n");
    fmt::print(fp, "  layout={};\n", layout);
    fmt::print(fp, "  node [shape=box, style=filled, fillcolor=azure2, height=0.5, width=1.5];\n");
    fmt::print(fp, "  edge [color=dimgray];\n");
    fmt::print(fp, "\n");

    // Write nodes
    for (const auto& node_id : node_indices) {
        fmt::print(fp, "  node{} [label=\"{}\"];\n", node_id, nodes[node_id]->formula);
    }
    fmt::print(fp, "\n");

    // Write edges, some of them are high-lighted because of bigger counts
    float penwidth = 1.0f;
    int max_highlights = edge_indices.size() / 4;
    int curr_highlights = 0;
    int reaction_weights = 0;

    // std::sort(sorted_edge_id_count.begin(), sorted_edge_id_count.end(),
    //           [](const auto &a, const auto &b) { return a.second > b.second; });
    std::vector<std::pair<int, int>> sorted_edge_id_weight;
    for (const auto& edge_id : edge_indices) {
        reaction_weights = edge_reaction_counts[edge_id];
        sorted_edge_id_weight.emplace_back(std::pair(edge_id, reaction_weights));
    }
    std::sort(sorted_edge_id_weight.begin(), sorted_edge_id_weight.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // Write edges based on sorted reaction counts
    std::string edge_line;
    for (const auto& id_weight : sorted_edge_id_weight) {
        penwidth = std::min(5.0, 2.0 + log(id_weight.second));

        if (write_atom_transfer) {
            edge_line = fmt::format(" node{} -> node{} [label=\"R={} AT={}\", penwidth={}",
                                    edges[id_weight.first].first, edges[id_weight.first].second, id_weight.second,
                                    edge_atom_transfer_counts[id_weight.first], penwidth);
        } else {
            edge_line = fmt::format(" node{} -> node{} [label=\"{}\", penwidth={}", edges[id_weight.first].first,
                                    edges[id_weight.first].second, id_weight.second, penwidth);
        }

        if (curr_highlights < max_highlights) {
            edge_line += ", color=goldenrod];\n";
        } else {
            edge_line += "];\n";
        }
        fmt::print(fp, edge_line);
        curr_highlights++;
    }

    fmt::print(fp, "}}\n");
    fclose(fp);
}

/**
 * @brief Saves the reaction flow graph and related visualizations
 * @param output_dir Directory to save the output files
 * @param max_reactions Maximum number of reactions to include
 * @param draw_molecules Whether to draw molecular structures
 * @param reduce_reactions Whether to reduce the reaction graph
 */
void ReaxFlow::save_graph() {
    if (!FLAG_NO_REDUCE_REACTIONS) {
        reduce_graph();
    }

    // Sort edges by reaction count
    std::vector<std::pair<int, int>> sorted_edge_id_count;
    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        if (edge_reaction_counts[edge_id] > 0) { // only consider reactions with positive counts
            sorted_edge_id_count.push_back({edge_id, edge_reaction_counts[edge_id]});
        }
    }
    std::sort(sorted_edge_id_count.begin(), sorted_edge_id_count.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    // Resize the graph by max reactions
    if (MAX_REACTIONS > 0 && MAX_REACTIONS < sorted_edge_id_count.size()) {
        sorted_edge_id_count.resize(MAX_REACTIONS);
    }

    std::vector<int> selected_edge_indices;
    for (const auto& pair : sorted_edge_id_count) {
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

    if (selected_edge_indices.size() < edges.size()) {
        std::vector<int> all_edge_indices;
        for (int tmp_id = 0; tmp_id < edges.size(); tmp_id++) {
            all_edge_indices.push_back(tmp_id);
        }
        write_dot_file("reactions.dot", selected_edge_indices, true);
        write_dot_file("reactions_full.dot", all_edge_indices, false);
        fmt::print("Note: Graphs too complex, write full graph (*full.dot) and default subgraph seperately.\n");
    } else {
        write_dot_file("reactions.dot", selected_edge_indices, true);
        fmt::print("Note: write reactions to reactions.dot.\n");
    }

    save_molecule_centered_subgraphs();
    // fmt::print("Reaction graphs centered on key moleulces saved to {}\n", output_dir);
}

/**
 * @brief Saves molecule-centered subgraphs of the reaction network
 * @param output_dir Directory to save the output files
 * @param write_atom_transfer Whether to write atom transfer counts in graph file and csv file
 * @details Generates subgraphs centered on key molecules, showing their reactions with other molecules
 */
void ReaxFlow::save_molecule_centered_subgraphs(bool write_atom_transfer, bool csv_only) {
    std::map<int, int> node_degrees;
    std::map<int, int> node_in_degrees;
    std::map<int, int> node_out_degrees;

    std::map<int, int> node_degrees_at;
    std::map<int, int> node_in_degrees_at;
    std::map<int, int> node_out_degrees_at;

    int max_key_molecules = 20;
    int max_molecule_size = 100;

    int source_id = -1;
    int target_id = -1;

    int curr_degree;
    int curr_degree_at;

    for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
        source_id = edges[edge_id].first;
        target_id = edges[edge_id].second;

        curr_degree = edge_reaction_counts[edge_id];
        curr_degree_at = edge_atom_transfer_counts[edge_id];

        node_degrees[source_id] += curr_degree;
        node_degrees[target_id] += curr_degree;
        node_in_degrees[target_id] += curr_degree;
        node_out_degrees[source_id] += curr_degree;

        node_degrees_at[source_id] += curr_degree_at;
        node_degrees_at[target_id] += curr_degree_at;
        node_in_degrees_at[target_id] += curr_degree_at;
        node_out_degrees_at[source_id] += curr_degree_at;
    }

    // Sort nodes by degree in descending order
    // Only output molecules smaller than max_molecule_size, we don't consider big polymers as a "main character"
    std::vector<std::pair<int, int>> sort_nodes;
    for (int node_id = 0; node_id < node_degrees.size(); node_id++) {
        if (nodes[node_id]->atom_ids.size() > max_molecule_size)
            continue;
        sort_nodes.emplace_back(std::pair(node_id, node_degrees[node_id]));
    }

    std::sort(sort_nodes.begin(), sort_nodes.end(),
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second > b.second; });

    if (sort_nodes.size() > max_key_molecules) {
        sort_nodes.resize(max_key_molecules);
    }

    // Write molecule-centered graph dot files
    FILE* fp_csv = create_file("key_molecules_reactions.csv");

    std::string csv_header = "";
    if (write_atom_transfer) {
        csv_header = "molecule,total reactions,in reaction,out reaction,total atom transfer,in atom transfer,out atom "
                     "transfer,from 1,from 2,from 3,from 4,from 5,to 1,to 2,to 3,to 4,to 5,\n";
    } else {
        csv_header = "molecule,total reactions,in reaction,out reaction,from 1,from 2,from 3,from 4,from 5,to 1,to "
                     "2,to 3,to 4,to 5,\n";
    }
    fmt::print(fp_csv, csv_header);

    for (const auto& pair : sort_nodes) {
        std::vector<int> subgraph_edge_indices;
        std::vector<std::pair<int, int>> from_nodes;
        std::vector<std::pair<int, int>> to_nodes;
        // Build molecule centered subgraph
        for (size_t edge_id = 0; edge_id < edges.size(); edge_id++) {
            if (edges[edge_id].first == pair.first) {
                to_nodes.emplace_back(std::pair(edges[edge_id].second, edge_reaction_counts[edge_id]));
                subgraph_edge_indices.emplace_back(edge_id);
            } else if (edges[edge_id].second == pair.first) {
                from_nodes.emplace_back(std::pair(edges[edge_id].first, edge_reaction_counts[edge_id]));
                subgraph_edge_indices.emplace_back(edge_id);
            }
        }

        std::string curr_formula = nodes[pair.first]->formula;
        // Write dot file.
        if (!csv_only) {
            std::string basename = fmt::format("reactions_centered_on_{}.dot", curr_formula);
            write_dot_file(basename, subgraph_edge_indices, write_atom_transfer, "dot");
        }
        // Sort from nodes and to nodes by reaction weight;
        std::sort(from_nodes.begin(), from_nodes.end(),
                  [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second > b.second; });
        std::sort(to_nodes.begin(), to_nodes.end(),
                  [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second > b.second; });

        size_t max_neigh_to_output = 5;
        std::string from_string;
        std::string to_string;
        for (size_t tmp_id = 0; tmp_id < max_neigh_to_output; tmp_id++) {
            if (tmp_id < from_nodes.size())
                from_string += fmt::format("{},", nodes[from_nodes[tmp_id].first]->formula);
            else
                from_string += ",";
        }
        for (size_t tmp_id = 0; tmp_id < max_neigh_to_output; tmp_id++) {
            if (tmp_id < to_nodes.size())
                to_string += fmt::format("{},", nodes[to_nodes[tmp_id].first]->formula);
            else
                to_string += ",";
        }

        std::string csv_record_string = "";

        if (write_atom_transfer) {
            csv_record_string =
                fmt::format("{},{},{},{},{},{},{},{}{}\n", curr_formula, node_degrees[pair.first],
                            node_in_degrees[pair.first], node_out_degrees[pair.first], node_degrees_at[pair.first],
                            node_in_degrees_at[pair.first], node_out_degrees_at[pair.first], from_string, to_string);
        } else {
            csv_record_string =
                fmt::format("{},{},{},{},{}{}\n", curr_formula, node_degrees[pair.first], node_in_degrees[pair.first],
                            node_out_degrees[pair.first], from_string, to_string);
        }

        fmt::print(fp_csv, csv_record_string);
    }
    fclose(fp_csv);
}

#ifndef WASM_MODE
void ReaxFlow::dump_smiles() {
    FILE* fp = create_file("molecules_smiles.csv");
    for (const auto& pair : nodes) {
        fmt::print(fp, "{},{}\n", pair.second->formula, rdkit_smiles(*pair.second));
    }

    fclose(fp);
}
#endif

#ifndef WASM_MODE
void ReaxFlow::draw_molecules() {
    for (const auto& pair : nodes) {
        rdkit_draw_molecule(*pair.second);
    }
}
#endif
