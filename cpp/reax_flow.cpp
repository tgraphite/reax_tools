#include "reax_flow.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <random>
#include <stack>
#include <unordered_map>

#include "argparser.h"
#include "fmt/format.h"
#include "string_tools.h"
#include "universe.h"

#ifndef WASM_MODE
#include "rdkit_utils.h"
#endif

/**
 * @brief Constructs a Node object representing a molecule in the reaction graph
 * @param mol Pointer to the molecule to be represented by this node
 * @throws ReaxFlowException if mol is nullptr
 * @note Creates a deep copy of the molecule to ensure the original can be
 * safely destroyed
 */
Node::Node(Molecule* mol) {
    // deep copy a molecule in system, thus the origin one can be detroyed safely.
    if (mol == nullptr) {
        throw ReaxFlowException("Invalid molecule when construct node.");
    }

    molecule = new Molecule(*mol);
    hash = molecule->hash;
}

/**
 * @brief Destructor for Node object
 * @note Safely deallocates the molecule pointer to prevent memory leaks
 */
Node::~Node() {
    if (molecule != nullptr) {
        delete molecule;
    }
}

/**
 * @brief Updates the degree statistics for this node based on reaction
 * participation
 * @param source_or_target True if this node is the source (outgoing), false if
 * target (incoming)
 * @param count Number of reactions involving this node
 * @param atom_transfer_count Total number of atoms transferred in reactions
 * involving this node
 * @note Updates topological degrees (edge counts) and weighted degrees (reaction counts)
 * @note degree = precursor_count + derivative_count (topological)
 * @note reaction_count = in_reaction_count + out_reaction_count (weighted)
 */
void Node::add_degrees(bool source_or_target, unsigned int count, unsigned int atom_transfer_count) {
    if (source_or_target) {
        // As reactant/source: forms derivatives
        derivative_count++;                     // topological
        derivative_reactions += count;          // weighted
        derivative_atom_transfer += atom_transfer_count;
    }
    else {
        // As product/target: has precursors
        precursor_count++;                      // topological
        precursor_reactions += count;           // weighted
        precursor_atom_transfer += atom_transfer_count;
    }

    // Update totals
    degree = precursor_count + derivative_count;
    reaction_count = precursor_reactions + derivative_reactions;
    atom_transfer = precursor_atom_transfer + derivative_atom_transfer;
}

/**
 * @brief Constructs an Edge object representing a reaction between two
 * molecules
 * @param from_node Pointer to the source node (reactant)
 * @param to_node Pointer to the target node (product)
 * @throws ReaxFlowException if either node is nullptr
 * @note Automatically generates a hash for efficient lookup
 */
Edge::Edge(Node* from_node, Node* to_node) {
    if (from_node == nullptr || to_node == nullptr) {
        throw ReaxFlowException("Invalid node when construct edge.");
    }

    source = from_node;
    target = to_node;
    hash = get_edge_hash(source, target);
}

/**
 * @brief Destructor for Edge object
 * @note No dynamic memory to deallocate, but provides virtual destructor for
 * inheritance
 */
Edge::~Edge() {}

/**
 * @brief Default constructor for ReaxFlow object
 * @note Initializes empty reaction graph with thread-safe mutex
 */
ReaxFlow::ReaxFlow() {}

/**
 * @brief Destructor for ReaxFlow object
 * @note Safely deallocates all nodes and edges to prevent memory leaks
 */
ReaxFlow::~ReaxFlow() {
    for (auto& edge : edges) {
        if (edge != nullptr) {
            delete edge;
        }
    }
    edges.clear();

    for (auto& node : nodes) {
        if (node != nullptr) {
            delete node;
        }
    }
    nodes.clear();

    // Clear the hash map
    molecule_hash_to_node.clear();
}

/**
 * @brief Retrieves a node by molecule pointer using hash-based lookup
 * @param mol Pointer to the molecule to search for
 * @return Node* Pointer to the node if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Node* ReaxFlow::get_node(Molecule* mol) {
    if (mol == nullptr) {
        return nullptr;
    }

    // Use hash map for O(1) lookup
    auto it = molecule_hash_to_node.find(mol->hash);
    if (it != molecule_hash_to_node.end()) {
        return it->second;
    }

    return nullptr;
}

/**
 * @brief Retrieves a node by hash value using hash-based lookup
 * @param hash Hash value of the molecule to search for
 * @return Node* Pointer to the node if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Node* ReaxFlow::get_node(unsigned int hash) {
    auto it = molecule_hash_to_node.find(hash);
    if (it != molecule_hash_to_node.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Retrieves an edge between two nodes using hash-based lookup
 * @param from_node Source node of the edge
 * @param to_node Target node of the edge
 * @return Edge* Pointer to the edge if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Edge* ReaxFlow::get_edge(Node* from_node, Node* to_node) {
    if (from_node == nullptr || to_node == nullptr) {
        return nullptr;
    }

    unsigned int hash = get_edge_hash(from_node, to_node);
    auto it = edge_hash_to_edge.find(hash);
    if (it != edge_hash_to_edge.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Retrieves an edge by hash value using hash-based lookup
 * @param hash Hash value of the edge to search for
 * @return Edge* Pointer to the edge if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Edge* ReaxFlow::get_edge(unsigned int hash) {
    auto it = edge_hash_to_edge.find(hash);
    if (it != edge_hash_to_edge.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Adds a molecule to the reaction flow system
 * @param mol Pointer to the molecule to be added
 * @return Node* Pointer to the node representing the molecule
 * @note If the molecule already exists, returns its existing node. Otherwise
 * creates a new node.
 * @note Thread-safe implementation with automatic node creation
 */
Node* ReaxFlow::add_molecule(Molecule* mol) {
    if (mol == nullptr) {
        return nullptr;
    }

    // Check if molecule already exists using hash map
    auto it = molecule_hash_to_node.find(mol->hash);
    if (it != molecule_hash_to_node.end()) {
        return it->second;  // Return existing node
    }

    // Create new node
    Node* new_node = new Node(mol);
    nodes.insert(new_node);
    molecule_hash_to_node[mol->hash] = new_node;

    return new_node;
}

/**
 * @brief Adds a reaction between two molecules to the system
 * @param frame_id Frame number of the reaction (currently unused but preserved
 * for future use)
 * @param atom_transfer_count Number of atoms transferred in the reaction
 * @param source Source molecule of the reaction (reactant)
 * @param target Target molecule of the reaction (product)
 * @return Edge* Pointer to the edge representing the reaction
 * @note Thread-safe implementation with mutex lock
 * @note If the reaction already exists, updates the existing edge counts
 * @note Automatically creates nodes for molecules if they don't exist
 */
Edge* ReaxFlow::add_reaction(const int& frame_id, const int& atom_transfer_count, Molecule* source, Molecule* target) {
    // lock for parallel, dynamic find & create in unordered_set / unordered_map
    // is thread unsafe.
    std::lock_guard<std::mutex> lock(reaxflow_mutex);

    Node* source_node = add_molecule(source);
    Node* target_node = add_molecule(target);

    if (source_node == nullptr || target_node == nullptr) {
        return nullptr;
    }

    // Check if edge already exists
    Edge* existing_edge = get_edge(source_node, target_node);
    if (existing_edge != nullptr) {
        // Update existing edge
        existing_edge->count++;
        existing_edge->atom_transfer += atom_transfer_count;

        return existing_edge;
    }

    // Create new edge
    Edge* new_edge = new Edge(source_node, target_node);
    new_edge->count++;
    new_edge->atom_transfer += atom_transfer_count;

    edges.insert(new_edge);
    edge_hash_to_edge[new_edge->hash] = new_edge;

    return new_edge;
}

/**
 * @brief Reduces the reaction graph by handling reversible reactions
 * @details Combines forward and reverse reactions, keeping only the net
 * reaction counts and atom transfer counts. Removes zero-count reactions from
 * the graph.
 */
void ReaxFlow::reduce_graph() {
    std::vector<Edge*> edges_to_remove;
    std::unordered_set<Edge*> processed_edges;

    for (Edge* edge1 : edges) {
        // Skip if already processed or marked for removal
        if (processed_edges.find(edge1) != processed_edges.end()) {
            continue;
        }

        Edge* edge2 = get_edge(edge1->target, edge1->source);

        if (edge2 == nullptr) {
            continue;
        }

        // Skip if edge2 is already processed or marked for removal
        if (processed_edges.find(edge2) != processed_edges.end()) {
            continue;
        }

        // Check if these are reverse edges
        if (edge1->source == edge2->target && edge1->target == edge2->source) {
            int reaction_count_diff = edge1->count - edge2->count;
            int atom_transfer_count_diff = edge1->atom_transfer - edge2->atom_transfer;

            if (reaction_count_diff > 0) {
                // Keep edge1, remove edge2
                edge1->count = reaction_count_diff;
                edge1->atom_transfer = atom_transfer_count_diff;
                edges_to_remove.push_back(edge2);
                processed_edges.insert(edge1);
                processed_edges.insert(edge2);
            }
            else if (reaction_count_diff < 0) {
                // Keep edge2, remove edge1
                edge2->count = -reaction_count_diff;
                edge2->atom_transfer = -atom_transfer_count_diff;
                edges_to_remove.push_back(edge1);
                processed_edges.insert(edge1);
                processed_edges.insert(edge2);
            }
            else {
                // Both cancel each other, remove both
                edges_to_remove.push_back(edge1);
                edges_to_remove.push_back(edge2);
                processed_edges.insert(edge1);
                processed_edges.insert(edge2);
            }
        }
    }

    // Remove edges marked for removal
    for (Edge* edge : edges_to_remove) {
        edges.erase(edge);
        edge_hash_to_edge.erase(edge->hash);
        delete edge;
    }
}

/**
 * @brief Sorts nodes and edges by their degree/count for efficient reporting
 * @note Calls calc_node_degrees() to ensure degree data is current
 * @note Clears and rebuilds sorted_nodes and sorted_edges vectors
 * @note Sorts in descending order (highest degree/count first)
 */
void ReaxFlow::update_graph() {
    // Clear all degree statistics and adjacency relationships
    for (auto& node : nodes) {
        // Topological degrees (precursors and derivatives)
        node->degree = 0;
        node->precursor_count = 0;
        node->derivative_count = 0;

        // Weighted degrees (reaction counts)
        node->reaction_count = 0;
        node->precursor_reactions = 0;
        node->derivative_reactions = 0;

        // Atom transfers
        node->atom_transfer = 0;
        node->precursor_atom_transfer = 0;
        node->derivative_atom_transfer = 0;

        // Adjacency (precursors and derivatives)
        node->from_nodes.clear();
        node->to_nodes.clear();
    }

    // Build adjacency relationships and calculate degrees
    for (const auto& edge : edges) {
        edge->source->add_degrees(true, edge->count, edge->atom_transfer);
        edge->target->add_degrees(false, edge->count, edge->atom_transfer);

        // Build adjacency relationships
        edge->source->to_nodes.insert(edge->target);
        edge->target->from_nodes.insert(edge->source);
    }

    sorted_nodes.clear();
    sorted_edges.clear();

    for (const auto& node : nodes) {
        sorted_nodes.emplace_back(std::pair(node, node->degree));
    }
    std::sort(sorted_nodes.begin(), sorted_nodes.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });

    for (const auto& edge : edges) {
        sorted_edges.emplace_back(std::pair(edge, edge->count));
    }
    std::sort(sorted_edges.begin(), sorted_edges.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });
}

/**
 * @brief Generates a brief report of the reaction flow statistics
 * @note Displays top 10 molecules by degree and top 20 reactions by count
 * @note Shows in-degree and out-degree for molecules, reaction counts and atom
 * transfers for edges
 * @note Automatically calls sort_nodes_and_edges() to ensure data is current
 */
void ReaxFlow::brief_report() {
    update_graph();

    unsigned int max_node_display = std::min(10, int(sorted_nodes.size()));
    unsigned int max_edge_display = std::min(20, int(sorted_edges.size()));

    fmt::print("=== Reaction Flow Report ===\n");
    fmt::print("Top {} key molecules:\n", max_node_display);
    fmt::print("{:<12s}{:<12s}{:<12s}\n", "molecule", "precursors", "derivatives");
    Node* tmp_node = nullptr;
    for (size_t i = 0; i < max_node_display; i++) {
        tmp_node = sorted_nodes[i].first;
        fmt::print("{:<12s}{:<12d}{:<12d}\n", tmp_node->molecule->formula, tmp_node->precursor_count, tmp_node->derivative_count);
    }

    fmt::print("\n");
    fmt::print("Top {} key reactions:\n", max_edge_display);
    fmt::print("{:<12s}{:<12s}{:<8s}{:<15s}\n", "from", "to", "count", "atom transfered");
    Node* tmp_source = nullptr;
    Node* tmp_target = nullptr;
    for (size_t i = 0; i < max_edge_display; i++) {
        tmp_source = sorted_edges[i].first->source;
        tmp_target = sorted_edges[i].first->target;
        fmt::print("{} -> {} = R:{} AT:{}\n", tmp_source->molecule->formula, tmp_target->molecule->formula,
            sorted_edges[i].second, sorted_edges[i].first->atom_transfer);
    }

    tmp_node = nullptr;
    tmp_source = nullptr;
    tmp_target = nullptr;
}

/**
 * @brief Writes the reaction flow graph to a DOT file for visualization
 * @param basename Path to the output DOT file
 * @param edges_to_write Set of edges to include in the graph
 * @param write_atom_transfer Whether to use atom transfer counts as edge
 * weights
 * @param layout Graph layout algorithm to use (e.g., "circo", "dot", "neato")
 * @note Creates a directed graph with nodes representing molecules and edges
 * representing reactions
 * @note Edge thickness is proportional to reaction count (logarithmic scale)
 * @note Highlights top 25% of edges in goldenrod color
 */
void ReaxFlow::write_dot_file(std::string basename, const std::vector<Edge*>& edges_to_write,
    bool write_atom_transfer, std::string layout) {
    FILE* fp = create_file(basename);

    std::unordered_set<Node*> nodes_to_write;
    for (const auto& edge : edges_to_write) {
        nodes_to_write.insert(edge->source);
        nodes_to_write.insert(edge->target);
    }

    // Dot file header.
    fmt::print(fp, "digraph ReactionFlow {{\n");
    fmt::print(fp, "  rankdir=LR;\n");
    fmt::print(fp, "  layout={};\n", layout);
    fmt::print(fp,
        "  node [shape=box, style=filled, fillcolor=azure2, height=0.5, "
        "width=1.5];\n");
    fmt::print(fp, "  edge [color=dimgray];\n");
    fmt::print(fp, "\n");

    // Write nodes
    Node* node = nullptr;
    for (const auto& node : nodes_to_write) {
        // Use formula_hash format for unique but readable IDs
        fmt::print(fp, "  \"{}_{}\" [label=\"{}\"];\n",
            node->molecule->formula, node->hash, node->molecule->formula);
    }
    fmt::print(fp, "\n");

    // Write edges, some of them are high-lighted because of bigger counts
    float penwidth = 1.0f;
    int max_highlights = edges_to_write.size() / 4;
    int curr_highlights = 0;

    // Write edges based on sorted reaction counts
    std::string edge_line;

    for (const auto& edge : edges_to_write) {
        penwidth = std::min(5.0, 2.0 + log(edge->count));

        // Use formula_hash format for node references
        std::string src_id = fmt::format("{}_{}", edge->source->molecule->formula, edge->source->hash);
        std::string tgt_id = fmt::format("{}_{}", edge->target->molecule->formula, edge->target->hash);

        if (write_atom_transfer) {
            edge_line = fmt::format(" \"{}\" -> \"{}\" [label=\"R={} AT={}\", penwidth={}",
                src_id, tgt_id, edge->count, edge->atom_transfer, penwidth);
        }
        else {
            edge_line = fmt::format(" \"{}\" -> \"{}\" [label=\"{}\", penwidth={}",
                src_id, tgt_id, edge->count, penwidth);
        }

        if (curr_highlights < max_highlights) {
            edge_line += ", color=goldenrod];\n";
        }
        else {
            edge_line += "];\n";
        }
        fmt::print(fp, edge_line);
        curr_highlights++;
    }

    fmt::print(fp, "}}\n");
    fclose(fp);
}

/**
 * @brief Writes a DOT file for the most significant nodes in the reaction
 * graph
 * @param basename Path to the output DOT file
 * @param max_nodes Maximum number of top nodes to include based on degree
 * @param write_atom_transfer Whether to use atom transfer counts as edge
 * weights
 * @param layout Graph layout algorithm to use (e.g., "circo", "dot", "neato")
 * @note Selects top nodes by degree and includes their connecting edges
 * @note Calls write_dot_file() internally to generate the DOT file
 */
void ReaxFlow::write_dot_file_significant_nodes(std::string basename, int max_nodes,
    bool write_atom_transfer, std::string layout) {

    max_nodes = std::min(max_nodes, int(sorted_nodes.size()));
    int max_depth = 2 * max_nodes;

    std::vector<Node*> selected_nodes;
    std::vector<std::pair<Edge*, int>> organized_connections;

    // Organized connections: find paths between node pairs in full graph, 
    // for each path, the weight is the min count along the path
    // for the final connection, sum all weights if multiple paths exist.
    for (int i = 0; i < max_nodes; i++) {
        selected_nodes.push_back(sorted_nodes[i].first);
    }

    std::vector<std::pair<Node*, Node*>> possible_node_pairs;
    for (const auto& left_node : selected_nodes) {
        for (const auto& right_node : selected_nodes) {
            if (left_node == right_node) {
                continue;
            }
            possible_node_pairs.push_back(std::make_pair(left_node, right_node));
        }
    }

    for (const auto& node_pair : possible_node_pairs) {
        Node* left_node = node_pair.first;
        Node* right_node = node_pair.second;
        if (left_node == nullptr || right_node == nullptr) {
            continue;
        }
        // Use DFS to find all paths from left_node to right_node
        std::vector<std::vector<Edge*>> all_paths;
        std::vector<Edge*> current_path;
        std::unordered_set<Node*> visited;

        std::function<void(Node*)> dfs = [&](Node* current_node) {
            if (current_node == right_node) {
                all_paths.push_back(current_path);
                return;
            }

            if (visited.find(current_node) != visited.end()) {
                return;
            }

            if (visited.size() >= max_depth || visited.find(current_node) != visited.end()) {
                return;
            }

            visited.insert(current_node);

            for (const auto& neighbor : current_node->to_nodes) {
                Edge* edge = get_edge(current_node, neighbor);
                if (edge != nullptr && visited.find(neighbor) == visited.end()) {
                    current_path.push_back(edge);
                    dfs(neighbor);
                    current_path.pop_back();
                }
            }
            visited.erase(current_node);
            };

        dfs(left_node);

        if (all_paths.empty()) {
            continue;
        }

        // Process all found paths to get edge weights
        int connection_weight = 0;
        for (const auto& path : all_paths) {
            int path_min = std::numeric_limits<int>::max();
            for (const auto& edge : path) {
                path_min = std::min(path_min, edge->count);
            }
            connection_weight = std::max(connection_weight, path_min);
        }

        // Create connection
        Edge* this_connection = new Edge(left_node, right_node);
        this_connection->count = connection_weight;
        this_connection->atom_transfer = connection_weight; // not used here, just placeholder

        organized_connections.push_back(std::make_pair(this_connection, connection_weight));
    }

    // Sort connections by weight
    std::sort(organized_connections.begin(), organized_connections.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });

    std::vector<Edge*> edges_to_write;
    for (const auto& pair : organized_connections) {
        edges_to_write.push_back(pair.first);
    }

    write_dot_file(basename, edges_to_write, write_atom_transfer, layout);

    for (const auto& conn : organized_connections) {
        delete conn.first;
    }
}

/**
 * @brief Saves the reaction flow graph and related visualizations
 * @note Conditionally calls reduce_graph() unless FLAG_NO_REDUCE_REACTIONS is
 * set
 * @note Generates both full and simplified graph files if graph is too complex
 * @note Creates molecule-centered subgraphs for detailed analysis
 * @note Automatically sorts nodes and edges before saving
 */
void ReaxFlow::save_graph() {
    // Integrated filtering workflow
    fmt::print("\n=== Integrated Reaction Network Filtering ===\n");

    // Step 1: Reduce reversible reactions
    if (!FLAG_NO_REDUCE_REACTIONS) {
        fmt::print("Step 1: Reducing reversible reactions...\n");
        reduce_graph();
    }

    // Step 2: Atom economy filtering
    fmt::print("Step 2: Atom economy filtering...\n");
    filter_by_atom_economy(0.3);

    // Step 3: Edge betweenness ranking if still too many edges
    if (edges.size() > MAX_REACTIONS) {
        fmt::print("Step 3: Edge betweenness ranking ({} > {})...\n", edges.size(), MAX_REACTIONS);
        filter_by_betweenness(MAX_REACTIONS);
    }

    // Step 4: Clean up isolated nodes and invalid edges
    fmt::print("Step 4: Cleaning up isolated nodes...\n");
    cleanup_isolated_nodes();

    update_graph();
    brief_report();

    std::vector<Edge*> all_edges;

    for (size_t i = 0; i < sorted_edges.size(); i++) {
        all_edges.push_back(sorted_edges[i].first);
    }

    if (all_edges.size() > MAX_REACTIONS) {
        std::vector<Edge*> selected_edges;
        for (size_t i = 0; i < MAX_REACTIONS; i++) {
            selected_edges.push_back(sorted_edges[i].first);
        }
        write_dot_file("reactions_simplified.dot", selected_edges, false);
        write_dot_file("reactions_full.dot", all_edges, false);
        fmt::print(
            "Note: Graphs too complex, write full graph (*full.dot) and default "
            "subgraph seperately.\n");
    }
    else {
        write_dot_file("reactions_full.dot", all_edges, false);
    }

    // write_dot_file_significant_nodes("reactions_main_nodes.dot", 20, true);

    save_molecule_centered_subgraphs(false, false, false);
    // save_molecule_centered_subgraphs(false, false, true);
}

/**
 * @brief Saves molecule-centered subgraphs of the reaction network
 * @param write_atom_transfer Whether to write atom transfer counts in graph
 * file and csv file
 * @param csv_only Whether to only generate CSV output (currently unused
 * parameter)
 * @details Generates subgraphs centered on key molecules, showing their
 * reactions with other molecules
 * @note Creates CSV file with detailed reaction statistics for top molecules
 * @note Automatically sorts nodes by degree before processing
 * @note Limits output to top 10 molecules by degree
 */
void ReaxFlow::save_molecule_centered_subgraphs(bool write_atom_transfer, bool csv_only, bool use_hash) {
    update_graph();

    if (sorted_nodes.size() > MAX_KEY_MOLECULES) {
        sorted_nodes.resize(MAX_KEY_MOLECULES);
    }

    // Output as Markdown instead of CSV
    std::string save_path = "key_molecules_reactions.md";
    FILE* fp_md = create_file(save_path);

    // Markdown header
    fmt::print(fp_md, "# Key Molecules Reaction Network\n\n");
    fmt::print(fp_md, "This report shows the most important molecules in the reaction network and their relationships.\n\n");
    fmt::print(fp_md, "- **Precursors**: Molecules that react to form this species (incoming)\n");
    fmt::print(fp_md, "- **Derivatives**: Molecules formed from this species (outgoing)\n\n");
    fmt::print(fp_md, "---\n\n");

    int mol_index = 1;
    for (const auto& [node, degree] : sorted_nodes) {
        std::string species_name = node->molecule->formula;

        // Species header
        fmt::print(fp_md, "## {}. Species: **{}**\n\n", mol_index++, species_name);

        // Topology statistics
        fmt::print(fp_md, "### Network Statistics\n\n");
        fmt::print(fp_md, "- **Total connections**: {}\n", node->degree);
        fmt::print(fp_md, "- **Precursors**: {}\n", node->precursor_count);
        fmt::print(fp_md, "- **Derivatives**: {}\n", node->derivative_count);

        if (write_atom_transfer) {
            fmt::print(fp_md, "- **Total reaction count**: {}\n", node->reaction_count);
            fmt::print(fp_md, "- **Precursor reactions**: {}\n", node->precursor_reactions);
            fmt::print(fp_md, "- **Derivative reactions**: {}\n", node->derivative_reactions);
            fmt::print(fp_md, "- **Total atom transfer**: {}\n", node->atom_transfer);
        }
        fmt::print(fp_md, "\n");

        // Precursors section
        fmt::print(fp_md, "### Precursors ({} total)\n\n", node->precursor_count);
        if (!node->from_nodes.empty()) {
            fmt::print(fp_md, "Top precursors (sorted by reaction count):\n\n");
            
            // Get edges from precursors to this node
            std::vector<std::pair<Node*, Edge*>> precursor_edges;
            for (auto* precursor : node->from_nodes) {
                Edge* edge = get_edge(precursor, node);
                if (edge) {
                    precursor_edges.emplace_back(precursor, edge);
                }
            }
            
            // Sort by reaction count (edge->count)
            std::sort(precursor_edges.begin(), precursor_edges.end(),
                [](const auto& a, const auto& b) { return a.second->count > b.second->count; });
            
            int count = 0;
            for (auto& [precursor, edge] : precursor_edges) {
                if (++count > 5) break;
                fmt::print(fp_md, "- **{}** → {} (reactions: {})\n", 
                    precursor->molecule->formula, 
                    node->molecule->formula,
                    edge->count);
            }
            fmt::print(fp_md, "\n");
        }
        else {
            fmt::print(fp_md, "_No precursors (initial species)_\n\n");
        }

        // Derivatives section
        fmt::print(fp_md, "### Derivatives ({} total)\n\n", node->derivative_count);
        if (!node->to_nodes.empty()) {
            fmt::print(fp_md, "Top derivatives (sorted by reaction count):\n\n");
            
            // Get edges from this node to derivatives
            std::vector<std::pair<Node*, Edge*>> derivative_edges;
            for (auto* derivative : node->to_nodes) {
                Edge* edge = get_edge(node, derivative);
                if (edge) {
                    derivative_edges.emplace_back(derivative, edge);
                }
            }
            
            // Sort by reaction count (edge->count)
            std::sort(derivative_edges.begin(), derivative_edges.end(),
                [](const auto& a, const auto& b) { return a.second->count > b.second->count; });

            int count = 0;
            for (auto& [derivative, edge] : derivative_edges) {
                if (++count > 5) break;
                fmt::print(fp_md, "- **{}** → {} (reactions: {})\n",
                    node->molecule->formula,
                    derivative->molecule->formula,
                    edge->count);
            }
            fmt::print(fp_md, "\n");
        }
        else {
            fmt::print(fp_md, "_No derivatives (terminal species)_\n\n");
        }
        fmt::print(fp_md, "---\n\n");
    }

    // Add reaction pathways section
    fmt::print(fp_md, "# Reaction Pathways\n\n");
    fmt::print(fp_md, "This section shows the longest continuous reaction pathways in the network.\n\n");

    find_and_print_pathways(fp_md);

    fclose(fp_md);
}

/**
 * @brief Merges molecules by element count ranges for simplified analysis
 * @param target_element Element symbol to group by (e.g., "C", "H", "O")
 * @param ranges Vector of range boundaries for grouping (e.g., {0, 5, 10}
 * creates groups 0-4, 5-9, 10+)
 * @note Creates grouped molecules with formula prefix "grp_" followed by
 * element and range
 * @note Molecules with target element count in each range are merged into a
 * single node
 * @note Useful for simplifying complex reaction networks by grouping similar
 * molecules
 */
void ReaxFlow::merge_by_element(std::string target_element, std::vector<int> ranges) {
    std::unordered_set<std::string> all_formulas;
    std::unordered_map<std::string, std::unordered_set<std::string>> formulas_map;
    std::map<std::string, int> elements_weights;

    for (const auto& node : nodes) {
        all_formulas.insert(node->molecule->formula);
    }

    for (size_t i = 0; i < ranges.size(); i++) {
        int start;
        int end;
        std::string new_formula;

        if (i < ranges.size() - 1) {
            start = ranges[i];
            end = ranges[i + 1] - 1;
            new_formula = fmt::format("grp_{}{}-{}", target_element, start, end);
        }
        else {
            start = ranges[i];
            end = 100000;  // If a molecule have > 10000 atoms, that's user's bad input,
            // will not provide anything makes sense from the beginning.
            new_formula = fmt::format("grp_{}{}-max", target_element, start);
        }

        formulas_map[new_formula] = {};

        for (const auto& formula : all_formulas) {
            if (starts_with("grp_", formula)) {
                continue;
            }

            elements_weights = parse_formula(formula);
            for (const auto& [elem, weight] : elements_weights) {
                if (target_element == elem && weight >= start && weight <= end) {
                    formulas_map[new_formula].insert(formula);
                }
            }
        }
    }

    for (const auto& [new_formula, old_formulas_set] : formulas_map) {
        merge_formulas(old_formulas_set, new_formula);
    }
}

/**
 * @brief Merges multiple molecules into a single representative node
 * @param formulas_set Set of formula strings to merge
 * @param new_formula New formula string for the merged node
 * @note Combines all nodes with formulas in formulas_set into a single node
 * @note Updates all edges to connect to the new merged node
 * @note Preserves reaction counts and atom transfer data
 * @note Deallocates old nodes and edges to prevent memory leaks
 */
void ReaxFlow::merge_formulas(const std::unordered_set<std::string>& formulas_set, const std::string& new_formula) {
    // Get target nodes to merge
    std::vector<Node*> nodes_to_merge;
    for (const auto& node : nodes) {
        if (formulas_set.count(node->molecule->formula)) {
            nodes_to_merge.push_back(node);
        }
    }

    // Nothing to merge
    if (nodes_to_merge.size() < 1) return;

    // Create a new merged node
    Node* merged_node = new Node(nodes_to_merge[0]->molecule);
    merged_node->molecule->formula = new_formula;
    merged_node->hash = std::hash<std::string>()(new_formula);

    // Add the merged node
    nodes.insert(merged_node);
    molecule_hash_to_node[merged_node->hash] = merged_node;

    // Collect all edges involving the nodes to be merged
    std::vector<Edge*> edges_to_remove;
    std::vector<Edge*> edges_to_add;

    for (const auto& node_to_merge : nodes_to_merge) {
        // Find all edges involving this node
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            Edge* edge = *it;

            if (edge->source == node_to_merge || edge->target == node_to_merge) {
                edges_to_remove.push_back(edge);

                // Create new edge with merged node
                Node* other_node = (edge->source == node_to_merge) ? edge->target : edge->source;
                Edge* new_edge = new Edge((edge->source == node_to_merge) ? merged_node : other_node,
                    (edge->target == node_to_merge) ? merged_node : other_node);
                new_edge->count = edge->count;
                new_edge->atom_transfer = edge->atom_transfer;
                edges_to_add.push_back(new_edge);
            }
        }

        // Remove the node
        nodes.erase(node_to_merge);
        molecule_hash_to_node.erase(node_to_merge->hash);
        delete node_to_merge;
    }

    // Remove old edges
    for (Edge* edge : edges_to_remove) {
        edges.erase(edges.find(edge));
        delete edge;
    }

    // Add new edges
    for (Edge* edge : edges_to_add) {
        edges.insert(edge);
    }
}

#ifndef WASM_MODE
/**
 * @brief Dumps SMILES representations of all molecules to a CSV file
 * @note Only processes molecules that don't start with "grp_" prefix
 * @note Creates molecules_smiles.csv with formula and SMILES columns
 * @note Requires RDKit library (only available when not in WASM mode)
 */
#ifndef WASM_MODE
 /**
  * @brief Dumps SMILES representations of all molecules to a CSV file
  * @note Only processes molecules that don't start with "grp_" prefix
  * @note Creates molecules_smiles.csv with formula and SMILES columns
  * @note Sorted by node degree (highest first)
  * @note Requires RDKit library (only available when not in WASM mode)
  */
void ReaxFlow::dump_smiles() {
    // Ensure degrees are calculated
    update_graph();

    // Sort nodes by degree (descending)
    std::vector<Node*> sorted_nodes_vec(nodes.begin(), nodes.end());
    std::sort(sorted_nodes_vec.begin(), sorted_nodes_vec.end(),
        [](Node* a, Node* b) { return a->degree > b->degree; });

    FILE* fp = create_file("molecules_smiles.csv");
    // Header with precursor/derivative naming
    fmt::print(fp, "hash,formula,smiles,total_connections,precursors,derivatives\n");

    for (const auto& node : sorted_nodes_vec) {
        if (starts_with(node->molecule->formula, "grp_")) continue;
        try {
            fmt::print(fp, "{},{},{},{},{},{}\n",
                node->hash,
                node->molecule->formula,
                rdkit_smiles(*node->molecule),
                node->degree,
                node->precursor_count,
                node->derivative_count);
        }
        catch (const std::exception& e) {
            fmt::print(fp, "Warning: SMILES of {} can not be solved.\n", node->molecule->formula);
        }
    }

    fclose(fp);
}
#endif
#endif

#ifndef WASM_MODE
/**
 * @brief Draws molecular structures for all molecules in the reaction network
 * @note Only processes molecules that don't start with "grp_" prefix
 * @note Uses RDKit to generate molecular structure images
 * @note Requires RDKit library (only available when not in WASM mode)
 * @note Creates individual image files for each molecule
 */
void ReaxFlow::draw_molecules() {
    for (const auto& node : nodes) {
        try {
            if (starts_with(node->molecule->formula, "grp_")) continue;
            rdkit_draw_molecule(*node->molecule);
        }
        catch (...) {
            fmt::print("Note: Molecule {} can not be drawn.\n", node->molecule->formula);
        }
    }
}
#endif

/**
 * @brief Imports molecules into the reaction network
 * @param initial_or_final True if importing initial molecules, false if final
 * @param mol_hashes Set of molecule hashes to import
 * @note Caches the hash values and delays actual node retrieval until
 * identify_candidates()
 * @note Only valid nodes are added to reactant_candidates or product_candidates
 */
void ReaxFlow::import_molecules(bool initial_or_final, const std::unordered_set<unsigned int>& mol_hashes) {
    // In this function just cache the hash,delay get_node() until the graph work
    // is done.

    if (initial_or_final) {
        reactant_candidates_hash = mol_hashes;
    }
    else {
        product_candidates_hash = mol_hashes;
    }
}

/**
 * @brief Identifies and selects candidate reactants and products
 * @note Sorts nodes and edges, then identifies valid candidates
 * @note Only valid nodes are added to reactant_candidates or product_candidates
 */
void ReaxFlow::identify_candidates() {
    update_graph();
    Node* node_to_insert = nullptr;
    for (const auto& hash : reactant_candidates_hash) {
        node_to_insert = get_node(hash);
        if (node_to_insert != nullptr) {
            reactant_candidates.insert(node_to_insert);
        }
    }

    for (const auto& hash : product_candidates_hash) {
        node_to_insert = get_node(hash);
        if (node_to_insert != nullptr) {
            product_candidates.insert(node_to_insert);
        }
    }

    /********** decide good candidates **********/
    float balance_factor = 0.0f;
    float reactant_threshold = 1.0f;
    float product_threshold = 1.0f;

    // use iterator to safe remove
    for (auto it = reactant_candidates.begin(); it != reactant_candidates.end();) {
        balance_factor = float((*it)->precursor_count) / float((*it)->derivative_count);
        if (balance_factor > reactant_threshold) {
            it = reactant_candidates.erase(it);
        }
        else {
            ++it;
        }
    }

    for (auto it = product_candidates.begin(); it != product_candidates.end();) {
        balance_factor = float((*it)->precursor_count) / float((*it)->derivative_count);
        if (balance_factor < product_threshold) {
            it = product_candidates.erase(it);
        }
        else {
            ++it;
        }
    }

    fmt::print("\n=== Network Flow Analysis ===\n");
    fmt::print("{} reactant candidates\n", reactant_candidates.size());
    for (const auto& node : reactant_candidates) {
        fmt::print("{} ", node->molecule->formula);
    }
    fmt::print("\n");

    fmt::print("{} product candidates\n", product_candidates.size());
    for (const auto& node : product_candidates) {
        fmt::print("{} ", node->molecule->formula);
    }
    fmt::print("\n");
}

/**
 * @brief Solves the network flow problem
 * @note Implements a network flow algorithm
 * @note Uses a graph representation for the network
 * @note Applies tarjan's algorithm to find strongly connected components and
 * calculate flow.
 */
void ReaxFlow::network_flow_solve() {
    /*
    algorithm tarjan is
        input: graph G = (V, E)
        output: set of strongly connected components (sets of vertices)

        index := 0
        S := empty stack
        for each v in V do
            if v.index is undefined then
                strongconnect(v)

        function strongconnect(v)
            // Set the depth index for v to the smallest unused index
            v.index := index
            v.lowlink := index
            index := index + 1
            S.push(v)
            v.onStack := true

            // Consider successors of v
            for each (v, w) in E do
                if w.index is undefined then
                    // Successor w has not yet been visited; recurse on it
                    strongconnect(w)
                    v.lowlink := min(v.lowlink, w.lowlink)
                else if w.onStack then
                    // Successor w is in stack S and hence in the current SCC
                    // If w is not on stack, then (v, w) is an edge pointing to an
    SCC already found and must be ignored
                    // See below regarding the next line
                    v.lowlink := min(v.lowlink, w.index)

            // If v is a root node, pop the stack and generate an SCC
            if v.lowlink = v.index then
                start a new strongly connected component
                repeat
                    w := S.pop()
                    w.onStack := false
                    add w to current strongly connected component
                while w ≠ v
                output the current strongly connected component
    */

    // Step 1: BFS to find all possible paths from reactant_candidates to
    // product_candidates

    struct Path {
        Node* reactant;
        Node* product;
        int significance = 0;
        unsigned int signature = 0;
        std::vector<Node*> nodes;
        int max_flow = 0;

        /**
         * @brief Checks if there is a path from reactant to product
         * @return true if there is a path, false otherwise
         * @note complexity: O(V + E)
         */
        bool connected() {
            std::queue<Node*> queue;
            std::unordered_set<Node*> visited;
            Node* current_node = reactant;

            queue.push(reactant);
            visited.insert(reactant);

            while (!queue.empty()) {
                current_node = queue.front();
                queue.pop();

                if (current_node == product) {
                    return true;
                }

                for (const auto& next_node : current_node->to_nodes) {
                    if (visited.find(next_node) == visited.end()) {
                        visited.insert(next_node);
                        queue.push(next_node);
                    }
                }
            }

            return false;
        }

        /**
         * @brief Finds the maximum flow using Dinic's algorithm
         * @return The maximum flow
         * @note complexity: O(V^2 * E)
         */
        int Dinic_max_flow() {
            // TODO: Implement Dinic's algorithm

            max_flow = 0;
            return max_flow;
        }
    };

    // Decide important reactants and products
    // std::vector<Node*> reactant_candidates_copy(reactant_candidates.begin(),
    // reactant_candidates.end()); std::vector<Node*>
    // product_candidates_copy(product_candidates.begin(),
    // product_candidates.end());

    // std::sort(reactant_candidates_copy.begin(), reactant_candidates_copy.end(),
    //           [](const Node* a, const Node* b) { return a->out_degree >
    //           b->out_degree; });
    // std::sort(product_candidates_copy.begin(), product_candidates_copy.end(),
    //           [](const Node* a, const Node* b) { return a->in_degree >
    //           b->in_degree; });

    // if (reactant_candidates_copy.size() > NETWORK_FLOW_MAX_REACTANTS) {
    //     reactant_candidates_copy.resize(NETWORK_FLOW_MAX_REACTANTS);
    // }
    // if (product_candidates_copy.size() > NETWORK_FLOW_MAX_PRODUCTS) {
    //     product_candidates_copy.resize(NETWORK_FLOW_MAX_PRODUCTS);
    // }

    // Find all possible paths
    std::vector<Path> all_paths;
    std::vector<Path> valid_paths;

    for (const auto& reactant : reactant_candidates) {
        for (const auto& product : product_candidates) {
            if (reactant->derivative_count == 0 || product->precursor_count == 0) {
                continue;
            }
            if (reactant->hash == product->hash) {
                continue;
            }

            Path path;
            path.reactant = reactant;
            path.product = product;
            path.significance = reactant->derivative_count * product->precursor_count;
            all_paths.push_back(path);
        }
    }

    // complexity: O(NETWORK_FLOW_MAX_REACTANTS * NETWORK_FLOW_MAX_PRODUCTS * (V +
    // E))
    for (auto& path : all_paths) {
        bool connected = path.connected();
        if (connected) {
            valid_paths.push_back(path);
        }
    }
    all_paths.clear();

    // std::sort(valid_paths.begin(), valid_paths.end(),
    //           [](const Path& a, const Path& b) { return a.significance >
    //           b.significance; });

    // for (const auto& path : valid_paths) {
    //     fmt::print("{} -> {} (significance: {})\n",
    //     path.reactant->molecule->formula, path.product->molecule->formula,
    //                path.significance);
    // }

    for (auto& path : valid_paths) {
        path.Dinic_max_flow();
    }

    // TODO: Step 2: Solve the max flow problem using the selected pairs as source
    // and sink, Dinic's algorithm
    // TODO: Step 3: Output results, like P -> T1 -> T2 -> ... -> Tn -> P (70% ->
    // 50% -> 30% -> 20% -> ...)
}

// ============================================================================
// Reaction Network Improvement Methods
// ============================================================================

/**
 * @brief Calculate atom economy for all reactions
 * @note Atom economy = atom_transfer / min(source_atoms, target_atoms)
 */
void ReaxFlow::calculate_atom_economy() {
    for (auto& edge : edges) {
        int src_atoms = edge->source->molecule->mol_atoms.size();
        int tgt_atoms = edge->target->molecule->mol_atoms.size();
        int min_atoms = std::min(src_atoms, tgt_atoms);

        if (min_atoms > 0) {
            edge->atom_economy = (double)edge->atom_transfer / min_atoms;
        }
        else {
            edge->atom_economy = 0.0;
        }
    }
}

/**
 * @brief Filter reactions by atom economy threshold
 * @param threshold Minimum atom economy ratio (default 0.3)
 * @note Preserves dissociation and association reactions with lower threshold
 */
void ReaxFlow::filter_by_atom_economy(double threshold) {
    calculate_atom_economy();

    std::vector<Edge*> edges_to_remove;
    int preserved_dissociation = 0;
    int preserved_association = 0;

    for (auto& edge : edges) {
        // Keep high atom economy reactions
        if (edge->atom_economy >= threshold) {
            continue;
        }

        // Special handling for dissociation/association
        int src_size = edge->source->molecule->mol_atoms.size();
        int tgt_size = edge->target->molecule->mol_atoms.size();

        bool is_dissociation = (src_size > tgt_size * 2);
        bool is_association = (tgt_size > src_size * 2);

        // Dissociation: large -> small (e.g., C10H20 -> C5H10 + C5H10)
        if (is_dissociation && edge->atom_economy > 0.15) {
            preserved_dissociation++;
            continue;
        }

        // Association: small -> large (e.g., 2CH3 -> C2H6)
        if (is_association && edge->atom_economy > 0.15) {
            preserved_association++;
            continue;
        }

        edges_to_remove.push_back(edge);
    }

    // Remove marked edges
    for (auto& edge : edges_to_remove) {
        edges.erase(edge);
        edge_hash_to_edge.erase(edge->hash);
        delete edge;
    }

    // Always print statistics for evaluation
    fmt::print("\n=== Atom Economy Filter ===\n");
    fmt::print("Total reactions: {}\n", edges.size() + edges_to_remove.size());
    fmt::print("Threshold: {:.2f}\n", threshold);
    fmt::print("Removed: {} reactions ({:.1f}%)\n",
        edges_to_remove.size(),
        100.0 * edges_to_remove.size() / (edges.size() + edges_to_remove.size()));
    fmt::print("Preserved dissociations (AE>0.15): {}\n", preserved_dissociation);
    fmt::print("Preserved associations (AE>0.15): {}\n", preserved_association);

    // Calculate atom economy distribution
    int high_ae = 0, mid_ae = 0, low_ae = 0;
    for (auto& edge : edges) {
        if (edge->atom_economy >= 0.5) high_ae++;
        else if (edge->atom_economy >= 0.3) mid_ae++;
        else low_ae++;
    }
    fmt::print("Distribution of kept reactions:\n");
    fmt::print("  High AE (>=0.5): {}\n", high_ae);
    fmt::print("  Mid AE (0.3-0.5): {}\n", mid_ae);
    fmt::print("  Low AE (<0.3, special): {}\n", low_ae);
}

/**
 * @brief Calculate edge betweenness centrality using Brandes algorithm
 * @note Uses sampling for large networks (max 100 source nodes)
 */
void ReaxFlow::calculate_edge_betweenness() {
    // Initialize
    for (auto& edge : edges) {
        edge->betweenness = 0.0;
    }

    if (nodes.empty()) return;

    // Convert to vector for sampling
    std::vector<Node*> node_list(nodes.begin(), nodes.end());

    // Sample at most 100 source nodes for efficiency
    int sample_size = std::min((int)node_list.size(), 100);
    std::shuffle(node_list.begin(), node_list.end(), std::mt19937(42));

    for (int i = 0; i < sample_size; i++) {
        Node* s = node_list[i];

        // BFS data structures
        std::stack<Node*> S;
        std::queue<Node*> Q;
        std::unordered_map<Node*, std::vector<Node*>> pred;
        std::unordered_map<Node*, int> sigma;
        std::unordered_map<Node*, int> dist;

        for (auto& w : nodes) {
            sigma[w] = 0;
            dist[w] = -1;
        }
        sigma[s] = 1;
        dist[s] = 0;
        Q.push(s);

        // BFS
        while (!Q.empty()) {
            Node* v = Q.front(); Q.pop();
            S.push(v);

            for (auto& w : v->to_nodes) {
                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    Q.push(w);
                }
                if (dist[w] == dist[v] + 1) {
                    sigma[w] += sigma[v];
                    pred[w].push_back(v);
                }
            }
        }

        // Back-propagation
        std::unordered_map<Node*, double> delta;
        for (auto& v : nodes) delta[v] = 0.0;

        while (!S.empty()) {
            Node* w = S.top(); S.pop();

            for (auto& v : pred[w]) {
                Edge* e = get_edge(v, w);
                if (!e) continue;

                double c = ((double)sigma[v] / sigma[w]) * (1.0 + delta[w]);
                e->betweenness += c;
                delta[v] += c;
            }
        }
    }

    // Normalize for sampling
    double norm_factor = (double)nodes.size() / sample_size;
    for (auto& edge : edges) {
        edge->betweenness *= norm_factor;
    }
}

/**
 * @brief Filter reactions by betweenness centrality
 * @param target_edge_count Target number of edges to keep
 * @note Combines betweenness, count, and atom_transfer criteria
 */
void ReaxFlow::filter_by_betweenness(int target_edge_count) {
    calculate_edge_betweenness();

    // Create set of edges to keep
    std::unordered_set<Edge*> edges_to_keep;

    // Sort by different criteria
    std::vector<Edge*> by_bc(edges.begin(), edges.end());
    std::sort(by_bc.begin(), by_bc.end(),
        [](Edge* a, Edge* b) { return a->betweenness > b->betweenness; });

    std::vector<Edge*> by_count(edges.begin(), edges.end());
    std::sort(by_count.begin(), by_count.end(),
        [](Edge* a, Edge* b) { return a->count > b->count; });

    std::vector<Edge*> by_transfer(edges.begin(), edges.end());
    std::sort(by_transfer.begin(), by_transfer.end(),
        [](Edge* a, Edge* b) { return a->atom_transfer > b->atom_transfer; });

    // Keep top from each category
    int n_bc = target_edge_count / 2;
    int n_count = target_edge_count / 3;
    int n_transfer = target_edge_count - n_bc - n_count;

    for (int i = 0; i < std::min(n_bc, (int)by_bc.size()); i++) {
        edges_to_keep.insert(by_bc[i]);
    }
    for (int i = 0; i < std::min(n_count, (int)by_count.size()); i++) {
        edges_to_keep.insert(by_count[i]);
    }
    for (int i = 0; i < std::min(n_transfer, (int)by_transfer.size()); i++) {
        edges_to_keep.insert(by_transfer[i]);
    }

    // Remove edges not in keep set
    std::vector<Edge*> edges_to_remove;
    for (auto& edge : edges) {
        if (edges_to_keep.find(edge) == edges_to_keep.end()) {
            edges_to_remove.push_back(edge);
        }
    }

    for (auto& edge : edges_to_remove) {
        edges.erase(edge);
        edge_hash_to_edge.erase(edge->hash);
        delete edge;
    }

    fmt::print("\n=== Edge Betweenness Filter ===\n");
    fmt::print("Target edges: {}\n", target_edge_count);
    fmt::print("Total before: {}\n", edges.size() + edges_to_remove.size());
    fmt::print("Kept: {} (BC top {} + count top {} + transfer top {})\n",
        edges_to_keep.size(), n_bc, n_count, n_transfer);
    fmt::print("Removed: {}\n", edges_to_remove.size());

    // Print top 5 highest betweenness edges
    fmt::print("\nTop 5 by betweenness centrality:\n");
    for (int i = 0; i < std::min(5, (int)by_bc.size()); i++) {
        fmt::print("  {} -> {} (BC={:.2f}, count={})\n",
            by_bc[i]->source->molecule->formula,
            by_bc[i]->target->molecule->formula,
            by_bc[i]->betweenness,
            by_bc[i]->count);
    }
}


/**
 * @brief Clean up isolated nodes (nodes with degree = 0)
 * @note Called after filtering to remove nodes that are no longer connected
 */
void ReaxFlow::cleanup_isolated_nodes() {
    std::vector<Node*> nodes_to_remove;

    for (auto& node : nodes) {
        // Check if node has any connected edges
        bool has_incoming = false;
        bool has_outgoing = false;

        for (auto& edge : edges) {
            if (edge->target == node) {
                has_incoming = true;
            }
            if (edge->source == node) {
                has_outgoing = true;
            }
            if (has_incoming && has_outgoing) break;
        }

        // Remove node with no connections
        if (!has_incoming && !has_outgoing) {
            nodes_to_remove.push_back(node);
        }
    }

    for (auto& node : nodes_to_remove) {
        nodes.erase(node);
        molecule_hash_to_node.erase(node->hash);
        delete node;
    }

    if (!nodes_to_remove.empty()) {
        fmt::print("Removed {} isolated nodes\n", nodes_to_remove.size());
    }
}


/**
 * @brief Find and print the longest reaction pathways in the network
 * @param fp File pointer for output
 * @note Uses DFS to find paths from sources (no precursors) to targets (no derivatives)
 * @note Handles cycles by tracking visited nodes
 */
void ReaxFlow::find_and_print_pathways(FILE* fp) {
    // Find source nodes (no precursors) and target nodes (no derivatives)
    std::vector<Node*> sources;
    std::vector<Node*> targets;

    for (auto& node : nodes) {
        if (node->precursor_count == 0 && node->derivative_count > 0) {
            sources.push_back(node);
        }
        if (node->derivative_count == 0 && node->precursor_count > 0) {
            targets.push_back(node);
        }
    }

    if (sources.empty() || targets.empty()) {
        fmt::print(fp, "No clear reaction pathways found (no distinct sources or targets).\n\n");
        return;
    }

    fmt::print(fp, "Found {} source species and {} target species.\n\n", sources.size(), targets.size());

    // Find paths from each source to each target
    std::vector<std::vector<Node*>> all_paths;

    for (auto* source : sources) {
        for (auto* target : targets) {
            if (source == target) continue;

            std::vector<Node*> path;
            std::unordered_set<Node*> visited;
            dfs_pathways(source, target, path, visited, all_paths, 10);
        }
    }

    // Sort paths by length (longest first)
    std::sort(all_paths.begin(), all_paths.end(),
        [](const std::vector<Node*>& a, const std::vector<Node*>& b) {
            return a.size() > b.size();
        });

    // Remove duplicate paths (same start and end)
    std::vector<std::vector<Node*>> unique_paths;
    std::set<std::pair<std::string, std::string>> seen_endpoints;

    for (auto& path : all_paths) {
        if (path.size() < 2) continue;

        std::string start = path.front()->molecule->formula;
        std::string end = path.back()->molecule->formula;
        auto key = std::make_pair(start, end);

        if (seen_endpoints.find(key) == seen_endpoints.end()) {
            seen_endpoints.insert(key);
            unique_paths.push_back(path);
        }
    }

    // Print top pathways
    fmt::print(fp, "## Top Reaction Pathways\n\n");

    int pathway_count = 0;
    for (auto& path : unique_paths) {
        if (++pathway_count > 10) break;  // Limit to top 10

        fmt::print(fp, "### Pathway {} ({} steps)\n\n", pathway_count, path.size() - 1);

        // Print path
        fmt::print(fp, "**Path**: ");
        for (size_t i = 0; i < path.size(); i++) {
            fmt::print(fp, "{}", path[i]->molecule->formula);
            if (i < path.size() - 1) {
                fmt::print(fp, " → ");
            }
        }
        fmt::print(fp, "\n\n");

        // Print details
        fmt::print(fp, "**Details**:\n\n");
        for (size_t i = 0; i < path.size() - 1; i++) {
            Node* from = path[i];
            Node* to = path[i + 1];

            // Find the edge
            Edge* edge = get_edge(from, to);
            if (edge) {
                fmt::print(fp, "- {} → {} (count: {})\n",
                    from->molecule->formula,
                    to->molecule->formula,
                    edge->count);
            }
            else {
                fmt::print(fp, "- {} → {}\n", from->molecule->formula, to->molecule->formula);
            }
        }
        fmt::print(fp, "\n---\n\n");
    }
}

/**
 * @brief DFS helper for finding pathways
 * @param current Current node
 * @param target Target node to reach
 * @param path Current path being explored
 * @param visited Set of visited nodes to avoid cycles
 * @param all_paths Vector to store all found paths
 * @param max_depth Maximum search depth to prevent infinite recursion
 */
void ReaxFlow::dfs_pathways(Node* current, Node* target, std::vector<Node*>& path,
    std::unordered_set<Node*>& visited, std::vector<std::vector<Node*>>& all_paths,
    int max_depth) {
    // Add current node to path
    path.push_back(current);
    visited.insert(current);

    // Check if we reached the target
    if (current == target) {
        all_paths.push_back(path);
        // Backtrack
        path.pop_back();
        visited.erase(current);
        return;
    }

    // Stop if max depth reached
    if (path.size() >= (size_t)max_depth) {
        path.pop_back();
        visited.erase(current);
        return;
    }

    // Explore derivatives (outgoing edges)
    for (auto* next_node : current->to_nodes) {
        if (visited.find(next_node) == visited.end()) {
            dfs_pathways(next_node, target, path, visited, all_paths, max_depth);
        }
    }

    // Backtrack
    path.pop_back();
    visited.erase(current);
}
