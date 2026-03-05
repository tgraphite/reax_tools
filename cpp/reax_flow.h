#pragma once
#include <fstream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "argparser.h"
#include "system.h"

struct Molecule;
struct Reaction;
struct Node;

class ReaxFlowException : public std::runtime_error {
public:
    /**
     * @brief Construct a new ReaxFlowException with a message.
     * @param msg The error message.
     */
    explicit ReaxFlowException(const std::string& msg) : std::runtime_error(msg) {}
};

struct Node {
    Molecule* molecule = nullptr;

    unsigned int hash = 0;

    int degree = 0;
    int in_degree = 0;
    int out_degree = 0;

    int degree_at = 0;
    int in_degree_at = 0;
    int out_degree_at = 0;

    std::unordered_set<Node*> from_nodes;
    std::unordered_set<Node*> to_nodes;

    Node(Molecule* mol);
    ~Node();

    void add_degrees(bool source_or_target, unsigned int count, unsigned int atom_transfer_count);
};

struct Edge {
    Node* source = nullptr;
    Node* target = nullptr;

    unsigned int hash = 0;
    int count = 0;
    int atom_transfer = 0;
    
    // New fields for reaction importance analysis
    double atom_economy = 0.0;  // Atom economy ratio
    double betweenness = 0.0;   // Edge betweenness centrality

    Edge(Node* from_node, Node* to_node);
    ~Edge();
};

inline unsigned int get_edge_hash(const Edge* edge) {
    return std::hash<unsigned int>()(edge->source->hash) ^ (std::hash<unsigned int>()(edge->target->hash) << 1);
}

inline unsigned int get_edge_hash(Node* from_node, Node* to_node) {
    return std::hash<unsigned int>()(from_node->hash) ^ (std::hash<unsigned int>()(to_node->hash) << 1);
}

struct node_hash {
    std::size_t operator()(const Node* p) const noexcept { return p->hash; }
};
struct edge_hash {
    std::size_t operator()(const Edge* p) const noexcept { return p->hash; }
};

struct node_equal {
    bool operator()(const Node* lhs, const Node* rhs) const noexcept { return lhs->hash == rhs->hash; }
};
struct edge_equal {
    bool operator()(const Edge* lhs, const Edge* rhs) const noexcept { return lhs->hash == rhs->hash; }
};

inline bool operator==(const Node& lhs, const Node& rhs) { return { lhs.hash == rhs.hash }; }
inline bool operator==(const Edge& lhs, const Edge& rhs) { return { lhs.hash == rhs.hash }; }

class ReaxFlow {
private:
    std::mutex reaxflow_mutex;

    std::unordered_set<Node*, node_hash, node_equal> nodes;
    std::unordered_set<Edge*, edge_hash, edge_equal> edges;

    // Fast lookup map from molecule hash to node pointer
    std::unordered_map<unsigned int, Node*> molecule_hash_to_node;
    std::unordered_map<unsigned int, Edge*> edge_hash_to_edge;

    std::vector<std::pair<Node*, unsigned int>> sorted_nodes;
    std::vector<std::pair<Edge*, unsigned int>> sorted_edges;

    std::unordered_set<unsigned int> reactant_candidates_hash;
    std::unordered_set<unsigned int> product_candidates_hash;
    std::unordered_set<Node*> reactant_candidates;
    std::unordered_set<Node*> product_candidates;

public:
    ReaxFlow();
    ~ReaxFlow();

    Node* get_node(Molecule* mol);
    Node* get_node(unsigned int hash);
    Edge* get_edge(Node* from_node, Node* to_node);
    Edge* get_edge(unsigned int hash);

    Node* add_molecule(Molecule* mol);
    Edge* add_reaction(const int& frame_id, const int& atom_transfer_count, Molecule* source, Molecule* target);

    void update_graph();
    void brief_report();

    void merge_formulas(const std::unordered_set<std::string>& formulas, const std::string& new_formula);
    void merge_by_element(std::string target_element, std::vector<int> ranges);

    void import_molecules(bool initial_or_final, const std::unordered_set<unsigned int>& mol_hashes);
    void identify_candidates();
    /**
     * @brief Solves the network flow problem using a multi-step approach
     * @note Step 1: BFS to find all possible paths from reactant_candidates to product_candidates
     * @note Step 2: Sort <reactant, product> pairs by min(reactant.out_degree, product.in_degree)
     * @note Step 3: Select top N pairs
     * @note Step 4: Solve the max flow problem using the selected pairs as source and sink, Dinic's algorithm
     * @note Step 5: Output results, like P -> T1 -> T2 -> ... -> Tn -> P (70% -> 50% -> 30% -> 20% -> ...)
     */
    void network_flow_solve();

    void write_dot_file(std::string basename, const std::vector<Edge*>& edges_to_write,
        bool write_atom_transfer = true, std::string layout = "dot");
    void write_dot_file_significant_nodes(std::string basename, int max_nodes,
        bool write_atom_transfer = true, std::string layout = "dot");
    void save_graph();
    void reduce_graph();
    void save_molecule_centered_subgraphs(bool write_atom_transfer = true, bool csv_only = true, bool use_hash = false);
    void draw_molecules();
    void dump_smiles();
    
    // New methods for reaction network improvement
    void calculate_atom_economy();
    void filter_by_atom_economy(double threshold = 0.3);
    void calculate_edge_betweenness();
    void filter_by_betweenness(int target_edge_count = 60);
};
