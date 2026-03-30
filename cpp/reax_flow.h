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

    // Topological degrees (edge counts) - renamed for chemical intuition
    int degree = 0;             // total connections (precursor_count + derivative_count)
    int precursor_count = 0;    // in_degree: number of precursors (incoming edges)
    int derivative_count = 0;   // out_degree: number of derivatives (outgoing edges)

    // Weighted degrees (reaction counts)
    int reaction_count = 0;         // total reaction count
    int precursor_reactions = 0;    // reactions as product (incoming)
    int derivative_reactions = 0;   // reactions as reactant (outgoing)

    // Atom transfer counts
    int atom_transfer = 0;          // total atom transfer
    int precursor_atom_transfer = 0;    // incoming atom transfer
    int derivative_atom_transfer = 0;   // outgoing atom transfer

    std::unordered_set<Node*> from_nodes;  // precursors
    std::unordered_set<Node*> to_nodes;    // derivatives

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
    void cleanup_isolated_nodes();
    
    // Pathway analysis
    void find_and_print_pathways(FILE* fp);
    void dfs_pathways(Node* current, Node* target, std::vector<Node*>& path, 
                      std::unordered_set<Node*>& visited, std::vector<std::vector<Node*>>& all_paths,
                      int max_depth = 10);
};
