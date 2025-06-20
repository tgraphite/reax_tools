#pragma once
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "defines.h"
#include "molecule.h"

class ReaxFlow {
   private:
    // Node: point to a specific molecule
    // Edge: reaction

    std::vector<Molecule*> nodes;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> edge_reaction_counts;
    std::vector<int> edge_atom_transfer_counts;

    // Hash maps for O(1) lookup
    std::unordered_map<std::string, int> molecule_map;
    std::mutex reaxflow_mutex;

   public:
    // return id of the node
    int add_molecule(Molecule* mol);
    std::pair<int, int> add_reaction(const int& frame, const int& atom_transfer_count, Molecule* source,
                                     Molecule* target);

    void brief_report();
    void write_dot_file(const std::string& output_file, const std::vector<int>& edge_indices,
                        bool write_atom_transfer_count = false, std::string layout = "circo");
    void write_molecule_centered_csv_file(const std::string& output_file, const std::vector<int>& edge_indices,
                                          bool write_atom_transfer = true);

    void save_graph(const std::string& output_dir, int& max_molecules, bool draw_molecules = false,
                    bool reduce_reactions = false);
    void reduce_graph();

    void save_molecule_centered_subgraphs(const std::string& output_dir, bool write_atom_transfer = true);

    Molecule* get_node_from_id(int id);

    void remove_zero_nodes();
    void remove_zero_edges();
    void remove_isolated_nodes();
};
