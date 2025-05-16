#pragma once
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#include "defines.h"
#include "molecule.h"

class ReaxFlow {
   private:
    // Node: molecule species
    // Edge: reaction

    std::vector<Molecule*> nodes;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> edge_reaction_counts;

    std::mutex reaxflow_mutex;

   public:
    // return id of the node
    int add_molecule(Molecule* mol);
    std::pair<int, int> add_reaction(int frame, Molecule* source, Molecule* target);

    std::vector<int> get_neighbors(int node_id);

    void brief_report();
    void save_graph(const std::string& output_dir, int& max_molecules, bool draw_molecules = false,
                    bool reduce_reactions = false);
    void reduce_graph();

    Molecule* get_node_from_id(int id);

    void remove_zero_nodes();
    void remove_zero_edges();
    void remove_isolated_nodes();
};
