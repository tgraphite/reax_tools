#pragma once
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "defines.h"
#include "molecule.h"

class ReaxFlow {
   private:
    // Node: molecule species
    // Edge: reaction
    struct Node {
        std::string formula;
        std::string smiles;
        std::vector<int> frames_appeared;
    };

    struct Edge {
        int source_node_id;
        int target_node_id;
        int reaction_count;
        std::vector<int> frames_occurred;
    };

    std::vector<Node*> nodes;
    std::vector<Edge*> edges;

    // Get node id from formula
    std::map<std::string, int> formula_to_node_id;

   public:
    void add_reaction(int frame, Molecule* source, Molecule* target);
    void brief_report();
    void save_graph(const std::string& raw_file_path, int& max_reactions);
};
