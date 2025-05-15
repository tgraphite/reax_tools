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
    struct Node {
        Molecule molecule;
        std::string formula;
        int reaction_count = 0;

        Node(const Molecule& mol) : molecule(mol) {}  // Copy constructor
    };

    struct Edge {
        Node* source_node;
        Node* target_node;
        int reaction_count = 0;
    };

    std::vector<Node*> nodes;
    std::vector<Edge*> edges;
    std::mutex reaxflow_mutex;

   public:
    void add_reaction(int frame, Molecule* source, Molecule* target);
    void brief_report();
    void save_graph(const std::string& raw_file_path, int& max_molecules, bool draw_molecules = false);
    void draw_molecules(const std::string& raw_file_path);
    void reduce_graph(int max_molecules);

    Node* get_node_from_molecule(Molecule* mol);
    Edge* get_edge_from_molecules(Molecule* source, Molecule* target);

    bool has_reaction(Molecule* source, Molecule* target);
};
