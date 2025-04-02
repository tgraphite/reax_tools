#pragma once
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "defines.h"
#include "molecule.h"

class ReaxFlow {
   private:
    // 节点表示分子物种
    struct Node {
        std::string formula;
        std::string smiles;
        std::vector<int> frames_appeared;
    };

    // 边表示反应关系
    struct Edge {
        int source_node_id;
        int target_node_id;
        int reaction_count;
        std::vector<int> frames_occurred;
    };

    std::vector<Node> nodes;
    std::vector<Edge> edges;

    // 用于快速查找节点
    std::map<std::string, int> formula_to_node_id;

   public:
    void add_reaction(int frame, const std::shared_ptr<Molecule>& source,
                      const std::shared_ptr<Molecule>& target);
    void brief_report();
    void save_graph(const std::string& raw_file_path);
};
