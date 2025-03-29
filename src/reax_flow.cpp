#include "reax_flow.h"
#include "universe.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "fmt/format.h"

// 添加反应关系
void ReaxFlow::add_reaction(int frame, const std::shared_ptr<Molecule>& source, const std::shared_ptr<Molecule>& target) {
    // 获取或创建源节点
    int source_id = -1;
    if (formula_to_node_id.find(source->formula) != formula_to_node_id.end()) {
        source_id = formula_to_node_id[source->formula];
        // 更新节点出现的帧
        if (std::find(nodes[source_id].frames_appeared.begin(), 
                     nodes[source_id].frames_appeared.end(), 
                     frame-1) == nodes[source_id].frames_appeared.end()) {
            nodes[source_id].frames_appeared.push_back(frame-1);
        }
    } else {
        // 创建新节点
        Node source_node;
        source_node.formula = source->formula;
        // TODO: 生成SMILES字符串
        source_node.smiles = ""; 
        source_node.frames_appeared.push_back(frame-1);
        
        nodes.push_back(source_node);
        source_id = nodes.size() - 1;
        formula_to_node_id[source->formula] = source_id;
    }
    
    // 获取或创建目标节点
    int target_id = -1;
    if (formula_to_node_id.find(target->formula) != formula_to_node_id.end()) {
        target_id = formula_to_node_id[target->formula];
        // 更新节点出现的帧
        if (std::find(nodes[target_id].frames_appeared.begin(), 
                     nodes[target_id].frames_appeared.end(), 
                     frame) == nodes[target_id].frames_appeared.end()) {
            nodes[target_id].frames_appeared.push_back(frame);
        }
    } else {
        // 创建新节点
        Node target_node;
        target_node.formula = target->formula;
        // TODO: 生成SMILES字符串
        target_node.smiles = ""; 
        target_node.frames_appeared.push_back(frame);
        
        nodes.push_back(target_node);
        target_id = nodes.size() - 1;
        formula_to_node_id[target->formula] = target_id;
    }
    
    // 检查是否已存在相同的边
    bool edge_exists = false;
    for (auto& edge : edges) {
        if (edge.source_node_id == source_id && edge.target_node_id == target_id) {
            // 边已存在，增加反应计数
            edge.reaction_count++;
            edge.frames_occurred.push_back(frame);
            edge_exists = true;
            break;
        }
    }
    
    // 如果边不存在，创建新边
    if (!edge_exists) {
        Edge new_edge;
        new_edge.source_node_id = source_id;
        new_edge.target_node_id = target_id;
        new_edge.reaction_count = 1;
        new_edge.frames_occurred.push_back(frame);
        edges.push_back(new_edge);
    }
}

// 生成反应流程报告
void ReaxFlow::brief_report() {
    std::cout << "=== Reaction Flow Report ===" << std::endl;
    std::cout << "Total nodes (species): " << nodes.size() << std::endl;
    std::cout << "Total edges (reactions): " << edges.size() << std::endl;
    
    // 按反应次数排序边
    std::vector<std::pair<int, int>> sorted_edges;
    for (size_t i = 0; i < edges.size(); i++) {
        sorted_edges.push_back({i, edges[i].reaction_count});
    }
    
    std::sort(sorted_edges.begin(), sorted_edges.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // 显示前10个最频繁的反应
    std::cout << "\nTop reactions:" << std::endl;
    int count = 0;
    for (const auto& [edge_idx, reaction_count] : sorted_edges) {
        if (count >= 10) break;
        
        const Edge& edge = edges[edge_idx];
        std::cout << fmt::format("{}: {} -> {} (count: {})", 
                                count+1,
                                nodes[edge.source_node_id].formula,
                                nodes[edge.target_node_id].formula,
                                reaction_count) << std::endl;
        count++;
    }
}

// 保存反应流程图为DOT格式
void ReaxFlow::save_graph(const std::string& raw_file_path) {
    std::string save_path = raw_file_path.substr(0, raw_file_path.find_last_of(".")) + ".dot";
    std::ofstream file(save_path);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << save_path << " for writing." << std::endl;
        return;
    }
    
    // 写入DOT文件头
    file << "digraph ReactionFlow {\n";
    file << "  rankdir=LR;\n";
    file << "  node [shape=box, style=filled, fillcolor=lightblue];\n\n";
    
    // 写入节点
    for (size_t i = 0; i < nodes.size(); i++) {
        file << "  node" << i << " [label=\"" << nodes[i].formula;
        if (!nodes[i].smiles.empty()) {
            file << "\\n" << nodes[i].smiles;
        }
        file << "\"];\n";
    }
    
    file << "\n";
    
    // 写入边
    for (const auto& edge : edges) {
        file << "  node" << edge.source_node_id << " -> node" << edge.target_node_id;
        file << " [label=\"" << edge.reaction_count << "\", penwidth=" << std::min(5.0, 1.0 + log(edge.reaction_count)) << "];\n";
    }
    
    file << "}\n";
    file.close();
    
    std::cout << "Reaction flow graph saved to " << save_path << std::endl;
}
