# ReaxTools 反应网络算法改进方案

## 目录
1. [当前算法分析与问题诊断](#1-当前算法分析与问题诊断)
2. [理论基础与启发式搜索](#2-理论基础与启发式搜索)
3. [改进方案一：多层次反应重要性评估](#3-改进方案一多层次反应重要性评估)
4. [改进方案二：基于社区发现的网络简化](#4-改进方案二基于社区发现的网络简化)
5. [改进方案三：动态主干路径提取](#5-改进方案三动态主干路径提取)
6. [改进方案四：时空关联反应聚类](#6-改进方案四时空关联反应聚类)
7. [实现路线图与优先级](#7-实现路线图与优先级)
8. [参考理论框架](#8-参考理论框架)

---

## 1. 当前算法分析与问题诊断

### 1.1 当前反应检测机制

```cpp
// system.cpp:907-936
for (const auto& curr_atom : this->atoms) {
    prev_atom = prev_sys->get_atom_by_id(curr_atom->id);
    if (!prev_atom) continue;
    prev_mol = prev_atom->belong_molecule;
    curr_mol = curr_atom->belong_molecule;
    if (prev_mol == nullptr || curr_mol == nullptr) continue;
    if (prev_mol->hash != curr_mol->hash) {
        mol_id_pair = { prev_mol->id, curr_mol->id };
        if (reaction_mol_id_pairs.find(mol_id_pair) == reaction_mol_id_pairs.end()) {
            reaction_mol_id_pairs.insert(mol_id_pair);
        }
    }
}

for (const auto& pair : reaction_mol_id_pairs) {
    intersection.clear();
    prev_mol = prev_sys->get_molecule_by_id(pair.first);
    curr_mol = this->get_molecule_by_id(pair.second);
    std::set_intersection(prev_mol->atom_ids.begin(), prev_mol->atom_ids.end(), 
                          curr_mol->atom_ids.begin(), curr_mol->atom_ids.end(), 
                          back_inserter(intersection));
    if (intersection.size() > 0) {
        reax_flow->add_reaction(this->frame_id, intersection.size(), prev_mol, curr_mol);
    }
}
```

**核心逻辑**：
- 帧间对比：比较相邻帧中每个原子所属分子是否变化
- 反应定义：只要有原子重叠（`intersection.size() > 0`）就认为是反应
- 反应权重：记录反应次数(`count`)和重叠原子数(`atom_transfer`)

### 1.2 当前图简化策略

| 功能 | 实现 | 问题 |
|------|------|------|
| 互逆反应抵消 | `reduce_graph()` - 简单相减 | 无法处理复杂循环反应 |
| 重要节点提取 | `write_dot_file_significant_nodes()` - DFS找路径 | 路径权重计算不科学，时间复杂度过高 |
| 网络流分析 | `network_flow_solve()` - 未完全实现 | Dinic算法为空(TODO) |
| 反应筛选 | 按`count`排序取前N个 | 忽略网络拓扑结构 |

### 1.3 核心问题诊断

#### 问题1：反应噪声过多
**现象**：燃烧体系产生数百个"反应"，大多数是碎片重组或瞬态碰撞

**根本原因**：
- 缺乏反应类型分类（断裂/结合/交换/重排）
- 没有区分热力学稳定产物与动力学瞬态
- 原子重叠阈值过低（只要有1个原子重叠就算反应）

**示例问题场景**：
```
C2H5OH + OH → C2H5O + H2O  (主干反应，应该保留)
CH3 + CH3 → C2H6           (重要反应)
C2H4 + H → C2H5 → C2H4 + H (瞬态吸附-脱附，应该过滤)
```

#### 问题2：缺乏反应重要性量化
**现象**：无法区分"燃料消耗→CO2"主干路径与"自由基→碎片"次要路径

**根本原因**：
- 只有局部统计（单个反应次数）
- 缺乏全局网络分析（PageRank、介数中心性等）
- 没有考虑物质转化率的时间累积效应

#### 问题3：图过于密集难以可视化
**现象**：即使取前60个反应，DOT图仍是一团乱麻

**根本原因**：
- 没有社区发现(Community Detection)来分组相关反应
- 缺乏层次化聚合（如按反应类型、按物质类别）
- 没有动态时间切片分析

---

## 2. 理论基础与启发式搜索

### 2.1 化学启发式规则

#### 规则1：稳定性启发式
```
稳定分子定义：
- 存在时间 > 系统总时间的 5%
- 最大分子数 > 系统总分子数的 1%
- 在产物侧出现次数明显多于反应物侧（净生成）

反应重要性 ∝ 产物稳定性 × 反应物稳定性
```

#### 规则2：原子经济启发式
```
有效反应：原子转移数 / min(反应物原子数, 产物原子数) > 0.3

高价值反应特征：
- 重排反应：原子转移率高(>0.8)，原子总数守恒
- 分解反应：大分子→小分子，原子转移率中等(0.3-0.7)
- 化合反应：小分子→大分子，原子转移率中等(0.3-0.7)
- 无效"反应"：原子转移率极低(<0.1)，通常是碰撞或识别误差
```

#### 规则3：路径连通性启发式
```
主干路径特征（借鉴代谢网络分析）：
- 从初始反应物到最终产物存在通路
- 路径上的反应具有时间连续性
- 路径上物质的"流量"（分子数×次数）单调递减不超过50%
```

### 2.2 图论算法工具箱

#### 算法A：带权介数中心性(Weighted Betweenness Centrality)
```
目标：识别网络中的"关键桥梁"反应

定义：
BC(e) = Σ (s≠v≠t) [σ_st(e) / σ_st]
其中：
- σ_st: 节点s到t的最短路径数
- σ_st(e): 经过边e的最短路径数
- 权重：边权重 = -log(反应概率)

优化（Brandes算法）：O(VE + V²logV)
```

#### 算法B：Leiden社区发现
```
目标：将反应分组为功能模块（如氧化模块、裂解模块）

优势（对比Louvain）：
- 保证社区连通性
- 更快的收敛速度
- 更适合有向加权图

适应度函数：
Q = (1/2m) Σ_ij [A_ij - γ(k_i^in × k_j^out)/2m] δ(c_i, c_j)
其中γ为分辨率参数，控制社区粒度
```

#### 算法C：最大流量-最小割(Max-Flow Min-Cut)
```
目标：识别限制整体转化速率的"瓶颈反应"

构建：
- 源点：初始高浓度反应物
- 汇点：最终稳定产物
- 边容量：反应速率（count/时间跨度）

应用：
- 瓶颈边 = 最小割集中的边
- 这些反应是调控整体反应网络的关键
```

#### 算法D：时间序列模式挖掘
```
目标：识别协同反应模式

方法：
1. 将每个反应的时序转化为二进制序列（发生/未发生）
2. 使用Apriori算法挖掘频繁项集
3. 支持度 > 阈值且置信度 > 阈值的反应组合构成"协同模块"
```

---

## 3. 改进方案一：多层次反应重要性评估

### 3.1 反应特征向量设计

```cpp
struct ReactionFeatures {
    // 基础统计
    int count;                          // 发生次数
    int atom_transfer;                  // 转移原子数
    double time_span;                   // 首次→末次发生的时间跨度
    
    // 拓扑特征
    double betweenness_centrality;      // 介数中心性
    double pagerank_source;             // 反应物PageRank
    double pagerank_target;             // 产物PageRank
    
    // 化学特征
    double atom_economy;                // 原子经济率
    bool is_rearrangement;              // 是否重排反应
    double reactant_stability;          // 反应物稳定性得分
    double product_stability;           // 产物稳定性得分
    
    // 时序特征
    double temporal_entropy;            // 时间分布熵（衡量是否均匀分布）
    int burst_count;                    // 爆发次数（短时间密集发生）
    
    // 综合重要性得分
    double importance_score;
};
```

### 3.2 稳定性评估算法

```cpp
// 基于物种时间序列计算稳定性
void calculate_stability_scores(Universe* uv) {
    auto& formulas_nums = uv->species_counter->formulas_nums;
    
    for (const auto& [formula, time_series] : formulas_nums) {
        // 1. 存在时间比例
        int exist_frames = 0;
        for (float count : time_series) {
            if (count > 0) exist_frames++;
        }
        double existence_ratio = (double)exist_frames / time_series.size();
        
        // 2. 最大浓度
        double max_conc = *std::max_element(time_series.begin(), time_series.end());
        double total_atoms = uv->system->atoms.size();  // 近似
        double peak_ratio = max_conc / total_atoms;
        
        // 3. 变化趋势（净生成 vs 净消耗）
        double first_half_avg = std::accumulate(time_series.begin(), 
            time_series.begin() + time_series.size()/2, 0.0) / (time_series.size()/2);
        double second_half_avg = std::accumulate(time_series.begin() + time_series.size()/2, 
            time_series.end(), 0.0) / (time_series.size()/2);
        double trend_ratio = second_half_avg / (first_half_avg + 1e-6);
        
        // 综合稳定性得分
        stability_score[formula] = 
            0.4 * existence_ratio + 
            0.3 * std::min(peak_ratio * 10, 1.0) +  // 归一化
            0.3 * (trend_ratio > 1 ? 0.5 : 0);      // 净生成加分
    }
}
```

### 3.3 综合重要性评分模型

```cpp
double calculate_importance(const ReactionFeatures& rf) {
    // 层次1：基础重要性（基于频次和原子经济）
    double base_importance = 
        0.3 * normalize(rf.count) +
        0.2 * rf.atom_economy +
        0.1 * (rf.time_span / total_simulation_time);
    
    // 层次2：拓扑重要性（基于网络位置）
    double topo_importance = 
        0.15 * normalize(rf.betweenness_centrality) +
        0.1 * (rf.pagerank_source + rf.pagerank_target) / 2;
    
    // 层次3：化学合理性（基于稳定性）
    double chem_importance =
        0.1 * (rf.reactant_stability + rf.product_stability) / 2 +
        0.05 * (rf.is_rearrangement ? 0.8 : 0.5);  // 重排反应通常更重要
    
    return base_importance + topo_importance + chem_importance;
}
```

---

## 4. 改进方案二：基于社区发现的网络简化

### 4.1 算法设计

```cpp
class ReactionCommunityDetector {
public:
    // Leiden算法适配有向加权图
    void detect_communities(ReaxFlow* rf) {
        // 步骤1：初始化每个节点为独立社区
        std::unordered_map<Node*, int> node_community;
        int next_comm_id = 0;
        for (auto& node : rf->nodes) {
            node_community[node] = next_comm_id++;
        }
        
        // 步骤2：局部移动阶段
        bool improved = true;
        while (improved) {
            improved = false;
            for (auto& node : rf->nodes) {
                int best_community = node_community[node];
                double best_gain = 0;
                
                // 尝试移动到邻居社区
                std::unordered_set<int> neighbor_communities;
                for (auto& neighbor : node->to_nodes) {
                    neighbor_communities.insert(node_community[neighbor]);
                }
                for (auto& neighbor : node->from_nodes) {
                    neighbor_communities.insert(node_community[neighbor]);
                }
                
                for (int comm : neighbor_communities) {
                    double gain = calculate_modularity_gain(node, comm, node_community, rf);
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_community = comm;
                    }
                }
                
                if (best_community != node_community[node]) {
                    node_community[node] = best_community;
                    improved = true;
                }
            }
        }
        
        // 步骤3：聚合社区为超节点
        aggregate_communities(node_community, rf);
    }
    
private:
    // 适应有向加权图的模块度增益计算
    double calculate_modularity_gain(Node* node, int target_comm, 
        const std::unordered_map<Node*, int>& node_community, ReaxFlow* rf) {
        
        double gain = 0;
        double m = total_edge_weight(rf);  // 总边权重
        
        // 移入社区的边权重
        double k_in = 0, k_out = 0;
        for (auto& neighbor : node->to_nodes) {
            if (node_community.at(neighbor) == target_comm) {
                Edge* e = rf->get_edge(node, neighbor);
                k_out += e ? e->count : 0;
            }
        }
        for (auto& neighbor : node->from_nodes) {
            if (node_community.at(neighbor) == target_comm) {
                Edge* e = rf->get_edge(neighbor, node);
                k_in += e ? e->count : 0;
            }
        }
        
        // 社区的度
        double sigma_in = community_in_degree(target_comm, node_community, rf);
        double sigma_out = community_out_degree(target_comm, node_community, rf);
        
        // 有向图模块度增益
        gain = (k_in + k_out) / m - 
               (node->in_degree * sigma_out + node->out_degree * sigma_in) / (m * m);
        
        return gain;
    }
};
```

### 4.2 社区聚合可视化

```cpp
void write_community_aggregated_dot(ReaxFlow* rf, 
    const std::unordered_map<Node*, int>& communities) {
    
    // 统计每个社区的特征
    std::map<int, CommunityStats> stats;
    for (auto& [node, comm_id] : communities) {
        stats[comm_id].node_count++;
        stats[comm_id].total_degree += node->degree;
        stats[comm_id].formulas.push_back(node->molecule->formula);
    }
    
    // 为社区命名（基于最频繁的分子类型）
    std::map<int, std::string> comm_names;
    for (auto& [comm_id, stat] : stats) {
        comm_names[comm_id] = classify_community(stat.formulas);
        // 例如："Oxidation_Pool", "C1_Cracking", "Aromatics_Formation"
    }
    
    // 构建社区间的边（聚合）
    std::map<std::pair<int, int>, int> comm_edges;
    for (auto& edge : rf->edges) {
        int src_comm = communities.at(edge->source);
        int tgt_comm = communities.at(edge->target);
        if (src_comm != tgt_comm) {
            comm_edges[{src_comm, tgt_comm}] += edge->count;
        }
    }
    
    // 输出简化的DOT文件
    FILE* fp = create_file("reactions_community_aggregated.dot");
    fmt::print(fp, "digraph CommunityReactionFlow {{\n");
    fmt::print(fp, "  rankdir=LR;\n");
    fmt::print(fp, "  node [shape=box, style=filled];\n\n");
    
    // 社区节点
    for (auto& [comm_id, stat] : stats) {
        std::string color = get_community_color(comm_id);
        fmt::print(fp, "  comm{} [label=\"{}\\n{} molecules\\nDeg:{}\", fillcolor={}];\n",
                   comm_id, comm_names[comm_id], stat.node_count, stat.total_degree, color);
    }
    
    // 社区间边
    for (auto& [pair, count] : comm_edges) {
        fmt::print(fp, "  comm{} -> comm{} [label=\"{}\", penwidth={}];\n",
                   pair.first, pair.second, count, 1 + log(count));
    }
    
    fmt::print(fp, "}}\n");
    fclose(fp);
}
```

---

## 5. 改进方案三：动态主干路径提取

### 5.1 时序感知的网络流分析

```cpp
class TemporalPathAnalyzer {
public:
    // 构建时间扩展网络
    struct TemporalEdge {
        int time_frame;
        Node* source;
        Node* target;
        int flow;
    };
    
    std::vector<std::vector<Node*>> extract_main_paths(ReaxFlow* rf, 
        Node* source, Node* target, int max_paths = 5) {
        
        std::vector<std::vector<Node*>> main_paths;
        
        // 使用Yen算法找K条最短路径（基于反应重要性）
        auto k_shortest = yen_k_shortest_paths(rf, source, target, max_paths);
        
        // 筛选符合时间顺序的路径
        for (auto& path : k_shortest) {
            if (is_temporally_valid(path)) {
                main_paths.push_back(path);
            }
        }
        
        return main_paths;
    }
    
private:
    // 检查路径是否符合时间顺序（反应物应先于产物出现）
    bool is_temporally_valid(const std::vector<Node*>& path) {
        // 基于物种浓度曲线的交叉相关性
        for (size_t i = 0; i < path.size() - 1; i++) {
            auto& formula_src = path[i]->molecule->formula;
            auto& formula_tgt = path[i+1]->molecule->formula;
            
            // 获取时间序列
            auto& ts_src = get_time_series(formula_src);
            auto& ts_tgt = get_time_series(formula_tgt);
            
            // 计算时间延迟互相关
            int lag = calculate_optimal_lag(ts_src, ts_tgt);
            
            // 如果产物先于反应物出现（lag < 0），可能是逆反应或数据噪声
            if (lag < -5) {  // 允许5帧的误差
                return false;
            }
        }
        return true;
    }
    
    // 计算两个时间序列的最优时间延迟
    int calculate_optimal_lag(const std::vector<float>& ts1, 
                              const std::vector<float>& ts2) {
        int max_lag = ts1.size() / 4;
        double best_corr = -1;
        int best_lag = 0;
        
        for (int lag = -max_lag; lag <= max_lag; lag++) {
            double corr = cross_correlation(ts1, ts2, lag);
            if (corr > best_corr) {
                best_corr = corr;
                best_lag = lag;
            }
        }
        return best_lag;
    }
};
```

### 5.2 瓶颈反应识别

```cpp
void identify_bottleneck_reactions(ReaxFlow* rf) {
    // 构建最大流网络
    // 源点：初始帧的高浓度物种
    // 汇点：最终帧的高浓度物种
    
    auto source_nodes = get_initial_abundant_species(rf);
    auto sink_nodes = get_final_abundant_species(rf);
    
    // 添加超级源点和超级汇点
    Node* super_source = new Node(nullptr);  // 虚拟节点
    Node* super_sink = new Node(nullptr);
    
    for (auto& src : source_nodes) {
        Edge* e = new Edge(super_source, src);
        e->count = src->degree;  // 初始供应量
        rf->edges.insert(e);
    }
    
    for (auto& sink : sink_nodes) {
        Edge* e = new Edge(sink, super_sink);
        e->count = sink->degree;  // 最终需求量
        rf->edges.insert(e);
    }
    
    // Dinic算法计算最大流
    int max_flow = dinic_max_flow(rf, super_source, super_sink);
    
    // 最小割 = 瓶颈反应
    auto min_cut = find_min_cut(rf, super_source);
    
    fmt::print("=== Bottleneck Reactions (Min-Cut) ===\n");
    fmt::print("Max Flow: {}\n", max_flow);
    fmt::print("Bottleneck count: {}\n", min_cut.size());
    
    for (auto& edge : min_cut) {
        fmt::print("BOTTLENECK: {} -> {} (count={})\n",
                   edge->source->molecule->formula,
                   edge->target->molecule->formula,
                   edge->count);
    }
}
```

---

## 6. 改进方案四：时空关联反应聚类

### 6.1 反应模式挖掘

```cpp
class ReactionPatternMiner {
public:
    // 挖掘频繁反应模式
    void mine_frequent_patterns(ReaxFlow* rf, double min_support = 0.1) {
        // 将反应网络转换为事务数据库
        // 每条"事务" = 一个时间窗口内发生的所有反应
        
        int window_size = std::max(10, rf->sorted_edges.size() / 20);
        std::vector<std::vector<std::string>> transactions;
        
        // 按时间窗口分组
        std::map<int, std::vector<std::string>> time_window_reactions;
        for (auto& edge : rf->edges) {
            // 假设Edge结构增加first_frame字段记录首次发生帧
            int window = edge->first_frame / window_size;
            std::string rxn_id = fmt::format("{}->{}",
                edge->source->molecule->formula,
                edge->target->molecule->formula);
            time_window_reactions[window].push_back(rxn_id);
        }
        
        // Apriori算法挖掘频繁项集
        auto frequent_itemsets = apriori(time_window_reactions, min_support);
        
        // 输出协同反应组
        fmt::print("=== Co-occurring Reaction Patterns ===\n");
        for (auto& [itemset, support] : frequent_itemsets) {
            if (itemset.size() >= 2) {
                fmt::print("Pattern (support={:.2f}):\n", support);
                for (auto& rxn : itemset) {
                    fmt::print("  - {}\n", rxn);
                }
            }
        }
    }
    
private:
    std::vector<std::pair<std::vector<std::string>, double>> apriori(
        const std::map<int, std::vector<std::string>>& transactions,
        double min_support) {
        
        // 计算单项支持度
        std::map<std::string, int> item_count;
        int total_trans = transactions.size();
        
        for (auto& [window, items] : transactions) {
            std::unordered_set<std::string> unique_items(items.begin(), items.end());
            for (auto& item : unique_items) {
                item_count[item]++;
            }
        }
        
        // 筛选频繁1-项集
        std::vector<std::string> frequent_1;
        for (auto& [item, count] : item_count) {
            if ((double)count / total_trans >= min_support) {
                frequent_1.push_back(item);
            }
        }
        
        // 迭代生成更高阶项集（简化实现）
        // ... 
        
        return result;
    }
};
```

### 6.2 反应爆发检测

```cpp
void detect_reaction_bursts(ReaxFlow* rf) {
    // 对每个反应，检测其时间序列中的爆发期
    
    std::map<Edge*, std::vector<std::pair<int, int>>> bursts;  // 反应 -> [(start,end),...]
    
    for (auto& edge : rf->edges) {
        // 假设Edge增加occurrences向量记录每次发生的帧号
        auto& occ = edge->occurrences;
        if (occ.size() < 10) continue;
        
        // 计算局部密度
        std::vector<double> densities = calculate_local_density(occ);
        
        // 基于阈值识别爆发区间
        double threshold = std::accumulate(densities.begin(), densities.end(), 0.0) / densities.size() * 2;
        
        int burst_start = -1;
        for (size_t i = 0; i < densities.size(); i++) {
            if (densities[i] > threshold && burst_start == -1) {
                burst_start = occ[i];
            } else if (densities[i] <= threshold && burst_start != -1) {
                bursts[edge].push_back({burst_start, occ[i]});
                burst_start = -1;
            }
        }
    }
    
    // 分析爆发同步性
    fmt::print("=== Synchronized Reaction Bursts ===\n");
    for (auto& [edge1, bursts1] : bursts) {
        for (auto& [edge2, bursts2] : bursts) {
            if (edge1 >= edge2) continue;
            
            // 计算爆发区间重叠
            int overlap = calculate_burst_overlap(bursts1, bursts2);
            if (overlap > 0) {
                fmt::print("Synchronized: {} -> {} and {} -> {} (overlap: {} frames)\n",
                    edge1->source->molecule->formula, edge1->target->molecule->formula,
                    edge2->source->molecule->formula, edge2->target->molecule->formula,
                    overlap);
            }
        }
    }
}
```

---

## 7. 实现路线图与优先级

### 7.1 阶段一：基础改进（1-2周）

| 优先级 | 任务 | 预期效果 |
|--------|------|----------|
| P0 | 实现反应特征向量计算 | 为所有后续算法提供基础数据 |
| P0 | 增加原子经济过滤（阈值0.3） | 立即过滤50%以上的噪声反应 |
| P1 | 实现物种稳定性评分 | 区分瞬态与稳定物种 |
| P1 | 基于稳定性的反应筛选 | 提高反应网络质量 |

### 7.2 阶段二：图算法增强（2-3周）

| 优先级 | 任务 | 预期效果 |
|--------|------|----------|
| P0 | 实现带权介数中心性 | 识别关键桥梁反应 |
| P0 | 完成Dinic最大流算法 | 识别瓶颈反应 |
| P1 | 实现Leiden社区发现 | 网络分组可视化 |
| P1 | 社区聚合DOT输出 | 大幅简化可视化效果 |

### 7.3 阶段三：时序分析（2-3周）

| 优先级 | 任务 | 预期效果 |
|--------|------|----------|
| P1 | 时序交叉相关分析 | 识别时间因果关系 |
| P1 | 时序有效路径提取 | 提取化学合理的反应路径 |
| P2 | 反应爆发检测 | 识别阶段性反应事件 |
| P2 | 频繁模式挖掘 | 发现协同反应模块 |

### 7.4 阶段四：整合优化（1-2周）

| 优先级 | 任务 | 预期效果 |
|--------|------|----------|
| P0 | 综合重要性评分整合 | 统一的反应排序标准 |
| P1 | 多层次DOT输出 | 用户可选择的简化级别 |
| P1 | 性能优化 | 确保大规模系统可用 |

---

## 8. 参考理论框架

### 8.1 化学信息学

1. **代谢网络分析** (Metabolic Network Analysis)
   - 参考：Barabási, A.-L., & Oltvai, Z. N. (2004). Network biology
   - 应用：将燃烧反应网络视为代谢网络，应用其分析方法

2. **反应力场理论** (ReaxFF Theory)
   - 参考：van Duin et al. (2001) ReaxFF
   - 应用：理解反应坐标和过渡态，指导反应分类

### 8.2 图论与网络科学

1. **社区发现算法**
   - Leiden算法: Traag et al. (2019) Nature
   - Louvain算法: Blondel et al. (2008)
   
2. **网络流算法**
   - Dinic算法: Dinic (1970)
   - 最小割理论: Ford & Fulkerson (1956)

3. **中心性度量**
   - Betweenness Centrality: Brandes (2001) Algorithm
   - PageRank: Page et al. (1999)

### 8.3 时间序列分析

1. **时序模式挖掘**
   - Apriori算法: Agrawal & Srikant (1994)
   - 爆发检测: Kleinberg (2003)

2. **交叉相关分析**
   - 时延相关性: Box & Jenkins (1970)

### 8.4 相关软件工具

| 工具 | 功能 | 参考价值 |
|------|------|----------|
| **Chemotion** | 反应网络可视化 | 界面设计理念 |
| **NetworKit** | 大规模图分析 | 算法实现参考 |
| **PathwayTools** | 代谢路径分析 | 路径提取算法 |
| **RDKit** | 化学信息学 | 分子特征计算 |

---

## 附录：关键公式汇总

### 反应重要性综合评分

$$I_{rxn} = \underbrace{0.3\frac{c}{c_{max}} + 0.2AE + 0.1\frac{\Delta t}{T}}_{\text{基础}} + \underbrace{0.15\frac{BC}{BC_{max}} + 0.1\frac{PR_s + PR_t}{2}}_{\text{拓扑}} + \underbrace{0.1\frac{S_s + S_t}{2}}_{\text{化学}}$$

其中：
- $c$: 反应次数, $c_{max}$: 最大反应次数
- $AE$: 原子经济率 = $\frac{\text{转移原子数}}{\min(N_{reactant}, N_{product})}$
- $\Delta t$: 反应时间跨度, $T$: 总模拟时间
- $BC$: 介数中心性, $PR$: PageRank值
- $S$: 物种稳定性得分

### 有向图模块度

$$Q = \frac{1}{m}\sum_{ij}\left[A_{ij} - \gamma\frac{k_i^{in}k_j^{out}}{m}\right]\delta(c_i,c_j)$$

其中：
- $m$: 总边权重
- $\gamma$: 分辨率参数（默认1.0）
- $\delta(c_i,c_j)$: 社区指示函数

---

**文档版本**: 1.0  
**作者**: ReaxTools Development Team  
**日期**: 2025-03-05
