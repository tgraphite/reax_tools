# ReaxTools 反应网络算法实现指南

本文档提供可直接实现的代码模板和逐步实施计划。

## 目录
1. [立即实施：噪声过滤](#1-立即实施噪声过滤)
2. [立即实施：反应分类](#2-立即实施反应分类)
3. [高优先级：介数中心性](#3-高优先级介数中心性)
4. [高优先级：最大流分析](#4-高优先级最大流分析)
5. [数据结构扩展](#5-数据结构扩展)
6. [性能优化建议](#6-性能优化建议)

---

## 1. 立即实施：噪声过滤

### 1.1 原子经济过滤

在 `reax_flow.cpp` 的 `save_graph()` 函数中添加过滤逻辑：

```cpp
void ReaxFlow::save_graph() {
    if (!FLAG_NO_REDUCE_REACTIONS) {
        reduce_graph();
    }
    
    // 新增：原子经济过滤
    filter_by_atom_economy(0.3);  // 阈值可配置
    
    update_graph();
    brief_report();
    // ... 后续代码
}

void ReaxFlow::filter_by_atom_economy(double threshold) {
    std::vector<Edge*> edges_to_remove;
    
    for (auto& edge : edges) {
        int reactant_atoms = edge->source->molecule->mol_atoms.size();
        int product_atoms = edge->target->molecule->mol_atoms.size();
        int min_atoms = std::min(reactant_atoms, product_atoms);
        
        if (min_atoms == 0) continue;  // 防止除零
        
        double atom_economy = (double)edge->atom_transfer / min_atoms;
        
        // 额外检查：排除极端情况
        bool is_extreme_case = (reactant_atoms > 50 && product_atoms < 3) ||
                               (reactant_atoms < 3 && product_atoms > 50);
        
        if (atom_economy < threshold && !is_extreme_case) {
            edges_to_remove.push_back(edge);
        }
    }
    
    // 移除低质量边
    for (auto& edge : edges_to_remove) {
        edges.erase(edge);
        edge_hash_to_edge.erase(edge->hash);
        delete edge;
    }
    
    if (!edges_to_remove.empty()) {
        fmt::print("Filtered {} low atom-economy reactions (threshold={})\n",
                   edges_to_remove.size(), threshold);
    }
}
```

### 1.2 瞬态物种过滤

```cpp
void ReaxFlow::filter_transient_species(SpeciesCounter* counter, double existence_threshold) {
    // 需要传入species_counter来获取时间序列数据
    std::unordered_set<unsigned int> transient_hashes;
    
    for (const auto& [formula, time_series] : counter->formulas_nums) {
        // 计算存在帧数
        int exist_frames = 0;
        for (float count : time_series) {
            if (count > 0.5) exist_frames++;  // 存在即计数
        }
        
        double existence_ratio = (double)exist_frames / time_series.size();
        
        // 计算最大浓度
        float max_count = *std::max_element(time_series.begin(), time_series.end());
        
        // 瞬态物种定义：存在时间短 且 最大浓度低
        if (existence_ratio < existence_threshold && max_count < 5.0) {
            // 找到对应的hash
            for (auto& node : nodes) {
                if (node->molecule->formula == formula) {
                    transient_hashes.insert(node->hash);
                    break;
                }
            }
        }
    }
    
    // 移除与瞬态物种相关的边（但保留其连接的重要反应）
    std::vector<Edge*> edges_to_remove;
    for (auto& edge : edges) {
        bool src_transient = transient_hashes.count(edge->source->hash);
        bool tgt_transient = transient_hashes.count(edge->target->hash);
        
        // 只有当两边都是瞬态时才移除
        if (src_transient && tgt_transient) {
            edges_to_remove.push_back(edge);
        }
    }
    
    for (auto& edge : edges_to_remove) {
        edges.erase(edge);
        edge_hash_to_edge.erase(edge->hash);
        delete edge;
    }
}
```

---

## 2. 立即实施：反应分类

### 2.1 反应类型识别

在 `reax_flow.h` 中添加：

```cpp
enum class ReactionType {
    UNKNOWN = 0,
    DISSOCIATION,      // 分解: A -> B + C
    ASSOCIATION,       // 化合: A + B -> C
    REARRANGEMENT,     // 重排: A -> B (原子数相同)
    SUBSTITUTION,      // 取代: A + B -> C + D (原子交换)
    TRANSFER,          // 转移: 原子/基团转移
    ELIMINATION        // 消去: 大分子失去小分子
};

struct ReactionClassifier {
    static ReactionType classify(Edge* edge) {
        int src_atoms = edge->source->molecule->mol_atoms.size();
        int tgt_atoms = edge->target->molecule->mol_atoms.size();
        int transferred = edge->atom_transfer;
        
        // 解析分子式获取元素组成
        auto src_formula = parse_formula_detailed(edge->source->molecule->formula);
        auto tgt_formula = parse_formula_detailed(edge->target->molecule->formula);
        
        // 判断逻辑
        if (src_atoms == tgt_atoms && transferred >= src_atoms * 0.8) {
            return ReactionType::REARRANGEMENT;
        }
        
        if (src_atoms > tgt_atoms && transferred >= tgt_atoms * 0.8) {
            // 检查是否是小分子脱落
            int diff = src_atoms - tgt_atoms;
            if (diff <= 3) {
                return ReactionType::ELIMINATION;
            }
            return ReactionType::DISSOCIATION;
        }
        
        if (tgt_atoms > src_atoms && transferred >= src_atoms * 0.8) {
            return ReactionType::ASSOCIATION;
        }
        
        if (transferred > 0 && transferred < std::min(src_atoms, tgt_atoms) * 0.5) {
            return ReactionType::SUBSTITUTION;
        }
        
        if (transferred > 0) {
            return ReactionType::TRANSFER;
        }
        
        return ReactionType::UNKNOWN;
    }
    
    static std::string type_to_string(ReactionType type) {
        switch (type) {
            case ReactionType::DISSOCIATION: return "DISSOCIATION";
            case ReactionType::ASSOCIATION: return "ASSOCIATION";
            case ReactionType::REARRANGEMENT: return "REARRANGEMENT";
            case ReactionType::SUBSTITUTION: return "SUBSTITUTION";
            case ReactionType::TRANSFER: return "TRANSFER";
            case ReactionType::ELIMINATION: return "ELIMINATION";
            default: return "UNKNOWN";
        }
    }
};
```

### 2.2 按类型输出统计

```cpp
void ReaxFlow::report_reaction_types() {
    std::map<ReactionType, int> type_counts;
    std::map<ReactionType, int> type_atom_transfers;
    
    for (auto& edge : edges) {
        auto type = ReactionClassifier::classify(edge);
        type_counts[type]++;
        type_atom_transfers[type] += edge->atom_transfer;
    }
    
    fmt::print("\n=== Reaction Type Statistics ===\n");
    fmt::print("{:<15s} {:<10s} {:<15s}\n", "Type", "Count", "Atom Transfers");
    
    for (auto& [type, count] : type_counts) {
        fmt::print("{:<15s} {:<10d} {:<15d}\n",
                   ReactionClassifier::type_to_string(type),
                   count,
                   type_atom_transfers[type]);
    }
}
```

---

## 3. 高优先级：介数中心性

### 3.1 Brandes算法实现

```cpp
void ReaxFlow::calculate_betweenness_centrality() {
    // 初始化
    for (auto& node : nodes) {
        node->betweenness = 0.0;
    }
    for (auto& edge : edges) {
        edge->betweenness = 0.0;
    }
    
    // 对每个节点作为源点执行BFS
    for (auto& s : nodes) {
        std::stack<Node*> S;
        std::queue<Node*> Q;
        std::unordered_map<Node*, std::vector<Node*>> P;  // 前驱
        std::unordered_map<Node*, int> sigma;              // 最短路径数
        std::unordered_map<Node*, int> dist;               // 距离
        std::unordered_map<Node*, double> delta;           // 依赖值
        
        // 初始化
        for (auto& w : nodes) {
            sigma[w] = 0;
            dist[w] = -1;
            delta[w] = 0.0;
        }
        
        sigma[s] = 1;
        dist[s] = 0;
        Q.push(s);
        
        // BFS遍历（带权图使用Dijkstra，这里简化为BFS）
        while (!Q.empty()) {
            Node* v = Q.front(); Q.pop();
            S.push(v);
            
            for (auto& w : v->to_nodes) {
                Edge* e = get_edge(v, w);
                if (!e) continue;
                
                // 使用反应次数作为距离权重（次数越多距离越短）
                int weight = std::max(1, 1000 / (e->count + 1));
                
                if (dist[w] < 0) {
                    dist[w] = dist[v] + weight;
                    Q.push(w);
                }
                
                if (dist[w] == dist[v] + weight) {
                    sigma[w] += sigma[v];
                    P[w].push_back(v);
                }
            }
        }
        
        // 回溯计算依赖值
        while (!S.empty()) {
            Node* w = S.top(); S.pop();
            
            for (auto& v : P[w]) {
                Edge* e = get_edge(v, w);
                if (!e) continue;
                
                double c = ((double)sigma[v] / sigma[w]) * (1.0 + delta[w]);
                delta[v] += c;
                e->betweenness += c;
            }
            
            if (w != s) {
                w->betweenness += delta[w];
            }
        }
    }
    
    // 归一化
    double max_bc = 0;
    for (auto& edge : edges) {
        max_bc = std::max(max_bc, edge->betweenness);
    }
    
    if (max_bc > 0) {
        for (auto& edge : edges) {
            edge->betweenness /= max_bc;
        }
    }
}
```

### 3.2 基于中心性的筛选

```cpp
void ReaxFlow::filter_by_betweenness(double percentile) {
    // 收集所有边的介数中心性
    std::vector<double> bc_values;
    for (auto& edge : edges) {
        bc_values.push_back(edge->betweenness);
    }
    
    std::sort(bc_values.begin(), bc_values.end());
    double threshold = bc_values[bc_values.size() * percentile];
    
    std::vector<Edge*> edges_to_remove;
    for (auto& edge : edges) {
        // 保留：高介数 或 高反应次数 或 高原子转移
        bool keep = (edge->betweenness >= threshold) ||
                    (edge->count >= 10) ||
                    (edge->atom_transfer >= 5);
        
        if (!keep) {
            edges_to_remove.push_back(edge);
        }
    }
    
    // 执行移除
    for (auto& edge : edges_to_remove) {
        edges.erase(edge);
        edge_hash_to_edge.erase(edge->hash);
        delete edge;
    }
    
    fmt::print("Filtered {} low-betweenness edges (threshold={:.3f})\n",
               edges_to_remove.size(), threshold);
}
```

---

## 4. 高优先级：最大流分析

### 4.1 Dinic算法完整实现

```cpp
class DinicMaxFlow {
public:
    struct FlowEdge {
        Node* to;
        int capacity;
        int flow;
        FlowEdge* reverse;
    };
    
    std::unordered_map<Node*, std::vector<FlowEdge*>> graph;
    std::unordered_map<Node*, int> level;
    std::unordered_map<Node*, int> ptr;
    
    void add_edge(Node* from, Node* to, int capacity) {
        FlowEdge* fwd = new FlowEdge{to, capacity, 0, nullptr};
        FlowEdge* rev = new FlowEdge{from, 0, 0, fwd};
        fwd->reverse = rev;
        
        graph[from].push_back(fwd);
        graph[to].push_back(rev);
    }
    
    bool bfs(Node* s, Node* t) {
        level.clear();
        std::queue<Node*> q;
        q.push(s);
        level[s] = 0;
        
        while (!q.empty()) {
            Node* v = q.front(); q.pop();
            
            for (auto& e : graph[v]) {
                if (e->flow < e->capacity && !level.count(e->to)) {
                    level[e->to] = level[v] + 1;
                    q.push(e->to);
                }
            }
        }
        
        return level.count(t);
    }
    
    int dfs(Node* v, Node* t, int pushed) {
        if (pushed == 0) return 0;
        if (v == t) return pushed;
        
        for (int& cid = ptr[v]; cid < graph[v].size(); cid++) {
            FlowEdge* e = graph[v][cid];
            
            if (e->flow >= e->capacity || level[e->to] != level[v] + 1)
                continue;
            
            int tr = dfs(e->to, t, std::min(pushed, e->capacity - e->flow));
            
            if (tr == 0) continue;
            
            e->flow += tr;
            e->reverse->flow -= tr;
            return tr;
        }
        
        return 0;
    }
    
    int max_flow(Node* s, Node* t) {
        int flow = 0;
        
        while (bfs(s, t)) {
            ptr.clear();
            while (int pushed = dfs(s, t, INT_MAX)) {
                flow += pushed;
            }
        }
        
        return flow;
    }
    
    // 找到最小割
    std::vector<FlowEdge*> min_cut(Node* s) {
        std::unordered_set<Node*> reachable;
        std::queue<Node*> q;
        q.push(s);
        reachable.insert(s);
        
        while (!q.empty()) {
            Node* v = q.front(); q.pop();
            
            for (auto& e : graph[v]) {
                if (e->flow < e->capacity && !reachable.count(e->to)) {
                    reachable.insert(e->to);
                    q.push(e->to);
                }
            }
        }
        
        // 割边 = 从reachable到非reachable的边
        std::vector<FlowEdge*> cut;
        for (auto& [node, edges] : graph) {
            if (!reachable.count(node)) continue;
            
            for (auto& e : edges) {
                if (!reachable.count(e->to) && e->capacity > 0) {
                    cut.push_back(e);
                }
            }
        }
        
        return cut;
    }
};
```

### 4.2 瓶颈分析集成

```cpp
void ReaxFlow::analyze_bottlenecks() {
    // 找出初始和最终的丰富物种
    std::vector<Node*> sources = get_abundant_species_at_frame(/*first*/ true);
    std::vector<Node*> sinks = get_abundant_species_at_frame(/*first*/ false);
    
    if (sources.empty() || sinks.empty()) {
        fmt::print("Warning: Cannot perform bottleneck analysis - no clear sources/sinks\n");
        return;
    }
    
    // 构建流网络
    DinicMaxFlow dinic;
    
    // 添加边（容量 = 反应次数）
    for (auto& edge : edges) {
        dinic.add_edge(edge->source, edge->target, edge->count);
    }
    
    // 创建超级源点和汇点
    Node* super_source = new Node(nullptr);
    Node* super_sink = new Node(nullptr);
    
    for (auto& src : sources) {
        dinic.add_edge(super_source, src, src->out_degree * 2);
    }
    
    for (auto& sink : sinks) {
        dinic.add_edge(sink, super_sink, sink->in_degree * 2);
    }
    
    // 计算最大流
    int max_flow = dinic.max_flow(super_source, super_sink);
    fmt::print("\n=== Network Flow Analysis ===\n");
    fmt::print("Max flow: {}\n", max_flow);
    
    // 找到瓶颈（最小割）
    auto cut = dinic.min_cut(super_source);
    
    fmt::print("Bottleneck reactions ({}):\n", cut.size());
    for (auto& e : cut) {
        if (e->to != super_sink && e->reverse->to != super_source) {
            fmt::print("  {} -> {} (capacity={}, flow={})\n",
                       e->reverse->to->molecule->formula,
                       e->to->molecule->formula,
                       e->capacity, e->flow);
        }
    }
    
    delete super_source;
    delete super_sink;
}
```

---

## 5. 数据结构扩展

### 5.1 Node和Edge结构扩展

```cpp
// reax_flow.h 中的修改

struct Node {
    Molecule* molecule = nullptr;
    unsigned int hash = 0;
    
    int degree = 0;
    int in_degree = 0;
    int out_degree = 0;
    int degree_at = 0;
    int in_degree_at = 0;
    int out_degree_at = 0;
    
    // 新增字段
    double betweenness = 0.0;           // 介数中心性
    double pagerank = 0.0;              // PageRank值
    double stability_score = 0.0;       // 稳定性得分
    
    // 时间信息
    int first_appearance = -1;
    int last_appearance = -1;
    std::vector<int> appearance_frames; // 出现的所有帧
    
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
    
    // 新增字段
    double betweenness = 0.0;           // 介数中心性
    int first_frame = -1;               // 首次发生帧
    std::vector<int> occurrences;       // 所有发生的帧号
    ReactionType type = ReactionType::UNKNOWN;  // 反应类型
    double importance_score = 0.0;      // 综合重要性
    
    Edge(Node* from_node, Node* to_node);
    ~Edge();
};
```

### 5.2 时间序列数据访问

```cpp
// 在SpeciesCounter中添加快速查找接口
class SpeciesCounter {
    // ... 现有代码
    
public:
    // 获取物种的时间序列（如果不存在返回空向量）
    std::vector<float> get_time_series(const std::string& formula) const {
        auto it = formulas_nums.find(formula);
        if (it != formulas_nums.end()) {
            return it->second;
        }
        return {};
    }
    
    // 检查物种是否存在于某帧
    bool exists_at_frame(const std::string& formula, int frame) const {
        auto ts = get_time_series(formula);
        if (frame >= 1 && frame <= ts.size()) {
            return ts[frame - 1] > 0.5;
        }
        return false;
    }
};
```

---

## 6. 性能优化建议

### 6.1 稀疏图优化

```cpp
// 对于大规模系统，使用稀疏表示
class SparseReactionGraph {
    // 使用压缩稀疏行(CSR)格式存储
    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    std::vector<double> values;
    
    // 预分配避免动态扩展
    void reserve(size_t n_nodes, size_t n_edges) {
        row_ptr.reserve(n_nodes + 1);
        col_idx.reserve(n_edges);
        values.reserve(n_edges);
    }
};
```

### 6.2 并行化策略

```cpp
// 并行计算介数中心性（近似算法）
void calculate_betweenness_parallel() {
    std::atomic<double> max_bc{0};
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < sorted_nodes.size(); i++) {
        // 采样部分节点作为源点（近似算法）
        if (i % 10 != 0) continue;  // 只取10%的节点
        
        Node* s = sorted_nodes[i].first;
        // ... 计算从s出发的贡献
        
        #pragma omp critical
        {
            // 合并结果
        }
    }
}
```

### 6.3 内存管理优化

```cpp
// 使用对象池避免频繁new/delete
class NodePool {
    std::vector<std::unique_ptr<Node>> pool;
    size_t next_available = 0;
    
public:
    Node* acquire(Molecule* mol) {
        if (next_available < pool.size()) {
            Node* node = pool[next_available++].get();
            // 重置状态
            node->molecule = new Molecule(*mol);
            node->hash = node->molecule->hash;
            return node;
        }
        // 分配新的
        pool.push_back(std::make_unique<Node>(mol));
        return pool.back().get();
    }
    
    void release_all() {
        next_available = 0;
    }
};
```

---

## 附录：快速启动检查清单

### 第1天：环境准备
- [ ] 备份现有代码
- [ ] 创建新分支 `feature/reaction-network-enhancement`
- [ ] 添加新的数据结构字段（Node/Edge扩展）

### 第2-3天：基础过滤
- [ ] 实现 `filter_by_atom_economy`
- [ ] 实现 `ReactionClassifier`
- [ ] 添加反应类型统计输出
- [ ] 测试：燃烧体系噪声减少 > 30%

### 第4-7天：中心性分析
- [ ] 实现 `calculate_betweenness_centrality`
- [ ] 实现 `filter_by_betweenness`
- [ ] 在DOT输出中标注高中心性边
- [ ] 测试：主干反应识别准确率

### 第8-14天：最大流分析
- [ ] 完整实现Dinic算法
- [ ] 实现 `analyze_bottlenecks`
- [ ] 集成到主流程
- [ ] 测试：瓶颈反应与化学直觉一致性

### 持续优化
- [ ] 性能测试（10万原子×1000帧）
- [ ] 内存使用优化
- [ ] 并行化加速

---

**本文档与 `Reaction_Network_Algorithm_Improvements.md` 配合使用**
