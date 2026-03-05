# ReaxTools 反应网络核心算法速查

## 当前算法流程图

```
┌─────────────────────────────────────────────────────────────────┐
│                        反应检测流程                               │
├─────────────────────────────────────────────────────────────────┤
│  1. 帧间分子对比                                                   │
│     for each atom:                                                │
│         prev_mol = atom.prev_frame.molecule                       │
│         curr_mol = atom.curr_frame.molecule                       │
│         if prev_mol.hash != curr_mol.hash:                        │
│             record_reaction(prev_mol, curr_mol)                   │
│                                                                   │
│  2. 原子重叠验证                                                   │
│     intersection = prev_mol.atoms ∩ curr_mol.atoms                │
│     if intersection.size > 0:                                     │
│         add_reaction(intersection.size)  // atom_transfer         │
│                                                                   │
│  3. 图简化                                                        │
│     reduce_graph()  // 互逆反应抵消                                │
│                                                                   │
│  4. 排序输出                                                      │
│     sort_edges_by_count()                                         │
│     output_top_n(n=60)                                            │
└─────────────────────────────────────────────────────────────────┘
```

## 核心问题与解决方案对照表

| 问题现象 | 根本原因 | 解决方案 | 实现难度 | 预期效果 |
|----------|----------|----------|----------|----------|
| 反应过多(数百个) | 缺乏原子经济过滤 | 添加AE阈值过滤(0.3) | ⭐ | 过滤50%噪声 |
| 主次不分 | 无全局网络分析 | 介数中心性计算 | ⭐⭐ | 识别关键桥梁 |
| 图太密集 | 无社区聚合 | Leiden社区发现 | ⭐⭐⭐ | 分组可视化 |
| 缺乏路径信息 | 无时序分析 | 交叉相关+K最短路径 | ⭐⭐⭐ | 提取主干路径 |
| 不知道瓶颈 | 无流量分析 | Dinic最大流+最小割 | ⭐⭐ | 识别限速步骤 |

## 关键算法伪代码

### 算法1：原子经济过滤（1小时实现）

```
function filter_by_atom_economy(threshold=0.3):
    for each edge in edges:
        src_size = edge.source.atoms.size
        tgt_size = edge.target.atoms.size
        min_size = min(src_size, tgt_size)
        
        atom_economy = edge.atom_transfer / min_size
        
        if atom_economy < threshold:
            mark_for_removal(edge)
    
    remove_marked_edges()
```

### 算法2：介数中心性（Brandes算法）（1天实现）

```
function calculate_betweenness():
    for each node in nodes:
        node.betweenness = 0
    
    for each source in nodes:          // 可以采样加速
        // BFS
        distance[source] = 0
        sigma[source] = 1
        queue = [source]
        
        while queue not empty:
            v = queue.pop()
            for each neighbor w of v:
                if distance[w] < 0:
                    distance[w] = distance[v] + 1
                    queue.push(w)
                if distance[w] == distance[v] + 1:
                    sigma[w] += sigma[v]
                    predecessor[w].push(v)
        
        // Back-propagation
        while stack not empty:
            w = stack.pop()
            for each v in predecessor[w]:
                delta[v] += (sigma[v]/sigma[w]) * (1 + delta[w])
                edge(v,w).betweenness += delta[v]
            if w != source:
                w.betweenness += delta[w]
```

### 算法3：社区发现（Leiden算法）（2天实现）

```
function detect_communities():
    // 初始化
    for each node:
        node.community = unique_id
    
    repeat:
        improved = false
        
        // 局部移动
        for each node:
            best_comm = node.community
            best_gain = 0
            
            for each neighbor_comm of node:
                gain = modularity_gain(node, neighbor_comm)
                if gain > best_gain:
                    best_gain = gain
                    best_comm = neighbor_comm
            
            if best_comm != node.community:
                node.community = best_comm
                improved = true
        
        // 聚合
        if improved:
            aggregate_communities()
    
    until not improved
```

### 算法4：最大流/最小割（Dinic算法）（1天实现）

```
function max_flow(source, sink):
    flow = 0
    
    while bfs_level_graph(source, sink):
        ptr[node] = 0 for all nodes
        
        while pushed = dfs(source, sink, INF) > 0:
            flow += pushed
    
    return flow

function bfs_level_graph(s, t):
    level = {}  // 重置
    queue = [s]
    level[s] = 0
    
    while queue not empty:
        v = queue.pop()
        for each edge e from v:
            if e.flow < e.cap and level[e.to] < 0:
                level[e.to] = level[v] + 1
                queue.push(e.to)
    
    return level[t] >= 0  // t是否可达

function dfs(v, t, pushed):
    if v == t or pushed == 0:
        return pushed
    
    for i from ptr[v] to graph[v].edges.size:
        e = graph[v].edges[i]
        if level[e.to] == level[v] + 1 and e.flow < e.cap:
            tr = dfs(e.to, t, min(pushed, e.cap - e.flow))
            if tr > 0:
                e.flow += tr
                e.reverse.flow -= tr
                return tr
        ptr[v]++
    
    return 0
```

### 算法5：时序路径验证（2小时实现）

```
function is_temporally_valid(path):
    for i from 0 to path.size - 2:
        src = path[i]
        tgt = path[i+1]
        
        // 获取浓度时间序列
        ts_src = get_concentration_series(src.formula)
        ts_tgt = get_concentration_series(tgt.formula)
        
        // 计算最优滞后
        lag = calculate_cross_correlation_lag(ts_src, ts_tgt)
        
        // 产物应滞后于反应物
        if lag < -5:  // 允许5帧误差
            return false
    
    return true

function calculate_cross_correlation_lag(ts1, ts2, max_lag=20):
    best_corr = -1
    best_lag = 0
    
    for lag from -max_lag to max_lag:
        corr = 0
        for i from max(0, -lag) to min(ts1.size, ts2.size - lag) - 1:
            corr += ts1[i] * ts2[i + lag]
        
        if corr > best_corr:
            best_corr = corr
            best_lag = lag
    
    return best_lag
```

## 数据结构修改清单

### Node结构新增字段
```cpp
struct Node {
    // 现有字段...
    
    // 中心性度量
    double betweenness = 0.0;
    double pagerank = 0.0;
    
    // 稳定性
    double stability_score = 0.0;
    int first_frame = -1;
    int last_frame = -1;
    
    // 社区
    int community_id = -1;
};
```

### Edge结构新增字段
```cpp
struct Edge {
    // 现有字段...
    
    // 中心性
    double betweenness = 0.0;
    
    // 时间
    int first_frame = -1;
    std::vector<int> occurrences;
    
    // 分类
    enum Type { UNKNOWN, DISSOCIATION, ASSOCIATION, REARRANGEMENT, ... } type;
    
    // 综合评分
    double importance_score = 0.0;
};
```

## 性能优化提示

### 对于10万原子×1000帧系统：

| 操作 | 当前复杂度 | 优化后 | 方法 |
|------|----------|--------|------|
| 介数中心性 | O(VE) ≈ 10^8 | O(0.1VE) | 10%采样近似 |
| 社区发现 | O(V^2) | O(E log V) | Leiden算法 |
| 最大流 | O(V^2 E) | O(E√V) | Dinic+当前弧优化 |
| 内存占用 | O(V+E) | O(V+E) | 对象池复用 |

### 并行化机会

```cpp
// 可并行化的计算
#pragma omp parallel for
- 介数中心性（源点间独立）
- 社区发现（节点移动）
- 稳定性评分（物种间独立）
- 交叉相关（反应对间独立）

// 需要串行
- 图修改操作（加边/删边）
- 文件输出
```

## 测试验证方案

### 测试用例1：乙醇燃烧
- **输入**: 已知燃烧机理的乙醇体系
- **验证**: 是否能识别出已知的C2H5OH→C2H4→CO→CO2主路径
- **指标**: 主路径恢复率 > 80%

### 测试用例2：聚乙烯裂解
- **输入**: 聚乙烯高温裂解模拟
- **验证**: 是否能识别出链断裂→小分子→气体的级联反应
- **指标**: 社区分组与化学类别一致（烷烃/烯烃/芳烃）

### 测试用例3：催化反应
- **输入**: 表面催化反应
- **验证**: 是否能识别催化剂表面的吸附-反应-脱附循环
- **指标**: 瓶颈反应识别准确率 > 70%

## 与现有代码集成点

```
Universe::process_traj()          // 现有流程
    ↓
System::process_reax_flow()       // 现有：检测反应
    ↓
ReaxFlow::add_reaction()          // 现有：添加边
    ↓
ReaxFlow::save_graph()            // 现有入口
    ├── NEW: filter_by_atom_economy()     // 原子经济过滤
    ├── NEW: calculate_betweenness()      // 中心性计算
    ├── NEW: detect_communities()         // 社区发现
    ├── NEW: analyze_bottlenecks()        // 瓶颈分析
    ├── EXISTING: reduce_graph()          // 互逆抵消
    └── EXISTING: write_dot_file()        // 输出
```

## 最小可行改进（MVP）

**如果只能实现一个改进，选择：原子经济过滤**

理由：
1. 实现简单（< 50行代码）
2. 效果显著（过滤30-50%噪声）
3. 无副作用（不影响有效反应）
4. 为基础（其他算法在干净数据上效果更好）

**如果只能实现两个改进，添加：介数中心性**

理由：
1. 识别网络关键桥梁
2. 与原子经济互补（AE过滤局部，BC评估全局）
3. 已有成熟算法（Brandes）
4. 可视化效果好（可在DOT中高亮显示）

---

**建议阅读顺序**：
1. 本速查文档（了解全貌）
2. `Reaction_Network_Implementation_Guide.md`（动手实现）
3. `Reaction_Network_Algorithm_Improvements.md`（深入理解）
