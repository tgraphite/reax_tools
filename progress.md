# 进度日志

## 会话: 2026-03-31

### 任务
分析 cpp/ 目录下项目时间复杂度，识别并优化非 O(N) 算法。

### 进度
- [x] 创建规划文件
- [x] Phase 1: 探索项目结构 (完成)
- [x] Phase 2: 逐个模块分析 (完成)
- [x] Phase 3: 时间复杂度评估 (完成)
- [x] Phase 4: 优化方案研究 (完成)

### 分析结果摘要

#### 识别的高复杂度算法（非 O(N)）

| 优先级 | 模块 | 函数 | 当前复杂度 | 问题描述 |
|--------|------|------|-----------|----------|
| 🔴 P0 | system.cpp | search_neigh_naive() | O(N²) | 无边界条件下的暴力搜索 |
| 🔴 P0 | system.cpp | find_rings_from_atom() | 指数级 | 环检测DFS+子集检查 |
| 🔴 P0 | reax_flow.cpp | write_dot_file_significant_nodes() | O(N⁴) | 双重循环+DFS枚举所有路径 |
| 🟡 P1 | reax_flow.cpp | find_and_print_pathways() | O(b^D) | DFS路径搜索 |
| 🟡 P1 | reax_flow.cpp | cleanup_isolated_nodes() | O(N×E) | 双重循环 |
| 🟡 P1 | reax_track.cpp | find_or_create_tracked() | O(T) | 线性搜索 |
| 🟢 P2 | reax_flow.cpp | calculate_edge_betweenness() | O(S×(V+E)) | 采样BFS |
| 🟢 P2 | reax_flow.cpp | update_graph() | O(N log N) | 排序操作 |

### 关键发现
1. **最严重问题**: `write_dot_file_significant_nodes` 使用DFS枚举所有路径，在稠密图中复杂度达 O(N⁴)
2. **环检测**: 当前DFS+子集检查方法对多环分子（如富勒烯）呈指数级复杂度
3. **邻居搜索**: 已有O(N)的cell list实现，但naive方法在无边界时仍被调用

### 生成的文档
- task_plan.md - 任务计划和进度跟踪
- findings.md - 详细的时间复杂度分析报告
- progress.md - 本进度日志

### 下一步行动（如需要实施优化）
1. 重构 write_dot_file_significant_nodes - 使用DP替代DFS
2. 重写环检测算法 - 使用BFS+并查集
3. 添加哈希索引优化 find_or_create_tracked

---

## 会话: 2026-03-31 下午

### 任务
用 git 回退到稳定版本，研究图论环检测算法。

### 执行的操作
1. ✅ 回退到 `20260331-morning` 提交 (d9640ca)
2. ✅ 查阅当前 system.cpp 环检测实现
3. ✅ 研究图论环检测算法：
   - Johnson 算法（枚举所有简单环）
   - Paton 算法（基本环基）
   - Vismara 算法（分子图 MCB）
   - RP-Path 方法（大分子 SSSR）
4. ✅ 分析 NetworkX cycle_basis 实现

### 关键发现

#### 现有实现问题
```cpp
// 当前实现：对每个原子都DFS + O(R²)子集检查
for (auto& start_atom : molecule->mol_atoms) {
    find_rings_from_atom(start_atom, start_atom, 0, ...);
}
```
- 对每个原子都启动 DFS，大量重复工作
- 子集检查是 O(R² × A) 的双重循环

#### 推荐算法：Paton 算法
- **来源**: Paton, K. Comm. ACM 12, 9 (1969)
- **复杂度**: O(V²)，分子图上接近线性
- **实现**: NetworkX 的 cycle_basis 函数使用此算法
- **优点**: 简单，产生基本环基，避免指数级枚举

#### 算法对比
| 算法 | 复杂度 | 适用场景 |
|-----|--------|---------|
| 当前DFS | 指数级 | - |
| Paton | O(V²) | 推荐实现 |
| Vismara MCB | O(V) | 最优但复杂 |
| Johnson | O(V+E+ec) | 枚举所有环 |

### 生成的文档
- findings.md 已追加图论算法研究附录

### 下一步建议
1. **短期**: 实现 Paton 算法替代当前 DFS
2. **中期**: 添加双连通分量分解
3. **验证**: 用 test/graphite_ring_detection_backup 验证结果一致性

---

## 会话: 2026-03-31 晚上 - Horton MCB 算法实现

### 任务
实现 Horton Minimum Cycle Basis 算法替代指数级 DFS 环检测。

### 执行的操作
1. ✅ 实现 Horton MCB 算法：
   - BFS 最短路径计算（O(V(V+E))）
   - 候选环生成（每条边 + 最短路径）
   - GF(2) 高斯消元提取独立环基
2. ✅ 编译测试
3. ✅ 验证结果正确性

### 关键发现

#### Paton vs Horton
- **Paton 算法**: 产生 FCB (Fundamental Cycle Basis)，不满足化学需要
- **Horton 算法**: 产生真正的 MCB (Minimum Cycle Basis)，结果与备份完全匹配

#### 性能对比
| 版本 | 时间 | 提升 |
|-----|------|------|
| 原始 DFS | 200-250 秒 | - |
| Horton MCB | 148 秒 | **~40%** |

#### 验证结果
```
Frame 1: R5=165, R6=2336, R7=142, R8=74  ✅ 与备份完全匹配
Frame 2: R5=177, R6=2365, R7=149, R8=59  ✅ 与备份完全匹配
```

### 实现细节
```cpp
// 1. BFS from each vertex for shortest paths
for (int start = 0; start < n; start++) {
    run BFS, record distances and parents
}

// 2. Generate candidate cycles from each edge
for each edge (u, v):
    if shortest_path(v, u) exists:
        cycle = edge(u,v) + path(v->u)
        add to candidates

// 3. GF(2) Gaussian elimination
sort candidates by size (prefer smaller cycles)
for each candidate:
    reduce against existing basis using XOR
    if independent, add to basis
```

### 下一步建议
1. ✅ **已完成**: 环检测算法优化
2. 可选：进一步优化 reax_flow 的 O(N⁴) 算法（用户指示暂时不处理）
3. 可选：为超大型分子（>1000原子）添加并行 BFS

