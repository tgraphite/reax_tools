# 任务计划：C++ 项目时间复杂度分析与优化

## 目标
分析 cpp/ 目录下主项目各个模块的时间复杂度，识别非 O(N) 算法，并研究改进措施。

## 项目结构
- **核心模块**: system.cpp/h, universe.cpp/h
- **反应分析模块**: reax_flow.cpp/h, reax_counter.cpp/h, reax_track.cpp/h
- **工具模块**: argparser.cpp, cell_list.cpp, rdkit_utils.cpp, string_tools.cpp

## 阶段

### Phase 1: 探索项目结构 ✓
- [x] 查找 cpp/ 目录及其子目录结构
- [x] 识别主要模块和核心算法文件
- [x] 列出所有源文件（.cpp/.h/.hpp）

### Phase 2: 逐个模块分析 ✓
- [x] system.cpp - 邻居搜索、成键、分子构建、环检测
- [x] reax_flow.cpp - 反应网络图构建、路径分析
- [x] reax_track.cpp - 反应追踪
- [x] reax_counter.cpp - 物种计数
- [x] cell_list.cpp - 邻居搜索优化数据结构

### Phase 3: 时间复杂度评估 ✓
- [x] 识别每个核心算法的复杂度
- [x] 标记非 O(N) 的算法
- [x] 分析复杂度瓶颈原因

### Phase 4: 优化方案研究 ✓
- [x] 针对每个非 O(N) 算法研究改进措施
- [x] 评估优化可行性和收益
- [x] 制定具体优化建议

## 主要发现总结

| 模块 | 算法/函数 | 当前复杂度 | 严重度 |
|------|----------|-----------|--------|
| system.cpp | search_neigh_naive() | O(N²) | 🔴 高 |
| system.cpp | find_rings_from_atom() | 指数级 | 🔴 高 |
| reax_flow.cpp | write_dot_file_significant_nodes() | O(N⁴) | 🔴 极高 |
| reax_flow.cpp | reduce_graph() | O(E²) | 🟡 中 |
| reax_flow.cpp | calculate_edge_betweenness() | O(S×(V+E)) | 🟡 中 |
| reax_track.cpp | find_or_create_tracked() | O(T) | 🟢 低 |
| reax_track.cpp | detect_reactions_from_changes() | O(C×M) | 🟡 中 |

## 错误记录
| 错误 | 尝试 | 解决方案 |
|------|------|----------|
| 无 | - | - |

## 决策记录
| 日期 | 决策 | 原因 |
|------|------|------|
| 2026-03-31 | 使用 hash 集合优化 O(N²) 循环 | 提升大规模数据性能 |
| 2026-03-31 | 回退到稳定版本，研究图论算法 | BFS实现有内存问题，先研究理论 |
| 2026-03-31 | 推荐 Paton 算法替代当前 DFS | O(V²)复杂度，NetworkX使用，实现简单 |

---

## 图论算法研究总结 (2026-03-31)

### 研究问题
环检测算法的数学基础与最优算法选择。

### 候选算法

| 算法 | 时间复杂度 | 空间复杂度 | 优点 | 缺点 |
|-----|-----------|-----------|------|------|
| **当前实现** | 指数级 | O(V) | 能找出所有环 | 稠密图上性能极差 |
| **Paton** | O(V²) | O(V²) | 简单，稳定 | 产生的是基而非MCB |
| **Vismara MCB** | O(V) | O(V) | 最优，化学专用 | 实现复杂 |
| **Johnson** | O(V+E+ec) | O(V) | 枚举所有环 | ec可能指数级 |

### 推荐方案
**Paton 算法**（短期）+ **Vismara MCB**（长期）

- Paton算法是NetworkX的默认实现，经过充分验证
- 对于分子图，O(V²)实际接近线性
- 可先解决性能问题，再考虑MCB的化学准确性
