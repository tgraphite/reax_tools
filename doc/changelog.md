# ReaxTools 开发日志

## 版本信息

**当前版本**: 1.98 (dev → main merge)
**发布日期**: 2025-03-06

---

## 版本 1.98 更新内容

### 核心改进

#### 1. 反应网络智能过滤 (Reaction Network Filtering)
- **原子经济过滤** (Atom Economy Filtering): 基于原子经济性移除低质量反应
  - 阈值: 0.30 (默认)
  - 解离/化合反应保护阈值: 0.15
  - 公式: AE = 转移原子数 / min(反应物原子数, 产物原子数)
  
- **边介数中心性排序** (Edge Betweenness Ranking): 识别网络关键"桥梁"反应
  - Brandes算法实现（采样优化）
  - 综合策略: 50% BC + 33% count + 17% atom_transfer
  - 当反应数 > 60 时自动触发

- **整合工作流**: 互逆对消 → AE过滤 → BC排序 → 空节点清理

#### 2. 用户友好的输出格式
- **化学命名**: 使用 `precursors` (前驱体) 和 `derivatives` (衍生物) 替代 `in_degree`/`out_degree`
- **Markdown报告**: `key_molecules_reactions.md` 替代 CSV
  - 清晰的层次结构
  - 网络统计信息
  - 前驱体/衍生物详细列表
- **反应路径分析**: 自动识别最长连续反应路径
  - DFS算法，处理环路
  - 显示前10条最长路径

#### 3. 输出文件优化
- `molecules_smiles.csv`: 按degree排序，包含拓扑信息
- `key_molecules_reactions.md`: 人类可读的Markdown格式
- DOT文件: 使用 `formula_hash` 格式，便于识别
- 移除: `species_count_hash.csv` (减少冗余)

### 技术实现

#### 代码变更
- `cpp/reax_flow.h`: Node结构新增字段
  - `precursor_count`, `derivative_count`
  - `precursor_reactions`, `derivative_reactions`
  - `precursor_atom_transfer`, `derivative_atom_transfer`

- `cpp/reax_flow.cpp`:
  - `filter_by_atom_economy()`: AE过滤
  - `filter_by_betweenness()`: BC排序
  - `cleanup_isolated_nodes()`: 空节点清理
  - `find_and_print_pathways()`: 路径分析
  - `save_molecule_centered_subgraphs()`: Markdown输出

### 性能指标

| 测试文件 | 过滤效果 | PNG大小(1.98) | PNG大小(1.0) | 压缩比 |
|----------|----------|---------------|--------------|--------|
| 2.5k.xyz | 407→48 | 214 KB | 12 MB | **1/56** |
| kerogen.xyz | 386→50 | 194 KB | 6.5 MB | **1/34** |

### 质量测试
- **内存测试**: Valgrind检测通过
  - definitely lost: 0 bytes
  - indirectly lost: 0 bytes
- **度数一致性**: degree = precursors + derivatives ✅

---

## 历史开发记录

### 2025-03-05 开发会话

#### 会话目标
1. 更新项目README文档
2. 研究并改进反应网络分析算法
3. 实现并测试新的反应过滤算法
4. 与稳定版进行对比评估

#### 完成的工作

##### 文档更新
- **文件**: `README_new.md` (7469 bytes)
- **内容**: 基于飞书文档重新组织项目介绍
- **改进**:
  - 添加详细的输入文件格式说明
  - 系统化命令行参数表格
  - 添加性能参考数据
  - 完善学习资源链接

##### 反应网络算法研究与设计
- **调研范围**: 
  - 代谢网络分析算法（Betweenness Centrality, PageRank）
  - 社区发现算法（Girvan-Newman, Leiden）
  - 反应力场理论（ReaxFF）
  - 原子经济性概念（Green Chemistry）

##### 关键提交
- `8911865`: Add atom economy filtering
- `d8cc86a`: Enable edge betweenness centrality ranking
- `6767fc5`: Fix degree calculation logic
- `55af12d`: UX improvements (precursor/derivative naming, Markdown output)
- `6f95d83`: Show reaction counts for precursors/derivatives

---

## 文件变更汇总

### 新增文件
```
doc/
├── changelog.md                                  # 本文件
├── README_Reaction_Network_Improvements.md
├── Reaction_Network_Algorithm_Improvements.md
├── Reaction_Network_Implementation_Guide.md
├── Reaction_Network_Quick_Reference.md

test/
├── V3_*/                                         # 测试输出
└── *_REPORT.md                                   # 测试报告
```

### 修改文件
```
cpp/reax_flow.h       # Node结构扩展，新增方法声明
cpp/reax_flow.cpp     # 核心过滤算法实现
cpp/main.cpp          # 移除species_count_hash.csv输出
```

---

## 后续开发计划

### 版本 1.99 (dev分支)
- [ ] 命令行参数支持自定义AE阈值
- [ ] 多级别过滤选项（宽松/标准/严格）
- [ ] Web版本动态阈值调整
- [ ] 更多测试案例验证

### 版本 2.0 (未来)
- [ ] 实时反应网络可视化
- [ ] 机器学习辅助反应分类
- [ ] 跨平台GUI界面

---

## 参考资源

### 核心算法论文
1. Brandes, U. (2001). A faster algorithm for betweenness centrality.
2. Girvan & Newman (2002). Community structure in social and biological networks.
3. Trost, B. M. (1995). Atom Economy: A Search for Synthetic Efficiency.

### 相关软件
- NetworkX, RDKit, Graphviz

---

*最后更新: 2025-03-06*  
*版本: 1.98*  
*维护者: ReaxTools Dev Team*
