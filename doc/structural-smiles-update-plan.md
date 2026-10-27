# 结构身份升级方案：canonical SMILES v1

状态：**已实施（2.2）**  
基线：`main` @ `48b1a86` (chore: checkpoint before structural identity update)  
目标版本：2.2（输出 schema 破坏性变更，但不是核心化学算法不兼容）

## 0. 已确认决议

1. 验收标准：**图等价 round-trip**（不要求与 RDKit 字符串逐字节一致）。
2. SMILES **不包含键级**，全部单键。
3. H 的写法：**全部显式 [H]**；规范化在 H 抑制图上做，写出时展开为 [H] 分支。
4. ID：**64-bit 十进制字符串**。
5. `-me/-mr/-rc` 等合并/缩放功能**从 core 删除**，职责留给未来的 Python 套件。
6. reaction_events 列名与 `;` 分隔符按方案；不输出 atom ID 列表。
7. `web/` 暂不改；只改 Python `reax_tools_viz`。
8. 超过 **200 原子**不生成 SMILES，退化为化学式 + 公式 hash；性能目标 2x 基线。
9. SMILES **自写**，用 RDKit 充分验证正确性。
10. 直接破坏旧输出；core + Python 可视化一起发布，端到端负责。

## 0.1 实施结果

- SMILES 引擎：`src/cpp/smiles.h` / `smiles.cpp`，无外部依赖；200 原子 cutoff。
- RDKit 验证：289 个图（curated 异构体、C7–C10 全部烷烃、3/4-正则笼、
  自由基、芳环/桥环/螺环、铵/氧鎓、200/203 原子边界），**288 pass + 1
  expected fallback + 0 fail**；随机原子重排 3 次，SMILES 完全一致。
- 对 N 价 4、O 价 3 的超价原子写出形式电荷（如 `[NH4+]`），否则 RDKit
  拒绝解析；形式电荷只用于合法 valence，不是电荷平衡模型。
- 输出迁移：`molecules.csv`、宽表 `species_count.csv`、带 SMILES 的
  `reaction_events.csv`；删除 `reaction_snapshots/`、`molecules.json`、
  `-me/-mr/-rc`。
- Python：`reax_tools_viz` 已迁移；破旧的 `snapshots.py`/CLI 已删除；
  `molecules` 命令改用 SMILES + 可选 RDKit。
- 性能（`test/energetic.xyz`, 2592 原子 x 40 帧, `-nt 1`, 5 次中位数）：
  基线 0.476 s，新版本 0.607 s，约 **1.27x**，低于 2x 目标。
- 未改：`web/`；`test/energetic_v3/v4` 是旧 schema 历史输出；
  `tools/prepublication_validate.py` 等论文工具未迁移。

## 0.2 2.3 追加

- 删除 `reax_tools_manifest.json`；`reax_tools.log` 记录完整命令、参数、
  frames/atoms、species/formulas/events/transfer_edges 和输出列表。
- 新增 `plot` 的 `key_molecules_page1..3.png`：全局 top-48 结构，
  按时间平均含量 + 反应图加权度排序，RDKit 2D 绘制；不生成 CSV。
- N 价 4、O 价 3 写出形式电荷，保证 RDKit 可解析；验证扩到 289 个图。
- `(a)(b)` 仍只在 Python 端按全局 ID 排序分配，core 输出原始数值 ID。

## 1. 目标与非目标

### 目标

1. 每个分子有一个**结构身份**，不再由分子式 hash 代替：同分异构体必须分开。
2. 身份以 **canonical SMILES** 为特征码，同时保留**数值型 ID**。
3. 输出精简：`molecules.csv`、宽表 `species_count.csv`、带 SMILES 的
   `reaction_events.csv`；取消 `reaction_snapshots/`。
4. C++ core **不新增外部依赖**，SMILES 代码自带；与 RDKit 做幂等性/图等价验证。
5. 明确处理非稳态自由基：SMILES 必须显式固定 H 数，第三方解析器不得自行加氢。
6. Python 绘图套件同步迁移；性能与当前版本对比。

### 非目标

- 不引入电荷/自旋感知（当前 C++ 没有 charge 字段）。SMILES 表达的是
  “元素 + H 数 + 连接性”的 radical 图。
- 不做立体化学（顺反/手性）感知，除非后续单独加开关。
- 不要求与 RDKit `MolToSmiles` 的字符串逐字节相同；验收是图等价、
  自身确定、跨原子顺序不变（已确认）。

## 2. 新输出 schema

### 2.1 molecules.csv

```csv
id,formula,smiles
1124119999,C2H6O,[CH3][CH2][OH]
...
```

- 一行一个**结构**（不是一行一个 formula）。
- `id`：数值型，默认 64-bit 十进制字符串。
- 排序：`id` 数值升序。
- 不再包含 `atom_counts`、`stats`、`example`。

### 2.2 species_count.csv

经典宽表，列头为数值 ID：

```csv
frame,873214...,1124119999,...
1,0,1,...
2,0,1,...
```

- 第一列 `frame`，其余列是出现过的所有结构 ID（数值升序）。
- 单元格是该帧该结构的分子数；缺失为 0。
- 不再有 `formula`/`count` 长表列；formula/SMILES 查 `molecules.csv`。
- 必须保持“每帧所有分子数之和 = 该帧分子总数”。

### 2.3 reaction_events.csv

建议列（分隔符见待确认问题）：

```csv
event_id,frame,n_reactants,n_products,reactant_ids,product_ids,reactant_formulas,product_formulas,reactant_smiles,product_smiles
0,2,1,2,"1124119999","...;...","C2H6O","C2H5O;H","[CH3][CH2][OH]","...;[H]"
```

- 多值分隔符建议用 `;`，**不能用 `+`**，因为 `+` 是 SMILES 合法电荷符号。
- 不再输出 `*_coords`、`*_bonds`、`broken_bond_ids`、`created_bond_ids`。
- 不输出 bond 变化：SMILES 的差异已经表达了连接性变化。
- 是否保留 atom ID 列表待确认；默认不保留（`reaction_event_pairs.csv` 仍有实例重叠）。

### 2.4 reaction_event_pairs.csv（保留）

建议：

```csv
event_id,frame,source_id,target_id,source_formula,target_formula,source_smiles,target_smiles,tracked_source_id,tracked_target_id,atom_overlap
```

### 2.5 transfer_flow.csv（保留）

建议在现有列基础上增加公式/SMILES：

```csv
source_id,target_id,source_formula,target_formula,source_smiles,target_smiles,count,atom_transfer,self_loop
```

### 2.6 取消的输出

- `reaction_snapshots/` 整个目录，包括 `event_*.rxtsnap.json` 和 `snapshots_manifest.csv`。
- C++ 中 `AtomSnapshot`/`BondSnapshot`、`TrackedMolecule::atoms/bonds`、
  `save_reaction_snapshots()`、所有 snapshot JSON 打印函数。
- Python 中 `snapshots.py`、`snapshots` CLI 子命令。
- `bin/requirements.txt` 中的 `pyvista`、`pillow`（仅被 snapshots 使用）。

## 3. 核心设计：SMILES

### 3.1 图模型

使用**氢抑制图**做规范化，但写 SMILES 时 H 全部显式：

- 重原子顶点：`(element, h_count, charge=0)`。
- 终端 H（只连一个重原子）不单独进入规范化图，计入该重原子的 `h_count`。
- H-H 组分（H2、孤立 H）保留为显式 `[H]` 顶点。
- 边：来自 `Bond`，**只表达连接性，不含键级**（已确认）。

规范化在 H 抑制图上做，避免把 H 的排列对称性带进 IR 分支；写出时把每个
`h_count` 展开成 `([H])` 分支，所以最终 SMILES 中每个 H 都是显式的。

### 3.2 原子写法与自由基

所有重原子用 bracket 形式（无隐式 H），H 全部写成 `[H]`：

- 甲烷：`[C]([H])([H])([H])([H])`
- 甲基自由基：`[C]([H])([H])[H]`（3 H，1 未成对电子）
- 水：`[O]([H])[H]`
- 羟基自由基：`[O]([H])`
- 氨：`[N]([H])([H])[H]`
- 氢原子：`[H]`；氢分子：`[H][H]`

bracket 重原子没有 `H` 计数时表示 0 个隐式 H。这样 RDKit/Open Babel 等
解析器不会额外加氢。已用 RDKit 验证：`[C]([H])([H])[H]` 解析为 CH3，
不产生 CH4。

### 3.3 规范化算法

不依赖外部库，建议实现：

1. **初始颜色**：元素、H 数、度（不含键级）。
2. **1-WL 细化**：颜色 = hash(旧颜色, 排序后的邻居颜色多重集)，直到稳定。
3. **individualization-refinement**：
   - 若所有颜色类大小为 1，得到规范顺序。
   - 否则取最小非单例类，逐个 individualize，再细化，递归。
   - 叶子处直接生成 SMILES，用**字典序最小的 SMILES 字符串**作为规范形式；
     不再构造 O(n^2) 的邻接证书。
   - 对同色分支做等价剪枝；分子图最大度通常 <= 4，分支可控。
4. **SMILES 写出**：
   - 按规范顺序做 DFS，邻居按规范 rank 排序。
   - 非树边作为环闭合，环号按规范顺序分配；支持 `%nn`。
   - bracket 原子 + 显式 `[H]` 分支；分支必须先于“延续边”写出，
     否则括号会挂到子树的最后一个原子上（已修复并加入回归用例）。
   - 图连通，不需要 `.`；单原子分子直接输出 `[H]`/`[C]` 等。

**规模守卫**：默认 `max_atoms = 200`（含 H）。超过 200 原子不生成 SMILES，
直接退化为 formula 身份；<=200 原子时仍设 `max_refine_rounds`、
`max_ir_leaves`、`max_cell` 上限，超限同样退化并写日志。

### 3.4 数值 ID

- 正常分子：`structure_id = FNV1a64("smiles|" + canonical_smiles)`。
- >200 原子或 SMILES 生成失败：`structure_id = FNV1a64("formula|" + formula)`，
  `smiles` 留空；这是 formula 级身份，可能合并同分异构体（按约定可接受）。
- 输出为十进制字符串，保证 WASM/native 一致。
- 64-bit 在 10^6 物种时碰撞概率约 2.7e-8；若检测到同 ID 不同 SMILES，
  记录警告并用二次探测/128-bit 备用方案。
- 不用 32-bit：10^4 物种碰撞约 1.2%，10^5 约 69%。

### 3.5 与 RDKit 的幂等性验证

“官方包”指 Python RDKit。验证脚本：
`analysis/smiles_validate/validate_smiles.py`。

1. **图等价**：`RDKit(MolFromSmiles(our_smiles))` 与原始图在
   element/adjacency 上同构（显式 H 展开后）。
2. **canonical round-trip**：把原始显式图按全单键建 RDKit mol，
   `MolToSmiles` 应与解析我们 SMILES 后的 canonical 结果一致。
3. **自身幂等**：随机打乱原子顺序 3 次，SMILES 完全一致。

当前语料 284 个图（curated 异构体、C7–C10 全部烷烃、3/4-正则笼、
自由基、芳环/桥环/螺环、200 与 203 原子边界）：**283 pass，1 expected
fallback，0 fail**。

## 4. C++ 代码改动

### 新增

- `src/cpp/smiles.h` / `src/cpp/smiles.cpp`
  - `Result canonical_smiles(const std::vector<Atom>&, const std::vector<Bond>&, const Options&)`
  - `uint64_t smiles_id(const std::string&)` / `uint64_t formula_id(const std::string&)`
  - 内部：H 抑制图构建、WL、IR、SMILES writer、200 原子守卫。
- `tests/cpp/smiles_selftest.cpp`（独立 CMake target，不进生产 binary）
  - 固定用例：甲烷/乙烷/乙醇/二甲醚、C6H6、苯、环己烷、立方烷、
    萘、C12 环、C66/C67 边界、自由基、H2、单原子。
  - 每个用例做 5 次原子顺序重排，要求 SMILES 完全一致。
- `tests/cpp/smiles_cli.cpp`：stdin 读显式图、stdout 输出 SMILES 的测试 CLI。
- `analysis/smiles_validate/validate_smiles.py`：RDKit 图等价/round-trip/
  重排验证；结果写 `validation.csv`。

### 修改

- `system.h` / `system.cpp`
  - `Molecule` 增加 `std::string smiles; uint64_t structure_id;`
  - `update_topology()` 计算 formula、smiles、structure_id。
  - 拷贝构造/`operator==` 同步。
  - `process_counters()` 传 `(id, formula, smiles)` 列表；删除 `hash_counter`。
  - 删除 `process_reax_flow()`（它只服务于 molecules.json 的 stats；
    反应事件由 ReactionTracker 负责）。
- `reax_counter.h` / `reax_counter.cpp`
  - `SpeciesCounter` 改为按 `structure_id` 计数，保存 formula/smiles。
  - `save_file()` 写宽表 `frame,id,id,...`。
  - `brief_report` 适配；删除 `merge_by_element`/`rescale_all_by_element`
    及 `-me/-mr/-rc` 参数。
- `reax_flow.h` / `reax_flow.cpp`
  - `save_molecules_json()` -> `save_molecules_csv()`，只写 id/formula/smiles。
  - 删除 Node/Edge/stats/example 相关代码；保留一个轻量 species registry。
- `reax_track.h` / `reax_track.cpp`
  - `TrackedMolecule` 只保留 id/start/end/atom_ids/formula/smiles/structure_id/
    predecessors/successors。
  - 删除 AtomSnapshot/BondSnapshot、坐标 BFS 展开、snapshot 输出。
  - `find_or_create_tracked()` 比较 atom_ids + structure_id。
  - `save_raw_events()`/`save_raw_event_pairs()`/`save_transfer_flow()`
    使用 structure_id 和 smiles；多值用 `;`。
  - 删除 `save_events`/`save_molecule_lifetimes`（当前无调用，死代码）。
- `universe.cpp`
  - 删除 `save_reaction_snapshots()` 调用；改名 molecules 保存；删除
    `process_reax_flow()` 调用（若确认删除）。
- `main.cpp`
  - manifest/log 文件列表、版本、`identity_model = canonical-smiles-v1`。
- `VERSION` -> 2.2；`main.cpp` 版本字符串同步。

## 5. Python 绘图套件改动

- `io.py`：读 `molecules.csv`；宽表 `species_count.csv`；
  reaction events 新列；验证所有列 ID 在 molecules.csv 中。
- `models.py`：`MoleculeRecord(id, formula, smiles)`；
  `ReactionEvent` 增加 `reactant_smiles/product_smiles`。
- `plots.py`：宽表直接使用；物种列标签用 formula，若同 formula 多 ID
  追加短 ID 后缀（如 `C2H6O#1124`），避免图例重复。
- `filters.py`：多值分隔符改为 `;`；标签从 molecules.csv 取。
- `network.py`/`sankey.py`/`export.py`：transfer_flow 新列。
- `molecule_draw.py`：从 SMILES 画图。RDKit 作为**可选** extra；
  未安装时明确报错/跳过，不进入核心依赖。
- `cli.py`：删除 `snapshots` 子命令；`molecules` 读 CSV；
  物种选择支持 id/formula/SMILES。
- 删除 `snapshots.py`；`requirements.txt` 删除 `pyvista`/`pillow`。
- `tools/check_reax_outputs.py`、`audit_event_network_consistency.py`、
  `prepublication_validate.py` 同步 schema。
- README、`doc/methods-and-limits.md`、`doc/community_update_v2.1.md`
  更新输出表、身份模型、自由基/键级限制。

## 6. 测试计划

### 单元/属性测试

- 规范形式：同一分子随机重排原子顺序 5 次，SMILES 完全一致。
- 图等价：ours -> RDKit -> 原始图同构。
- 公式一致：RDKit formula == Molecule.formula（常见元素）。
- 自由基：CH3/CH2/CH/OH/NH2/NO/NO2/O2/N2/H/O/N 不额外加 H。
- 环：苯、环己烷、萘、立方烷、桥环、螺环、>9 环闭合。
- 单原子/双原子/H2/未知元素 X。
- 随机图：沿用 `analysis/hash_bench/datasets.py` 的 curated/alkanes/regular 图，
  转成显式 H 图后跑 ours -> RDKit 对比。

### 集成测试

- `test/energetic.xyz`：总运行时间、输出 schema、每帧分子数守恒。
- `test/isomer_probe/`：乙醇/二甲醚必须不同 ID；混合物 count=1+1；
  1->1 异构化必须产生事件。
- `tools/check_reax_outputs.py` 更新后通过。
- `reax_tools plot/network/events/molecules` 在新输出目录上运行通过。
- 旧输出不再兼容（版本 2.2），在 README 明确。

## 7. 性能计划

基线（当前 commit，`test/energetic.xyz`，2592 原子 x 40 帧，-nt 1）：
- 总时间约 **0.66 s**，~500 分子/帧，最大分子 32 原子（含 H）。
- 当前版本为 snapshot 做 BFS 展开/坐标/键拷贝，本身有成本；删除 snapshot
  可抵消一部分 SMILES 开销。

测量：

- 总墙钟时间、峰值 RSS、SMILES 生成时间、分子/秒。
- 分子大小分桶：1-5、6-10、11-20、21-50、51-200、>200 原子。
- 与 RDKit canonical SMILES（仅基准，不进生产）对比。
- 验收建议：`test/energetic.xyz` 总时间 <= 2x 基线（~1.3 s）；
  若超过，加入 species cache（formula + 度序列 + WL bucket + 精确同构校验）。

## 8. 分阶段

1. **Phase 1（风险最高）**：SMILES engine + 独立 selftest + RDKit 验证。
2. **Phase 2**：Molecule/SpeciesCounter/TrackedMolecule/ID 集成。
3. **Phase 3**：输出 schema（molecules.csv、宽 species_count、reaction_events、
   删除 snapshots）。
4. **Phase 4**：Python 绘图套件 + tools + README。
5. **Phase 5**：集成测试 + 性能对比 + 版本发布。

## 9. 风险

| 风险 | 影响 | 缓解 |
| --- | --- | --- |
| canonical labeling 实现错误 | 同结构不同 SMILES / 不同结构同 ID | RDKit 图等价 + 重排测试；selftest |
| 高对称大分子 IR 分支爆炸 | 性能/卡死 | 分支上限、等价剪枝、按分子大小分级 |
| 键级贪心不确定 | 若包含键级，ID 不稳定 | 不包含键级（推荐）或先规范化键级 |
| SMILES 过长 | 输出膨胀 | >200 原子直接不生成 SMILES，退化为 formula ID |
| 64-bit ID 在 JS 失精度 | Web 端错误 | 输出十进制字符串；Web 类型改 string |
| 旧输出/旧脚本不兼容 | 用户迁移成本 | 版本 2.2 + 迁移文档；不双 schema |

## 10. 默认实现假设

已确认项见第 0 节；其余默认：

1. 验收用**图等价 round-trip**，不要求与 RDKit 字符串完全一致。
2. SMILES **不包含键级**（全单键），先保证连接性/自由基；键级后续加开关。
3. 所有重原子 bracket；H 全部显式 `[H]`，不依赖隐式 H。
4. ID 用 **64-bit 十进制字符串**。
5. species_count 列按 ID 数值升序，缺失为 0。
6. `-me/-mr/-rc` 从 core 删除；未来由 Python 套件实现（已确认）。
7. reaction_events 多值用 `;`；保留 reaction_event_pairs.csv；
   默认不输出 atom IDs。
8. 绘图只改 Python `reax_tools_viz`；`web/` 不在 git 内，暂不改。
9. 删除 snapshots + pyvista/pillow；molecule 画图用可选 RDKit。
10. 性能目标 <= 2x 基线；超过则加 cache。
