# ReaxTools

ReaxTools 是面向 ReaxFF 与反应分子动力学轨迹的本地后处理工具：从 XYZ 或 LAMMPS dump 中识别键和分子，按结构统计物种、反应事件和原子转移，输出可审计的 CSV 与图表。

## 核心功能

| 功能 | 说明 |
| --- | --- |
| 轨迹读取 | XYZ 与 LAMMPS dump，支持周期性正交盒子 |
| 键与分子 | 按距离和元素价态成键，连通分量定义分子，统计 5–8 元环 |
| 结构身份 | canonical SMILES 与 64-bit ID，区分构造异构体 |
| 反应分析 | 结构 ID 变化触发事件，按原子重叠生成 transfer flow |
| 数据与绘图 | core 写 CSV 和 log；Python 负责标签、筛选、分页和绘图 |

一次分析会生成以下文件：

| 文件 | 内容 |
| --- | --- |
| `species_count.csv` | 宽表，第一列 `frame`，其余每列一个结构 ID，值为该帧计数；缺失为 0 |
| `molecules.csv` | 结构身份表：`id,formula,smiles` |
| `reaction_events.csv` | 每个事件的 frame、反应物/产物 ID、分子式和 SMILES |
| `reaction_event_pairs.csv` | 每个事件内反应物-产物实例的原子重叠对，用于审计 transfer flow |
| `transfer_flow.csv` | 结构间聚合转移边：source/target、formula、SMILES、count、atom_transfer、self_loop |
| `bond_count.csv` | 每帧各元素对的键数 |
| `atom_bonded_num_count.csv` | 每帧各元素的配位数分布 |
| `ring_count.csv` | 每帧 5–8 元环的数量 |
| `reax_tools.log` | 运行记录：版本、完整命令、参数、帧数/原子数、结构数、事件数、输出列表 |

结构 ID 是 canonical SMILES 的 64-bit FNV-1a hash。SMILES 只表达连接性，全部按单键写出，H 显式列出；N 价 4、O 价 3 时补形式电荷，保证 SMILES 可被 RDKit 解析。超过 200 个原子的结构退化为分子式身份，`molecules.csv` 中 `smiles` 留空。

`species_count.csv` 的列头是数值 ID。Python 读取后按 formula 分组、组内按 ID 升序分配 `C2H6O (a)`、`C2H6O (b)` 这样的显示标签；标签只用于图、事件聚合和 transfer 标签，不写回 core 生成的原始 CSV。

## 安装与快速使用

```bash
git clone https://github.com/tgraphite/reax_tools
cd reax_tools
bash install_reax_tools.sh
```

安装脚本会编译 `bin/reax_tools_core`，创建项目内的 Python 虚拟环境，安装绘图依赖，并把 `bin/` 加入 shell 启动文件。重新打开终端或 source 对应文件后即可使用 `reax_tools`。

依赖：

- CMake 3.14+，C++17 编译器。
- Python 3.8+，以及 `numpy`、`pandas`、`matplotlib`、`PyYAML`。
- 网络图需要 Graphviz 和 PyGraphviz。
- RDKit 可选，用于 `molecules` 命令和 `plot` 生成的 2D 关键结构图。

也可以直接编译 C++ core：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./bin/reax_tools_core -f test/energetic.xyz -o output
```

最常用的两步：

```bash
reax_tools -f test/energetic.xyz -o output
reax_tools plot -f output
```

LAMMPS dump 需要指定类型映射：

```bash
reax_tools -f dump.lammpstrj -t C,H,O,N -o output
```

常用参数：

| 参数 | 作用 |
| --- | --- |
| `-f, --traj` | 输入轨迹，`.xyz` 或 `.lammpstrj` |
| `-o, --output` | 输出目录 |
| `-t, --types` | LAMMPS dump 的类型顺序，如 `C,H,O,N` |
| `-r, --rescale-vdw` | 全局 vdW 半径缩放，默认 1.2 |
| `-tr, --type-radius` | 单元素半径覆盖，如 `N:1.5 O:1.6` |
| `-tv, --type-valence` | 单元素最大价态覆盖，如 `N:4 O:3` |
| `-nt, --threads` | 帧处理线程数 |
| `--no-rings` | 关闭环检测 |
| `--no-reactions` | 关闭 transfer-flow 分析 |
| `--no-track-reactions` | 关闭反应事件跟踪 |

绘图命令：

```bash
reax_tools plot -f output
reax_tools counts -f output
reax_tools events -f output
reax_tools network -f output
reax_tools flow -f output
reax_tools focus -f output
reax_tools molecules -f output
reax_tools summary -f output
```

`plot` 会生成物种计数分页图、键/配位/环统计图、反应事件分页图、分子转移图、反应转移图，以及 `key_molecules_page1.png` 到 `key_molecules_page3.png`。关键结构图取全局含量与反应图度数综合排名前 48 的结构，每页 16 个；RDKit 缺失时跳过该组图，其他图不受影响。绘图用的键级是显示用的 residual-valence 贪心结果，不进入结构 ID。

本地 20k 原子测例（未随仓库发布）为 20,736 原子 × 40 帧；4 线程耗时约 9.4 s，得到 192 个结构、10,410 个事件和 1,203 条 transfer 边。该数字仅用于说明规模。

## 原理简述

1. **键与分子**：对每对原子使用 `r_ij = scale * (R_i + R_j) / 2` 作为截断，周期性体系用 cell list 搜索近邻；按归一化距离排序后，在各自最大价态约束内贪心接受键。分子是键图的连通分量。
2. **结构身份**：在去掉 H 的图上做 1-WL 颜色细化与 individualization–refinement，枚举规范化标号，取字典序最小的 SMILES。ID 来自最终 SMILES，因此同分子式的不同连接方式会得到不同 ID。
3. **物种统计**：逐帧把分子按结构 ID 计数，写成宽表；同一结构在不同帧复用同一列。
4. **反应事件**：相邻帧中结构 ID 发生变化的原子作为种子，按旧/新分子划分求闭包；事件记录反应物与产物的实例集合。纯键级变化不改变结构 ID，因此不会单独触发事件；同分子式的连接性变化会触发。
5. **原子转移**：对事件内每个反应物/产物实例计算原子 ID 交集；非零交集形成 pair 记录，再按结构对聚合为 `transfer_flow.csv`。
6. **显示层**：Python 负责标签、筛选、分页和绘图。core 的 CSV 保持数值 ID 和原始字段，便于审计与二次处理。

当前限制：只处理正交盒子，不支持三斜晶胞和逐轴周期性；键级是几何/价态启发式，不代表电子结构；结构身份不区分键级、自旋和立体化学；超过 200 原子的结构按分子式合并。

## 社区发表论文

以下为部分已发表工作。体系一栏只写主要对象，方便了解 ReaxTools 的使用范围。

| 体系 | 论文 | 链接 |
| --- | --- | --- |
| CO₂ 捕集胺液热降解 | Thermal degradation and inhibition of long-chain triamines for CO2 capture, `Chem. Eng. J.` 2026 | [DOI](https://doi.org/10.1016/j.cej.2026.174329) |
| 等离子体石墨烯纳米片 | Experimental and ReaxFF MD simulation study on plasma gas-phase synthesis of graphene nanoflakes, `Diamond Relat. Mater.` 2025 | [DOI](https://doi.org/10.1016/j.diamond.2025.113141) |
| 三相凝胶泡沫灭火材料 | Optimization and performance assessment of an eco-friendly three-phase gel foam for liquid fires, `Colloids Surf. A` 2026 | [DOI](https://doi.org/10.1016/j.colsurfa.2026.141146) |
| 乙醇/Al(111) 表面反应 | Adsorption and thermal decomposition mechanisms of ethanol on Al(111), `J. Anal. Appl. Pyrolysis` 2026 | [DOI](https://doi.org/10.1016/j.jaap.2026.107932) |
| 铜表面脉冲电位重构 | Atomistic Origin of Non-equilibrium Surface Restructuring on Copper under Pulsed Potentials, `JACS` 2026 | [DOI](https://doi.org/10.1021/jacs.6c09411) |
| 桦甸油页岩干酪根热解 | A ReaxFF molecular dynamics study on the bond cleavage priority in the pyrolysis process of Huadian oil shale kerogen, `Fuel` 2026 | [DOI](https://doi.org/10.1016/j.fuel.2026.141081) |
| 纤维素蒸汽气化 | Fe2O3-assisted steam gasification of cellulose: a ReaxFF molecular dynamics study, `Process Saf. Environ. Prot.` 2026 | [DOI](https://doi.org/10.1016/j.psep.2026.109259) |

## 引用

当前版本为 2.3。软件尚无单独的方法论文；引用时请使用仓库地址、版本号和实际使用的 commit。测试输出和运行日志会记录版本、命令与参数，便于在论文中复现。
