# ReaxTools

ReaxTools 是面向 ReaxFF/机器学习势/AIMD 等任意反应分子动力学轨迹的本地后处理工具。ReaxTools从 XYZ 或 LAMMPS dump 中识别键和分子，按结构统计物种、反应事件和原子转移，输出可审计的 CSV 与图表。

## 核心功能

下表为 `test/energetic_20k` 的实际运行结果；左侧是功能，右侧是 `reax_tools plot` 生成的图。

| 功能 | energetic_20k 示例 |
| --- | --- |
| 物质含量 | <img src="test/figures/species_count.png" width="620" alt="物质含量"> |
| 键统计 | <img src="test/figures/bond_count.png" width="620" alt="键统计"> |
| 配位数 | <img src="test/figures/atom_bonded_num_count.png" width="620" alt="配位数"> |
| 环统计 | <img src="test/figures/ring_count.png" width="620" alt="环统计"> |
| 反应事件 | <img src="test/figures/reaction_events.png" width="620" alt="反应事件"> |
| 分子转移 | <img src="test/figures/molecule_transfer.png" width="620" alt="分子转移"> |
| 反应转移 | <img src="test/figures/reaction_transfer.png" width="620" alt="反应转移"> |
| 关键结构 | <img src="test/figures/key_molecules.png" width="620" alt="关键结构"> |

一次分析会生成以下文件：

| 文件                        | 内容                                                                              |
| --------------------------- | --------------------------------------------------------------------------------- |
| `species_count.csv`         | 宽表，第一列 `frame`，其余每列一个结构 ID，值为该帧计数；缺失为 0                 |
| `molecules.csv`             | 结构身份表：`id,formula,smiles`                                                   |
| `reaction_events.csv`       | 每个事件的 frame、反应物/产物 ID、分子式和 SMILES                                 |
| `reaction_event_pairs.csv`  | 每个事件内反应物-产物实例的原子重叠对，用于审计 transfer flow                     |
| `transfer_flow.csv`         | 结构间聚合转移边：source/target、formula、SMILES、count、atom_transfer、self_loop |
| `bond_count.csv`            | 每帧各元素对的键数                                                                |
| `atom_bonded_num_count.csv` | 每帧各元素的配位数分布                                                            |
| `ring_count.csv`            | 每帧 5–8 元环的数量                                                               |
| `reax_tools.log`            | 运行记录：版本、完整命令、参数、帧数/原子数、结构数、事件数、输出列表             |

## 安装与快速使用

```bash
git clone https://github.com/tgraphite/reax_tools
cd reax_tools
bash install_reax_tools.sh
```

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

典型使用方法：

```bash
reax_tools -f test/energetic.xyz -o output
reax_tools plot -f output
```

LAMMPS dump 需要指定类型映射：

```bash
reax_tools -f dump.lammpstrj -t C,H,O,N -o output
```

常用参数：

| 参数                   | 作用                                 |
| ---------------------- | ------------------------------------ |
| `-f, --traj`           | 输入轨迹，`.xyz` 或 `.lammpstrj`     |
| `-o, --output`         | 输出目录                             |
| `-t, --types`          | LAMMPS dump 的类型顺序，如 `C,H,O,N` |
| `-r, --rescale-vdw`    | 全局 vdW 半径缩放，默认 1.2          |
| `-tr, --type-radius`   | 单元素半径覆盖，如 `N:1.5 O:1.6`     |
| `-tv, --type-valence`  | 单元素最大价态覆盖，如 `N:4 O:3`     |
| `-nt, --threads`       | 帧处理线程数                         |
| `--no-rings`           | 关闭环检测                           |
| `--no-reactions`       | 关闭 transfer-flow 分析              |
| `--no-track-reactions` | 关闭反应事件跟踪                     |

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

`plot` 默认一次性生成所有推荐的图表。

## 当前限制

1. 不支持三斜晶胞和逐轴周期性；
2. 与电子结构无关；
3. 仅区分中小分子构造异构体。

## 社区发表论文

截至2026年10月，社区已发表20+篇相关工作，以下为部分内容。

| 体系                  | 论文                                                                                                                               | 链接                                                  |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| CO₂ 捕集胺液热降解    | Thermal degradation and inhibition of long-chain triamines for CO2 capture, `Chem. Eng. J.` 2026                                   | [DOI](https://doi.org/10.1016/j.cej.2026.174329)      |
| 等离子体石墨烯纳米片  | Experimental and ReaxFF MD simulation study on plasma gas-phase synthesis of graphene nanoflakes, `Diamond Relat. Mater.` 2025     | [DOI](https://doi.org/10.1016/j.diamond.2025.113141)  |
| 三相凝胶泡沫灭火材料  | Optimization and performance assessment of an eco-friendly three-phase gel foam for liquid fires, `Colloids Surf. A` 2026          | [DOI](https://doi.org/10.1016/j.colsurfa.2026.141146) |
| 乙醇/Al(111) 表面反应 | Adsorption and thermal decomposition mechanisms of ethanol on Al(111), `J. Anal. Appl. Pyrolysis` 2026                             | [DOI](https://doi.org/10.1016/j.jaap.2026.107932)     |
| 铜表面脉冲电位重构    | Atomistic Origin of Non-equilibrium Surface Restructuring on Copper under Pulsed Potentials, `JACS` 2026                           | [DOI](https://doi.org/10.1021/jacs.6c09411)           |
| 桦甸油页岩干酪根热解  | A ReaxFF molecular dynamics study on the bond cleavage priority in the pyrolysis process of Huadian oil shale kerogen, `Fuel` 2026 | [DOI](https://doi.org/10.1016/j.fuel.2026.141081)     |
| 纤维素蒸汽气化        | Fe2O3-assisted steam gasification of cellulose: a ReaxFF molecular dynamics study, `Process Saf. Environ. Prot.` 2026              | [DOI](https://doi.org/10.1016/j.psep.2026.109259)     |

## 引用

Hanxiang Chen, ReaxTools: A high performace ReaxFF/AIMD trajectory analysis tool based on graph theory. https://github.com/tgraphite/reax_tools [Computer Software], retrived on [date]
或相似表述。