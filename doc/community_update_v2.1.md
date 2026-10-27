# ReaxTools v2.1 更新说明

ReaxTools v2.1 是一次比较大的结构重构版本，核心目标是让程序输出更稳定、更干净、更可审计。

## 主要变化

1. 统一入口

现在用户只需要使用一个命令：

```bash
reax_tools -f traj.xyz -o output_dir
reax_tools plot -f output_dir
reax_tools network -f output_dir
reax_tools events -f output_dir
reax_tools counts -f output_dir
```

直接 `reax_tools -f ...` 等价于分析模式；`plot/network/events/counts/flow/focus` 等子命令负责后处理和画图。

2. C++ 和 Python 职责分离

C++ 只负责高性能分析和可审计 raw 输出：

- `species_count.csv`
- `bond_count.csv`
- `atom_bonded_num_count.csv`
- `ring_count.csv`
- `reaction_events.csv`
- `reaction_event_pairs.csv`
- `transfer_flow.csv`
- `molecules.json`
- `reax_tools.log`
- `reax_tools_manifest.json`

Python 负责过滤、画图和展示。这样以后大多数用户体验和绘图改动都可以在 Python 层快速更新，不需要频繁改 C++ 核心。

3. 移除不必要外部依赖

C++ 核心构建不再依赖 RDKit、Boost 或 Graphviz。完整 Python 绘图环境仍安装 PyGraphviz 和 PyVista/VTK，网络布局需要 Graphviz，不能把核心依赖与完整安装依赖混为一谈。

4. 反应网络审计增强

v2.1 新增 `reaction_event_pairs.csv`，它记录每个 raw reaction event 内 reactant-product 的原子重叠 pair。

`transfer_flow.csv` 是 `reaction_event_pairs.csv` 的逐边聚合结果。也就是说，完整 network 不再只是“画出来的图”，而是可以从 raw reaction events 审计出来的网络。

我们增加了审计工具：

```bash
python3 tools/check_reax_outputs.py output_dir
python3 tools/audit_event_network_consistency.py output_dir
```

2026-09-11 复核：`energetic_v3` 和本地 `gpumd_v3` 的 pair/flow 聚合审计均通过；`energetic_v3` 的旧事件表头不兼容当前 `check_reax_outputs.py`，结构检查请使用 `energetic_v4` 或新生成结果。

5. network 和 flow 的定位更清楚

- `transfer_flow.csv`：从已检测事件的原子重叠 pair 聚合得到的完整转移表。
- `network`：该转移表的过滤视图；默认 top-N 图不包含所有边。
- `flow`：实验性 Sankey 风格叙事图，用来帮助快速看主线，但会压缩和丢弃部分循环/同层信息。

正式统计以 raw CSV 为准，并记录过滤参数。聚合一致不等于反应识别准确：当前物种按分子式区分，不能区分异构体；`--stable-time` 尚未接入原始事件筛选。详见[方法与限制](methods-and-limits.md)。

6. README 和输出风格更新

新版 README 以 `test/energetic_v3` 为样例，展示 species、reaction events、transfer network 和 experimental flow。命令行输出也比旧版更克制，详细信息放到 log 和 manifest 中。

## 推荐安装

```bash
git clone https://github.com/tgraphite/reax_tools
cd reax_tools
bash install_reax_tools.sh
```

## 推荐使用

```bash
reax_tools -f trajectory.xyz -o reax_tools_output
reax_tools plot -f reax_tools_output
reax_tools network -f reax_tools_output
```

LAMMPS dump 文件通常需要指定元素类型：

```bash
reax_tools -f dump.lammpstrj -t C,H,O,N -o reax_tools_output
```

## 当前适用范围（2026-09-11 核查）

ReaxTools v2.1 提供事件—原子重叠—物种转移表的审计链。它适合检查统计来源、探索轨迹中的组分变化；化学准确性还需对成键参数、采样间隔和代表性事件分别验证。本说明保留版本发布背景；当前输入限制和安装要求以 README 为准。
