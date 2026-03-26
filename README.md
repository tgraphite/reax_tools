# ReaxTools: 高性能反应分子动力学后处理工具

<p align="center">
  <b>高性能 | 多格式支持 | 反应网络分析 | 在线即用</b>
</p>

<p align="center">
  <a href="https://github.com/tgraphite/reax_tools/releases"><img src="https://img.shields.io/badge/version-1.98-blue" alt="Version"></a>
  <a href="http://cc-portal.xyz/reax_tools_pro"><img src="https://img.shields.io/badge/web-在线试用-green" alt="Web"></a>
</p>

ReaxTools 是一款专为反应分子动力学（ReaxFF/AIMD/MLP-MD）模拟数据设计的高性能后处理分析工具，支持 LAMMPS、GPUMD、CP2K 等主流分子动力学软件产生的轨迹文件。

---

## 核心功能

### 物质与结构分析
- **物种统计** - 追踪分子种类及数量随时间的演化
- **化学键统计** - 分析各类化学键的形成与断裂
- **原子成键数统计** - 统计原子配位数变化（如 C(3) sp²、C(2) sp 等）
- **环结构检测** - 识别并统计环状分子的数量与类型

### 反应网络分析
- **物质转化网络** - 构建物质转化关系图（Graph/Network）
- **反应路径追踪** - 分析分子间的反应路径与转化关系
- **原子转移分析** - 追踪特定原子/官能团在反应中的流动
- **关键分子识别** - 基于反应通量识别关键中间体

### 输出与可视化
- **SMILES 表示** - 输出分子的 SMILES 字符串
- **分子结构绘制** - 自动生成分子结构图像
- **反应网络可视化** - 生成 DOT 格式网络图（可用 Graphviz 渲染）
- **CSV 数据导出** - 所有统计结果导出为 CSV 格式便于后续处理

---

## 使用方式

### 方式一：Web 在线版 (体验)

**访问地址**: http://cc-portal.xyz/reax_tools

在线版本无需安装，即开即用，支持：
- 完整的分析功能
- 自动数据可视化
- 交互式反应网络查看
- 一键导出所有结果

**推荐配置**: Chrome 浏览器，PC 端，轨迹大小 < 10,000 原子 × 200 帧

### 方式二：本地版本（正式计算工作）

#### 安装步骤（无需编译）

```bash
# 1. 克隆仓库
git clone https://github.com/tgraphite/reax_tools
cd reax_tools/bin

# 2. 运行安装脚本
bash install_reax_tools.sh
```

**推荐配置**: Linux, 8 核, 轨迹大小 < 4,000,000 原子 × 10,000 帧

### 方式三：手动编译（非必须）

#### 编译依赖

- **编译器**: GCC (支持 C++17) 或 Emscripten (Web 版)
- **依赖库**: 
  - RDKit (分子处理与绘图)
  - Boost 1.88.0+
  - fmtlib (格式化输出)

#### 编译命令

```bash
# 本地版本
bash build_local.sh

# Web 版本 (需要 Emscripten)
bash build_wasm.sh
```

---

## 输入文件格式

ReaxTools 支持以下三种格式的轨迹文件：

### 1. LAMMPS dump 文件 (`.lammpstrj`)
```
ITEM: TIMESTEP
1000
ITEM: NUMBER OF ATOMS
10000
ITEM: BOX BOUNDS pp pp pp
0.0 40.0 0.0 40.0 0.0 40.0
ITEM: ATOMS id type x y z
1 1 0.000 0.000 0.000
2 2 1.000 1.000 1.000
...
```
**注意**: 需配合 `-t C,H,O,N` 参数指定 type_id 对应的元素符号。

### 2. Extended XYZ (推荐)
```
10000
Lattice="40.000 0.0 0.0 0.0 40.000 0.0 0.0 0.0 40.000" Properties=species:S:1:pos:R:3
1 C 0.000 0.000 0.000
2 H 1.000 1.000 1.000
...
```
GPUMD、CP2K 直接输出，或 OVITO 转格式获得。已包含元素符号，无需 `-t` 参数。

### 3. 普通 XYZ
```
10000
Any comment
C 0.000 0.000 0.000
H 1.000 1.000 1.000
...
```
适用于非周期性体系（如 ORCA、CP2K 的 AIMD）。

**格式要求**: 各帧原子顺序需保持一致（用于追踪原子）。

---

## 命令行参数

### 基础参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `-t` | 定义元素符号（用于 lammpstrj） | `-t C,H,O,N` |
| `-r` | vdW 半径缩放因子（默认 1.2） | `-r 1.0` |
| `-tr` | 设置特定元素半径 | `-tr N:1.5` |
| `-tv` | 设置特定元素化合价 | `-tv N:4,S:6` |
| `--order` | 化学式元素排序 | `--order C,H,O,N` |

### 高级分析参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `-me` | 按元素含量归纳 | `-me C` |
| `-mr` | 归纳范围（与 `-me` 配合） | `-mr 1,4,8,16` |
| `-rc` | 输出原子数而非分子数 | `-rc` |
| `-norr` | 保留互逆反应（不抵消） | `-norr` |
| `--max-reactions` | 最大反应数（默认 60） | `--max-reactions 100` |

### 自定义元素类型（原子追踪）

通过同时定义 `-tr` 和 `-tv` 可创建新元素类型，实现原子标记追踪：
```bash
# 创建氘(D)和氚(T)追踪氢同位素
-tr D:1.05 T:1.05 -tv D:1 T:1
```

---

## 使用示例

### 基础分析
```bash
reax_tools -t C,H,O,N trajectory.lammpstrj
```

### 严格成键判据（适合燃烧体系）
```bash
reax_tools -t C,H,O -r 1.0 dump.xyz
```

### 碳数分组统计
```bash
reax_tools -t C,H,O,N -me C -mr 1,4,8,16 -rc trajectory.xyz
```
输出: C1-C3、C4-C7、C8-C15、C16+ 四组的碳原子数随时间变化。

### 网络分析（保留所有反应）
```bash
reax_tools -t C,H,O,N --max-reactions 200 -norr trajectory.xyz
```

---

## 输出文件说明

运行后会在 `reax_tools_output/` 目录生成：

| 文件名 | 内容 |
|--------|------|
| `species_count.csv` | 各物种分子数随时间变化 |
| `species_count_hash.csv` | 物种哈希值映射（唯一标识） |
| `bond_count.csv` | 各类化学键数量统计 |
| `atom_bonded_num_count.csv` | 原子成键数统计 |
| `ring_count.csv` | 各类环结构数量统计 |
| `reactions.dot` | 主要反应网络（DOT 格式） |
| `key_molecule_reactions.md` |  |

**可视化 DOT 文件**:
```bash
dot -Tpng reactions.dot -o reactions.png
```

---

## 性能参考

### 本地版本（WSL, i5-13500H, 4 线程）

| 体系规模 | 原子数 | 帧数 | 文件大小 | 耗时 |
|----------|--------|------|----------|------|
| 2.5K | 2,592 | 40 | 2.7 MB | 1.3 s |
| 10K | 10,368 | 40 | 13 MB | 2.5 s |
| 50K | 51,840 | 40 | 76 MB | 9.4 s |
| 160K | 165,888 | 40 | 246 MB | 24.5 s |
| 1M | 995,328 | 40 | 1.6 GB | 147.6 s |

---

## 学习资源

- **详细教程**: http://cc-portal.xyz/reax_tools/tutorial.html
- **参数说明**: http://cc-portal.xyz/reax_tools/inputargs.html
- **文件格式**: http://cc-portal.xyz/reax_tools/fileformat.html
- **常见问题**: http://cc-portal.xyz/reax_tools/faqandhint.html

---

## 引用

如果您在研究中使用了 ReaxTools，请引用：

```
Hanxiang Chen. ReaxTools: A high performance ReaxFF/AIMD/MLP-MD 
post-process code [Computer software]. 
https://github.com/tgraphite/reax_tools
```

---

## 项目统计

- ✅ 30+ 个验证项目和论文
- ✅ 10+ 篇社区发表论文
- ✅ 200+ 社区成员（持续增长中）

---

## 联系方式

- **作者**: Graphite (Hanxiang Chen)
- **邮箱**: chenhx_cpu@163.com
- **GitHub**: https://github.com/tgraphite/reax_tools

---

## 示例图像

### 物种统计
![物种统计](https://github.com/user-attachments/assets/01b5aa96-2839-469b-9daa-fb76e4515301)

### 环检测
![环检测](https://github.com/user-attachments/assets/e06cf000-80a2-4181-9ab8-492cd9663b97)

### 原子成键数统计
![原子成键数](https://github.com/user-attachments/assets/99dcab57-4596-4977-824b-6588975111ce)

### 化学键统计
![化学键统计](https://github.com/user-attachments/assets/0a033448-6cb6-4d3e-bccc-7cf62a125646)

### 原子转移分析
![原子转移](https://github.com/user-attachments/assets/0260a588-4a39-4efc-b375-bad7ded1bb20)

### 反应网络
![反应网络](https://github.com/user-attachments/assets/ed571f01-c24d-403c-aac2-86876443ffaa)

### 环数量统计
![环数量](https://github.com/user-attachments/assets/892d0f0c-2ea6-464f-b467-23ba884dc984)
