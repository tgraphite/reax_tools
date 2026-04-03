# ReaxTools: 高性能反应分子动力学后处理工具

> **ReaxTools** is a high-performance post-processing tool for reactive molecular dynamics (ReaxFF / AIMD / MLP-MD). It reads trajectory files from LAMMPS, GPUMD, CP2K, and other MD engines, and produces species statistics, bond analysis, ring detection, reaction network graphs, and reaction event tracking. Available as a CLI binary, a WebAssembly web app, and Python visualization scripts.

<p align="center">
  <b>高性能 | 多格式支持 | 反应网络分析 | 在线即用</b>
</p>

<p align="center">
  <a href="https://github.com/tgraphite/reax_tools/releases"><img src="https://img.shields.io/badge/version-2.0-blue" alt="Version"></a>
  <a href="http://cc-portal.xyz/reax_tools_pro"><img src="https://img.shields.io/badge/web-在线试用-green" alt="Web"></a>
</p>

ReaxTools 是一款专为反应分子动力学（ReaxFF/AIMD/MLP-MD）模拟数据设计的高性能后处理分析工具，支持 LAMMPS、GPUMD、CP2K 等主流分子动力学软件产生的轨迹文件。

---

## 核心功能

### 物质与结构分析
- **物种统计** - 追踪分子种类及数量随时间的演化
- **化学键统计** - 分析各类化学键的形成与断裂
- **原子成键数统计** - 统计原子配位数变化（如 C(3) sp2、C(2) sp 等）
- **环结构检测** - 识别并统计环状分子的数量与类型

### 物质转化网络分析
- **物质转化网络** - 构建物质转化关系图（Graph/Network）
- **智能过滤** - 基于原子经济性（Atom Economy）自动过滤低质量反应，基于边介数中心性（Betweenness Centrality）排序保留关键反应
- **转化路径追踪** - 自动识别并输出最长连续反应路径
- **原子转移分析** - 追踪特定原子/官能团在反应中的流动
- **关键分子识别** - 基于拓扑度识别关键中间体，输出前驱体/衍生物关系

### 反应事件跟踪分析
- **反应事件检测** - 基于原子连接关系变化自动识别反应事件
- **反应频率统计** - 统计各类型反应的发生频次（Top 20 反应可视化）
- **分子生命周期追踪** - 追踪单个分子的生成、转化与消亡过程
- **前驱体/后继体分析** - 建立分子间的反应前驱与后继关系

### 输出与可视化
- **SMILES 表示** - 输出分子的 SMILES 字符串（基于 RDKit）
- **分子结构绘制** - 自动生成分子结构 SVG 图像
- **反应网络可视化** - 生成 DOT 格式网络图（可用 Graphviz 渲染）
- **CSV/Markdown 数据导出** - 统计结果导出为 CSV，关键分子报告导出为 Markdown
- **Python 可视化脚本** - 附带 `reax_plot.py` 和 `reax_network.py` 自动生成出版级图表

---

## 使用方式

### 方式一：Web 在线版（体验）

**访问地址**: http://cc-portal.xyz/reax_tools

在线版本无需安装，即开即用，支持：
- 完整的分析功能
- 自动数据可视化
- 交互式反应网络查看
- 一键导出所有结果

**推荐配置**: Chrome 浏览器，PC 端，轨迹大小 < 10,000 原子 x 200 帧

### 方式二：本地版本（正式计算工作）

#### 安装步骤（无需编译）

```bash
# 1. 克隆仓库
git clone https://github.com/tgraphite/reax_tools
cd reax_tools/bin

# 2. 运行安装脚本
bash install_reax_tools.sh
```

**推荐配置**: Linux, 8 核, 轨迹大小 < 4,000,000 原子 x 10,000 帧

### 方式三：自建 Web 服务器 (Web 2.0 新版)

使用服务器端执行 + 任务队列的版本，支持更大文件和并发任务：

```bash
# 1. 克隆仓库
git clone https://github.com/tgraphite/reax_tools
cd reax_tools

# 2. Docker Compose 一键启动
docker-compose up -d

# 3. 访问 http://localhost
```

**特性**:
- 任务队列系统，支持多任务并发
- 实时进度推送 (SSE)
- 支持更大文件 (最高 500MB)
- 现代化数据大屏界面 (Tailwind CSS)
- 响应式设计，支持深色模式

详见 [server/README.md](server/README.md)

### 方式四：手动编译（非必须）

#### 编译依赖

- **编译器**: GCC (支持 C++17) 或 Emscripten (Web 版)
- **依赖库**:
  - [RDKit](https://www.rdkit.org/) (分子处理与绘图)
  - [Boost](https://www.boost.org/) 1.88.0+
  - [fmtlib](https://fmt.dev/) (格式化输出，已内置)

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
0.0 40.0
0.0 40.0
0.0 40.0
ITEM: ATOMS id type x y z
1 1 0.000 0.000 0.000
2 2 1.000 1.000 1.000
...
```
**注意**: 需配合 `-t C,H,O,N` 参数指定 type_id 对应的元素符号。

### 2. Extended XYZ（推荐）
```
10000
Lattice="40.000 0.0 0.0 0.0 40.000 0.0 0.0 0.0 40.000" Properties=species:S:1:pos:R:3
C 0.000 0.000 0.000
H 1.000 1.000 1.000
...
```
GPUMD、CP2K 直接输出，或**通过 OVITO 转格式获得**。已包含元素符号，无需 `-t` 参数。
**重要**如果因为非标准的输出而不能读取其他格式，优先尝试OVITO转格式。

**LAMMPSTRJ转标准Extended XYZ的方法**
```
1. 在OVITO中正确载入Lammpstrj文件。
2. 在Particle Type中，手动为type 1...type N的原子设置正确的元素名。
3. 如果有PBC问题，可以用Affine Transformation重新定义盒子，但是要注意原子坐标是否需要被重新缩放。
4. 导出为Extended XYZ文件。
```

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
| `-f` | 文件名 | `-f traj.xyz/traj.lammpstrj` |
| `-t` | 定义元素符号（用于 lammpstrj） | `-t C,H,O,N` |
| `-r` | vdW 半径缩放因子（默认 1.2） | `-r 1.0` |
| `-tr` | 设置特定元素半径 | `-tr N:1.5 O:1.5` |
| `-tv` | 设置特定元素化合价 | `-tv N:4 S:6` |
| `--order` | 化学式元素排序 | `--order C,H,O,N` |

### 高级分析参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `-me` | 按元素含量归纳 | `-me C` |
| `-mr` | 归纳范围（与 `-me` 配合） | `-mr 1,4,8,16` |
| `-rc` | 输出原子数而非分子数 | `-rc` |
| `-norr` | 保留互逆反应（不抵消） | `-norr` |
| `--max-reactions` | 最大反应数（默认 60） | `--max-reactions 100` |
| `--no-track-reactions` | 禁用反应跟踪（默认启用） | `--no-track-reactions` |
| `--stable-time` | 反应稳定时间阈值（帧，默认 3） | `--stable-time 5` |

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
reax_tools  -f trajectory.lammpstrj -t C,H,O,N
```

### 严格成键判据
```bash
reax_tools -f dump.xyz -t C,H,O -r 1.0 
```

### 碳数分组统计
```bash
reax_tools -f trajectory.xyz -t C,H,O,N -me C -mr 1,4,8,16 -rc 
```
输出: C1-C3、C4-C7、C8-C15、C16+ 四组的碳原子数随时间变化。

---

## 输出文件说明

运行后会在 `reax_tools_output/` 目录生成：

| 文件名 | 内容 |
|--------|------|
| `species_count.csv` | 各物种分子数随时间变化 |
| `bond_count.csv` | 各类化学键数量统计 |
| `atom_bonded_num_count.csv` | 原子成键数统计 |
| `ring_count.csv` | 各类环结构数量统计 |
| `reaction_flow_full.dot` | 完整反应流网络图（DOT 格式） |
| `reaction_flow_simplified.dot` | 简化版反应流网络图（仅在反应数过多时生成） |
| `key_molecules_flow.md` | 关键分子反应流汇总报告（Markdown） |
| `molecules_smiles.csv` | 分子 SMILES、拓扑度信息 |
| `molecule_pictures/*.svg` | 分子结构图（SVG 格式） |
| `reaction_track_events.csv` | 反应事件频率统计 |
| `reaction_track_molecules.csv` | 分子生命周期追踪数据 |

**可视化 DOT 文件**:
```bash
dot -Tpng reaction_flow_full.dot -o reaction_flow.png
```

### Python 可视化工具

`bin/` 目录下附带 Python 脚本，可自动生成出版级图表：

```bash
# 安装依赖
pip install -r bin/requirements.txt

# 绘制物种、键、环等统计图
python bin/reax_plot.py -d reax_tools_output/

# 生成高清反应网络图（需要安装graphviz以获取dot命令）
dot -Tpng reax_tools_output/reaction_flow_full.dot -o reaction_flow_full.png
```

---

## 项目结构

```
reax_tools/
├── cpp/                    # C++ 源代码
│   ├── main.cpp            # 入口
│   ├── universe.cpp/h      # 轨迹处理主循环
│   ├── system.cpp/h        # 单帧体系：成键、找分子、反应检测
│   ├── reax_flow.cpp/h     # 反应网络图（含过滤算法）
│   ├── reax_counter.cpp/h  # 物种统计与归纳
│   ├── reax_track.cpp/h    # 反应事件跟踪
│   ├── rdkit_utils.cpp/h   # RDKit SMILES/绘图接口
│   ├── argparser.cpp/h     # 命令行参数解析
│   ├── cell_list.cpp/h     # 近邻搜索（Cell List）
│   ├── string_tools.cpp/h  # 化学式解析工具
│   └── fmt/                # fmtlib（内置）
├── bin/                    # 预编译二进制 & Python 工具
│   ├── reax_tools           # Linux 二进制
│   ├── reax_plot.py         # 统计图绘制
│   ├── reax_network.py      # 反应网络可视化
│   └── install_reax_tools.sh
├── web/                    # WebAssembly 在线版 (旧版)
├── web-new/                # Web 2.0 - 服务器端版本 (新版)
│   ├── index.html          # 主页面
│   ├── src/                # JS 源码
│   └── dist/               # 编译后的 CSS
├── server/                 # Web 2.0 后端服务器
│   ├── src/                # Node.js 源码
│   └── package.json
├── doc/                    # 开发文档 & 变更记录
├── test/                   # 测试数据
├── CMakeLists.txt
├── build_local.sh
├── build_wasm.sh
└── docker-compose.yml      # Docker 部署配置
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
- **变更记录**: [doc/changelog.md](doc/changelog.md)

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

- 30+ 个验证项目和论文
- 10+ 篇社区发表论文
- 200+ 社区成员（持续增长中）

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

### 环数量统计
![环数量](https://github.com/user-attachments/assets/892d0f0c-2ea6-464f-b467-23ba884dc984)

### 物质转化网络
![物质转化网络](https://github.com/user-attachments/assets/ed571f01-c24d-403c-aac2-86876443ffaa)

### 反应跟踪分析
![反应跟踪分析](https://github.com/user-attachments/assets/e1959e5e-8b18-41bd-8d45-75a9156b050b)
