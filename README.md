### 使用reax_tools全自动分析反应轨迹  

[English software description (old)](README_en.md)

----

#### reax_tools介绍
reax_tools是C++写的一个反应MD后分析工具，能够分析lammpstrj和xyz格式的轨迹，无论轨迹是ReaxFF、AIMD、NEP产生的。  

其算法实现简要描述如下：   
1. 构建邻居列表
2. 根据vdW半径构建原子间化学键（支持周期性），虽然很LOW但是意外的很好用
3. 深度优先搜索构建分子拓扑
4. 根据图连通性计算3-8元独立环的数量
5. 根据前后帧分子间Jaccard相似性构建物质转化的图（Graph）

能够自动输出的结果包括：   
1. 各帧物种化学式和数量
2. 各帧键的类型和数量
3. 各帧环的数量
4. 物种转化的图
5. 各分子的结构图

软件能够并行运行，并且对内存管理进行了优化。不需要额外设置，能够适用数百万原子的体系，即使以低配手机或平板电脑的性能也可以运行（另一个帖子已进行过测试）。每1GB lammpstrj轨迹约耗时2-3分钟。  
附有绘图python脚本，易于全自动产生结果。 

#### 基本使用

1. 从Github页面(https://github.com/tgraphite/reax_tools/releases)下载最新linux版本。(对于Windows用户，更新会更慢，使用方法一致，只不过是使用cmd/powershell而不是linux shell，环境变量请自行Google。或者使用WSL。)
2. 解压，假设解压到`/your/reax_tools/path`
3. 设置PATH和LD_LIBRARY_PATH。
    ```
    export PATH=${PATH}:/your/reax_tools/path
    export LD_LIBRARAY_PATH=${LD_LIBRARAY_PATH}:/your/reax_tools/path/lib
    ```
4. 按以下方式使用
    ```
    reax_tools -f <.xyz/.lammpstrj file> -t <element1,element2...>

    # 例如：reax_tools -f traj.lammpstrj -t C,H,O,N
    # 对于有元素名标注的.xyz文件，-t选项不是必要的
    ```
    ```
    其他选项：
    -s <lammps species file> 替代-f，不是读取轨迹，只是清洗lammps reaxff/species的输出。
    
    -r <vdw scaling factor> vdW半径的缩放因子，默认1.2（与OVITO一致）。一般而言1.2左右通用。调小这个值会使得判断更严格，容易判断出更小的分子、分子碎片，反之亦然。具体看体系尝试。

    --dump 每一帧输出加了bond的lammps data文件，方便用OVITO或者VMD读取绘图。默认不使用。

    -nt 线程数，增大不一定会变快，因为有部分算法无法并行。默认4。

    -me <element> 按照某种元素的含量合并物种群，例如合并成C1-C4、C5-C8等。默认不使用，如果只设置了-mr而忘记设置-me，默认为C。

    -mr <mranges,range,range...> 上一个选项对应的合并的上下界，左闭右开区间，默认1,4,8,16。

    -rc 使用合并功能后对物种群的数量（权重）进行重算，输出原子数而不是分子数。默认不使用。

    --order <element,element,...> 输出的化学式元素顺序，例如把H2CO重写成CH2O，默认C,H,O,N,S,F,P。
    ```

5. 结果输出到```name_reax_tools```新目录中，文件包括：
    ```
    species_count.csv # 物种（物种群）数量
    bond_count.csv    # 键数量
    ring_count.csv    # 环数量
    reaction_flow.dot # 主要物质（默认转化频率最高的30个）直接的物质转化图
    molecule_*_*.svg  # 分子结构图（默认转化频率最高的30个）
    ```

6. 三个csv文件可以通过提供的python脚本自动绘图(使用matplotlib和pandas)：`python autoplot -d <dir> -m <max_species>`

7. dot文件为图标注文件，先安装graphviz(`apt install graphviz`)，然后用`dot -Tpng reaction_flow.dot -o reaction_flow.png`自动绘图

#### 结果示例

以下结果全部为自动绘图，没有人工修改。

1. 分子结构图

  1 | 2 | 3 (自由基)  
---- | ---- | ----  
<img src=materials/molecule_1.svg> | <img src=materials/molecule_2.svg> | <img src=materials/molecule_3.svg>  

2. 物种数量

<img src=materials/species_count.png width=600>

3. 键数量

<img src=materials/bond_count.png width=600>

4. 反应的图

<img src=materials/reaction_flow.png width=600>

#### 提交需求和报告bug

1. Github(https://github.com/tgraphite/reax_tools/releases)
2. QQ群（561184358）

#### 使用许可证和题外话

MIT许可证，想用就用。  

对于代算的朋友：拿去做单子随意，别直接卖软件（特别是闲鱼的，没出息）。