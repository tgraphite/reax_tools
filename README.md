### 使用reax_tools全自动分析反应轨迹  

----

### reax_tools介绍

reax_tools是C++写的一个反应MD后分析工具，能够分析lammpstrj和xyz格式的轨迹，无论轨迹是ReaxFF、AIMD、NEP产生的。  

其算法实现简要描述如下：   

1. 构建邻居列表
2. 根据vdW半径构建原子间化学键（支持三维周期性或无周期性）
3. 搜索并构建分子拓扑
4. 根据图连通性计算3-8元独立环的数量、原子键连数、物质数量变化
5. 根据前后帧分子间相关性、原子转移关系构建物质转化的网络

能够自动输出的结果包括：   

1. 各帧物种化学式和数量(species_count.csv)
2. 各帧键的类型和数量(bond_count.csv)
3. 各帧环的数量(ring_count.csv)
4. 物质转化网络(reactions.dot)
5. 以关键分子为中心的物质转化子网络(reactions_cenetered_on_*_.dot)和关键分子的流入/流出信息(key_molecule_reactions.csv)

软件能够并行运行，并且对内存管理进行了优化。不需要额外设置，能够适用数百万原子的体系。每1GB lammpstrj轨迹约耗时2-3分钟。推荐的最低配置为2-4核CPU、4GB内存。  

附有绘图python脚本reax_plot.py，易于全自动产生结果。 

----

### 基本使用

#### 1. 从Github页面 (https://github.com/tgraphite/reax_tools/releases) 下载最新版本。无后缀名的是Linux版，后缀名为exe的是Windows版。**作者推荐使用Linux版二进制程序**。
#### 2. 解压，假设解压到`/your/reax_tools/path`
#### 3. 设置PATH和LD_LIBRARY_PATH。
    export PATH=${PATH}:/your/reax_tools/path

#### 4. 最简使用
    reax_tools -f <.xyz/.lammpstrj file> -t <element1,element2...> 

    # 例如：reax_tools -f traj.lammpstrj -t C,H,O,N
    # 对于有元素名标注的.xyz文件，-t选项不是必要的

#### 5. 关于文件格式  

推荐使用以下两种方式产生的输入文件:   

- LAMMPS的lammpstrj文件格式，后缀名为lammpstrj，原子记录中需要包括id type x y z或者id type xs ys zs字段。并且原子类型对应的实际元素名不能重复，例如type 1-4 = C,H,O,N，但是不可以type 1-5 = C,H,H,O,N。使用这种lammpstrj格式的文件需要在使用时加上如-t C,H,O,N标定元素符号。  
- 扩展的xyz格式(extended xyz)，可以通过OVITO导入文件，再以xyz (extended)的方式导出，字段为元素名 X Y Z。CP2K和GPUMD输出的轨迹如traj.pos-1.xyz和dump.xyz可以直接使用。  

**如果初次使用不能产生结果，大部分情况下是文件格式或读取问题，请特别注意文件格式和类型标定**  
**对于如type 1-5 = C,H,H,O,N这种有重复元素名对应的lammpstrj文件，建议用OVITO导入，然后在界面中设定Particle Type对应的Name，这样导出xyz文件时就会自动合并映射到唯一元素名**  

#### 6. 全部选项

    -s <lammps species file> 替代-f，不是读取轨迹，只是清洗lammps reaxff/species的输出。（旧reax_species.py的功能）
    
    -r <vdw scaling factor> vdW半径的缩放因子，默认1.2（与OVITO一致）。一般而言1.2左右通用。调小这个值会使得判断更严格，容易判断出更小的分子、分子碎片，反之亦然。具体看体系尝试。

    --dump 每一帧输出加了bond的lammps data文件，方便用OVITO或者VMD读取绘图。默认不使用。

    -nt 线程数，增大不一定会变快，因为有部分算法无法并行。默认4。

    -me <element> 按照某种元素的含量合并物种群，例如合并成C1-C4、C5-C8等。默认不使用，如果只设置了-mr而忘记设置-me，默认为C。

    -mr <mranges,range,range...> 上一个选项对应的合并的上下界，左闭右开区间，默认1,4,8,16。

    -rc 使用合并功能后对物种群的数量（权重）进行重算，输出原子数而不是分子数。默认不使用。

    --order <element,element,...> 输出的化学式元素顺序，例如把H2CO重写成CH2O，默认C,H,O,N,S,F,P。

    -rr 对于一对互逆的反应，只保留净值。对过于繁杂的反应网络很有用。

----

### 图表型结果文件
reax_tools是纯计算程序，默认的数据绘图功能使用附带的reax_plot.py实现，reax_plot.py的用法如下：  
    
    reax_plot.py -d <目录名>            # 自动绘制reax_tools的各类图表型结果
    reax_plot.py -c <文件1，文件2...>   # 对于类型相同、可用于比对的各个csv文件，自动合并同类项并绘制对比图，例如有system1_bond_count.csv system2_bond_count.csv ...等6个平行体系，抽取可供对比的键类型。species_count等同理。   
    
**reax_plot会根据文件类型和具体数值，自动决定不输出一些不重要的结果，例如全程几乎为0的物质、数量几乎不变的键等等。reax_plot对于字段过多的情况还会自动分页。**  

reax_plot对图表型数据的绘制结果：

- species_count.csv

<img src="doc/species_count.png" width=600> 

95% | Avg | 5% 分别表示95百分位值（相对高值）、平均值、5百分位值（相对低值）

- bond_count.csv

<img src="doc/bond_count.png" width=600>

- ring_count.csv

<img src="doc/ring_count.png" width=600>

R3-R8表示3元-8元环，注意有些体系如含有金属块体的，环数可能没有意义。如果希望忽略金属块体，在初始输入时，-t中把对应元素名改为X，如此这种元素就全程都不会纳入考虑。

- atom_bonded_num_count.csv  

<img src="doc/atom_bonded_num_count.png" width=600>

C(3)表示连接有3个原子的碳，假如在研究煤、干酪根、聚合物等有机体系，那么近似等于sp2碳，其他如O(2)即连了2个原子的氧，H(0)即孤立氢原子。

- compare_*.png

<img src="doc/compare_species_count_NO2.png" width=600>

当你有多个平行体系时，使用如```python reax_plot.py -c system1_species_count.csv ...```产生的结果，本质上就是合并同类项，输入的若干文件类型必须相同，reax_plot.py会自动识别有比较价值、并且在每个文件中都有出现的字段。

#### 网络型结果文件

像工作流程图一类的顶点和边构成的结构，为图论意义上的图（Graph），为了消歧义也称为网络（Network）。描述一个网络本质上只需要给出顶点的列表和边的列表，dot文件就是一种常用的记录网络的格式。

输出的dot文件推荐使用graphviz可视化，graphviz可以直接用apt或yum安装(如`apt install graphviz`)，然后用`dot -Tpng reactions.dot -o reactions.png`绘制。

批处理时，可以使用如以下命令自动绘制当前目录下所有dot文件。

    for file in *.dot;do
        dot -Tpng $file -o ${file/.dot/.png}
    done
 
需要注意对于特别复杂的反应网络，dot文件过大时绘制时间可能会长达几十分钟，这时reax_tools会给出reactions.dot和reactions_full.dot两个文件，前者是后者的子网络（反应数最多的那部分），你可以只绘制前者。

- reactions.dot

<img src="doc/reactions.png" width=600>

其中如5:10表示这个路径有效反应了5次，向下游转移了10个原子。最频繁的一些反应会自动用黄色高亮。

额外说明：

1. 模拟时往往并不会像我们预想的那样，如C4H8 -> C2H4 + C2H4这样能够完美配平，从轨迹分析而言，我们只能知道C4H8有时向C2H4转移了2个原子、有时转移了4个原子等。所以这里将反应次数和转移的总原子量都显示出来。

2. 很多时候都是互逆的转化关系，如H2O <-> OH，当使用-rr选项时，其中一个数值较小的会被消除，剩下那个减去它的反应数和转移的总原子量，是否要加这个选项，根据情况自己把握。

3. 题外话：如上图，初始物质C4H6N4O11可以通过4条渠道转化成终端物质NO2，其中有最为便捷的直接转化，也有通过两种中间体的渠道。可以手动计算初始物质到终端物质的**最大网络流和反应比例**，以后可能加入自动识别和计算的功能。

- reactions_centered_on_*.dot

<img src="doc/reactions_centered_on_C4H6N2O7.png" width=600>

自动识别的关键分子的上下游关系，相当于一个子网络。标注方式同上。

此外，key_molecules_reactions.csv还会输出当前体系下关键分子的流入（流出）反应数、流入（流出）原子数、5个最大来源和5个最大去路。相当于以下结果。

<img src="doc/compare_HO_atom_transfer.png" height=400>

### 提交需求和报告bug

QQ群（561184358）

### 使用许可证和题外话

MIT许可证，想用就用。  

对于代算的朋友：拿去做单子随意，别直接卖软件（特别是闲鱼的，没出息）。