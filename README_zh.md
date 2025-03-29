# ReaxTools 1.1版本

一个高性能的C++分子动力学轨迹分析工具，专注于反应系统分析。

[English Version](README.md)

## 1.1版本更新说明
---------
- 新增反应流程图功能，可以追踪分子间的反应关系
- 优化了分子识别算法，提高了准确性
- 支持.xyz和.lammpstrj文件格式
- 修复了多个稳定性问题

## 特性
-----
- 长期迭代更新，持续完善功能
- 采用C++完全重写，性能极大优化，资源占用低
- 直接读取轨迹功能，支持.lammpstrj文件和.xyz文件
- 保留了清洗Lammps的.species.out输出的功能
- 依赖极少，仅需C++标准库和fmt库

## 使用方法
---------
（输入-h或--help可查看完整帮助） 

-f .xyz/.lammpstrj file -> [读取轨迹] 采用van der Waals半径决定键、建立分子拓扑并分析
-s lammps reaxff/species file (species.out) -> [读取物种文件] 读取species.out文件清洗
-r value 调整van der Waals半径的缩放因子，默认1.2（OVITO的默认设置）
-t element,element... 逗号分割的元素名称, 如C,H,O,N,S,F（读取lammpstrj时必须）
-me element 通过元素的原子数量合并物种类型，例如C1~C4合并为group_C1-C4
-mr range1,range2... 通过原子数量合并物种类型的范围，例如1,4,8,16
-rc 按指定原子的数量重算物种的权重，如C2H4被认为是权重2，而不是1个分子
--order element,element... 整理输出化学式中元素的顺序 (默认: C,H,O,N,S,F,P)

### 示例

处理包含碳、氢、氧和氮原子的XYZ轨迹：

./reax_tools -f trajectory.xyz -t C,H,O,N

使用自定义van der Waals半径缩放处理LAMMPS轨迹：

./reax_tools -f trajectory.lammpstrj -t C,H,O,N -r 1.3

## 输出文件
---------
- .species.csv - 包含分子种类统计和演化数据
- .dot - 反应流程图，可用Graphviz等工具可视化

## 测试输出
```
Frame: 1 Atoms: 12902 Bonds: 10559 Mols: 5173
Frame: 2 Atoms: 12902 Bonds: 10223 Mols: 5180
Frame: 3 Atoms: 12902 Bonds: 10181 Mols: 5243
Frame: 4 Atoms: 12902 Bonds: 10073 Mols: 5284
Frame: 5 Atoms: 12902 Bonds: 10031 Mols: 5326
C1 : 86 74 95 108 104 
C102H8O44 : 0 0 0 1 0 
C104H8O43 : 0 0 1 0 0 
C10H1O3 : 0 0 0 0 1 
C10H1O4 : 0 0 1 0 1 
C10H1O5 : 0 0 0 0 1 
(...)

Save file ../examples/FeCHON_5frames.csv successfully.
=== Reaction Flow Report ===
Total nodes (species): 40
Total edges (reactions): 49

Top reactions:
1: C2H1 -> C1H1 (count: 17)
2: C1H1 -> C2H1 (count: 13)
3: C3H2 -> C2H1 (count: 9)
4: C2H1 -> C3H2 (count: 6)
5: C2O1 -> C1O1 (count: 5)
6: C2O1 -> C2 (count: 4)
7: C1H1O1 -> C2H2O1 (count: 4)
8: H2N1 -> H3N1 (count: 4)
9: C2H2 -> C3H3 (count: 3)
10: C3H3 -> C2H2 (count: 3)
Reaction flow graph saved to ../examples/FeCHON_5frames.dot
```

## 实现细节
---------
- Tick-Tock交替式读取轨迹，节省内存并支持帧间分析
- 使用K-D树算法高效搜索近邻原子
- 基于van der Waals半径判断化学键
- 使用深度优先搜索算法构建分子图
- 基于原子ID交集/并集计算分子相似度，追踪反应关系

## 性能
-----
- 测试用例：350MB的lammpstrj轨迹，1.3万原子，1000帧
- 性能表现：单核处理器用时约160秒，内存占用不到100MB
- 物种文件模式处理速度更快，通常只需几秒钟

## 后期开发计划
-----------
- 反应链：追踪分子从生成到销毁的完整生命周期
- 基团转移：分析分子片段在不同分子间的迁移
- 反应动力学：计算反应速率和反应级数
- 可视化：提供更丰富的图形输出和报表功能

## 从源代码构建

mkdir build
cd build
cmake ..
make -j4

## 许可证

没有许可证，爱咋用咋用。