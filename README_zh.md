# ReaxTools 1.1版本
[English Version](README.md)  
用于ReaxFF/AIMD模拟轨迹的后处理，产生清洗过的物种演变和反应关系信息。  
永久更新地址：https://github.com/tgraphite/reax_tools
下载地址：https://github.com/tgraphite/reax_tools/releases
如果对你有用，不妨在github上点个star，这对我搞项目很重要。谢谢！


## 1.1版本更新说明
---------
- 代码全部开源
- 新增反应流程图功能，可以追踪分子间的反应关系
- 在Cursor的帮助下，懒惰的作者已经恢复更新

## 特性
-----
- 极简工具，不需要任何依赖或安装过程
- 更好的性能，能够处理百万原子体系，低压笔记本1 GB轨迹耗时5-10分钟以内，内存消耗200 MB以内
- 直接读取轨迹功能，支持.lammpstrj文件和.xyz文件（如CP2K产生的pos-1.xyz等轨迹也可以用）
- 清洗Lammps的.species.out输出的功能
- 自定义键合半径缩放因子、自定义重算物种权重、整理分子组、元素顺序等实用功能

## 使用方法
---------
（输入-h可显示此帮助） 
```
-f <.xyz/.lammpstrj file>
// [读取轨迹] 采用van der Waals半径决定键、建立分子拓扑并分析

-s <lammps species.out file> 
// (Speices) [读取物种文件] 读取文件并清洗

-r <value>
// 调整van der Waals半径的缩放因子，默认1.2（同OVITO）
// 增大这个数值会使使得原子间更容易被判定成键、分子更大。反之亦然。

-t <element,element...>
// 逗号分割的元素名称, 如C,H,O,N,S,F（读取lammpstrj时必须）
// 如果某个元素你不想统计，例如固定的Fe基底，可以把它的元素符号设为X。

-me <element>
// 通过元素的原子数量合并物种类型，例如C1~C4合并为group_C1-C4。选项名是merge-element的缩写。

-mr <range1,range2...>
// 通过原子数量合并物种类型的范围，例如1,4,8,16。选项名是merge-range的缩写。
// e.g. -mr 1,4,8,16

-rc 
// 按-me指定原子的数量重算物种的权重，如C2H4被认为是权重2，而不是1个分子。选项名是recalc的缩写。

--order <element,element...>
// 整理输出化学式中元素的顺序 (默认: C,H,O,N,S,F,P)
// e.g. --order C,H,O,N
```
### 示例

处理包含碳、氢、氧和氮原子的XYZ轨迹：

```
./reax_tools -f trajectory.xyz -me C -rc -mr 1,4,6,8
./reax_tools -f trajectory.lammpstrj -t C,H,O,N -r 1.3
```

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

## 后期开发计划（投票）
-----------
- 反应网络：完整的反应网络展示。
- 基团转移：分析分子片段在不同分子间的迁移。（做出准确度比较困难）  
- 可视化：可以附带示例python一键绘图代码。

## 从源代码构建
------
build目录下已经有了编译好的文件，也可以直接从release中下载。  
如果需要从源码构建，和其他所有cmake编译的程序一样。  
```
mkdir build
cd build
cmake ..
cmake --build . -j8
```

## 许可证

MIT许可证，或者说爱咋用咋用，作者不会上门找你，也别指望我修bug多勤快。

## 旧版本下载次数统计区  

1.0版本[attach]98527[/attach]
0.3版本[attach]89685[/attach]
