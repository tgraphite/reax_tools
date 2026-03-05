## ReaxTools输出文件说明

### ReaxTools介绍
ReaxTools是一个高性能的C++反应动力学模拟后处理软件，输出各类纯文本分析结果。
此说明用于指导编写其对应的可视化组件，基于Python输出单一HTML文件（嵌入js和css），从而为用户提供可交互、易读的报告。

### 术语介绍
- hash: 指体系中分子的唯一hash标识符，一个化学式或者一个物质组可能对应多种分子（即同分异构现象），但hash是唯一、最小的对象单元。
- smiles: 一种分子结构的线性表示，hash和smiles是对应的，但是注意有些结构比较特殊的分子，无法直接用smiles绘制出图片，对于这种分子，如果需要使用图片，只标注其化学式即可。
- species_count、hash_count、bond_count、atom_bonded_num_count、ring_count: 就是对象在不同时间点的计数，用于绘制各类曲线图和进行比较。
- reactions*.dot: 模拟过程中的反应网络，dot(Graphviz)格式。
- key_molecules_reactions(kmr)信息: 对于高权重（流入流出数）的一些分子，其流入流出数信息、关键上下游分子信息。

### 文件树
basedir/ (主输出目录可更改，默认为reax_tools_output/)
├── atom_bonded_num_count.csv       (原子类型-数量随时间变化)
├── bond_count.csv                  (键类型-数量随时间变化)
├── key_molecules_reactions.csv          (关键物质上下游文件，以化学式标注)
├── key_molecules_reactions_hash.csv     (关键物质上下游文件，以hash标注)
├── molecule_pictures
│   ├── 1100911419.svg              (分子图片，以hash命名)
│   ├── 1276603345.svg
│   ├── 1292406354.svg
│   ├── ....svg
├── molecules_smiles.csv            (分子hash-化学式-SMILES对照表)
├── reactions.dot                   (反应网络默认图)
├── reactions_full.dot              (反应网络全图，如果有这个文件，那么reactions.dot是简图，否则reactions.dot就是唯一全图)
├── ring_count.csv                  (环尺寸-数量随时间变化)
├── species_count.csv               (物质含量随时间变化，以化学式标注)
└── species_count_hash.csv          (物质含量随时间变化，以hash标注)

### 文件头示例
==> atom_bonded_num_count.csv <==
H-0,H-1,N-0,N-1,N-2,N-3,N-4,O-0,O-1,O-2
0,864,0,3,229,508,124,0,781,83
1,863,0,40,265,519,40,0,704,160
1,863,0,55,268,493,48,0,695,169

==> bond_count.csv <==
H-H,N-H,N-N,O-H,O-N,O-O
0,772,430,92,849,3
0,659,407,204,814,3
0,652,397,211,816,3

==> key_molecules_reactions.csv <==
molecule,total reactions,in reaction,out reaction,total atom transfer,in atom transfer,out atom transfer,from 1,from 2,from 3,from 4,from 5,to 1,to 2,to 3,to 4,to 5,
H3N,195,11,30,331,228,103,,H4N,H2N,H3N4O3,HN3O2,HN3O6,HN3O4,H4NO,H2O,HN3O5,H3NO
H4N,152,0,18,305,0,305,,,,,,,H3N,HN3O4,H4NO,H2O,HN3O3
N3O4,133,1,35,464,4,460,,N5O4,,,,,N3O3,NO2,HN3O4,H2O,HN3O3

==> key_molecules_reactions_hash.csv <==
molecule,total reactions,in reaction,out reaction,total atom transfer,in atom transfer,out atom transfer,from 1,from 2,from 3,from 4,from 5,to 1,to 2,to 3,to 4,to 5,
3974656353,195,11,30,331,228,103,,1415123626,849149194,2844822199,1315470353,1728850689,2874893556,2828709746,329115944,3494029282,299135856
1415123626,152,0,18,305,0,305,,,,,,,3974656353,2874893556,2828709746,329115944,4242784642
4044036521,133,1,35,464,4,460,,1110129394,,,,,1116761848,3630255785,2874893556,329115944,4242784642

==> molecules_smiles.csv <==
2881671817,N2O3,NON(O)=O
2103245073,H2N3O4,[H]ON(O)(O[H])ON#N
2057507959,N3O4,O=NONON=O
4259180543,H2N4O5,[H]ON([H])N(N(O)=O)=N(O)O

==> ring_count.csv <==
ring size 5,ring size 6,ring size 7,ring size 8
0,0,0,0
0,0,0,0
0,0,0,0

==> species_count.csv <==
H,H2N,H2N2O,H2N2O2,H2N2O4,H2N3O
0,0,0,0,0,0
1,3,0,1,0,0
1,4,0,1,0,1

==> species_count_hash.csv <==
1063985685,1086850095,1093753013,1100911419,1110129394,1116761848
0,2,0,3,0,10
1,0,0,4,0,38
0,0,0,1,0,38

==> reactions.dot <==
digraph ReactionFlow {
  rankdir=LR;
  layout=circo;
  node [shape=box, style=filled, fillcolor=azure2, height=0.5, width=1.5];
  edge [color=dimgray];

  node315248181 [label="H2N4O3"];
  node1116761848 [label="N3O3"];
  node329115944 [label="H2O"];

  node2874893556 -> node4242784642 [label="R=5 AT=20", penwidth=3.609438, color=goldenrod];
  node696609698 -> node3630255785 [label="R=13 AT=23", penwidth=4.5649495, color=goldenrod];
  node1415123626 -> node299135856 [label="R=6 AT=8", penwidth=3.7917595, color=goldenrod];
  node3974656353 -> node329115944 [label="R=7 AT=7", penwidth=3.9459102, color=goldenrod];
}

**说明: 如species_count_hash.csv中，1063985685在molecules_smiles.csv中、reactions.dot中必有正确对应，可以查到其化学式和SMILES**