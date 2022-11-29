# nanotube

## 模型

### 模型来源

利用梯度下降法模拟三角形纳米管在位错下的稳定形状，复刻论文《可弯曲管状晶体在移动的位错的相互作用下的形状多稳定性》（[Shape multistability in ﬂexible tubular crystals through interactions of mobile dislocations](https://www.pnas.org/doi/pdf/10.1073/pnas.2115423119)）

> Zakharov A, Beller D A. Shape multistability in flexible tubular crystals through interactions of mobile dislocations[J]. Proceedings of the National Academy of Sciences, 2022, 119(6): e2115423119.

![模型直观图](https://gitee.com/Bovera/nanotube/raw/master/fig/model.png)  
<p align="center"> 图 1. 模型直观图 </p>

### 模型表示

如图 2，边长为 $m$，$n$ 的平行四边形沿 $A$、$B$ 点重合后，将数轴下方平移拼接至上方，即构成了叶序数为 $(m,n)$ 的圆柱。其半径、角度等参数很容易得出，如（其他详见论文）

$$R=\frac{|AB|}{2\pi}=\frac{a\sqrt{m^2+n^2-mn}}{2\pi}.$$

![模型示意](https://gitee.com/Bovera/nanotube/raw/master/fig/tube.png)  
<p align="center"> 图 2. 模型示意 </p>

当然，为了获得更长的管子，我们通过上图所示基本单元沿 $m$ 方向的平移延长管子。

### 位错表示

下面是错位的表示，其中绿线代表生成的新键，蓝点代表增加或减少的原子，红点为位错原点（初始化时位错发生的起点）。

![位错生成示意](https://gitee.com/Bovera/nanotube/raw/master/fig/dislocation.png)  
<p align="center"> 图 3. 位错生成示意 </p>

滑移相对与整体有六个方向，如图 4 所示，direction 每加 1，逆时针旋转 120 度；变为负号则沿着滑移方向镜像。

![滑移方向约定](https://gitee.com/Bovera/nanotube/raw/master/fig/direction.png)  
<p align="center"> 图 4. 滑移方向约定 </p>

### 动力学

通过在引入 `molucule.h` 前定义 `DYNAMICS` 宏选择动力学方法，`DYNAMICS` 值分别为
| 值 | 模型           | 说明 |
| -- | ---           | --- |
|0  | 梯度下降         | 等价于质量温度均为 0，阻尼为 1 的朗之万动力学 |
|1  | 朗之万动力学      |  |
|2  | 过阻尼朗之万动力学 | 等价于质量为 0 的朗之万动力学 |
注意对于梯度下降和过阻尼朗之万动力学，部分参数无效，详见后文。

## 函数及其参数

### 模型参数设置方法

模型参数可以在两处设置：

1. para 结构体，初始化时设置。主要调控模型的拓扑结构，设置能量参数、步长和朗之万动力学参数；
2. set 函数，运行时动态修改。设置能量参数、步长精度和朗之万动力学参数。

> 注意：
> 1. 初始化之前必须设置好拓扑参数，初始化后不可更改；
> 2. 初始化后再修改 para 结构体无效，请用 set 函数；
> 3. 精度一般不用设置，要设置只能用 set 函数。

使用方法示例如下：

```cpp
para ppara = default_para;      // 设置 para 为默认值
ppara.repeat = 10;              // 修改 para 成员参数值
model this_model(ppara);        // 以 ppara 参数初始化
this_model.set_precision(1e-9); // 动态更新参数
```

对于能量参数，设置时没有名称，需使用 `set_paras` 函数，使用如下
```cpp
this_model.set_paras(0.1, 2); // 设置 tau
```
第一个参数是值，第二个参数指定参数，0 为 rest_len，1 为 k，2 为 tau。

### 模型参数说明

1. 标 √ 的表示理论上会影响稳定时能量和形状的参数，其它参数理论上不应该影响结果。第三组参数呈比例增加理论上也不该影响结果。
2. 标 ⋯ 的表示可以使用 set 函数修改，其余只能初始化时修改

**第一组：完美晶体的拓扑描述**

|| 参数                       | 解释              |
|--| ---                        | ---              |
|√|`int m = 13, n = 13`        | 叶序螺旋数         |
||`int repeat = 3`            | 模型重复的次数（增长模型）|

**第二组：晶体的缺陷表示**

|| 参数                       | 解释              |
|--| ---                        | ---              |
|√|`int glide = 3`             | 位错滑移的步数 |
|√|`int climb = 3`             | 位错攀移的步数 |
|√|`int direction = 1`         | 位错发生的方位    |
||`int bn = 0`                | 位错初始位置（以中点为 0）|

**第三组：原长及能量系数**

|| 参数                       | 解释              |
|--| ---                        | ---             |
|⋯ √|`double rest_len = 0.1`     | 键的平衡长度      |
|⋯ √|`double k = 1`             | 键能系数 |
|⋯ √|`double tau = 3e-7`        | 曲率能量系数 |

**第四组：步长及精度**

|| 参数                      | 解释          |
|--| ---                      | ---         |
|⋯|`double step = 1e-6`      | 步长，每一步的时间间隔 |
|⋯|`double precision = 1e-9` | 求导精度（差分增加的值）|

**第五组：朗之万动力学参数**

|| 参数                      | 解释          |
|--| ---                      | ---         |
|⋯|`double mass = 1`        | 粒子质量（对过阻尼朗之万与梯度下降无效） |
|⋯|`double damp = 1`        | 朗之万阻尼系数（对梯度下降无效） |
|⋯|`double tempr = 273`     | 系统温度（对梯度下降无效） |

### 输出函数参数

`dump` 函数输出所需值，其原型为
```cpp
molecule::dump(std::string fname, nano::dump_t dump_type);
```

第一个参数表示输出的文件名（注意不带后缀名 .data 和 .dump），第二个为调控输出内容的参数，使用方法是将所需 `nano::dump_t` 常量逻辑与，例如
```cpp
this_model.dump("test", nano::EMPHASIS | nano::LAN_FORCE);
```

`dump_t` 常量及其含义如下（省略了命名空间 `nano`）

| 常量名        | 作用         |
| ---          | ---         |
|`NOTHING`     | 啥也不干  |
|`DATA_FILE`   | 输出 .data 文件  |
|`EMPHASIS`    | 强调位错粒子（类型设为 2）    |
|`VELOCITY`    | 输出粒子速度信息      |
|`DIV_FORCE`   | 输出梯度力（能量梯度）     |
|`LAN_FORCE`   | 输出朗之万力（包括了阻尼和随机项）|
|`K_ENERGY`    | 输出动能信息  |
|`P_ENERGY`    | 输出局部势能（键能取一半） |
|`GAUSS_CURVE` | 输出高斯曲率 |
|`MEAN_CURVE`  | 输出平均曲率 |

> 注意：
> 1. 一定会输出 .dump 文件（不含键），使用了 `DATA_FILE` 才会输出 .data 文件（包含键信息）；
> 2. .data 多次输出会在不同文件，文件名自动加 `_n.data`，n 为当前运行步数；.dump 会在同一个文件，文件名加 `_0.dump`；
> 3. 命令行可以使用 `ovito name*` （不包含自动加的后缀）同时导入 `.data` 和 `.dump`，也可以导入 `.data` 文件后使用 load trajactory；
> 4. `LAN_FORCE` 和 `DIV_FORCE` 选项对于梯度下降和过阻尼朗之万无意义。

### 其他函数

其他可调用的函数原型有：
```cpp
// 位错对间距，返回值第一个为垂轴，第二个为沿轴
nano::pair<double> model::dis_pair_distance();
// 总能量，第一个参数为下界，第二个参数为上界
double molecule::total_energy(nano::vector bound1, nano::vector bound2);
void molecule::update(); // 运行一次
```

### 储存与恢复

将当前状态储存为二进制文件：
```cpp
this_model.store("restore.bin"); // 参数为文件名
```
恢复时不需要 model 类，直接调用 molecule 初始化函数即可
```cpp
molecule this_model("restart.bin", my_node_energy, my_bond_energy);
```

### 文件概览

1. `nmath.h/nmath.cpp` 定义了数学相关类，如三维矢量、数组、随机函数以及一些常数；
2. `molecule.h/molecule.cpp` 定义了 molecule 类，这是一个可以根据自定义能量函数、拓扑结构，进而通过朗之万动力学进行演化的类；
3. `energy.cpp` 定义了两个能量函数；
4. `model.h/model.cpp` 定义了 model 类，这是 molecule 类的子类，初始化时进行建模，继承了 molecule 的函数；
5. `main.cpp` 主函数，调用其他类进行演化，根据需要选择运行步骤次数以及动态更改参数。