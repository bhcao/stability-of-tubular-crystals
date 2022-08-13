# nanotube

### 介绍
利用梯度下降法模拟三角形纳米管在位错下的稳定形状，复刻论文《可弯曲管状晶体在移动的位错的相互作用下的形状多稳定性》（Shape multistability in ﬂexible tubular crystals through interactions of mobile dislocations） [^论文引用]

![模型直观图](https://gitee.com/Bovera/nanotube/raw/master/model.png)  
图 1. 模型直观图

### 模型表示

如图 2，边长为 $m$，$n$ 的平行四边形沿 $A$、$B$ 点重合后，将数轴下方平移拼接至上方，即构成了叶序数为 $(m,n)$ 的圆柱。其半径、角度等参数很容易得出，如（其他详见论文）

$$R=\frac{|AB|}{2\pi}=\frac{a\sqrt{m^2+n^2-mn}}{2\pi}.$$

![模型示意](https://gitee.com/Bovera/nanotube/raw/master/tube.png)  
图 2. 模型示意

当然，为了获得更长的管子，我们通过上图所示基本单元沿 $m$ 方向的平移延长管子。

### 位错表示

下面是错位的表示，其中绿线代表生成的新键，蓝点代表增加或减少的原子，红点为位错原点（初始化时位错发生的起点）。

![位错生成示意](https://gitee.com/Bovera/nanotube/raw/master/dislocation.png)  
图 3. 位错生成示意

滑移相对与整体有三个方向，如图 4 所示，右下方键的移动方法与已有重复，故不再考虑。生成滑移时，滑移先从原点产生，接着位错远点根据滑移步数正负号移动。蓝色箭头即移动方向。

![滑移方向约定](https://gitee.com/Bovera/nanotube/raw/master/direction.png)  
图 4. 滑移方向约定

### 主函数示例
```cpp
para ppara = default_para; // 设置参数为模型默认参数[1]
ppara.repeat = 10;         // 修改模型参数
model this_model(ppara);   // 以 ppara 参数生成模型
// 生成模型后再修改模型参数没用

// 求位错间的距离，x 为水平距离（沿管面），z 为轴向距离
double x, z;
this_model.dis_pair_distance(&x, &z);
std::cout << x << z;

this_model.precision = 1e-7;     // 更改运行时参数[2]
this_model.disorganize();        // 退火
std::ofstream fout("test.dump"); // dump 文件
figure energy_change;            // 声明图片类

for (int k=0; k<100; k++) {
    this_model.update();         // 模型计算更新
    // 输出当前模型至 dump 文件，第二个参数指定文件类型[3]
    this_model.dump(fout, DUMP_FILE);       
    // 计算总能量并增加至图片
    energy_change << this_model.total_energy();
}

energy_change.draw("test.png"); // 画图
fout.close();                   // 保存 dump 文件
```

### 模型参数
模型构建参数：
> **说明：** 对应主函数示例注释 [1]，等号后面是参数默认值（下同）

| 参数                       | 解释          |
| ---                        | ---         |
|`int m = 13, n = 13`        | 叶序螺旋数    |
|`double rest_len = 0.1`     | 键的平衡长度  |
|`int direction = 1`         | 位错发生的方位 |
|`int glide = 3, climb = 3`  | 位错滑移和攀移的步数 |
|`int repeat = 3`            | 模型重复的次数（增长模型）|
|`double k = 1, tau = 3e-7`  | 能量参数（第一项键能、第二项曲率）|

运行时参数：
> **说明：** 对应主函数示例注释 [2]

| 参数                      | 解释          |
| ---                      | ---         |
|`double step = 2e-5`      |  更新的步长，向梯度方向前进多少|
|`double precision = 1e-6` | 求导精度（求导时增加的 $\Delta$ 值）
|`double range = 1.0`      | 退火时粒子随机移动的高斯分布的方差 |

文件类型：
> **说明：** 对应主函数示例注释 [3]，与 LAMMPS 文件类型一致。对文件的详细描述见 [LAMMPS 官网](https://www.lammps.org/)，可视化读取文件可使用 [Ovito](https://www.ovito.org/)。我们对位错原子采取了特别标记，在 Ovito 中会展示不同颜色。

| 类型                      | 解释          |
| ---                      | ---         |
|`DATA_FILE`               | .data 文件，包含键信息，同一个文件只能输出一次 |
|`DUMP_FILE`               | .dump 文件，不包含键信息，同一个文件可以输出多次 |

### 并行测试
1. 粒子数小于 1000 时，单线程或极少线程不建议开启 Cuda，反而拖慢速度；若开启 Cuda，线程数可以设置为电脑 CPU 核数的 2-3 倍；
2. 粒子数大于 1000 小于 10000 时，Cuda 加速明显，可以开启；线程数可以超过 CPU 核数一些；
3. 粒子数大于 10000 时，速度限制因素是 CPU，线程数应小于或等于 CPU 核数。

> **提醒：** CUDA 函数执行错误会以数字代码的形式返回，在不进行错误处理时会自动忽略，不会报错（错误处理暂时还未完成）；因而请仔细检查显卡驱动是否为 Nvidia 驱动等，避免无法调用却不报错。

### 编译选项
> **提醒：** Windows 手动编译 CPU 版本时需将 update.cu、adjacent.cu 重命名为 update.cpp、adjacent.cpp

1. `make gpu` 使用 Cuda 加速，需要安装 nvcc。手动编译命令如下
```
nvcc -D USE_CUDA main.cpp energy.cpp model.cpp molecule.cpp update.cu adjacent.cu -o main.out --expt-relaxed-constexpr
```

2. `make cpu` 不使用任何加速，默认编译器 g++，如有需要自行更换。手动编译命令如下
```
g++ main.cpp energy.cpp model.cpp molecule.cpp update.cpp adjacent.cpp -o main.out
```

3. `make gpu-mpi` 使用 Cuda 加速，同时使用 mpich 并行，需要安装 nvcc 和 mpich。手动编译命令如下（两步）
```
nvcc -D USE_CUDA -lib energy.cpp model.cpp molecule.cpp update.cu adjacent.cu --expt-relaxed-constexpr -o nano.a
mpic++ -D USE_MPI main.cpp nano.a -L/usr/lib/cuda/lib64 -lcudart -o main.out
```

> **提醒：** 第二步链接时 -L 后路径替换为为 Cuda 运行时 cudart.dll / libcudart.so 的路径。

4. `make cpu-mpi` 使用 mpich 并行，需要安装 mpich。手动编译命令如下
```
mpic++ -D USE_MPI main.cpp energy.cpp model.cpp molecule.cpp update.cpp adjacent.cpp -o main.out
```

5. 其他选项有 `make clean` 清理所有编译出来的文件；以及 `make debug`（默认值）用默认 C++ 编译器编译后在 gdb 中打开调试。

### 目标
1. 优化相关算法，充分利用显卡，提高运行速度，减少不必要的内存拷贝；
2. 重构代码，使得格式更加清新，模块化更加突出；
3. 测试参数，修改 bug，完成对论文的复刻；
4. 尝试利用计算机代数系统（Python 的 sympy 等）进行梯度计算，以提高速度同时保证准确性；
5. 完成自己的项目，进行拓展，考虑四方格子的情况等。

[^论文引用]: Zakharov A, Beller D A. [Shape multistability in flexible tubular crystals through interactions of mobile dislocations](https://www.pnas.org/doi/pdf/10.1073/pnas.2115423119)[J]. Proceedings of the National Academy of Sciences, 2022, 119(6): e2115423119.