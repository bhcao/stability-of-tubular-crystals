# nanotube

### 介绍
利用梯度下降法模拟三角形纳米管在位错下的稳定形状，复刻论文

[Shape multistability in ﬂexible tubular crystals through interactions of mobile dislocations](https://www.pnas.org/doi/pdf/10.1073/pnas.2115423119)

### 参数
模型构建参数（等号后面是参数默认值）：
1. int m = 13, n = 13;         叶序螺旋数
2. double rest_len = 0.1;      键的平衡长度
3. int direction = 1;          位错发生的方位
4. int glide = 3, climb = 3;   位错滑移和攀移的步数
5. int repeat = 3;             模型重复的次数（增长模型）
6. double k = 1, tau = 3e-7;   能量参数（第一项键能、第二项曲率）

运行时参数：
1. double step = 2e-5;         更新的步长，就是向梯度方向前进多少
2. double precision = 1e-6;    求导精度（求导时增加的 delta 长度）
3. double range = 1.0;         打乱的范围，高斯分布的方差（退火时粒子随机移动的范围）

### 并行测试说明
1. 粒子数小于 1000 时，单线程或极少线程不建议开启 Cuda，反而拖慢速度；若开启 Cuda，线程数可以设置为电脑 CPU 核数的 2-3 倍；
2. 粒子数大于 1000 小于 10000 时，Cuda 加速明显，可以开启；线程数可以超过 CPU 核数一些；
3. 粒子数大于 10000 时，速度限制因素是 CPU，线程数应小于或等于 CPU 核数。

### 编译选项
Windows 手动编译 CPU 版本时需将 update.cu 重命名为 update.cpp

1. `make gpu` 使用 Cuda 加速，需要安装 nvcc。手动编译命令如下
> `nvcc -D USE_CUDA main.cpp energy.cpp model.cpp molecule.cpp update.cu -o main.out --expt-relaxed-constexpr`

2. `make cpu` 不使用任何加速，默认编译器 g++，如有需要自行更换。手动编译命令如下
> `g++ main.cpp energy.cpp model.cpp molecule.cpp update.cpp -o main.out`

3. `make gpu-mpi` 使用 Cuda 加速，同时使用 mpich 并行，需要安装 nvcc 和 mpich。手动编译命令如下（两步）
> `nvcc -D USE_CUDA -lib energy.cpp model.cpp molecule.cpp update.cu --expt-relaxed-constexpr -o nano.a`  
> `mpic++ -D USE_MPI main.cpp nano.a -L/usr/lib/cuda/lib64 -lcudart -o main.out`

注意第二步链接时 -L 后路径替换为为 Cuda 运行时 cudart.dll 的路径

4. `make cpu-mpi` 使用 mpich 并行，需要安装 mpich。手动编译命令如下
> `mpic++ -D USE_MPI main.cpp energy.cpp model.cpp molecule.cpp update.cpp -o main.out`
