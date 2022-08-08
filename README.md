# nanotube

### 介绍
纳米管模拟 C++ 程序

### 编译选项
`make gpu` 使用 Cuda 加速，需要安装 nvcc

`make cpu` 不使用任何加速，默认编译器 g++，如有需要自行更换

`make gpu-mpi` 使用 Cuda 加速，同时使用 mpich 并行，需要安装 nvcc 和 mpich

`make cpu-mpi` 使用 mpich 并行，需要安装 mpich