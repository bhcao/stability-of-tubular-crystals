/* -*- c++ -*- ----------------------------------------------------------
   AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)

   分子基类，子类提供能量函数和初始模型，molecule 进行演化，包括
   1. 更新函数（演化）、输出函数
   2. 求导精度步长、朗之万动力学物理参数及其设置函数
------------------------------------------------------------------------- */

#ifndef NANO_MOLECULE_H_
#define NANO_MOLECULE_H_

#include <string>

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include "nmath.h"
#include "narray.h"

namespace nano {

// dump 参数，const C++ 会优化掉的，大胆用
typedef int64_t dump_t;
const dump_t NOTHING  = 0; // 空
// 输出 DATA 文件，包含键信息，使用后一切其他选项失效
const dump_t DATA_FILE  = 1 << 0;
const dump_t EMPHASIS   = 1 << 1; // 强调粒子
// dump 文件时输出其他信息
const dump_t VELOCITY   = 1 << 2; // 输出动能信息
const dump_t DIV_FORCE  = 1 << 3; // 输出梯度力
const dump_t LAN_FORCE  = 1 << 4; // 输出朗之万力
const dump_t K_ENERGY   = 1 << 5; // 输出动能信息
const dump_t P_ENERGY   = 1 << 6; // 输出局部势能
const dump_t GAUSS_CURVE = 1 << 6; // 输出高斯曲率
const dump_t MEAN_CURVE = 1 << 7; // 输出平均曲率

// 检查 i 中是否有宏 name
#define DUMP_CHECK(name, i) ((name & i) != 0)

}

// 采用的动力学模型，0 为梯度下降，1 为朗之万动力学，2 为过阻尼朗之万动力学
#ifndef DYNAMICS
    #define DYNAMICS 2
#endif

// 集合了节点、键、相邻的类，必须子类包装
class molecule {
public:
    // 初始化，求导精度有需要自行用 set_precision 调整，其他参数输入。bond_div_node_n 键数比原子数
    inline molecule(double step, double mass, double damp, double tempr, int node_n,
        int bond_div_node_n, int argc, char* argv[]): time(0), step(step), precision(1e-10), mass(mass), 
        damp(damp), tempr(tempr), nodes(node_n), velocities(node_n), adjacents(node_n),
        bonds(bond_div_node_n * node_n) {
    #ifdef USE_KOKKOS
        Kokkos::initialize(argc, argv);
    #endif
    }
    
#ifdef USE_KOKKOS
    inline ~molecule() { Kokkos::finalize(); }
#endif
    
    // restart 方法，从文件中读取，储存到文件中（加两个能量函数设置）
    molecule(std::string fname, double(*node_energy)(nano::vector, nano::sarray<nano::vector>,
        nano::sarray<double>), double(*bond_energy)(nano::vector, nano::vector, nano::sarray<double>),
        int argc, char* argv[]);
    void store(std::string fname);

    double total_energy(nano::vector range_l, nano::vector range_r);  // 系统整体的能量
    void update();          // 更新函数

    // 输出到以 fname 为名的文件，注意 .dump 会追加而不是覆盖，尽量删除再运行（默认 dump 什么也不输出）
    void dump(std::string fname, nano::dump_t dump_type = 0);

    // 设置参数，如果小于零只输出当前参数，funcname 函数名，paraname 参数名
    #define SET_PARA_FUNC(name) \
        inline double set_ ## name(double para) { \
            return (para < 0) ? this->name : (this->name = para); \
        }
    
    // 包装的目的是防止乱动
    SET_PARA_FUNC(step)
    SET_PARA_FUNC(precision)
    SET_PARA_FUNC(mass)
    SET_PARA_FUNC(damp)
    SET_PARA_FUNC(tempr)
    
    // 如果大于数组长度扩充长度，para < 0 返回 -1
    inline double set_paras(double para, int i) {
        if (i < this->paras.size())
            return (para < 0) ? this->paras[i] : (this->paras[i] = para);
        if (para >= 0) {
            this->paras.set_size(i+1);
            return this->paras[i] = para;
        }
        return -1;
    }

protected:
    // 坐标、速度、相邻粒子（代替键）
    nano::darray<nano::vector> nodes;
    nano::darray<nano::vector> velocities;
    nano::darray<nano::sarray<int>> adjacents;

    // 键，只在 dump 中有用，计算时无用
    nano::darray<nano::pair<int>> bonds;
    nano::sarray<int> emphasis; // 强调的粒子，输出时会用不同 type

    // 函数指针，以粒子为中心的能量，其余顺序已排好
    double (*node_energy)(nano::vector center, nano::sarray<nano::vector> others, nano::sarray<double> paras);
    // 函数指针，以键为中心的能量
    double (*bond_energy)(nano::vector a, nano::vector b, nano::sarray<double> paras);
    
    // 更新速度，主要方便初始化
    void update_velocity(int i);

    // 随机函数池
    nano::rand_pool prand_pool;

private:
    // 运行参数
    double step;      // 更新的步长
    int time;         // 当前时间（运行次数）
    double precision; // 求导精度（增加的 delta 长度）

    // 能量参数，会传递至 node_energy、bond_energy 函数
    nano::sarray<double> paras;

    // 物理参数
    double mass;      // 粒子质量
    double damp;      // 朗之万阻尼系数
    double tempr;     // 系统的温度
    
    // 局部能量，第一个键能平分，算局部和总体能量，第二个键能不平分
    double local_energy(nano::vector center, nano::sarray<nano::vector> others);
    double local_energy_for_update(nano::vector center, nano::sarray<nano::vector> others);

    // 梯度
    inline nano::vector div(nano::vector center, nano::sarray<nano::vector> others) {
        nano::vector temp[3] = {center, center, center};
        temp[0][0] += precision;
        temp[1][1] += precision;
        temp[2][2] += precision;

        #define DIV(i) (local_energy_for_update(temp[i], others) - \
            local_energy_for_update(center, others)) / this->precision
        return { DIV(0), DIV(1), DIV(2) };
    }    
};

#endif // NANO_MOLECULE_H_
