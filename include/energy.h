/* -*- c++ -*- ----------------------------------------------------------
   AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)

   各种能量函数：
   1. 键能
   2. 曲率能，高斯曲率与平均曲率
------------------------------------------------------------------------- */

#ifndef NADO_ENERGY_H_
#define NANO_ENERGY_H_

#include <nmath.h>
#include <narray.h>

namespace energy_func {

// 周围区域面积与所有角
double size_around(nano::vector center, nano::sarray<nano::vector> others);
nano::sarray<double> angles_around(nano::vector center, nano::sarray<nano::vector> others);

// 高斯曲率与平均曲率
double gauss_curvature(nano::vector center, nano::sarray<nano::vector> others, 
        nano::sarray<double> angle, double size);
double mean_curvature(nano::vector center, nano::sarray<nano::vector> others, 
        nano::sarray<double> angle, double size);

// 两类能量函数，paras 顺序为 rest_len, k, tau
double node_energy_func(nano::vector center, nano::sarray<nano::vector> others, nano::sarray<double> paras);
double bond_energy_func(nano::vector p1, nano::vector p2, nano::sarray<double> paras);

}

#endif // NANO_ENERGY_H_
