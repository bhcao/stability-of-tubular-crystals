/*********************************************************************//**
 * AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)
 * 
 * @file include/energy.h
 * @brief 能量函数
 * @version 0.2.0
 * 
 * 各种能量函数：
 * 1. 键能
 * 2. 曲率能，高斯曲率与平均曲率
 **************************************************************************/

#ifndef NADO_ENERGY_H_
#define NANO_ENERGY_H_

#include <nmath.h>
#include <narray.h>

namespace energy_func {

//! 周围区域总面积除 3
double size_around(nano::vector center, nano::sarray<nano::vector> others);

//! 高斯曲率
double gauss_curvature(nano::vector center, nano::sarray<nano::vector> others, double size);
//! 平均曲率
double mean_curvature(nano::vector center, nano::sarray<nano::vector> others, double size);

//! 曲率能量函数，paras 顺序为 rest_len, k, tau
double node_energy_func(nano::vector center, nano::sarray<nano::vector> others, nano::sarray<double> paras);
//! 键能函数，paras 顺序为 rest_len, k, tau
double bond_energy_func(nano::vector p1, nano::vector p2, nano::sarray<double> paras);

}

#endif // NANO_ENERGY_H_
