/*********************************************************************//**
 * AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)
 * 
 * @file src/energy.cpp
 * @brief 能量函数的实现
 * @version 0.2.0
 **************************************************************************/

#include <cmath>

#include "energy.h"
#include "nmath.h"
#include "narray.h"

namespace energy_func {

double bond_energy_func(nano::vector p1, nano::vector p2, nano::sarray<double> paras) {
    double rest_len = paras[0], k = paras[1];
    return k/2 * (nano::mod(p1-p2) - rest_len) * (nano::mod(p1-p2) - rest_len);
}

// 计算周围区域面积
double size_around(nano::vector center, nano::sarray<nano::vector> others) {
    double size = 0;
    for (int i=0; i<others.size()-1; i++)
        // 叉乘的模即三角形面积
        size += nano::mod(nano::cross(others[i]-center, others[i+1]-center))/2;
    size += nano::mod(nano::cross(others[others.size()-1]-center, others[0]-center))/2;
    return size/3;
}

// 一个点周围所有角
nano::sarray<double> angles_around(nano::vector center, nano::sarray<nano::vector> others) {
    nano::sarray<double> angle;
    for (int i=0; i<others.size()-1; i++)
        // 计算任意三点之间的所成角
        angle.push_back(std::acos(((others[i+1]-center)*(others[i]-center))/
            nano::mod(center - others[i])/nano::mod(center - others[i+1])));
    angle.push_back(std::acos(((others[others.size()-1]-center)*(others[0]-center))/
        nano::mod(center - others[0])/nano::mod(center - others[others.size()-1])));

    return angle;
}

double gauss_curvature(nano::vector center, nano::sarray<nano::vector> others, 
        nano::sarray<double> angle, double size) {    
    // 计算高斯曲率
    double K = 2*PI;
    for (int i=0; i<angle.size(); i++)
        K -= angle[i];
    return K/size;
}

double mean_curvature(nano::vector center, nano::sarray<nano::vector> others, 
        nano::sarray<double> angle, double size) {

    nano::sarray<double> angle1;
    for (int i=0; i<others.size()-1; i++)
        // 计算任意三点之间的所成角
        angle1.push_back(std::acos(((center-others[i+1])*(others[i]-others[i+1]))/
            nano::mod(center - others[i+1])/nano::mod(others[i] - others[i+1])));
    angle1.push_back(std::acos(((center-others[0])*(others[others.size()-1]-others[0]))/
        nano::mod(center - others[0])/nano::mod(others[others.size()-1] - others[0])));
    
    nano::sarray<double> angle2;
    angle2.push_back(std::acos(((center-others[others.size()-1])*(others[0]-others[others.size()-1]))/
        nano::mod(center - others[others.size()-1])/nano::mod(others[0] - others[others.size()-1])));
    for (int i=0; i<others.size()-1; i++)
        // 计算任意三点之间的所成角
        angle2.push_back(std::acos(((center-others[i])*(others[i+1]-others[i]))/
            nano::mod(center - others[i])/nano::mod(others[i+1] - others[i])));

    // 计算平均曲率中的向量
    nano::vector H_vec;
    for (int i=0; i<others.size(); i++)
        H_vec += (1/std::tan(angle1[i]) + 1/std::tan(angle2[i])) * (center-others[i]);

    return nano::mod(H_vec)/4/size;
}

double node_energy_func(nano::vector center, nano::sarray<nano::vector> others, nano::sarray<double> paras) {
    double tau = paras[2];
    nano::sarray<double> angle = angles_around(center, others);
    double size = size_around(center, others);
    double K = gauss_curvature(center, others, angle, size);    
    double H = mean_curvature(center, others, angle, size);
    return tau * (2*H*H - K);
}

}