#include <cmath>
#include <iostream>

#include "model.h"
#include "nmath.h"

double my_bond_energy(nano::vector p1, nano::vector p2, double k, double rest_len) {
    return k/2 * (nano::mod(p1-p2) - rest_len) * (nano::mod(p1-p2) - rest_len);
}

double my_node_energy(nano::vector center, nano::sarray<nano::vector> others, double tau) {
    // 一个点周围所有角，避免重复计算
    nano::sarray<double> angle;
    for (int i=0; i<others.size()-1; i++) {
        // 计算任意三点之间的所成角，n1, n3 为两边，n2 为顶点
        double l_1 = nano::mod(center - others[i]);
        double l_2 = nano::mod(center - others[i+1]);
        double product = (others[i+1]-center)*(others[i]-center);
        angle.push_back(std::acos(product/l_2/l_1));
    }
    
    double l_1 = nano::mod(center - others[0]);
    double l_2 = nano::mod(center - others[others.size()-1]);
    double product = (others[others.size()-1]-center)*(others[0]-center);
    angle.push_back(std::acos(product/l_2/l_1));

    // 计算周围区域面积
    double size = 0;
    for (int i=0; i<others.size()-1; i++) {
        // 叉乘的模即三角形面积
        size += 1.0/2 * nano::mod(nano::cross(others[i]-center, others[i+1]-center));
    }
    size += 1.0/2 * nano::mod(nano::cross(others[others.size()-1]-center, others[0]-center));
    
    // 计算高斯曲率
    double K = 2*PI;
    for (int i=0; i<angle.size(); i++) {
        K -= angle[i];
    }
    K = K * 3 / size;
    
    // 计算平均曲率中的向量
    nano::vector H_vec;
    for (int i=1; i<others.size(); i++) {
        H_vec = H_vec + (1/std::tan(angle[i-1])/std::tan(angle[i-1]) + 
            1/std::tan(angle[i])/tan(angle[i])) * (center-others[i]);
    }
    H_vec = H_vec + (1/std::tan(angle[angle.size()-1])/std::tan(angle.size()-1) +
        1/std::tan(angle[0])/std::tan(angle[0])) * (center-others[0]);
    
    double H = nano::mod(H_vec)/4/size;
    return tau * (2*H*H - K);
}