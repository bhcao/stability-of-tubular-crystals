#include <cmath>

#include "vector.h"
#include "molecule.h"
#include "model.h"

double cpu_bond_energy(node p1, node p2, double k, double rest_len) {
    return k/2 * (dist(p1, p2)-rest_len) * (dist(p1, p2)-rest_len);
}

double cpu_curve_energy(node center, nano::s_vector<node> others, double tau) {
    // 一个点周围所有角，避免重复计算
    nano::s_vector<double> angle = {0};
    for (int i=0; i<others.size()-1; i++) {
        // 计算任意三点之间的所成角，n1, n3 为两边，n2 为顶点
        double l_1 = dist(center, others[i]);
        double l_2 = dist(center, others[i+1]);
        double product = (others[i+1]-center)*(others[i]-center);
        angle.push_back(std::acos(product/l_2/l_1));
    }
    
    double l_1 = dist(center, others[0]);
    double l_2 = dist(center, others[others.size()-1]);
    double product = (others[others.size()-1]-center)*(others[0]-center);
    angle.push_back(std::acos(product/l_2/l_1));

    // 计算周围区域面积
    double size = 0;
    for (int i=0; i<others.size()-1; i++) {
        // 叉乘的模即三角形面积
        size += 1.0/2 * dist(cross(others[i]-center, others[i+1]-center));
    }
    size += 1.0/2 * dist(cross(others[others.size()-1]-center, others[0]-center));
    
    // 计算高斯曲率
    double K = 2*PI;
    for (int i=0; i<angle.size(); i++) {
        K -= angle[i];
    }
    K = K * 3 / size;
    
    // 计算平均曲率中的向量
    node H_vec = {0};
    for (int i=1; i<others.size(); i++) {
        H_vec = H_vec + (1/std::tan(angle[i-1])/std::tan(angle[i-1]) + 
            1/std::tan(angle[i])/tan(angle[i])) * (center-others[i]);
    }
    H_vec = H_vec + (1/std::tan(angle[angle.size()-1])/std::tan(angle.size()-1) +
        1/std::tan(angle[0])/std::tan(angle[0])) * (center-others[0]);
    
    double H = dist(H_vec)/4/size;
    return tau * (2*H*H - K);
}

// 总能量
double model::total_energy() {
    double count = 0;
    for (int i=0; i<this->nodes.size(); i++) {
        count += cpu_curve_energy(this->nodes[i], this->adjacents[i], ppara.tau);
    }
    for (int i=0; i<this->bonds.size(); i++) {
        count += cpu_bond_energy(this->nodes[this->bonds[i].a],
            this->nodes[this->bonds[i].b], this->ppara.k, this->ppara.rest_len);
    }
    return count;
}
