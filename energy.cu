#include <cmath>
#include <vector>

#include "molecule.h"
#include "model.h"

using std::vector, std::acos, std::tan;

// 虚函数的实现，生成一个点周围的能量，便于更新
double model::energy(node center) {
    vector<node> others = this->adjacents[center.id].others;
    double count = 0;
    for (int i=0; i<others.size(); i++) {
        count += bond_energy(center, others[i]);
    }
    count += curve_energy(center, others);
    return count;
}

double model::bond_energy(node p1, node p2) {
    return this->ppara.k/2 * (dist(p1, p2) - this->ppara.rest_len)*
        (dist(p1, p2) - this->ppara.rest_len);
}

// 一个点周围所有角，避免重复计算
vector<double> angle_around(node center, vector<node> others) {
    vector<double> angle;
    for (int i=0; i<others.size()-1; i++) {
        // 计算任意三点之间的所成角，n1, n3 为两边，n2 为顶点
        double l_12 = dist(center, others[i]);
        double l_23 = dist(center, others[i+1]);
        double l_13 = dist(others[i], others[i+1]);
        angle.push_back(acos(l_12/2/l_23 + l_23/2/l_12 - l_13*l_13/2/l_23/l_12));
    }
    
    double l_12 = dist(center, others[0]);
    double l_23 = dist(center, others[others.size()-1]);
    double l_13 = dist(others[0], others[others.size()-1]);
    angle.push_back(acos(l_12/2/l_23 + l_23/2/l_12 - l_13*l_13/2/l_23/l_12));
    return angle;
}

double model::curve_energy(node center, vector<node> others) {
    
    vector<double> angle = angle_around(center, others);
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
        H_vec = H_vec + (1/tan(angle[i-1])/tan(angle[i-1]) + 1/tan(angle[i])/tan(angle[i])) * (center-others[i]);
    }
    H_vec = H_vec + (1/tan(angle[angle.size()-1])/tan(angle.size()-1) + 1/tan(angle[0])/tan(angle[0])) * (center-others[0]);
    
    double H = dist(H_vec)/4/size;
    return this->ppara.tau * (2*H*H - K);
}

// 总能量
double model::total_energy() {
    double count = 0;
    for (int i=0; i<this->nodes.size(); i++) {
        count += curve_energy(this->nodes[i], this->adjacents[i].others);
    }
    for (int i=0; i<this->bonds.size(); i++) {
        count += bond_energy(this->nodes[this->bonds[i].a],
            this->nodes[this->bonds[i].b]);
    }
    return count;
}