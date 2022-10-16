#include <cmath>
#include <iostream>

#include "vector.h"
#include "molecule.h"
#include "model.h"

__device__ double bond_energy(node p1, node p2, double k, double rest_len) {
    return k/2 * (dist(p1, p2)-rest_len) * (dist(p1, p2)-rest_len);
}

__device__ double curve_energy(node center, nano::s_vector<node> others, double tau) {
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

// 虚函数的实现，生成一个点周围的能量，便于更新
__device__ double energy(node center, nano::s_vector<node> others, double k, double rest_len, double tau) {
    double count = 0;
    for (int i=0; i<others.size(); i++) {
        count += bond_energy(center, others[i], k, rest_len);
    }
    count += curve_energy(center, others, tau);
    return count;
}

// 所有力总和
#define OTHER_PARA others, ppara.k, ppara.rest_len, ppara.tau
__device__ node accelerate(node center, node speed, nano::s_vector<node> others, para ppara, double precision) {
    node temp = center;
    temp.i += precision;
    double i = (energy(temp, OTHER_PARA) - energy(center, OTHER_PARA)) / precision;

    temp = center;
    temp.j += precision;
    double j = (energy(temp, OTHER_PARA) - energy(center, OTHER_PARA)) / precision;

    temp = center;
    temp.k += precision;
    double k = (energy(temp, OTHER_PARA) - energy(center, OTHER_PARA)) / precision;
    
    node f_c = {i, j, k};
    return (-1)*f_c/ppara.mass - ppara.damp*speed + std::sqrt(2*ppara.damp*ppara.tempr*K_B/ppara.mass)*randnode();
}

__global__ void cudaUpdate(node *center, node *speed, nano::s_vector<node> *others, nano::s_vector<int> *cuda_others_id, 
        para ppara, double step, double precision) {
    // CUDA 线程编号
    #ifdef USE_CUDA
    int a = blockIdx.x;
    #else
    int a = 0;
    #endif

    speed[a] = speed[a] + accelerate(center[a], speed[a], others[a], ppara, precision)*step;
    center[a] = center[a] + speed[a]*step;
    
    #ifdef USE_CUDA
    // 更新邻域，准备下一次计算，CPU 上分立以保证正确性
    for (int j=0; j<others[a].size(); j++) {
        others[a][j] = center[cuda_others_id[a][j]];
    }
    #endif
}

#ifdef USE_CUDA
void model::update() {
    node *cuda_nodes = this->nodes.get_gpu_data();
    node *cuda_speeds = this->speeds.get_gpu_data();
    nano::s_vector<node> *cuda_others = this->adjacents.get_gpu_data();
    nano::s_vector<int> *cuda_others_id = this->adjacents_id.get_gpu_data();
    
    cudaUpdate<<<this->nodes.size(), 1>>>(cuda_nodes, cuda_speeds, cuda_others, cuda_others_id, this->ppara, 
        this->step, this->precision);
    
    this->speeds.cpu_synchro();
    this->nodes.cpu_synchro();
    this->adjacents.cpu_synchro();
}
#endif

#ifndef USE_CUDA
void model::update() {
    for (int a=0; a<this->nodes.size(); a++) {
        cudaUpdate(&this->nodes[a], &this->speeds[a], &this->adjacents[a], &this->adjacents_id[a],
            this->ppara, this->step, this->precision);
    }
    for (int a=0; a<this->nodes.size(); a++) {
        for (int j=0; j<this->adjacents[a].size(); j++) {
            this->adjacents[a][j] = this->nodes[this->adjacents_id[a][j]];
        }
    }
}
#endif
