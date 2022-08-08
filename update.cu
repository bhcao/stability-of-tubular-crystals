#include <cmath>

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
        double l_12 = dist(center, others[i]);
        double l_23 = dist(center, others[i+1]);
        double l_13 = dist(others[i], others[i+1]);
        angle.push_back(std::acos(l_12/2/l_23 + l_23/2/l_12 - l_13*l_13/2/l_23/l_12));
    }
    
    double l_12 = dist(center, others[0]);
    double l_23 = dist(center, others[others.size()-1]);
    double l_13 = dist(others[0], others[others.size()-1]);
    angle.push_back(std::acos(l_12/2/l_23 + l_23/2/l_12 - l_13*l_13/2/l_23/l_12));

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

#define OTHER_PARA others[a], ppara.k, ppara.rest_len, ppara.tau
__global__ void cudaUpdate(node *center, nano::s_vector<node> *others, nano::s_vector<int> *cuda_others_id, 
        para ppara, double step, double precision) {
    // CUDA 线程编号
    #ifdef USE_CUDA
    int a = blockIdx.x;
    #else
    int a = 0;
    #endif

    node temp = center[a];
    temp.i += precision;
    // 偏导数乘以步长
    double i = step * (energy(temp, OTHER_PARA) - energy(center[a], OTHER_PARA)) / precision;
    temp = center[a];
    temp.j += precision;
    double j = step * (energy(temp, OTHER_PARA) - energy(center[a], OTHER_PARA)) / precision;
    temp = center[a];
    temp.k += precision;
    double k = step * (energy(temp, OTHER_PARA) - energy(center[a], OTHER_PARA)) / precision;
    node div = {i, j, k};
    center[a] = center[a] - div;

    // 更新邻域，准备下一次计算，CPU 上有不同步风险，只能尽可能忽略其影响
    for (int j=0; j<others[a].size(); j++) {
        others[a][j] = center[cuda_others_id[a][j]];
    }
}

void model::update() {
    #ifdef USE_CUDA
    node *cuda_nodes = this->nodes.get_gpu_data();
    nano::s_vector<node> *cuda_others = this->adjacents.get_gpu_data();
    nano::s_vector<int> *cuda_others_id = this->adjacents_id.get_gpu_data();

    cudaUpdate<<<1,this->nodes.size()>>>(cuda_nodes, cuda_others, cuda_others_id, this->ppara, 
        this->step, this->precision);
	
    this->nodes.cpu_synchro();
    this->adjacents.cpu_synchro();
    #else
    for (int a=0; a<this->nodes.size(); a++) {
        node others;
        cudaUpdate(&this->nodes[a], &this->adjacents[this->nodes[a].id], &this->adjacents_id[a],
            this->ppara, this->step, this->precision);
    }
    #endif
}
