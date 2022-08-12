#include <cmath>

#include "vector.h"
#include "molecule.h"
#include "model.h"

__device__ void rearrange(nano::s_vector<int> *padjacent_id, node *nodes, node *pnode, node *plane_vec) {
    if ((*padjacent_id).size() == 0 || (*padjacent_id).size() == 1 ) return;
    if ((*padjacent_id).size() == 2) {
        *plane_vec = cross(nodes[(*padjacent_id)[0]]-*pnode, nodes[(*padjacent_id)[1]]-*pnode);
        *plane_vec = *plane_vec / dist(*plane_vec);
        return;
    }
    double angle = 0;
    node last = nodes[(*padjacent_id)[(*padjacent_id).size()-1]] - (*pnode);
    last = last / dist(last);
    for (int i=0; i<(*padjacent_id).size()-1; i++) {
        node temp = nodes[(*padjacent_id)[0]] - *pnode;
        temp = temp / dist(temp);
    }
}

__global__ void cudaGenerate_adjacent(node *pnode, bond *bonds, int bonds_len,
        nano::s_vector<node> *padjacent, nano::s_vector<int> *padjacent_id, node *nodes) {
#ifdef USE_CUDA
    int a = blockIdx.x;
#else
    int a = 0;
#endif
    
    node plane_vec;
    padjacent_id[a] = {0};
    padjacent[a] = {0};
    for (int i = 0; i < bonds_len; i++) {
        // 如果键的一端是 center
        if (bonds[i].a == a) {
            // 记录位置
            padjacent_id[a].push_back(bonds[i].b);
        } else if (bonds[i].b == a) {
            padjacent_id[a].push_back(bonds[i].a);
        } else continue;
        rearrange(padjacent_id, nodes, pnode + a, &plane_vec);
    }

    for (int i=0; i<padjacent_id[a].size(); i++) {
        padjacent[a].push_back(nodes[i]);
    }
}

// 找到所有原子与之相邻的点
#ifdef USE_CUDA
void model::generate_adjacent() {
    node *nodes = this->nodes.get_gpu_data();
    bond *bonds = this->bonds.get_gpu_data();
    nano::s_vector<node> *adjacents = this->adjacents.get_gpu_data();
    nano::s_vector<int> *adjacents_id = this->adjacents_id.get_gpu_data();

    this->nodes.gpu_synchro();
    this->adjacents_id.gpu_synchro();
    this->adjacents.size() = this->adjacents_id.size() = this->nodes.size();

    cudaGenerate_adjacent<<<this->nodes.size(), 1>>>(nodes, bonds,
        this->bonds.size(), adjacents, adjacents_id, nodes);

    this->adjacents.cpu_synchro();
    this->adjacents_id.cpu_synchro();
}
#endif

#ifndef USE_CUDA
void model::generate_adjacent() {
    this->adjacents.size() = this->adjacents_id.size() = this->nodes.size();

    for (int a = 0; a < this->nodes.size(); a++) {
        cudaGenerate_adjacent(&this->nodes[a], &this->bonds[0], this->bonds.size(),
            &this->adjacents[a], &this->adjacents_id[a], &this->nodes[0]);
    }
}
#endif