/********************************************************
 *  命名约定：类型实例前面加 p，类型矢量实例后面加 s
 *    中心点为 center，相邻点为 other，所有相邻点为 others
 ********************************************************/

#ifndef NADO_MODEL_H_
#define NANO_MODEL_H_

#include <cmath>
#include <fstream>

#include "vector.h"
#include "molecule.h"

// 参数默认值
#define default_para {13,13,0.1,1,3,3,3,1,3e-7,1,1,273}
#define K_B 1.380649e-23 // 玻尔兹曼常数

class pos2d {
    friend pos2d operator+(pos2d &p1, pos2d &p2) {
        return {p1.i+p2.i, p1.j+p2.j};
    }
    friend pos2d operator-(pos2d &p1, pos2d &p2) {
        return {p1.i-p2.i, p1.j-p2.j};
    }
public:
    int i, j;
};

// 参数
typedef struct {
    int m, n;
    double rest_len;
    int direction;
    int glide, climb;
    int repeat;
    double k, tau;
    double damp;
    double mass;
    double tempr;
} para;

inline node average(node p1, node p2) {
    return {(p1.i+p2.i)/2, (p1.j+p2.j)/2, (p1.k+p2.k)/2};
}

// 集合了节点、键、相邻的类
class model: public molecule {
public:
    // 初始化生成 nodes、bonds、adjacents
    model(para ppara);
    
    para ppara;

    // 通过数组储存相邻结构
    nano::vector<nano::s_vector<int>> adjacents_id;
    nano::vector<nano::s_vector<node>> adjacents;

    // 由节点的二维坐标计算标识符
    inline int flat(int i, int j) {
        int i2 = i;
        i %= this->ppara.n;
        j += this->ppara.m * (i2/this->ppara.n);
        if (i < 0) {
            i += this->ppara.n;
            j -= this->ppara.m;
        }
        j %= this->ppara.m * this->ppara.repeat;
        if (j < 0) {
            j += this->ppara.m * this->ppara.repeat;
        }
        return i + this->ppara.n*j;
    }
    inline int flat(pos2d i) { return flat(i.i, i.j); }
    inline int flat(int i[2]) { return flat(i[0], i[1]); }

    // f* nvidia... Cuda 函数没法用
    inline void dis_pair_distance(double *x, double *z) {
        node begin = average(this->nodes[this->pdis_pair.begin[0]],
            this->nodes[this->pdis_pair.begin[1]]);
        node end = average(this->nodes[this->pdis_pair.end[0]],
            this->nodes[this->pdis_pair.end[1]]);
        node temp = {begin.i - end.i, begin.j - end.j, begin.k - end.k};
        double m = (double)this->ppara.m;
        double n = (double)this->ppara.n;
        double r = this->ppara.rest_len * 1/PI/2 * std::sqrt(m*m+n*n-m*n);
        double half_dist = std::sqrt(temp.i*temp.i+temp.j*temp.j)/2;
        *x = r * std::asin(half_dist/r);
        *z = temp.k;
    }

    // 更新函数
    void update();

    // 系统整体的能量
    double total_energy();
    
private:
    // 生成所有原子与之相邻的点
    void generate_adjacent();
    // 生成点
    void generate_nodes();
    // 生成键
    void generate_bonds();
    // 滑移和攀移
    void glide_climb();
};

#endif // NANO_MODEL_H_
