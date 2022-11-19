/* -*- c++ -*- ----------------------------------------------------------
   AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)

   molecule 子类，建模函数，调用 molecule update

   命名约定：类型实例前面加 p，类型矢量实例后面加 s
      中心点为 center，相邻点为 other，所有相邻点为 others
------------------------------------------------------------------------- */

#ifndef NADO_MODEL_H_
#define NANO_MODEL_H_

#include <cmath>

#include "nmath.h"
#include "molecule.h"

// 参数默认值
#define default_para {13,13,0.1,1,3,3,3,1,3e-7,1,1,273}

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
    // 模型构建参数
    int m, n;
    int direction;
    int glide, climb;
    int repeat;

    // 能量参数（键原长、系数）
    double rest_len, k, tau;

    // 

} para;

inline nano::vector average(nano::vector p1, nano::vector p2) {
    return (p1 + p2)/2;
}

typedef struct {
    int begin[2], end[2];
} dis_pair;

// 集合了节点、键、相邻的类
class model: public molecule {
public:
    // 初始化生成 nodes、bonds、adjacents
    model(para ppara);

private:
    dis_pair pdis_pair;

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
        nano::vector begin = average(this->nodes[this->pdis_pair.begin[0]],
            this->nodes[this->pdis_pair.begin[1]]);
        nano::vector end = average(this->nodes[this->pdis_pair.end[0]],
            this->nodes[this->pdis_pair.end[1]]);
        nano::vector temp = begin - end;
        double m = (double)this->ppara.m;
        double n = (double)this->ppara.n;
        double r = this->ppara.rest_len * 1/PI/2 * std::sqrt(m*m+n*n-m*n);
        double half_dist = std::sqrt(temp[0]*temp[0]+temp[1]*temp[1])/2;
        *x = r * std::asin(half_dist/r);
        *z = temp[2];
    }

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
