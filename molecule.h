#ifndef NANO_MOLECULE_H_
#define NANO_MOLECULE_H_

#include <fstream>
#include <cmath>
#include <vector>

#include "vector.h"

// 符号常量 PI
#define PI 3.141592653589793

// 节点笛卡尔坐标
class node {
    // 重载相减、相加
    __device__ friend inline node operator-(node p1, node p2) {
        return {p1.i-p2.i, p1.j-p2.j, p1.k-p2.k};
    }

    __device__ friend inline node operator+(node p1, node p2) {
        return {p1.i+p2.i, p1.j+p2.j, p1.k+p2.k};
    }

    // 数量积两种写法
    __device__ friend inline node operator*(double n1, node p2) {
        return {n1*p2.i, n1*p2.j, n1*p2.k};
    }

    __device__ friend inline node operator*(node p2, double n1) {
        return {n1*p2.i, n1*p2.j, n1*p2.k};
    }
    
    // 除法，只能除数
    __device__ friend inline node operator/(node p2, double n1) {
        return {p2.i/n1, p2.j/n1, p2.k/n1};
    }

    // 内积
    __device__ friend inline double operator*(node p1, node p2) {
        return {p1.i*p2.i + p1.j*p2.j + p1.k*p2.k};
    }

    // 两个向量之间叉乘
    __device__ friend inline node cross(node p1, node p2) {
        return { p1.j*p2.k - p1.k*p2.j, p1.k*p2.i - p1.i*p2.k,
            p1.i*p2.j - p1.j*p2.i};
    }

    // node 间的距离
    __device__ friend inline double dist(node p1, node p2 = {0}) {
        using std::sqrt;
        node temp = p1 - p2;
        return sqrt(temp.i*temp.i + temp.j*temp.j + temp.k*temp.k);
    }
    
public:
    node() = default;
    double i, j, k;
    int id;
};

// 两个节点之间的键
class bond {
    // 判断两个键是否相等
    friend inline bool operator==(bond b1, bond b2) {
        return (b1.a == b2.a && b1.b == b2.b) || (b1.a == b2.b && b1.b == b2.a);
    }

public:
    bond() = default;
    int a, b;
};

typedef struct {
    int begin[2], end[2];
} dis_pair;

// 集合了节点、键、相邻的类
class molecule {
public:
    // 初始化生成 nodes、bonds、adjacents
    inline molecule(int nodes_len, int bonds_len): time(0), step(2e-5),
        precision(1e-6), range(1.0), nodes(nodes_len), bonds(bonds_len) {}

    nano::vector<node> nodes;
    nano::vector<bond> bonds;
    
    dis_pair pdis_pair;

    double step;      // 更新的步长
    int time;         // 当前时间（运行次数）
    double precision; // 求导精度（增加的 delta 长度）
    double range;     // 打乱的范围
    
    // 打乱所有坐标，进行退火计算
    void disorganize();
    
    // 输出至文件
    void dump(std::ofstream &file);
};

class figure {
    // 重载输入运算符
    inline friend figure& operator<<(figure &p, double n) {
        p.data.push_back(n);
        return p;
    }
    
public:
    figure() = default;
    
    // 绘制函数
    void draw(const char* name);

private:
    std::vector<double> data;
};

#endif // NANO_MOLECULE_H_