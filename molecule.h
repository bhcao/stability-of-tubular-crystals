#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include <vector>
#include <fstream>
#include <cmath>
using std::vector;

// 符号常量 PI
#define PI 3.141592653589793

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

// 节点笛卡尔坐标
class node {
    // 重载相减、相加
    friend inline node operator-(node p1, node p2) {
        return {p1.i-p2.i, p1.j-p2.j, p1.k-p2.k};
    }

    friend inline node operator+(node p1, node p2) {
        return {p1.i+p2.i, p1.j+p2.j, p1.k+p2.k};
    }

    // 数量积两种写法
    friend inline node operator*(double n1, node p2) {
        return {n1*p2.i, n1*p2.j, n1*p2.k};
    }

    friend inline node operator*(node p2, double n1) {
        return {n1*p2.i, n1*p2.j, n1*p2.k};
    }

    // 内积
    friend inline double operator*(node p1, node p2) {
        return {p1.i*p2.i + p1.j*p2.j + p1.k*p2.k};
    }

    // 两个向量之间叉乘
    friend inline node cross(node p1, node p2) {
        return { p1.j*p2.k - p1.k*p2.j, p1.k*p2.i - p1.i*p2.k,
            p1.i*p2.j - p1.j*p2.i};
    }

    // node 间的距离
    friend inline double dist(node p1, node p2 = {0}) {
        using std::sqrt;
        node temp = p1 - p2;
        return sqrt(temp.i*temp.i + temp.j*temp.j + temp.k*temp.k);
    }

public:
    node() = default;
    double i, j, k;
    int id;
};

// 集合了节点、键、相邻的类
class molecule {
public:
    // 初始化生成 nodes、bonds、adjacents
    inline molecule(): time(0), step(2e-5),
        precision(1e-6), range(1.0) {}

    vector<node> nodes;
    vector<bond> bonds;
    
    double step;      // 更新的步长
    int time;         // 当前时间（运行次数）
    double precision; // 求导精度（增加的 delta 长度）
    double range;     // 打乱的范围
    
    // 打乱所有坐标，进行退火计算
    void disorganize();
    
    // 系统更新一次
    void update();
    
    // 输出至文件
    void dump(std::ofstream &file);
    
    // 通过能量求导，纯虚函数
    virtual double energy(node center) = 0;
};

#endif // _MOLECULE_H_