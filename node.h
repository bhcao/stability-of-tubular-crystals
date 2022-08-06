#ifndef _NODE_H_
#define _NODE_H_

#include <cmath>

// 两个节点之间的键
typedef struct {
    int a, b;
} bond;

// 节点笛卡尔坐标
typedef struct {
    double i, j, k;
} node;

// 重载相减、相加
inline node operator-(node p1, node p2) {
    return {p1.i-p2.i, p1.j-p2.j, p1.k-p2.k};
}

inline node operator+(node p1, node p2) {
    return {p1.i+p2.i, p1.j+p2.j, p1.k+p2.k};
}

// 数量积两种写法
inline node operator*(double n1, node p2) {
    return {n1*p2.i, n1*p2.j, n1*p2.k};
}

inline node operator*(node p2, double n1) {
    return {n1*p2.i, n1*p2.j, n1*p2.k};
}

// 内积
inline double operator*(node p1, node p2) {
    return {p1.i*p2.i + p1.j*p2.j + p1.k*p2.k};
}

// 两个向量之间叉乘
inline node cross(node p1, node p2) {
    return { p1.j*p2.k - p1.k*p2.j, p1.k*p2.i - p1.i*p2.k,
        p1.i*p2.j - p1.j*p2.i};
}

// node 间的距离
inline double dist(node p1, node p2 = {0}) {
    using std::sqrt;
    node temp = p1 - p2;
    return sqrt(temp.i*temp.i + temp.j*temp.j + temp.k*temp.k);
}

// 判断两个键是否相等
inline bool operator==(bond b1, bond b2) {
    return (b1.a == b2.a && b1.b == b2.b) || (b1.a == b2.b && b1.b == b2.a);
}

#endif // _NODE_H_
