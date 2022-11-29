/* -*- c++ -*- ----------------------------------------------------------
   AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)

   定义了基础数学类，包括
   1. 三维矢量 vector 及其运算（加减、数乘、内外积、模）；
   2. pair（pair<int> 做键、二维点阵坐标，pair<double> 做位错距输出）；
   3. 数学常数 PI、K_B，高斯随机数 rand_vector。
------------------------------------------------------------------------- */

#ifndef NANO_NMATH_H_
#define NANO_NMATH_H_

#include <cmath>
#include <mutex>

#include "narray.h"

namespace nano {

// 三维矢量，表示位置、速度等
class vector {
public:
    inline vector(): data{0} {}
    inline vector(double i, double j, double k): data{i, j, k} {}
    inline double& operator[](int i) { return this->data[i];}

private:
    double data[3];
};

// 重载相减、相加、自增
inline vector operator+(vector p1, vector p2) {
    return vector(p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]);
}
inline vector operator-(vector p) { return vector(-p[0], -p[1], -p[2]); }
inline vector operator-(vector p1, vector p2) { return p1 + (-p2); }
inline void operator+=(vector &p, vector p_a) { p = p + p_a; }

// 数量积两种写法
inline vector operator*(double n, vector p) {
    return vector(n * p[0], n * p[1], n * p[2]);
}
inline vector operator*(vector p, double n) { return n * p; }
    
// 除法，只能除数
inline vector operator/(vector p, double n) {
    return vector(p[0] / n, p[1] / n, p[2] / n);
}

// 内积
inline double operator*(vector p1, vector p2) {
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

// 两个向量之间叉乘
inline vector cross(vector p1, vector p2) {
    return vector( p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2],
        p1[0]*p2[1] - p1[1]*p2[0]);
}

// 两个向量比较
inline bool operator==(vector p1, vector p2) {
    return p1[0]==p2[0] && p1[1]==p2[1] && p1[2]==p2[2];
}

// vector 间的距离
inline double mod(vector p) { return std::sqrt(p * p); }

// 模板仅仅是为了位错对（注意 T 为非数会报错）
template <typename T> class pair {
public:
    inline pair(): data{0} {}
    inline pair(T a, T b): data{a, b} {}
    
    // 找键中是否连接粒子 a，无返回 -1，有返回另一个粒子
    inline T find(T a) {
        if (this->data[0] == a) return this->data[1];
        if (this->data[1] == a) return this->data[0];
        return -1;
    }

    // 替换 a -> b，失败 false，成功 true
    inline T replace(T a, T b) {
        if (this->data[0] == a) {
            this->data[0] = b; return true;
        }
        if (this->data[1] == a) {
            this->data[1] = b; return true;
        }
        return false;
    }

    inline T& operator[](int i) { return this->data[i]; }

private:
    T data[2];
};

// 比较、加减法（加减法是为了代替二维坐标）
template <typename T> inline bool operator==(pair<T> b1, pair<T> b2) {
    return (b1[0] == b2[0] && b1[1] == b2[1]) || (b1[1] == b2[0] && b1[0] == b2[1]);
}
template <typename T> inline pair<T> operator-(pair<T> b1, pair<T> b2) {
    return pair<T>(b1[0] - b2[0], b1[1] - b2[1]);
}
template <typename T> inline pair<T> operator+(pair<T> b1, pair<T> b2) {
    return pair<T>(b1[0] + b2[0], b1[1] + b2[1]);
}

// 常数：圆周率、玻尔兹曼常数
#define PI 3.141592653589793
#define K_B 1.380649e-23

// 随机数种子
struct rand_seed {
    std::mutex lock;  // 互斥锁
    uint64_t seed;    // 线性同余法随机种子
    double V1, V2, S; // 高斯随机数种子
    bool phase;       // 高斯随机数相
};

#ifndef RAND_POOL_SIZE
#define RAND_POOL_SIZE 32
#endif

// 高斯分布的随机位置
class rand_pool {
public:
    rand_pool();
    nano::vector gen_vector();
private:
    double rand(int i);            // 使用第 i 个随机数种子生成线性随机数（0-1）
    double gauss_rand(int i);      // 使用第 i 个随机数种子生成高斯随机数
    nano::darray<rand_seed> pool;  // 随机池，大小以后设成动态
};

}

#endif // NANO_NMATH_H_