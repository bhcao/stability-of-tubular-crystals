/*********************************************************************//**
 * AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)
 * 
 * @file include/model.h
 * @brief 建模类，molecule 子类
 * @version 0.2.0
 * 
 * molecule 子类，建模函数，调用 molecule update
 *
 * 命名约定：类型实例前面加 p，类型矢量实例后面加 s
 *    中心点为 center，相邻点为 other，所有相邻点为 others
 **************************************************************************/

#ifndef NADO_MODEL_H_
#define NANO_MODEL_H_

#include <cmath>

#include "nmath.h"
#include "narray.h"
#include "molecule.h"

//! 参数默认值
#define default_para {13,13,3,1,3,3,0,0.1,1,3e-7,1e-6,1,1,273}

//! 参数结构体，初始化参数
typedef struct {
    // 模型参数
    int m;  //!< 叶序数 m
    int n;  //!< 叶序数 n
    int repeat;  //!< 重复次数
    int direction; //!< 位错发生方向，正负号表示手性
    int glide;  //!< 滑移步数
    int climb;  //!< 攀移步数
    int bn; //!< 位错在管上的初始位置

    // 能量参数
    double rest_len; //!< 键原长
    double k;   //!< 键能系数
    double tau; //!< 曲率能系数
    
    // 运行参数
    double step; //!< 更新的步长（要调节用 set 函数）
    double mass; //!< 粒子质量
    double damp; //!< 朗之万阻尼系数
    double tempr; //!< 系统的温度
} para;

//! 初始化建模类
class model: public molecule {
public:
    //! 初始化生成 nodes、bonds、adjacents
    model(para ppara, int argc, char* argv[]);
    
    //! 位错对间距
    inline nano::pair<double> dis_pair_distance() {
        nano::vector temp = (this->nodes[this->begin[0]] + this->nodes[this->begin[1]] -
            this->nodes[this->end[0]] - this->nodes[this->end[1]]) / 2;
        double m = (double)this->ppara.m, n = (double)this->ppara.n;
        double r = this->ppara.rest_len * 1/PI/2 * std::sqrt(m*m+n*n-m*n);
        double half_dist = std::sqrt(temp[0]*temp[0]+temp[1]*temp[1])/2;
        return nano::pair<double>(r * std::asin(half_dist/r), temp[2]);
    }

private:
    //! 位错对
    int begin[2], end[2];
    int checkpoint; //!< 检查点，在滑移攀移中有作用

    para ppara; //!< 参数，重复很多，仅仅是为了方便使用

    //! 由节点的二维坐标计算标识符
    inline int flat(int i, int j) {
        j += this->ppara.m * (i/this->ppara.n);
        i %= this->ppara.n;
        if (i < 0) {
            i += this->ppara.n;
            j -= this->ppara.m;
        }
        j %= this->ppara.m * this->ppara.repeat;
        if (j < 0)
            j += this->ppara.m * this->ppara.repeat;
        return i + this->ppara.n*j;
    }
    inline int flat(nano::pair<int> i) { return flat(i[0], i[1]); }
    inline int flat(int i[2]) { return flat(i[0], i[1]); }

    //! ck 关于 a,b 的对偶点，就是既与 a 又与 b 相连的两个点中的另一个
    inline int dual(int ck, int a, int b) {
        int i = this->adjacents[a].find(b);
        if (i == -1)
            return -1; // 错误（提示，其他错误不会检查）
        // i 左右邻居
        int il = (i == 0) ? this->adjacents[a].size()-1 : i-1;
        int ir = (i == this->adjacents[a].size()-1) ? 0 : i+1;
        return (this->adjacents[a][il] == ck) ?
            this->adjacents[a][ir] : this->adjacents[a][il];
    }
    
    //! 移除键
    inline void remove(nano::pair<int> bond) {
        this->adjacents[bond[0]].remove(bond[1]);
        this->adjacents[bond[1]].remove(bond[0]);
        this->bonds.remove(bond);
    }

    void remove(int node); //!< 移除粒子及其相连键
    //! 增添键，n1 是 bond[0] 插入位置，n2 是 bond[1] 插入位置
    inline void insert(nano::pair<int> bond, int n1, int n2) {
        this->bonds.push_back(bond);
        this->adjacents[bond[0]].insert(bond[1], n1);
        this->adjacents[bond[1]].insert(bond[0], n2);
    }
    
    //! 在位置 p 插入粒子，返回粒子序号
    inline int insert(nano::vector p) {
        this->nodes.push_back(p);
        return this->nodes.size() - 1;
    }
    
    //! 在键中间的位置
    inline int between(int i, int j, int k) {
        // 键插入位置
        int m = this->adjacents[i].find(j);
        int n = this->adjacents[i].find(k);
        // m=0 在 n 后，n=0 在 m 后，否则在中间
        return (m == 0) ? n+1 : ((n == 0) ? m+1 : ((n>m) ? n : m));
    }

    //! 生成点位置
    void perfect_model_position();
    //! 键、邻接（拓扑结构）
    void perfect_model_topology();
    //! 滑移和攀移（根据拓扑结构）
    void add_dislocation();
};

#endif // NANO_MODEL_H_
