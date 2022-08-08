/********************************************************
 *  命名约定：类型实例前面加 p，类型矢量实例后面加 s
 *    中心点为 center，相邻点为 other，所有相邻点为 others
 ********************************************************/

#ifndef _MODEL_H_
#define _MODEL_H_

#include <fstream>

#include "vector.h"
#include "molecule.h"

// 参数默认值
#define default_para {13,13,0.1,1,3,3,3,1,3e-7}

// 参数
typedef struct {
    int m, n;
    double rest_len;
    int direction;
    int glide, climb;
    int repeat;
    double k, tau;
} para;

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
        i %= this->ppara.n;
        j += this->ppara.m * (i/this->ppara.n);
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
    
    // 更新函数
    void update();

    // 系统整体的能量
    // double total_energy();
    
private:
    // 生成所有原子与之相邻的点
    void generate_adjacent();
    // 生成点
    void generate_nodes();
    // 生成键
    void generate_bonds();
    // 滑移
    void glide_bond();
    // 攀移
    void climb_bond();
    // 交换两根键
    void replace_bond(bond b1, bond b2);
};

#endif // _MODEL_H_
