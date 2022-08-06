/********************************************************
 *  命名约定：类型实例前面加 p，类型矢量实例后面加 s
 *    中心点为 center，相邻点为 other，所有相邻点为 others
 ********************************************************/

#ifndef _MODEL_H_
#define _MODEL_H_

#include <fstream>
#include <vector>

#include "node.h"
using std::vector;

// 符号常量 PI
#define PI 3.141592653589793

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

// 通过数组储存相邻结构
typedef struct {
    vector<int> others_id;
    vector<node> others;
} adjacent;

// 集合了节点、键、相邻的类
class model {
public:
    // 初始化生成 nodes、bonds、adjacents
    inline model(para ppara): ppara(ppara),
        time(0), step(2e-5), precision(1e-6), range(1.0) {
        generate_nodes();
        generate_bonds();
        glide_bond();
        climb_bond();
        generate_adjacent();
    }
    
    para ppara;
    vector<node> nodes;
    vector<bond> bonds;
    vector<adjacent> adjacents;
    
    double step;      // 更新的步长
    int time;         // 当前时间（运行次数）
    double precision; // 求导精度（增加的 delta 长度）
    double range;     // 打乱的范围

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
    
    // 打乱所有坐标，进行退火计算
    void disorganize();
    
    // 系统整体的能量
    double total_energy();
    // 系统更新一次
    void update();
    
    // 输出至文件
    void dump(std::ofstream &file);
    
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
    
    // 求两个键之间的键能
    double bond_energy(node p1, node p2);
    // 求一个点周围曲率能量，即能量函数第二项
    double curve_energy(node center, vector<node> others);
    // 生成一个点周围的能量
    double round_energy(node center, vector<node> others);
};

#endif // _MODEL_H_
