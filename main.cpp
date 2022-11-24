#include <iostream>
#include <string>

#include "molecule.h"
#include "nmath.h"

// #define RESTART

// 如果声明了重启，只导入能量函数；否则导入模型类
#ifdef RESTART
double my_bond_energy(nano::vector, nano::vector, nano::sarray<double>);
double my_node_energy(nano::vector, nano::sarray<nano::vector>, nano::sarray<double>);
#else
#include "model.h"
#endif

#define PARA_rest_len 0
#define PARA_k 1
#define PARA_tau 2

int main(int argc, char* argv[]) {

// 如果声明了重启，调用直接 molecule 初始化，否则调用 model 初始化
#ifdef RESTART
    molecule this_model("restart.bin", &my_node_energy, &my_bond_energy);
#else
    // 默认模型参数 default_para
    para ppara = default_para;
    // 修改模型参数，比如 m,n
    ppara.glide = 3;
    ppara.climb = 2;
    // 以 ppara 参数生成模型
    model this_model(ppara);
#endif

    // 更改运行时参数
    this_model.set_precision(1e-9);
    this_model.set_step(2e-5);
    this_model.set_paras(1, PARA_rest_len);

#ifndef RESTART
    // 输出位错距离
    nano::pair<double> x = this_model.dis_pair_distance();
    std::cout << "位错距离为 x=" << x[0] << ", z=" << x[1] << std::endl;

    // 输出初始位型
    this_model.dump("test", nano::DATA_FILE | nano::EMPHASIS);
    this_model.dump("test", nano::LAN_FORCE | nano::P_ENERGY | nano::EMPHASIS);
#endif
    
    for (int k=0; k<1000; k++) {
        // 模型计算更新
        this_model.update();
        if (k%100==0) {
            this_model.dump("test", nano::LAN_FORCE | nano::P_ENERGY | nano::EMPHASIS);
            std::cout << this_model.total_energy() << std::endl; // 输出总能量
        }
    }

    // 储存文件以便继续运行
    this_model.store("restart.bin");

    return 0;
}
