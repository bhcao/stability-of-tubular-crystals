/*********************************************************************//**
 * AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)
 * 
 * @file main.cpp
 * @brief 主函数
 * @version 0.2.0
 * 
 * 主函数示例，请根据需要自行更改
 **************************************************************************/

#include <iostream>
#include <string>
#include <cmath>

#define DYNAMICS 2 // 过阻尼朗之万

#include "molecule.h"
#include "nmath.h"
#include "energy.h"

// #define RESTART

// 如果声明了重启，只导入能量函数；否则导入模型类
#ifndef RESTART
#include "model.h"
#endif

#define PARA_rest_len 0
#define PARA_k 1
#define PARA_tau 2

int main(int argc, char* argv[]) {

// 如果声明了重启，调用直接 molecule 初始化，否则调用 model 初始化
#ifdef RESTART
    molecule this_model("restart.bin", energy_func::node_energy_func, energy_func::bond_energy_func, argc, argv);
#else
    // 默认模型参数 default_para
    para p = default_para;
    // 修改模型参数，比如 m,n
    p.m = 13, p.n = 13, p.repeat = 4;
    p.direction = 3, p.glide = 2, p.climb = 2;
    p.rest_len = 1, p.k = std::sqrt(3)/2, p.tau = 0.1;
    p.damp = 1, p.tempr = 0;
    // 以 p 参数生成模型
    model this_model(p, argc, argv);
#endif
    
    // 输出总能量，设置 z 方向计入的范围，x、y 设无穷即可
    #define OUT_ENERGY this_model.total_energy(nano::vector(-INFINITY, -INFINITY, 7), \
        nano::vector(INFINITY, INFINITY, 30))
    #define OUT_DUMP this_model.dump("test", nano::DATA_FILE | nano::VELOCITY | nano::P_ENERGY | nano::MEAN_CURVE \
        | nano::GAUSS_CURVE | nano::EMPHASIS);

    // 更改运行时参数
    this_model.set_step(1e-6);

#ifndef RESTART
    // 输出位错距离
    nano::pair<double> x = this_model.dis_pair_distance();
    std::cout << "位错距离为 x=" << x[0] << ", z=" << x[1] << std::endl;
#endif

    for (int k=200000; k>0; k--) {
        if (k%10000==0) {
            std::cout << OUT_ENERGY << std::endl;
            OUT_DUMP
        }
        // this_model.set_tempr((k+5000)*1e19);
        // 模型计算更新
        this_model.update();
    }

    // 储存文件以便继续运行
    this_model.store("restart.bin");

    return 0;
}
