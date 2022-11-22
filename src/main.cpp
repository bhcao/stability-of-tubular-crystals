#include <iostream>
#include <string>

#include "molecule.h"
#include "model.h"
#include "nmath.h"

int main(int argc, char* argv[]) {
    for (int climb=-3; climb<4; climb++)
    for (int glide=0; glide<5; glide++) {

        // 默认模型参数 default_para
        para ppara = default_para;
        // 修改模型参数，比如 m,n
        ppara.n = 9;
        ppara.repeat = 4;
        ppara.glide = glide;
        ppara.climb = climb;
        ppara.tau = 1e-5;
        // 以 ppara 参数生成模型
        // 模型已经生成，之后再修改模型参数也没用
        model this_model(ppara);
    
        // 更改运行时参数
        this_model.set_precision(1e-9);
        this_model.set_step(1e-7);
    
        nano::pair<double> x = this_model.dis_pair_distance();
        std::cout << "位错距离为 x=" << x[0] << ", z=" << x[1] << std::endl;

        // 文件名
        std::string name = "test_" + std::to_string(glide) + "_" + std::to_string(climb);

        this_model.dump(name, nano::DATA_FILE | nano::EMPHASIS);
        this_model.dump(name, nano::LAN_FORCE | nano::P_ENERGY | nano::EMPHASIS);
    
        for (int k=0; k<100; k++) {
            // 模型计算更新
            this_model.update();
            if (k%10==0) {
                // 输出当前模型至 dump 文件
                this_model.dump(name, nano::LAN_FORCE | nano::P_ENERGY | nano::EMPHASIS);
                // 计算总能量并增加至图片
                std::cout << this_model.total_energy();
            }
        }
    }

    return 0;
}
