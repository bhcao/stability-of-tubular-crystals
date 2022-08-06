#include <fstream>

#include "node.h"
#include "model.h"
#include "figure.h"

/***********************************************************
下面是参数说明，等号后面是参数默认值

模型构建参数：
    int m = 13, n = 13;         // 叶序螺旋数
    double rest_len = 0.1;      // 键的平衡长度
    int direction = 1;          // 位错发生的方位
    int glide = 3, climb = 3;   // 位错滑移和攀移的步数
    int repeat = 3;             // 模型重复的次数（增长模型）
    double k = 1, tau = 3e-7;   // 能量参数（第一项键能、第二项曲率）

运行时参数：
    double step = 2e-5;      // 更新的步长，就是向梯度方向前进多少
    double precision = 1e-6; // 求导精度（求导时增加的 delta 长度）
    double range = 1.0;     // 打乱的范围（退火时粒子随机移动的范围）
************************************************************/

int main() {
    
    // 默认模型参数 default_para
    para ppara = default_para;
    // 修改模型参数，比如 m,n
    ppara.m = 12;
    
    // 以 ppara 参数生成模型
    // 模型已经生成，之后再修改模型参数也没用
    model this_model(ppara);
    
    // 更改运行时参数
    this_model.precision = 1e-7;
    this_model.step = 1e-8;
    
    // 退火
    this_model.disorganize();
    
    // dump 文件
    std::ofstream fout("test.dump");
    
    // 声明图片类
    figure energy_change;
    for (int k=0; k<1000; k++) {
        // 模型计算更新
        this_model.update();
        if (k%10==0) {
            // 输出当前模型至 dump 文件
            this_model.dump(fout);
            // 计算总能量并增加至图片
            energy_change << this_model.total_energy();
        }
    }
    // 画图
    energy_change.draw("\"energy.png\"");
    
    // 保存 dump 文件
    fout.close();
    
    return 0;
}