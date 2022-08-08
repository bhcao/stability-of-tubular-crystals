#include <fstream>
#include <cstdio>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "molecule.h"
#include "model.h"

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

void thread(int glide, int n) {
    // 默认模型参数 default_para
    para ppara = default_para;
    // 修改模型参数，比如 m,n
    ppara.n = n;
    ppara.repeat = 100;
    ppara.glide = glide;
    // 以 ppara 参数生成模型
    // 模型已经生成，之后再修改模型参数也没用
    model this_model(ppara);
    
    // 更改运行时参数
    this_model.precision = 1e-7;
    this_model.step = 1e-8;
    
    // 退火
    this_model.disorganize();
    
    // dump 文件
    char name[20];
    std::sprintf(name, "test_%d_%d.dump", glide, n);
    std::ofstream fout(name);
    
    // 声明图片类
    figure energy_change;
    for (int k=0; k<100; k++) {
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
    std::sprintf(name, "energy_%d_%d.png", glide, n);
    energy_change.draw(name);
    
    // 保存 dump 文件
    fout.close();
}

int main(int argc, char* argv[]) {

    #ifdef USE_MPI
    int size;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int task = 35 / size + 1;
    int count = 0;
    #endif

    for (int n=7; n<14; n++)
    for (int glide=0; glide<5; glide++) {

        #ifdef USE_MPI
        if (count % task == rank) {
        #endif

        thread(glide, n);

        #ifdef USE_MPI
        }
        #endif
    
    }

    #ifdef USE_MPI
    MPI_Finalize();
    #endif

    return 0;
}
