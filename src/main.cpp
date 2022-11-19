#include <fstream>
#include <cstdio>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "molecule.h"
#include "model.h"

void thread(int glide, int n) {
    // 默认模型参数 default_para
    para ppara = default_para;
    // 修改模型参数，比如 m,n
    ppara.n = n;
    ppara.repeat = 4;
    ppara.glide = -glide;
    ppara.tau = 1e-5;
    // 以 ppara 参数生成模型
    // 模型已经生成，之后再修改模型参数也没用
    model this_model(ppara);
    
    // 更改运行时参数
    this_model.precision = 1e-9;
    this_model.step = 1e-7;
    
    double x, z;
    this_model.dis_pair_distance(&x, &z);
    std::cout << "位错距离为 x=" << x << ", z=" << z << std::endl;

    // 退火
    // this_model.disorganize();
    
    // dump 文件
    char name[20];
    std::sprintf(name, "test_%d_%d.dump", glide, n);
    std::ofstream fout(name);
    // data 文件
    std::sprintf(name, "test_%d_%d.data", glide, n);
    std::ofstream datafile(name);
    // data 文件只能输出一次
    this_model.dump(datafile, DATA_FILE);
    this_model.dump(fout, DUMP_FILE);
    
    // 声明图片类
    figure energy_change;
    for (int k=0; k<10000; k++) {
        // 模型计算更新
        this_model.update();
        if (k%100==0) {
            // 输出当前模型至 dump 文件
            this_model.dump(fout, DUMP_FILE);
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
