/*********************************************************************//**
 * AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)
 * 
 * @file src/molecule.cpp
 * @brief 运行类的实现
 * @version 0.2.0
 **************************************************************************/

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <stdio.h>

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include "nmath.h"
#include "narray.h"
#include "molecule.h"
#include "energy.h"

double molecule::local_energy(nano::vector center, nano::sarray<nano::vector> others) {
    double count = 0;
    for (int i=0; i<others.size(); i++) {
        count += this->bond_energy(center, others[i], this->paras)/2;
    }
    count += this->node_energy(center, others, this->paras);
    return count;
}

// 借助 local_energy 键要算两次，效率低，但可以省去键结构，逻辑也更清晰
double molecule::total_energy(nano::vector range_l, nano::vector range_r) {
    double count = 0;
    // 以粒子为中心的能量
    for (int i=0; i<this->nodes.size(); i++) {
        // 相邻粒子转化成坐标
        nano::sarray<nano::vector> adjacent;
        for (int j=0; j<this->adjacents[i].size(); j++)
            adjacent.push_back(this->nodes[this->adjacents[i][j]]);
        if (this->nodes[i][0] > range_l[0] && this->nodes[i][0] < range_r[0] &&
            this->nodes[i][1] > range_l[1] && this->nodes[i][1] < range_r[1] &&
            this->nodes[i][2] > range_l[2] && this->nodes[i][2] < range_r[2])
            count += local_energy(this->nodes[i], adjacent);
    }
    return count;
}

// 粒子数
int molecule::total_particle(nano::vector range_l, nano::vector range_r) {
    int count = 0;
    for (int i=0; i<this->nodes.size(); i++) {
        if (this->nodes[i][0] > range_l[0] && this->nodes[i][0] < range_r[0] &&
            this->nodes[i][1] > range_l[1] && this->nodes[i][1] < range_r[1] &&
            this->nodes[i][2] > range_l[2] && this->nodes[i][2] < range_r[2])
            count++;
    }
    return count;
}

double molecule::local_energy_for_update(nano::vector center, nano::sarray<nano::vector> others) {
    double count = 0;
    for (int i=0; i<others.size(); i++) {
        count += this->bond_energy(center, others[i], this->paras);
    }
    count += this->node_energy(center, others, this->paras);
    return count;
}

void molecule::update_velocity(int i) {
    nano::sarray<nano::vector> adjacent;
    for (int j=0; j<this->adjacents[i].size(); j++) 
        adjacent.push_back(this->nodes[this->adjacents[i][j]]);

    #if DYNAMICS == 0 // 梯度下降
    this->velocities[i] = -div(this->nodes[i], adjacent);
    #elif DYNAMICS == 1  // 求加速度（朗之万），朗之万方程三项
    nano::vector accelerate = -div(this->nodes[i], adjacent)/this->mass - this->damp*
        this->velocities[i] + std::sqrt(2*this->damp*this->tempr*K_B/this->mass)*
        this->prand_pool.gen_vector();
    this->velocities[i] += accelerate * this->step;
    #elif DYNAMICS == 2  // 过阻尼朗之万
    this->velocities[i] = -div(this->nodes[i], adjacent)/this->damp;// + std::sqrt(2*this->damp*
        // this->tempr*K_B/this->mass)*this->prand_pool.gen_vector()/this->damp;
    #endif
}

void molecule::update() {
#ifdef USE_KOKKOS
    Kokkos::parallel_for(this->nodes.size(), KOKKOS_LAMBDA(int i) {
#else
    for (int i=0; i<this->nodes.size(); i++) {
#endif
        update_velocity(i);
        if (this->adjacents[i].size() == 4) {
            nano::vector velocity;
            for (int j=0; j<4; j++) {
                if (this->adjacents[this->adjacents[i][j]].size() != 4)
                    velocity += this->velocities[this->adjacents[i][j]];
            }
            this->nodes[i] += velocity/2 * this->step;
        } else
            this->nodes[i] += this->velocities[i] * this->step;
    }
#ifdef USE_KOKKOS
    );
#endif
    this->time++;
}

void molecule::border_update(nano::darray<nano::sarray<int>> *rigid) {
    // 质心、惯性张量（主轴近似）、合力、合力矩、角加速度
    nano::vector up_center, up_tensor, up_force, up_moment;
    for (int i=0; i<(*rigid).size(); i++) {
        up_center += this->nodes[(*rigid)[i][2]];
    }
    up_center = up_center / (*rigid).size();
    for (int i=0; i<(*rigid).size(); i++) {
        double rest_len = this->paras[0], k = this->paras[1];
        nano::vector x1 = this->nodes[(*rigid)[i][2]] - this->nodes[(*rigid)[i][1]];
        nano::vector x2 = this->nodes[(*rigid)[i][2]] - this->nodes[(*rigid)[i][0]];
        nano::vector force = k * (rest_len/nano::mod(x1) - 1) * x1 +
            k * (rest_len/nano::mod(x2) - 1) * x2;
        up_force += force;
        up_moment += nano::cross(this->nodes[(*rigid)[i][2]]-up_center, force);
        nano::vector t = this->nodes[(*rigid)[i][2]]-up_center;
        up_tensor += nano::vector(t[1]*t[1]+t[2]*t[2], t[0]*t[0]+t[2]*t[2], t[0]*t[0]+t[1]*t[1]);
    }
    up_tensor = (*rigid).size() * nano::vector(up_moment[0]/up_tensor[0], 
        up_moment[1]/up_tensor[1], up_moment[2]/up_tensor[2]);
    for (int i=0; i<(*rigid).size(); i++) {
        nano::vector t = this->nodes[(*rigid)[i][2]]-up_center;
        this->nodes[(*rigid)[i][2]] += (up_force + nano::cross(up_tensor, t)) *
            this->step / this->damp;
    }
}

// 边缘刚体模型更新
void molecule::update_rigid() {
    // 确定粒子序号
    nano::darray<nano::sarray<int>> up_rigid(50), down_rigid(50);
    for (int i=0; i<this->nodes.size(); i++)
    if (this->adjacents[i].size() == 4) {
        nano::sarray<int> adj;
        for (int j=0; j<4; j++) {
            if (this->adjacents[this->adjacents[i][j]].size() != 4)
                adj.push_back(this->adjacents[i][j]);
        }
        adj.push_back(i);
        if (this->nodes[i][2]<20)
            up_rigid.push_back(adj);
        else
            down_rigid.push_back(adj);
    }
    border_update(&up_rigid);
    border_update(&down_rigid);

#ifdef USE_KOKKOS
    Kokkos::parallel_for(this->nodes.size(), KOKKOS_LAMBDA(int i) {
#else
    for (int i=0; i<this->nodes.size(); i++) {
#endif
        update_velocity(i);
        if (this->adjacents[i].size() != 4)
            this->nodes[i] += this->velocities[i] * this->step;
    }
#ifdef USE_KOKKOS
    );
#endif
    this->time++;
}

// 输出
void molecule::dump(std::string fname, nano::dump_t dump_type) {
    // 设置边界
    double boundary[6] = {0, 30*this->paras[0], 0, 30*this->paras[0], 0, 30*this->paras[0]};
    
    // 强调变色（类型变为 2）
    nano::darray<int> types(this->nodes.size());
    for (int i=0; i<this->nodes.size(); i++) 
        types.push_back(1);
    if (DUMP_CHECK(nano::EMPHASIS, dump_type))
    for (int i=0; i<this->emphasis.size(); i++)
        types[this->emphasis[i]] = 2;

// 如果要求输出 data 文件同时文件还不存在，如存在则跳过
if (DUMP_CHECK(nano::DATA_FILE, dump_type)) {
if (FILE *file = std::fopen((fname + ".data").c_str(), "r")) {
    std::fclose(file); // 判断文件是否存在
} else {
    std::ofstream fout(fname + ".data"); // 以 fname_n.data 打开文件以写入
    // 文件头
    fout << "# Model for nanotube. AUTO generated, DO NOT EDIT\n\n"
        << this->nodes.size() << "\tatoms\n" // 原子、键数
        << this->bonds.size() << "\tbonds\n\n"
        << (DUMP_CHECK(nano::EMPHASIS, dump_type) ? 2 : 1) // 原子、键类型数
        << "\tatom types\n1\tbond types\n\n"
        << boundary[0] << "\t" << boundary[1] << "\txlo xhi\n" // 边界
        << boundary[2] << "\t" << boundary[3] << "\tylo yhi\n"
        << boundary[4] << "\t" << boundary[5] << "\tzlo zhi\n\n"
        << "Masses\n\n" << "1\t" << this->mass << "\n";  // 质量
    if (DUMP_CHECK(nano::EMPHASIS, dump_type))
        fout << "2\t" << this->mass << "\n";

    fout << "\nAtoms\n\n"; // 原子
    for (int i=0; i<this->nodes.size(); i++)
        fout << i << "\t " << types[i] << "\t" // 基础 id type x y z
            << this->nodes[i][0] << '\t' << this->nodes[i][1] << '\t' << this->nodes[i][2] << '\n';

    fout << "\nBonds\n\n"; // 键
    for (int i=0; i<this->bonds.size(); i++)
        fout << i << "\t1\t" << this->bonds[i][0] << '\t' << this->bonds[i][1] << '\n';

    fout.close();
}}

    // 打开文件以写入，文件名为 fname.dump（追加）或 fname.1.data
    std::ofstream fout(fname + ".dump", (this->time != 0) ? std::ios::app : (std::ios::out |
        std::ios::trunc));
    
    // 文件头
    fout << "ITEM: TIMESTEP\n" << this->time << "\nITEM: NUMBER OF ATOMS\n"
        << this->nodes.size() << "\nITEM: BOX BOUNDS ss ss ss\n" // 边界
        << boundary[0] << " " << boundary[1] << "\n"
        << boundary[2] << " " << boundary[3] << "\n"
        << boundary[4] << " " << boundary[5] << "\n"
        << "ITEM: ATOMS id type x y z";
    // 检查是否要输出这些内容
    if (DUMP_CHECK(nano::VELOCITY, dump_type)) fout << " vx vy vz";
    if (DUMP_CHECK(nano::DIV_FORCE, dump_type)) fout << " dfx dfy dfz";
    if (DUMP_CHECK(nano::LAN_FORCE, dump_type)) fout << " fx fy fz";
    if (DUMP_CHECK(nano::K_ENERGY, dump_type)) fout << " ke";
    if (DUMP_CHECK(nano::P_ENERGY, dump_type)) fout << " pe";
    if (DUMP_CHECK(nano::GAUSS_CURVE, dump_type)) fout << " Gaussian_Curvature";
    if (DUMP_CHECK(nano::MEAN_CURVE, dump_type)) fout << " Mean_Curvature";
    fout << "\n";
    
    // 正文数据
    for (int i=0; i<this->nodes.size(); i++) {
        nano::sarray<nano::vector> adjacent;
        for (int j=0; j<this->adjacents[i].size(); j++)
            adjacent.push_back(this->nodes[this->adjacents[i][j]]);
        double size = energy_func::size_around(this->nodes[i], adjacent);
        nano::sarray<double> angle = energy_func::angles_around(this->nodes[i], adjacent);

        fout << i << "\t " << types[i] << "\t" // 基础 id type xs ys zs
            << this->nodes[i][0] << '\t' << this->nodes[i][1] << '\t' << this->nodes[i][2];
        if (DUMP_CHECK(nano::VELOCITY, dump_type))
            fout << "\t" << this->velocities[i][0] << "\t" << this->velocities[i][1]
                << "\t" << this->velocities[i][2];
        if (DUMP_CHECK(nano::DIV_FORCE, dump_type)) {
            nano::vector temp = div(this->nodes[i], adjacent);
            fout << "\t" << temp[0] << "\t" << temp[1] << "\t" << temp[2];
        }
        if (DUMP_CHECK(nano::LAN_FORCE, dump_type)) {
            nano::vector accelerate = -div(this->nodes[i], adjacent)/this->mass - 
                this->damp*this->velocities[i] + std::sqrt(2*this->damp*this->tempr*
                K_B/this->mass)*this->prand_pool.gen_vector();
            nano::vector temp = this->mass * accelerate;
            fout << "\t" << temp[0] << "\t" << temp[1] << "\t" << temp[2];
        }
        if (DUMP_CHECK(nano::K_ENERGY, dump_type))
            fout << "\t" << this->mass/2*nano::mod(this->velocities[i]);
        if (DUMP_CHECK(nano::P_ENERGY, dump_type))
            fout << "\t" << local_energy(this->nodes[i], adjacent);
        if (DUMP_CHECK(nano::GAUSS_CURVE, dump_type))
            fout << "\t" << energy_func::gauss_curvature(this->nodes[i], adjacent, angle, size);
        if (DUMP_CHECK(nano::MEAN_CURVE, dump_type))
            fout << "\t" << energy_func::mean_curvature(this->nodes[i], adjacent, angle, size);
        fout << '\n';
    }
    fout.close();

    return;
} 

molecule::molecule(std::string fname, double(*in_node_energy)(nano::vector, nano::sarray<nano::vector>,
    nano::sarray<double>), double(*in_bond_energy)(nano::vector, nano::vector, nano::sarray<double>),
    int argc, char* argv[]): node_energy(in_node_energy), bond_energy(in_bond_energy) {
    std::ifstream fin(fname, std::ios::in | std::ios::binary);

#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#endif

    #define READ(name) fin.read((char*)&this->name, sizeof(this->name));
    READ(step) READ(precision) READ(mass) READ(damp) 
    READ(tempr) READ(time) READ(paras) READ(emphasis)

    this->nodes.deserialize(fin);
    this->velocities.deserialize(fin);
    this->adjacents.deserialize(fin);
    this->bonds.deserialize(fin);

    fin.close();
}

void molecule::store(std::string fname) {
    std::ofstream fout(fname, std::ios::out | std::ios::binary);

    #define WRITE(name) fout.write((char*)&this->name, sizeof(this->name));
    WRITE(step) WRITE(precision) WRITE(mass) WRITE(damp) 
    WRITE(tempr) WRITE(time) WRITE(paras) WRITE(emphasis)
    
    this->nodes.serialize(fout);
    this->velocities.serialize(fout);
    this->adjacents.serialize(fout);
    this->bonds.serialize(fout);

    fout.close();
}