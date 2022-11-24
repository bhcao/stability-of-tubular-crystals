#include <fstream>
#include <iostream>
#include <cmath>
#include <string>

#include "nmath.h"
#include "narray.h"
#include "molecule.h"

double molecule::local_energy(nano::vector center, nano::sarray<int> others) {
    // 相邻粒子转化成坐标
    nano::sarray<nano::vector> adjacent;
    for (int i=0; i<others.size(); i++)
        adjacent.push_back(this->nodes[others[i]]);

    double count = 0;
    for (int i=0; i<others.size(); i++) {
        count += this->bond_energy(center, adjacent[i], this->paras)/2;
    }
    count += this->node_energy(center, adjacent, this->paras);
    return count;
}

// 借助 local_energy 键要算两次，效率低，但可以省去键结构，逻辑也更清晰
double molecule::total_energy() {
    double count = 0;
    // 以粒子为中心的能量
    for (int i=0; i<this->nodes.size(); i++) {
        count += local_energy(this->nodes[i], this->adjacents[i]);
    }
    return count;
}

double molecule::local_energy_for_update(nano::vector center, nano::sarray<int> others) {
    // 相邻粒子转化成坐标
    nano::sarray<nano::vector> adjacent;
    for (int i=0; i<others.size(); i++)
        adjacent.push_back(this->nodes[others[i]]);

    double count = 0;
    for (int i=0; i<others.size(); i++) {
        count += this->bond_energy(center, adjacent[i], this->paras);
    }
    count += this->node_energy(center, adjacent, this->paras);
    return count;
}

void molecule::update() {
    for (int i=0; i<this->nodes.size(); i++) {
        this->velocities[i] += accelerate(this->nodes[i], this->velocities[i], this->adjacents[i]) * this->step;
        this->nodes[i] += this->velocities[i] * this->step;
        this->time++;
    }
}

// 输出
void molecule::dump(std::string fname, nano::dump_t dump_type) {
    
    // 打开文件以写入，文件名为 fname.dump（追加）或 fname.1.data
    std::ios::openmode file_mode = (!DUMP_CHECK(nano::DATA_FILE, dump_type) && (this->time != 0)) ?
        std::ios::app : (std::ios::out | std::ios::trunc);
    fname += !DUMP_CHECK(nano::DATA_FILE, dump_type) ? ".dump" : 
        ("." + std::to_string(this->time) + ".data");
    std::ofstream fout;
    fout.open(fname, file_mode);

    // 设置边界
    double boundary[6] = {0, 1, 0, 1, 0, 1};
    
    // 文件头，data 文件的
    if (DUMP_CHECK(nano::DATA_FILE, dump_type)) {
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
        fout << "\nAtoms\n\n";

    } else { // dump 文件的
        fout << "ITEM: TIMESTEP\n" << this->time << "\nITEM: NUMBER OF ATOMS\n"
            << this->nodes.size() << "\nITEM: BOX BOUNDS ss ss ss\n" // 边界
            << boundary[0] << " " << boundary[1] << "\n"
            << boundary[2] << " " << boundary[3] << "\n"
            << boundary[4] << " " << boundary[5] << "\n"
            << "ITEM: ATOMS id type xs ys zs";

        // 检查是否要输出这些内容
        if (DUMP_CHECK(nano::VELOCITY, dump_type)) fout << " vx vy vz";
        if (DUMP_CHECK(nano::DIV_FORCE, dump_type)) fout << " dfx dfy dfz";
        if (DUMP_CHECK(nano::LAN_FORCE, dump_type)) fout << " fx fy fz";
        if (DUMP_CHECK(nano::K_ENERGY, dump_type)) fout << " ke";
        if (DUMP_CHECK(nano::P_ENERGY, dump_type)) fout << " pe";
        fout << "\n";
    }
    
    // 强调变色（类型变为 2）
    nano::darray<int> types(this->nodes.size());
    for (int i=0; i<this->nodes.size(); i++) 
        types.push_back(1);
    if (DUMP_CHECK(nano::EMPHASIS, dump_type))
    for (int i=0; i<this->emphasis.size(); i++)
        types[this->emphasis[i]] = 2;
    
    // 正文数据
    for (int i=0; i<this->nodes.size(); i++) {
        fout << i << "\t " << types[i] << "\t" // 基础 id type xs ys zs
            << this->nodes[i][0] << '\t' << this->nodes[i][1] << '\t' << this->nodes[i][2];
        if (!DUMP_CHECK(nano::DATA_FILE, dump_type)) {
            if (DUMP_CHECK(nano::VELOCITY, dump_type))
                fout << "\t" << this->velocities[i][0] << "\t" << this->velocities[i][1]
                    << "\t" << this->velocities[i][2];
            if (DUMP_CHECK(nano::DIV_FORCE, dump_type)) {
                nano::vector temp = div(this->nodes[i], this->adjacents[i]);
                fout << "\t" << temp[0] << "\t" << temp[1] << "\t" << temp[2];
            }
            if (DUMP_CHECK(nano::LAN_FORCE, dump_type)) {
                nano::vector temp = this->mass * accelerate(this->nodes[i], this->velocities[i],
                    this->adjacents[i]);
                fout << "\t" << temp[0] << "\t" << temp[1] << "\t" << temp[2];
            }
            if (DUMP_CHECK(nano::K_ENERGY, dump_type))
                fout << "\t" << this->mass/2*nano::mod(this->velocities[i]);
            if (DUMP_CHECK(nano::P_ENERGY, dump_type))
                fout << "\t" << local_energy(this->nodes[i], this->adjacents[i]);
        }
        fout << '\n';
    }

    // 原子间的键，data 文件才输出
    if (DUMP_CHECK(nano::DATA_FILE, dump_type)) {
        fout << "\nBonds\n\n";
        for (int i=0; i<this->bonds.size(); i++) {
            fout << i << "\t1\t" << this->bonds[i][0] << '\t' << this->bonds[i][1] << '\n';
        }
    }

    fout.close();
    return;
} 

molecule::molecule(std::string fname, double(*in_node_energy)(nano::vector, nano::sarray<nano::vector>,
    nano::sarray<double>), double(*in_bond_energy)(nano::vector, nano::vector, nano::sarray<double>)):
    node_energy(in_node_energy), bond_energy(in_bond_energy) {
    std::ifstream fin(fname, std::ios::in | std::ios::binary);

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