#include <fstream>
#include <cstdlib>
#include <ctime>

#include "molecule.h"

using std::rand;

// 打乱
void molecule::disorganize() {
    std::srand(std::time(0));
    for (int i=0; i<this->nodes.size(); i++) {
        // 随机偏差
        node rand_deviation = {double(rand()), double(rand()), double(rand())};
        this->nodes[i] = this->nodes[i] + this->range/RAND_MAX*rand_deviation;
    }
}

void molecule::update() {
    for (int a=0; a<this->nodes.size(); a++) {
        node temp = this->nodes[a];
        temp.i += this->precision;
        // 偏导数乘以步长
        double i = this->step * (energy(temp) - energy(this->nodes[a])) / this->precision;
        temp = this->nodes[a];
        temp.j += precision;
        double j = this->step * (energy(temp) - energy(this->nodes[a])) / this->precision;
        temp = this->nodes[a];
        temp.k += precision;
        double k = this->step * (energy(temp) - energy(this->nodes[a])) / this->precision;
        node div = {i, j, k};
        this->nodes[a] = this->nodes[a] - div;
    }
}

// 输出
void molecule::dump(std::ofstream &fout) {
    fout << "ITEM: TIMESTEP\n" << this->time << "\nITEM: NUMBER OF ATOMS\n"
        << this->nodes.size() << "\nITEM: BOX BOUNDS ss ss ss\n";
    
    double boundary[6] = {0, 1, 0, 1, 0, 1};
    for (int i=0; i<6; i++) {
        if (i%2==0) {
            fout << boundary[i] << ' ';
        } else {
            fout << boundary[i] << '\n';
        }
    }
    
    fout << "ITEM: ATOMS id type xs ys zs\n";
    for (int i=0; i<this->nodes.size(); i++) {
        fout << i << "\t1\t" << this->nodes[i].i << ' ' << this->nodes[i].j
            << ' ' << this->nodes[i].k << '\n';
    }
}