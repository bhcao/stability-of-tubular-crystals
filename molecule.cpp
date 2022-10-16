#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cstdio>

#include "molecule.h"

using std::rand;

double gaussrand() {
    static double V1, V2, S;
    static int phase = 0;
    double X;
    if (phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
            V1 = 2*U1 - 1;
            V2 = 2*U2 - 1;
            S = V1*V1 + V2*V2;
        } while (S >= 1 || S == 0);
        X = V1 * std::sqrt(-2 * std::log(S) / 5);
    } else {
        X = V2 * std::sqrt(-2 * std::log(S) / 5);
    }
    phase = 1 - phase;
    return X;
}

node randnode() {
    std::srand(std::time(0));
    double r = gaussrand();
    double phi = 2*PI * (double)rand() / RAND_MAX;
    double theta = PI/2 * ((double)rand() / RAND_MAX - 0.5);
    node temp = {std::cos(phi)*std::sin(theta), std::sin(phi)*std::sin(theta),
            std::cos(theta)};
    return r*temp;
}

// 随机函数随距离成高斯分布
void molecule::disorganize() {
    
    for (int i=0; i<this->nodes.size(); i++) {
        // 随机偏差
        this->nodes[i] = this->nodes[i] + this->range*randnode();
    }
}

// 输出
void molecule::dump(std::ofstream &fout, enum dump_type dtype) {
    
    double boundary[6] = {0, 1, 0, 1, 0, 1};
    
    // 文件头
    if (dtype == DUMP_FILE) {
        fout << "ITEM: TIMESTEP\n" << this->time << "\nITEM: NUMBER OF ATOMS\n"
        << this->nodes.size() << "\nITEM: BOX BOUNDS ss ss ss\n" // 边界
        << boundary[0] << " " << boundary[1] << "\n"
        << boundary[2] << " " << boundary[3] << "\n"
        << boundary[4] << " " << boundary[5] << "\n"
        << "ITEM: ATOMS id type xs ys zs\n";
    } else if (dtype == DATA_FILE) {
        fout << "# Model for nanotube. AUTO generated, DO NOT EDIT\n\n"
        << this->nodes.size() << "\tatoms\n" // 内容
        << this->bonds.size() << "\tbonds\n\n"
        << "2\tatom types\n1\tbond types\n\n"
        << boundary[0] << "\t" << boundary[1] << "\txlo xhi\n" // 边界
        << boundary[2] << "\t" << boundary[3] << "\tylo yhi\n"
        << boundary[4] << "\t" << boundary[5] << "\tzlo zhi\n\n"
        << "Masses\n\n" << "1\t0.01\n2\t0.01\n"
        << "\nAtoms\n\n";
    }
    
    for (int i=0; i<this->nodes.size(); i++) {
        if (i == this->pdis_pair.begin[0] || i == this->pdis_pair.begin[1] ||
            i == this->pdis_pair.end[0] || i == this->pdis_pair.end[1]) {
            fout << i << "\t2\t" << this->nodes[i].i << '\t' << this->nodes[i].j
                << '\t' << this->nodes[i].k << '\n';
        } else {
            fout << i << "\t1\t" << this->nodes[i].i << '\t' << this->nodes[i].j
                << '\t' << this->nodes[i].k << '\n';
        }
    }

    // 原子间的键
    if (dtype == DATA_FILE) {
        fout << "\nBonds\n\n";
        for (int i=0; i<this->bonds.size(); i++) {
            fout << i << "\t1\t" << this->bonds[i].a << '\t' << this->bonds[i].b << '\n';
        }
    }
    return;
}

void figure::draw(const char* name){
    std::ofstream fout("__nano__temp__.py");
    fout << "from matplotlib import pyplot as plt\n"
             "import numpy as np\nx = np.linspace(1, " 
         << this->data.size() << ", "
         << this->data.size() << ")\ny = np.array([";
    for (int i=0; i<this->data.size()-1; i++) {
        fout << this->data[i] << ", ";
    }
    fout << this->data[this->data.size()-1];
    fout << "])\nplt.figure()\n"
            "plt.plot(x, y)\nplt.savefig(\""
         << name << "\")";
    fout.close();
    std::system("python __nano__temp__.py");
    std::remove("__nano__temp__.py");
}
