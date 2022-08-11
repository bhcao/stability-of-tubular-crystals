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

// 随机函数随距离成高斯分布
void molecule::disorganize() {
    std::srand(std::time(0));
    for (int i=0; i<this->nodes.size(); i++) {
        // 随机偏差
        double r = gaussrand();
        double phi = 2*PI * (double)rand() / RAND_MAX;
        double theta = PI/2 * ((double)rand() / RAND_MAX - 0.5);
        node rand_deviation = {std::cos(phi)*std::sin(theta), std::sin(phi)*std::sin(theta),
            std::cos(theta)};
        this->nodes[i] = this->nodes[i] + this->range*r*rand_deviation;
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
        if (i == this->pdis_pair.begin[0] || i == this->pdis_pair.begin[1] ||
            i == this->pdis_pair.end[0] || i == this->pdis_pair.end[1]) {
            fout << i << "\t2\t" << this->nodes[i].i << ' ' << this->nodes[i].j
                << ' ' << this->nodes[i].k << '\n';
        } else {
            fout << i << "\t1\t" << this->nodes[i].i << ' ' << this->nodes[i].j
                << ' ' << this->nodes[i].k << '\n';
        }
        
    }
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
    fout << this->data[this->data.size()];
    fout << "])\nplt.figure()\n"
            "plt.plot(x, y)\nplt.savefig(\""
         << name << "\")";
    fout.close();
    std::system("python __nano__temp__.py");
    std::remove("__nano__temp__.py");
}