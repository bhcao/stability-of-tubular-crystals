#include <cmath>
#include <vector>
#include <fstream>

#include "model.h"

using std::vector;

// 输出
void model::dump(std::ofstream &fout) {
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

// 找到所有原子与之相邻的点
void model::generate_adjacent() {
    for (int center = 0; center < this->nodes.size(); center++) {
        adjacent padjacent;
        for (int i = 0; i < this->bonds.size(); i++) {
            // 如果键的一端是 center
            if (this->bonds[i].a == center) {
                node other = this->nodes[this->bonds[i].b];
                // 将坐标添加到末尾
                padjacent.others.push_back(other);
                // 同时记录位置
                padjacent.others_id.push_back(this->bonds[i].b);
            } else if (this->bonds[i].b == center) {
                node other = this->nodes[this->bonds[i].a];
                padjacent.others.push_back(other);
                padjacent.others_id.push_back(this->bonds[i].a);
            }
        }
        this->adjacents.push_back(padjacent);
    }
}

// 生成点
void model::generate_nodes() {
    using std::atan, std::sin, std::cos, std::sqrt;
    
    // 转化为浮点数避免整数舍入
    double m = (double)this->ppara.m;
    double n = (double)this->ppara.n;
    
    // 角度 A、B
    double A, B;
    // C 语言分母为零必须单独处理
    if (2*n == m) {
        A = PI / 2;
    } else {
        A = atan(sqrt(3)/2 * m / (n - m/2));
    }

    if (2*m == n) {
        B = PI / 2;
    } else {
        B = atan(sqrt(3)/2 * n / (m - n/2));
    }
    
    // 半径 r（归一化）
    double r = 1/PI/2 * sqrt(m*m+n*n-m*n);

    for (int j = 0; j < this->ppara.repeat*this->ppara.m; j++) {
        for (int i = 0; i < this->ppara.n; i++) {
            double x, y, z;
            x = -(i-n)*cos(A) + j*cos(B);
            z = (i-n)*sin(A) + j*sin(B);

            // 如果在数轴下方
            if (z < 0) {
                x += this->ppara.repeat*m*cos(B);
                z += this->ppara.repeat*m*sin(B);
            }
                       
            // 实际坐标
            y = r * (sin(x/r) + 1);
            x = r * (cos(x/r) + 1);
            
            // 增加到队列
            this->nodes.push_back({this->ppara.rest_len*x,
                this->ppara.rest_len*y, this->ppara.rest_len*z});
        }
    }
}

// 生成键
void model::generate_bonds() {
    using std::abs;
    
    for (int j = 0; j < this->ppara.repeat*this->ppara.m; j++) {
        for (int i = 0; i < this->ppara.n; i++) {
            if (abs(this->nodes[flat(i, j)].k - this->nodes[flat(i+1, j)].k) <
                2*this->ppara.rest_len) {
                // 增加到队列
                this->bonds.push_back({flat(i, j), flat(i+1, j)});
            }
            if (abs(this->nodes[flat(i, j)].k - this->nodes[flat(i, j+1)].k) <
                2*this->ppara.rest_len) {
                this->bonds.push_back({flat(i, j), flat(i, j+1)});
            }
            if (abs(this->nodes[flat(i, j)].k - this->nodes[flat(i-1, j+1)].k) <
                2*this->ppara.rest_len) {
                this->bonds.push_back({flat(i, j), flat(i-1, j+1)});
            }
        }
    }
}

void model::replace_bond(bond from, bond to) {
    for (int i=0; i<this->bonds.size(); i++) {
        if (this->bonds[i] == from) {
            this->bonds[i] = to;
            return;
        }
    }
}

// TODO
void model::glide_bond() {
    // 初始原子 j 位置 begin number
    double bn = (this->ppara.repeat + 1) * double(this->ppara.m) / 2;
    switch (this->ppara.direction) {
    case 1:
        if (this->ppara.glide < 0) {
            for (int i=0; i < -this->ppara.glide; i++) {
                //replace_bond({flat(0, bn+i+1), flat(1, bn+i)}, {flat(0, bn+i), flat(1, bn+i+1)});
            }
        } else if (this->ppara.glide > 0) {
            for (int i=0; i < this->ppara.glide; i++) {
                replace_bond({flat(0, bn+i+1), flat(1, bn+i)}, {flat(0, bn+i), flat(1, bn+i+1)});
            }
        }
        break;
    case 2:
        if (this->ppara.glide < 0) {
            for (int i=0; i < -this->ppara.glide; i++) {
                //replace_bond({flat(0, bn+i+1), flat(1, bn+i)}, {flat(0, bn+i), flat(1, bn+i+1)});
            }
        } else if (this->ppara.glide > 0) {
            for (int i=0; i < this->ppara.glide; i++) {
                //replace_bond({flat(0, bn+i+1), flat(1, bn+i)}, {flat(0, bn+i), flat(1, bn+i+1)});
            }
        }
        break;
    case 3:
        if (this->ppara.glide < 0) {
            for (int i=0; i < -this->ppara.glide; i++) {
                //replace_bond({flat(0, bn+i+1), flat(1, bn+i)}, {flat(0, bn+i), flat(1, bn+i+1)});
            }
        } else if (this->ppara.glide > 0) {
            for (int i=0; i < this->ppara.glide; i++) {
                //replace_bond({flat(0, bn+i+1), flat(1, bn+i)}, {flat(0, bn+i), flat(1, bn+i+1)});
            }
        }
        break;
    }
}

// TODO
void model::climb_bond() {
    /*if (climb > 0) {

    } else if (climb < 0) {

    }*/
}


