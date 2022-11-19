#include <cmath>

#include "model.h"
#include "nmath.h"

inline int relu(int i) {
    if (i<0) return 0;
    else return i;
}

// 保证增原子攀移时不会越界，同时防止减原子攀移初始时越界
#define NUM ppara_in.m*ppara_in.n*ppara_in.repeat + relu(ppara_in.climb)
model::model(para ppara_in): time(0), step(2e-5),
        precision(1e-6), range(1.0), nodes(NUM), bonds(3*NUM), speeds(NUM),
        ppara(ppara_in), adjacents(NUM), adjacents_id(NUM) {
    generate_nodes();
    generate_bonds();
    glide_climb();
    // 与显存结构共用，无法初始化，只能手动初始化
    double v = std::sqrt(3*K_B*this->ppara.tempr*this->ppara.mass);
    for (int i=0; i<this->nodes.size(); i++) {
        this->nodes[i].id = i;
        this->speeds.push_back(v*randnode());
        this->adjacents[i].len = 0;
        this->adjacents_id[i].len = 0;
    }
    generate_adjacent();
}

// 生成点
void model::generate_nodes() {
    
    // 转化为浮点数避免整数舍入
    double m = (double)this->ppara.m;
    double n = (double)this->ppara.n;
    
    // 角度 A、B
    double A, B;
    // C 语言分母为零必须单独处理
    if (2*n == m) {
        A = PI / 2;
    } else {
        A = std::atan(sqrt(3)/2 * m / (n - m/2));
    }

    if (2*m == n) {
        B = PI / 2;
    } else {
        B = std::atan(sqrt(3)/2 * n / (m - n/2));
    }
    
    // 半径 r（归一化）
    double r = 1/PI/2 * std::sqrt(m*m+n*n-m*n);

    for (int j = 0; j < this->ppara.repeat*this->ppara.m; j++) {
        for (int i = 0; i < this->ppara.n; i++) {
            double x, y, z;
            x = -(i-n)*std::cos(A) + j*std::cos(B);
            z = (i-n)*std::sin(A) + j*std::sin(B);

            // 如果在数轴下方
            if (z < 0) {
                x += this->ppara.repeat*m*std::cos(B);
                z += this->ppara.repeat*m*std::sin(B);
            }
                       
            // 实际坐标
            y = r * (std::sin(x/r) + 1);
            x = r * (std::cos(x/r) + 1);
            
            // 增加到队列
            this->nodes.push_back({this->ppara.rest_len*x,
                this->ppara.rest_len*y, this->ppara.rest_len*z});
        }
    }
}

// 生成键
void model::generate_bonds() {
    
    for (int j = 0; j < this->ppara.repeat*this->ppara.m; j++) {
        for (int i = 0; i < this->ppara.n; i++) {
            if (std::abs(this->nodes[flat(i, j)].k - this->nodes[flat(i+1, j)].k) <
                2*this->ppara.rest_len) {
                // 增加到队列
                this->bonds.push_back({flat(i, j), flat(i+1, j)});
            }
            if (std::abs(this->nodes[flat(i, j)].k - this->nodes[flat(i, j+1)].k) <
                2*this->ppara.rest_len) {
                this->bonds.push_back({flat(i, j), flat(i, j+1)});
            }
            if (std::abs(this->nodes[flat(i, j)].k - this->nodes[flat(i-1, j+1)].k) <
                2*this->ppara.rest_len) {
                this->bonds.push_back({flat(i, j), flat(i-1, j+1)});
            }
        }
    }
}

void model::glide_climb() {
    // 初始原子 j 位置 begin number
    pos2d bn = {0, (this->ppara.repeat + 1) * this->ppara.m / 2};
    this->pdis_pair.begin[1] = flat(bn);

    pos2d left, right, center;
    switch (this->ppara.direction) {
    case 0:
        left = {0,-1}; right = {1,0};
        break;
    case 1:
        left = {1,-1}; right = {0,1};
        break;
    case 2:
        left = {1,0}; right = {-1,1};
    }
    
    center = left + right;
    pos2d go = this->ppara.glide > 0 ? left: right;
    pos2d other = this->ppara.glide < 0 ? left: right;
    this->pdis_pair.begin[0] = flat(bn + other);
    
    // 滑移
    if (this->ppara.glide != 0) {
        bond from, to;
        for (int i=0; i < std::abs(this->ppara.glide); i++) {
            from = {flat(bn), flat(bn + center)};
            to = {flat(bn + left), flat(bn + right)};
            this->bonds.replace(from, to);
            bn = bn + go;
        }
    }
    
    // 攀移，尽量增大 repeat（大于等于 3），避免与最后几个原子重合造成隐患
    if (this->ppara.climb < 0) {
        // 减原子攀移
        for (int i=0; i<-this->ppara.climb; i++) {
            this->nodes[flat(bn + other)] = this->nodes.pop_back();
            bn = bn + other;
            this->bonds.remove({flat(bn), flat(bn + left)});
            this->bonds.remove({flat(bn), flat(bn - left)});
            this->bonds.remove({flat(bn), flat(bn + right)});
            this->bonds.remove({flat(bn), flat(bn - right)});
            this->bonds.replace({flat(bn), flat(bn + center)},
                {flat(bn + go), flat(bn - go)});
            this->bonds.push_back({flat(bn+go), flat(bn+other)});
            for (int i=0; i<this->bonds.size(); i++) {
                if (this->bonds[i].a == this->nodes.size()) {
                    this->bonds[i].a = flat(bn);
                }
                if (this->bonds[i].b == this->nodes.size()) {
                    this->bonds[i].b = flat(bn);
                }
            }
            bn = bn + go;
        }

    } else if (this->ppara.climb > 0) {
        // 增原子攀移
        int theright = flat(bn + other);
        // 不知为什么不能连写
        pos2d temp = bn + other;
        int theleft = flat(temp - go);
        for (int i=0; i<this->ppara.climb; i++) {
            // 更新
            this->nodes.push_back(average(this->nodes[flat(bn)],
                this->nodes[flat(bn-go)]));
            int thethis = this->nodes.size()-1;
            this->bonds.push_back({thethis, theright});
            this->bonds.push_back({thethis, theleft});
            this->bonds.push_back({thethis, flat(bn - center)});
            this->bonds.replace({flat(bn), flat(bn - go)},
                {thethis, flat(bn - go)});
            this->bonds.replace({flat(bn), theleft},
                {thethis, flat(bn)});
            bn = bn - center;
            // 下一次更新
            theright = thethis;
            theleft = flat(bn + other);
        }
        this->pdis_pair.end[0] = flat(bn);
        this->pdis_pair.end[1] = theright;
        return;
    }
    
    this->pdis_pair.end[0] = flat(bn);
    this->pdis_pair.end[1] = flat(bn + other);
}
