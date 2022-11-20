#include <cmath>

#include "model.h"
#include "nmath.h"

// 原子数和键数保证增原子攀移时不会越界，同时防止减原子攀移初始时越界
model::model(para p): molecule(p.step, p.mass, p.damp, p.tempr, 
        p.m * p.n * p.repeat + ((p.climb < 0) ? 0 : p.climb), 3) {

    // 能量参数传递给 molecule，参数储存
    this->paras[0] = p.rest_len; this->paras[1] = p.k; this->paras[2] = p.tau;
    this->ppara = p;
    
    // 能量函数绑定，迫不得已使用 wrap 的参数，因为 lambda 表达式不能捕获局部变量
    this->bond_energy = my_bond_energy;
    this->node_energy = my_node_energy;
    
    // 模型构建过程
    perfect_model_position();
    this->adjacents.set_size(this->nodes.size()); // 邻接增加
    perfect_model_topology();
    add_dislocation();

    // 随机初始速度（平均速度与温度有关）
    double v = std::sqrt(3*K_B*this->ppara.tempr*this->ppara.mass);
    for (int i=0; i<this->nodes.size(); i++)
        this->velocities.push_back(v*nano::rand_vector());
}

// 生成点
void model::perfect_model_position() {
    
    // 转化为浮点数避免整数舍入
    double m = (double)this->ppara.m, n = (double)this->ppara.n;
    
    // 角度 A、B，分母为零必须单独处理
    double A = (2*n == m) ? PI / 2 : std::atan(sqrt(3)/2 * m / (n - m/2));
    double B = (2*m == n) ? PI / 2 : std::atan(sqrt(3)/2 * n / (m - n/2));
    
    // 半径 r（归一化）
    double r = 1/PI/2 * std::sqrt(m*m+n*n-m*n);

    for (int j = 0; j < this->ppara.repeat*this->ppara.m; j++)
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

void model::perfect_model_topology() {
    // 生成键和邻接结构
    #define NEAR_THEN_PUSH_ALL(a) \
        if (std::abs(this->nodes[flat(i, j)][2] - this->nodes[a][2]) < 2*this->ppara.rest_len) { \
            this->bonds.push_back({flat(i, j), a});   \
            this->adjacents[flat(i, j)].push_back(a); \
        }
    #define NEAR_THEN_PUSH(a) \
        if (std::abs(this->nodes[flat(i, j)][2] - this->nodes[a][2]) < 2*this->ppara.rest_len) \
            this->adjacents[flat(i, j)].push_back(a);
    
    // 键只三个方向生成，因为剩余三个方向已经生成了；邻接留个方向均生成
    for (int i = 0; i < this->ppara.n; i++) 
    for (int j = 0; j < this->ppara.repeat*this->ppara.m; j++) {
        NEAR_THEN_PUSH_ALL(flat(i+1, j))
        NEAR_THEN_PUSH_ALL(flat(i, j+1))
        NEAR_THEN_PUSH_ALL(flat(i-1, j+1))
        NEAR_THEN_PUSH(flat(i-1, j))
        NEAR_THEN_PUSH(flat(i, j-1))
        NEAR_THEN_PUSH(flat(i+1, j-1))
    }
}

void model::add_dislocation() {
    // 初始原子 j 位置 begin number，pos2d
    nano::pair<int> bn = {0, (this->ppara.repeat + 1) * this->ppara.m / 2};
    this->begin[1] = flat(bn);

    nano::pair<int> left, right, center;
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
    nano::pair<int> go = this->ppara.glide > 0 ? left: right;
    nano::pair<int> other = this->ppara.glide < 0 ? left: right;
    this->begin[0] = flat(bn + other);
    
    // 滑移
    if (this->ppara.glide != 0) {
        nano::pair<int> from, to; // 键
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
                if (this->bonds[i].replace(this->nodes.size(), flat(bn))) 
                    break;
            }
            bn = bn + go;
        }

    } else if (this->ppara.climb > 0) {
        // 增原子攀移
        int theright = flat(bn + other);
        int theleft = flat(bn + other - go);
        for (int i=0; i<this->ppara.climb; i++) {
            // 更新
            this->nodes.push_back((this->nodes[flat(bn)] + this->nodes[flat(bn-go)])/2);
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
        this->end[0] = flat(bn);
        this->end[1] = theright;
        return;
    }
    
    this->end[0] = flat(bn);
    this->end[1] = flat(bn + other);
}
