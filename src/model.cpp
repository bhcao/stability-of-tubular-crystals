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
    
    // 初始化位错，设置检查点等
    this->begin[1] = this->end[0] = flat(0, (p.repeat + 1) * p.m / 2);
    nano::sarray<int>& center = this->adjacents[this->end[0]];
    if (p.direction > 0)
        this->begin[0] = this->end[1] = center[std::abs(p.direction)];
    if (p.direction < 0)
        this->begin[0] = this->end[1] = 
            center[p.direction == -1 ? 5 : (std::abs(p.direction) - 2)];
    this->checkpoint = dual(center[std::abs(p.direction) - 1], this->end[0], this->end[1]);

    add_dislocation();
    
    // 如果发生了位错，加入强调的粒子中
    if (this->begin[0] != this->end[1]) {
        this->emphasis.push_back(this->begin[0]);
        this->emphasis.push_back(this->begin[1]);
        this->emphasis.push_back(this->end[0]);
        this->emphasis.push_back(this->end[1]);
    }

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
    
    // 键只三个方向生成，因为剩余三个方向已经生成了；邻接六个方向均生成
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

// 移除粒子及其相连键
void model::remove(int node) {
    this->nodes[node] = this->nodes[this->nodes.size()-1];
    this->nodes.pop_back();
    // 移除相邻键
    for (int i = 0; i<this->adjacents[node].size(); i++)
        this->remove(nano::pair<int>(node, this->adjacents[node][i]));
    this->adjacents[node] = this->adjacents[this->nodes.size()-1];
    this->adjacents.pop_back();
}

void model::add_dislocation() {    
    // 滑移
    for (int i=0; i < this->ppara.glide; i++) {
        int new_end0 = dual(this->checkpoint, this->end[0], this->end[1]);
        int new_end1 = dual(this->end[0], new_end0, this->end[1]);

        // 键插入位置
        int n1 = between(this->end[0], this->end[1], new_end0);
        int n2 = between(new_end1, this->end[1], new_end0);

        remove(nano::pair<int>(this->end[1], new_end0));
        insert(nano::pair<int>(this->end[0], new_end1), n1, n2);
        this->checkpoint = this->end[0];
        this->end[0] = new_end0;
        this->end[1] = new_end1;
    }
    
    /* // 攀移，尽量增大 repeat（大于等于 3），避免与最后几个原子重合造成隐患
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
    this->end[1] = flat(bn + other);*/
}
