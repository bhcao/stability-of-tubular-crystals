#include <ctime>
#include <cmath>

#include "nmath.h"

namespace nano {

rand_pool::rand_pool(): pool(RAND_POOL_SIZE) {
    this->pool.set_size(RAND_POOL_SIZE);
    int time = std::time(0);
    // 初始化随机种子
    for (int i=0; i<this->pool.size(); i++)
        this->pool[i].seed = time + 3*i;
}

// 不够时会频繁调用第一个池，一定要设置足够池
nano::vector rand_pool::gen_vector() {
// 不使用 KOKKOS 时不需要上锁，但已有的锁结构保留（懒得改了）
#ifdef USE_KOKKOS
    int num_of_block = -1;
    // 一个个尝试获取锁
    for (int i=0; i<this->pool.size(); i++)
    if (this->pool[i].lock.try_lock()) {
        num_of_block = i;
        break;
    }
    // 没有获取到就在默认位置上等待
    if (num_of_block == -1) {
        num_of_block = 0;
        this->pool[num_of_block].lock.lock();
    }
#else
    int num_of_block = 1;
#endif
    // 球坐标参数
    double r = gauss_rand(num_of_block);
    double phi = 2*PI * rand(num_of_block);
    double theta = PI/2 * (rand(num_of_block) - 0.5);

#ifdef USE_KOKKOS
    this->pool[num_of_block].lock.unlock();
#endif

    return {r*std::cos(phi)*std::sin(theta),
        r*std::sin(phi)*std::sin(theta), r*std::cos(theta)};
}

double rand_pool::rand(int i) {
    // 参数选取参考 MMIX by Donald Knuth
    this->pool[i].seed = (6364136223846793005*this->pool[i].seed +
       1442695040888963407) % UINT64_MAX;
    return (double)this->pool[i].seed / UINT64_MAX;
}

double rand_pool::gauss_rand(int i) {
    rand_seed &p = this->pool[i];
    p.phase = !p.phase;
    if (p.phase) {
        do {
            p.V1 = 2 * rand(i) - 1;
            p.V2 = 2 * rand(i) - 1;
            p.S = p.V1*p.V1 + p.V2*p.V2;
        } while (p.S >= 1 || p.S == 0);
        return p.V1 * std::sqrt(-2 * std::log(p.S) / 5);
    }
    return p.V2 * std::sqrt(-2 * std::log(p.S) / 5);
}

}