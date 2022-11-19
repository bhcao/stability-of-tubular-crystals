#include <cstdlib>
#include <ctime>
#include <cmath>

#include "nmath.h"

namespace nano {

// 输出高斯分布的数
double gauss_rand() {
    static double V1, V2, S;
    static bool phase = false;
    phase = !phase;
    if (phase) {
        do {
            V1 = 2 * (double)rand() / RAND_MAX - 1;
            V2 = 2 * (double)rand() / RAND_MAX - 1;
            S = V1*V1 + V2*V2;
        } while (S >= 1 || S == 0);
        return V1 * std::sqrt(-2 * std::log(S) / 5);
    }
    return V2 * std::sqrt(-2 * std::log(S) / 5);
}

vector rand_vector() {
    std::srand(std::time(0));
    double r = gauss_rand();
    double phi = 2*PI * (double)rand() / RAND_MAX;
    double theta = PI/2 * ((double)rand() / RAND_MAX - 0.5);
    return {r*std::cos(phi)*std::sin(theta),
        r*std::sin(phi)*std::sin(theta), r*std::cos(theta)};
}

}