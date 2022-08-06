#ifndef _FIGURE_H_
#define _FIGURE_H_

#include <vector>
using std::vector;

class figure {
    // 重载输入运算符
    inline friend figure& operator<<(figure &p, double n) {
        p.data.push_back(n);
        return p;
    }
    
public:
    figure() = default;
    
    // 绘制函数
    void draw(const char* name);

private:
    vector<double> data;
};

#endif // _FIGURE_H_
