/* -*- c++ -*- ----------------------------------------------------------
   AVAILABLE https://gitee.com/Bovera/nanotube, COPYRIGHT Bovera (2022)

   定义了数组，包括
   1. 数组 sarray（栈）、darray（堆），SARRAY_SIZE 可以改大小；
   2. 序列化和反序列化函数。
------------------------------------------------------------------------- */

#ifndef NANO_NARRAY_H_
#define NANO_NARRAY_H_

#include <iostream>
#include <fstream>

namespace nano {

#ifndef SARRAY_SIZE
#define SARRAY_SIZE 8
#endif

// 静态数组，编译时分配大小，通过 SARRAY_SIZE 改变
template <typename T> class sarray {
public:
    inline sarray(): len(0) {}
    // size 数组长度，cap 容量大小
    inline int size() { return this->len; }
    inline T& operator[](int i) { return this->data[i]; }
    
    // 增加减少元素
    inline void push_back(T in) { this->data[this->len++] = in;}
    inline T pop_back() { return this->data[--this->len];}

    // 增加长度，仅仅为了方便 adjacent 使用
    inline void set_size(int i) { this->len = i; }

    int find(T from) { // 返回 from 所在位置
        for (int i = 0; i < this->len; i++)
            if ((*this)[i] == from)
                return i;
        return -1;
    }
    
    // 交换两个，from -> to
    bool replace(T from, T to) {
        for (int i=0; i<this->len; i++)
        if (this->data[i] == from) {
            this->data[i] = to;
            return true;
        }
        return false;
    }
    
    // 在点 n 插入，顺序保留
    void insert(T from, int n) {
        for (int i=this->len; i>n; i--)
            this->data[i] = this->data[i-1];
        this->data[n] = from;
        this->len++;
    }

    // 移除 from，顺序保留
    bool remove(T from) {
        for (int i=0; i<this->len; i++)
        if (this->data[i] == from) {
            for (int j=i; j<this->len-1; j++)
                this->data[j] = this->data[j+1];
            this->len--;
            return true;
        }
        return false;
    }
    
    // 比较（有纯虚函数的类不能做参数）
    bool operator==(sarray<T> p) {
        if (this->size() != p.size())
            return false;
        for (int i=0; i<this->size(); i++)
            if (this->data[i] != p.data[i])
                return false;
        return true;
    }
    
    // 重载赋值函数，避免指针被拷贝了
    inline void operator=(sarray<T> p) {
        for (int i=0; i<p.len; i++)
            this->data[i] = p.data[i];
        this->len = p.len;
    }

private:
    int len;
    T data[SARRAY_SIZE];
};

// 动态数组，运行时分配大小，为了效率不会检查越界，注意不会清零
template <typename T> class darray{
public:
    // 除了 cap，其余与 std::vector 兼容
    inline darray(int size): capacity(size), len(0) {
        this->data = new T[size];
    }

    // size 数组长度，cap 容量大小
    inline int size() { return this->len; }
    
    // 增加减少元素
    inline void push_back(T in) { this->data[this->len++] = in;}
    inline T pop_back() { return this->data[--this->len];}

    // 增加长度，仅仅为了方便 adjacent 使用
    inline void set_size(int i) { this->len = i; }

    int find(T from) { // 返回 from 所在位置
        for (int i = 0; i < this->len; i++)
            if (this->data[i] == from)
                return i;
        return -1;
    }
    
    // 交换两个，from -> to
    bool replace(T from, T to) {
        for (int i=0; i<this->len; i++)
        if (this->data[i] == from) {
            this->data[i] = to;
            return true;
        }
        return false;
    }

    inline T& operator[](int i) { return this->data[i]; }

    inline ~darray() { delete []this->data; }
    inline int cap() { return this->capacity; }

    // 移除 from，成功返回 true，否则 false（顺序不保存）
    bool remove(T from) {
        for (int i=0; i<this->len; i++)
        if (this->data[i] == from) {
            this->data[i] = this->data[--this->len];
            return true;
        }
        return false;
    }

    // 比较（有纯虚函数的类不能做参数）
    bool operator==(darray<T> p) {
        if (this->size() != p.size())
            return false;
        for (int i=0; i<this->size(); i++)
            if (this->data[i] != p.data[i])
                return false;
        return true;
    }

    // 序列化与反序列化（sarray 不需要因为没有指针）
    void serialize(std::ofstream &f) {
        f.write((char*)&this->len, sizeof(int));
        f.write((char*)this->data, this->len*sizeof(T));
    }
    
    inline darray() = default; // 默认构造，之后必须调用 deserialize
    void deserialize(std::ifstream &f) {
        f.read((char*)&this->len, sizeof(int));
        this->capacity = this->len;
        this->data = new T[this->len];
        f.read((char*)this->data, this->len*sizeof(T));
    }

private:
    int len;
    T* data;
    int capacity; // 容量大小
};

}

#endif // NANO_NARRAY_H_