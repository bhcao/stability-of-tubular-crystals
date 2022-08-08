#ifndef _VECTOR_H_
#define _VECTOR_H_

#ifdef USE_CUDA
#include <cuda_runtime.h>
#else // not define USE_CUDA
#define __global__
#define __device__
#endif // USE_CUDA

namespace nano {
    // 存储小向量的类，最大为 7，运行于 GPU
    template <typename T> class s_vector {
    public:
        // 非常倒霉，构造函数我无法使其为 __global__
        __device__ inline void push_back(T in) {
            this->data[len++] = in;
        }
        __device__ inline T& operator[](int i) {
            return data[i];
        }
        __device__ inline int size() {
            return this->len;
        }

        T data[7];
        int len;
    };

    // 储存大向量的类，运行于 CPU，为了效率，没有设置 private
    template <typename T> class vector {
    public:

        // 初始化定义了最大容量，而不是这么大的数组，超出最大容量不会管理！
        inline vector(int size): len(0) {
            this->data = new T[size];
        }
        
        inline ~vector() {
            free(this->data);
        }

        inline void push_back(T in) {
            this->data[len++] = in;
        }

        inline T& operator[](int i) {
            return data[i];
        }

        inline int size() {
            return this->len;
        }
        
        // 拷贝至显存
        #ifdef USE_CUDA
        inline T* gpu_copy() {
            T* temp;
            cudaMalloc((void**)&temp, this->size()*sizeof(T));
            cudaMemcpy(temp, this->data, this->size()*sizeof(T), cudaMemcpyHostToDevice);
            return temp;
        }

        inline void gpu_copy_back(T* temp) {
            cudaMemcpy(this->data, temp, this->size()*sizeof(T), cudaMemcpyDeviceToHost);
            cudaFree(temp);
        }
        #endif

    private:
        T* data;
        int len;
    };
}

#endif // _VECTOR_H_