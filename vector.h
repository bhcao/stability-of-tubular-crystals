#ifndef NANO_VECTOR_H_
#define NANO_VECTOR_H_

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

    // 储存大向量的类，内存初始化时确定，显存第一次同步时确定，析构时同时释放
    template <typename T> class vector {
    public:

        // 初始化定义了最大容量，而不是这么大的数组，超出最大容量不会管理！
        inline vector(int size): len(0) {
            this->data = new T[size];
            #ifdef USE_CUDA
            this->gpu_data = NULL;
            #endif
        }

        inline ~vector() {
            free(this->data);
            #ifdef USE_CUDA
            if (this->gpu_data != NULL) {
                cudaFree(this->gpu_data);
            }
            #endif
        }

        // std::vector 兼容函数
        inline void push_back(T in) { this->data[len++] = in;}
        inline T pop_back() { return this->data[--len];}
        inline T& operator[](int i) { return data[i]; }
        inline int size() { return this->len; }
        
        
        #ifdef USE_CUDA
        // 拷贝至显存
        inline void gpu_synchro() {
            if (this->gpu_data == NULL) {
                cudaMalloc((void**)&this->gpu_data, this->len*sizeof(T));
            }
            cudaMemcpy(this->gpu_data, this->data, this->len*sizeof(T),
                cudaMemcpyHostToDevice);
        }

        // 拷贝回内存
        inline void cpu_synchro() {
            cudaMemcpy(this->data, this->gpu_data, this->len*sizeof(T),
                cudaMemcpyDeviceToHost);
        }

        // 返回显存指针
        inline T* get_gpu_data() { return this->gpu_data; }
        #endif

    private:
        // 显存数据
        #ifdef USE_CUDA
        T* gpu_data;
        #endif

        T* data;
        int len;
    };
}

#endif // NANO_VECTOR_H_