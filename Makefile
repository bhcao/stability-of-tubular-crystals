TARGET = $(shell ls | grep .cpp$ | grep -v main.cpp)

gpu:
	nvcc -D USE_CUDA main.cpp $(TARGET) energy.cu -o main.out --expt-relaxed-constexpr

build:
	mv energy.cu energy.cpp
	g++ main.cpp $(TARGET) energy.cpp -o main.out
	mv energy.cpp energy.cu
