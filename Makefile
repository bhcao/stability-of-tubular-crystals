MAIN = main.cpp
SOURCE = figure.cpp model.cpp molecule.cpp
CPPENERGY = energy.cpp
CUDAENERGY = energy.cu
LIB = nano.a
LDFLAGS = -L/usr/lib/cuda/lib64 -lcudart

gpu:
	nvcc -D USE_CUDA $(MAIN) $(SOURCE) $(CUDAENERGY) -o main.out --expt-relaxed-constexpr

cpu:
	mv energy.cu energy.cpp
	g++ $(MAIN) $(SOURCE) $(CPPENERGY) -o main.out
	mv energy.cpp energy.cu

cpu-mpi:
	mv energy.cu energy.cpp
	mpic++ -D USE_MPI $(MAIN) $(SOURCE) $(CPPENERGY) -o main.out -g2
	mv energy.cpp energy.cu

gpu-mpi:
	nvcc -D USE_CUDA -lib $(SOURCE) $(CUDAENERGY) --expt-relaxed-constexpr -o $(LIB)
	mpic++ -D USE_MPI $(MAIN) $(LIB) $(LDFLAGS) -o main.out
	rm $(LIB)
