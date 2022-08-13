MAIN = main.cpp
SOURCE = energy.cpp model.cpp molecule.cpp
CPPENERGY = update.cpp adjacent.cpp
CUDAENERGY = update.cu adjacent.cu
LIB = nano.a
LDFLAGS = -L/usr/lib/cuda/lib64 -lcudart

NVCC = nvcc -D USE_CUDA --expt-relaxed-constexpr
MPICXX = mpic++ -D USE_MPI

debug: update.cpp adjacent.cpp
	$(CXX) $(MAIN) $(SOURCE) $(CPPENERGY) -o main.out -g2
	gdb ./main.out

gpu: adjacent.cu update.cu
	$(NVCC) $(MAIN) $(SOURCE) $(CUDAENERGY) -o main.out 

cpu: update.cpp adjacent.cpp
	$(CXX) $(MAIN) $(SOURCE) $(CPPENERGY) -o main.out

cpu-mpi: update.cpp adjacent.cpp
	$(MPICXX) $(MAIN) $(SOURCE) $(CPPENERGY) -o main.out

gpu-mpi: adjacent.cu update.cu
	$(NVCC) -lib $(SOURCE) $(CUDAENERGY) -o $(LIB)
	$(MPICXX) $(MAIN) $(LIB) $(LDFLAGS) -o main.out

clean: adjacent.cu update.cu
	rm -f *.dump *.data
	rm -rf build
	rm -f *.out *.a
	rm -f energy*.png
	rm -f __*__temp__*

adjacent.cu update.cu:
	mv adjacent.cpp adjacent.cu
	mv update.cpp update.cu

adjacent.cpp update.cpp:
	mv adjacent.cu adjacent.cpp
	mv update.cu update.cpp