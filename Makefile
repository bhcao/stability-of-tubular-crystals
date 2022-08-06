TARGET = $(shell ls | grep .cpp$ | grep -v main.cpp)

debug:
	g++ main.cpp $(TARGET) -g2 -o main.out
	./main.out

release:
	g++ main.c $(TARGET) -o main.out
