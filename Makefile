TARGET = $(shell ls | grep .c$ | grep -v main.c)

debug:
	g++ main.c $(TARGET) -g2 -o main.out
	./main.out

release:
	g++ main.c $(TARGET) -o main.out
