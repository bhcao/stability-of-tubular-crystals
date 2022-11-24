all:
	$(CXX) $(wildcard ./*.cpp) -o main.out -g2

debug: all
	gdb ./main.out

clean:
	rm -f *.dump *.data
	rm -rf build
	rm -f *.out *.a
	rm -f energy*.png
	rm -f __*__temp__*
	rm -f restart*.bin
