all: clean matmult

matmult: matmult.cpp
	g++ -std=c++17 matmult.cpp -o matmult

clean:
	rm -f matmult matmult.o