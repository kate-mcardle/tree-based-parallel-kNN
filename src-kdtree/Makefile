kdlib = ../libkdtree.a
CPP = gcc
CFLAGS = -Wall -Wextra -std=c++0x -lstdc++ -O2 -g3 -fopenmp -I..
test: 	test.cpp 
	$(CPP) $(CFLAGS) test.cpp  ${kdlib} -o test

clean:	
	rm -f test

run: 
	./test 20640 8 /scratch/02336/nrajani/data/cadata 32
