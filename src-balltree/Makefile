CPP = g++
CFLAGS = -O2 -fopenmp -std=c++11 
LIBS = 

test: balltree.cpp
	$(CPP) $(CFLAGS)  balltree.cpp ${LIBS} -o balltree

clean: 
	rm balltree

run:
	./balltree 10000 784 /Users/nrajani/kdtree-knn/data/mnist-test.dat 16 30
