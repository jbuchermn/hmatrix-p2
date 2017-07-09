CXX=clang++
LIBS=-L./lib -lstdc++ -llapacke -llapack -lcblas -lblas -lgfortran
INCS=-I./lib/lapacke
OPTS=-Wall -O3 -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

default: main

%.o: %.cpp
	$(CXX) -o $@ -c $< $(INCS) $(OPTS)

main: main.o hMatrixP2.o
	$(CXX) -o main main.o hMatrixP2.o $(LIBS)

clean:
	rm -f main
	rm -f *.o
