CXX=g++
LIBS=-L./lib -lstdc++ -llapacke -llapack -lcblas -lblas -lgfortran
INCS=-I./lib/lapacke -I./lib
OPTS=-Wall -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP 

#OPTS+=-DEIGEN_DONT_VECTORIZE # For comparison
OPTS+=-Ofast # Significant speedup compared to -O3

default: main

%.o: %.cpp
	$(CXX) -o $@ -c $< $(INCS) $(OPTS)

main: main.o hMatrixP2.o
	$(CXX) -o main main.o hMatrixP2.o $(LIBS)

clean:
	rm -f main
	rm -f *.o
