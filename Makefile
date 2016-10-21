default:
	make array
	make gamete_test

array:
	g++ -std=c++1y -O3 -o lib/array.o -c lib/array/src/array.cpp

gamete_test:
	g++ -std=c++1y lib/array.o -fopenmp -O3 src/gamete_test.cpp -o bin/gamete_test_omp
	g++ -std=c++1y lib/array.o -O3 src/gamete_test.cpp -o bin/gamete_test 

clean:
	rm -f bin/*

