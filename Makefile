default:
	make gamete_test


gamete_test:
	g++ -std=c++1y -fopenmp -O3 src/CSVRow.cpp src/metrics.cpp src/gamete_test.cpp -o bin/gamete_test_omp -pedantic
	g++ -std=c++1y -O3 src/CSVRow.cpp src/gamete_test.cpp src/metrics.cpp -o bin/gamete_test -pedantic

clean:
	rm -f bin/*

