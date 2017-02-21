all: driver




all: driver

driver: driver.o
	g++ -g -std=c++11 driver.o -o nsdg -llapacke -lblas

driver.o: driver.cpp
	g++ -g -std=c++11 -c driver.cpp -llapacke -lblas

clean:
	rm -rf *.o test *.gnu *.temp *.jpg load *.png *.dat *.vtk *~ *.gch