CC      = gcc
CXX     = g++
MPICXX  = mpicxx
CFLAGS  = -g -I../
CXXFLAGS = $(CFLAGS)
LDFLAGS = -g
AWK     = awk
GNUPLOT = gnuplot

default: Generator_test Generator_test_c Lattice_test LU_test LesSolver_test Graphite_lattice

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<


Lattice_test: Lattice_test.o ../Lattice.o
	$(CXX) $^ $(LDFLAGS) -o $@

Generator_test: Generator_test.o ../Generator.o ../Lattice.o ../Basis.o ../LU.o ../LA.o ../LesSolver.o
	$(CXX) $^ $(LDFLAGS) -o $@

Generator_test_c: Generator_test_c.o ../Generator.o ../Lattice.o ../Basis.o ../LU.o ../LA.o ../LesSolver.o
	$(CXX) $^ $(LDFLAGS) -o $@

Parallel_Generator_test.o: Parallel_Generator_test.cpp
	$(MPICXX) $^ $(CXXFLAGS) -c -o $@
Parallel_Generator_test: Parallel_Generator_test.o ../Generator.o ../Lattice.o ../Basis.o ../LU.o ../LA.o ../LesSolver.o
	$(MPICXX) $^ $(LDFLAGS) -o $@

LU_test: LU_test.o ../LU.o ../LA.o
	$(CXX) $^ $(LDFLAGS) -o $@

LesSolver_test: LesSolver_test.o ../LU.o ../LA.o ../LesSolver.o
	$(CXX) $^ $(LDFLAGS) -o $@

Graphite_lattice: Graphite_lattice.o ../Generator.o ../Lattice.o ../Basis.o ../LU.o ../LA.o ../LesSolver.o
	$(CXX) $^ $(LDFLAGS) -o $@

graphite-wall.dat: Graphite_lattice
	$(shell ./Graphite_lattice | $(AWK) '$$1 == 0 {print $2" "$3" "$4;}' > graphite-wall.dat)
Graphite_plots: graphite-wall.dat
	$(GNUPLOT) hexplot.gp

.PHONY: clean

clean:
	$(RM) *.o
	$(RM) Generator_test Generator_test_c Lattice_test LU_test LesSolver_test Graphite_lattice
	$(RM) *.dat
	$(RM) *.svg
