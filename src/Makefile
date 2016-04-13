#
# Copyright (c) 2012-2013 Christoph Niethammer <christoph.niethammer@gmail.com>
#
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER
# 

COMPILER ?= gcc
TARGET ?= DEBUG
MPI  ?= 0
OMP  ?= 0
ARCH ?= CPU
NVCC ?=
APP  ?= CMD

ifeq ($(MPI),1)
COMPILER = mpi
CFLAGS   = -DMPI
endif
ifeq ($(OMP),1)
CFLAGS   += $(OMPFLAGS) 
LDFLAGS   = $(OMPFLAGS)
else
ifeq ($(COMPILER),cray)
CFLAGS   += -h noomp
endif
endif

COMPILER_CONFIG = makefiles/$(COMPILER).mk
include $(COMPILER_CONFIG)

GENERAL_CFLAGS  = -I.
GENERAL_LDFLAGS = -lm


CFLAGS  += $(GENERAL_CFLAGS)  $(CFLAGS_$(TARGET))
CXXFLAGS = $(CFLAGS)
LDFLAGS += $(GENERAL_LDFLAGS) $(LDFLAGS_$(TARGET))

SOURCES  = $(shell find . -name "*.c" | grep -v test)
SOURCES_CXX = $(shell find . -name "*.cpp" | grep -v test)
HEADERS  = $(shell find . -name "*.h")
ifeq ($(MPI),0)
SOURCES := $(filter-out ./mpi_utils.c,$(SOURCES))
HEADERS := $(filter-out ./mpi_utils.h,$(HEADERS))
endif

PAPI_BASE   =$(shell which papi_avail|sed -e 's/\/bin\/papi\_avail//')
CUDA_DIR  =$(shell which nvcc|sed -e 's/\/bin\/nvcc//')
ifneq ($(PAPI_BASE),)
        LDFLAGS += -L$(PAPI_BASE)/lib -lpapi
endif

ifeq ($(ARCH),GPU)
NVCC = nvcc
LDFLAGS += -L$(CUDA_DIR)/lib64 -lcudart
CFLAGS += -DGPU
GPUFLAGS_R = -I. -O3 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_13,code=sm_13
GPUFLAGS_D = -I. -O0 -g -G  -gencode arch=compute_35,code=sm_35 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_13,code=sm_13
SOURCES  += $(shell find . -name "*.cu")
HEADERS  += $(shell find . -name "*.cuh")

ifeq ($(TARGET),DEBUG)
LDFLAGS += -O0
%.o:%.cu
	$(NVCC) $(GPUFLAGS_D) -c  -m64 $< -o $@
else
ifeq ($(TARGET),RELEASE)
LDFLAGS += -O3
%.o:%.cu
	$(NVCC) $(GPUFLAGS_R) -c  -m64 $< -o $@
endif
endif
endif

BINARIES = main
OBJECTS  := $(SOURCES:.c=.o) $(SOURCES_CXX:.cpp=.o)
OBJECTS  := $(OBJECTS:.cu=.o)

default: main
test: test_basicn2

main: $(OBJECTS)
	$(LD) -o $(APP)_$(ARCH) $(OBJECTS) $(LDFLAGS)
	
%.o:%.c
	$(CC) -c $(CFLAGS) $< -o $@ 
%.o:%.cpp                           
	$(CXX) -c $(CXXFLAGS) $< -o $@
	
uncrustify:
	uncrustify -c contrib/uncrustify.cfg $(SOURCES) $(HEADERS)
	./contrib/replace
        
.PHONY: clean clean_output distlcean

clean:
	rm -f $(BINARIES)
	find . -name "*.o" -delete
	find . -name "*~"  -delete

clean_output:
	rm -f *.dat
	rm -f *pov

distclean: clean clean_output


%.obj:%
	objdump -D $< > $@

.PHONY: linecount
linecount:
	find . -name "*.c" -or -name "*.h" -or -name Makefile | xargs wc -l	

.PHONY: testrun1 testrun2 testrun3 ljtest test_sc test_basicn2
test_basicn2: main
	./main --verbose=5 -N $(shell perl -e "print 2*8**3") -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 --ascii-output --povray-output --simulation-end-time=0.1 --molecule-container=BASICN2 --thermostat=velocity-scaling

test_sc: main
	./main --verbose=5 -N $(shell perl -e "print 8**3") -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 --ascii-output --povray-output --simulation-end-time=1
testrun1: main
	./main -v -N $(shell perl -e "print 4 * 4**3") -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=1 -m 1 --ascii-output --povray-output
testrun2: main
	./main -v -N $(shell perl -e "print 4 * 4**3") -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=3 -m 1 --ascii-output --povray-output
testrun3: main
	./main -v -N $(shell perl -e "print 4 * 4**3") -n 1 -T 0.85 --domain-type=cuboid --cutoff-radius=5 -m 1 --ascii-output --povray-output

# Scenario used in http://www.cstl.nist.gov/srs/LJ_PURE/md.htm
ljtest: main
	$(info Result schould be U_pot=-6.234(1), p=2.544(6), D=0.027)
	./main -v -N 500 -n 0.9 -T 0.85 --domain-type=cuboid --timestep-length=0.005 --cutoff-radius=3 -m 1 --simulation-end-time=150 --simulation-equlibration=50 --molecule-container=BASICN2 --thermostat=velocity-scaling

.PHONY: genpov
genpov:
	povray *.pov
