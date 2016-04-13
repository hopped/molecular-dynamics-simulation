CC = nvcc
OMPFLAGS = -fopenmp
CFLAGS_RELEASE = -arch=sm_13 --compiler-options '-O3'
CFLAGS_DEBUG   = -arch=sm_13 -G --compiler-options="-O0,  -pg"  
#CFLAGS = -O3 -msse4.1 -ftree-vectorizer-verbose=9

LD = gcc
LDFLAGS_RELEASE = -arch=sm_13 -O3
LDFLAGS_DEBUG   = arch=sm_13 -g -ggdb

