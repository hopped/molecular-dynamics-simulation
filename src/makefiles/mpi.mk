CC              = mpicc
CXX             = mpiCC
OMPFLAGS        = -fopenmp
CFLAGS_RELEASE  = -DNDEBUG -O3 -Wall -DNDEBUG
CFLAGS_PROF     = -DNDEBUG -O3 -Wall -DNDEBUG -pg
CFLAGS_DEBUG    = -O0 -g -ggdb -Wall

LD              = mpiCC
LDFLAGS_RELEASE = -O3
LDFLAGS_DEBUG   = -g -ggdb
LDFLAGS_PROF    = -O3 -pg

