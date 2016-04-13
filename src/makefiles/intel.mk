CC              = icc
CXX             = icpc
OMPFLAGS        = -openmp
CFLAGS_RELEASE  = -DNDEBUG -O3 -vec-report3
CFLAGS_DEBUG    = -O0 -g -vec-report3

LD              = icpc
LDFLAGS_RELEASE = -O3 -vec-report3
LDFLAGS_DEBUG   = -O0 -g -vec-report3

