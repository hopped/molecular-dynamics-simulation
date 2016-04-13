CC              = pgcc
CXX             = pgCC
OMPFLAGS        = -mp
CFLAGS_RELEASE  = -DNDEBUG -O3 -Minfo=opt -Mneginfo=opt -Mvect=levels:10,sse
CFLAGS_DEBUG    = -O0 -g -Minfo=opt -Mneginfo=opt

LD              = pgCC
LDFLAGS_RELEASE  = -O3
LDFLAGS_DEBUG    = -O0 -g

