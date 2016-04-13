#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "phasespace.h"

#ifdef MPI_UTILS_C
int Sx, Sy, Sz, NSD;
#else
extern int Sx, Sy, Sz, NSD;
#endif

void MPI_Create_dims(int myrank);
void MPI_Map(phasespace_t *psp, int myrank);

#endif

