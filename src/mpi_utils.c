#define MPI_UTILS_C
#include "mpi_utils.h"
#include "moleculecontainer/basicN2/basicN2.h"

void MPI_Create_dims(int myrank){
    int size;
    int dims[3] = {0, 0, 0};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Dims_create(size, 3, dims);
    Sx = dims[0];
    Sy = dims[1];
    Sz = dims[2];
    NSD = dims[0] * dims[1] * dims[2];
//    if (myrank == 0) {
//        printf("The number of subdomains are as follow: Sx=%d, Sy=%d, Sz=%d and Number of Subdomains are %d\n\n",
//               Sx,
//               Sy,
//               Sz,
//               NSD);
//    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_Map(phasespace_t *psp, int myrank){
    basicN2 * mc = (basicN2 *)psp->mc->container;
    mc->myrank = myrank;
    int layer = Sx * Sy;
    mc->Nx[0] = ( myrank / Sx ) * Sx + ( myrank + 1 ) % Sx;
    mc->Ny[0] = ( myrank + Sx ) % ( layer ) + ( myrank / layer ) * layer;
    mc->Nz[0] = ( myrank + layer ) % NSD;
    mc->Nx[1] = ( myrank / Sx ) * Sx + (( myrank - 1 ) % Sx + Sx ) % Sx;
    mc->Ny[1] = (( myrank - Sx ) % ( layer ) + layer ) % layer + ( myrank / layer ) * layer;
    mc->Nz[1] = (( myrank - layer ) % NSD + NSD ) % NSD;
}

