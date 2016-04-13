#include <iostream>

#include <mpi.h>

#include "Lattice.h"
#include "Basis.h"
#include "Generator.h"
#include "LA.h"


using namespace std;

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_size = 1;
    int comm_rank = 0;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);

    int dims[3] = {0,0,0};
    MPI_Dims_create(comm_size, 3, &dims[0]);

    double a[3] = {1.0, 0.0, 0.0};
    double b[3] = {1.0, 1.0, 0.0};
    double c[3] = {0.0, 0.0, 1.0};
    Lattice lat;
    lat.init(cubic, face, a, b, c);

    Basis base;

    molecule_t molecule;

    molecule.cid = 0;
    molecule.r[0] = 0.0;
    molecule.r[1] = 0.0;
    molecule.r[2] = 0.0;
    base.addMolecule(molecule);

    molecule.cid = 1;
    molecule.r[0] = 0.5*a[0];
    molecule.r[1] = 0.5*a[1];
    molecule.r[2] = 0.5*a[2];
    base.addMolecule(molecule);

    double origin[3] = {10.0, 0.0, 0.55};

    Generator gen;

    double bBoxMin[3] = {0, 0, 0};
    double bBoxMax[3] = {2, 1, 1};

    double sbBoxMin[3];
    double sbBoxMax[3];

    for(int i = 0; i < dims[0]; i++) {
        sbBoxMin[0] = bBoxMin[0] + (i + 0) * (bBoxMax[0] - bBoxMin[0]) / dims[0];
        sbBoxMax[0] = bBoxMin[0] + (i + 1) * (bBoxMax[0] - bBoxMin[0]) / dims[0];

        for(int j = 0; j < dims[1]; j++) {
            sbBoxMin[1] = bBoxMin[1] + (j + 0) * (bBoxMax[1] - bBoxMin[1]) / dims[1];
            sbBoxMax[1] = bBoxMin[1] + (j + 1) * (bBoxMax[1] - bBoxMin[1]) / dims[1];

            for(int k = 0; k < dims[2]; k++) {
                sbBoxMin[2] = bBoxMin[2] + (k + 0) * (bBoxMax[2] - bBoxMin[2]) / dims[2];
                sbBoxMax[2] = bBoxMin[2] + (k + 1) * (bBoxMax[2] - bBoxMin[2]) / dims[2];

                gen.init(lat, base, origin, sbBoxMin, sbBoxMax);
                if(comm_rank == 0) {
                    while(gen.getMolecule(&molecule)) {
                        cout << molecule.cid << "\t" << molecule.r[0] << ", " << molecule.r[1] << ", " << molecule.r[2] << endl;
                    }
                }
            }
        }
    }

    MPI_Finalize();
    return 0;
}

