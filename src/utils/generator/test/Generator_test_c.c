#include <stdio.h>

#include "Lattice.h"
#include "Basis.h"
#include "Generator.h"


int main(int argc, char *argv[]) {
    double a[3] = {1.0, 0.0, 0.0};
    double b[3] = {0.0, 1.0, 0.0};
    double c[3] = {0.0, 0.0, 1.0};
    long startdims[3] = {0, 0, 0};
    long dims[3] = {2, 2, 2};
    lattice_t *lattice;
    lattice = lattice_create();
    lattice_initDims(lattice, cubic, face, a, b, c, startdims, dims);

    basis_t *base;
    base = basis_create();

    molecule_t molecule;

    molecule.cid = 0;
    molecule.r[0] = 0.0;
    molecule.r[1] = 0.0;
    molecule.r[2] = 0.0;
    basis_addMolecule(base, molecule);

    molecule.cid = 1;
    molecule.r[0] = 0.5*a[0];
    molecule.r[1] = 0.5*a[1];
    molecule.r[2] = 0.5*a[2];
    basis_addMolecule(base, molecule);

    double origin[3] = {0.1, 0.2, 0.3};

    generator_t *gen;
    gen = generator_create();
    generator_init(gen, lattice, base, origin);

    while(generator_getMolecule(gen, &molecule)) {
        printf("%ld\t%lf, %lf, %lf\n", molecule.cid, molecule.r[0], molecule.r[1], molecule.r[2]);
    }
    return 0;
}

