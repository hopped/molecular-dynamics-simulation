#include <iostream>

#include "Lattice.h"
#include "Basis.h"
#include "Generator.h"
#include "LA.h"


using namespace std;

int main(int argc, char *argv[]) {
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
    //double origin[3] = {0., 0., 0.};

    Generator gen;
    //gen.init(lat, base, origin);

    //double bBoxMin[3] = {0.5, -1.0, 0.0};
    //double bBoxMax[3] = {4.0, 2.5, 3.5};
    double bBoxMin[3] = {0, 0, 0};
    double bBoxMax[3] = {2, 1, 1};
    gen.init(lat, base, origin, bBoxMin, bBoxMax);

    while(gen.getMolecule(&molecule)) {
        cout << molecule.cid << "\t" << molecule.r[0] << ", " << molecule.r[1] << ", " << molecule.r[2] << endl;
    }
    return 0;
}

