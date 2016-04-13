#include <iostream>

#include "Lattice.h"
#include "Basis.h"
#include "Generator.h"
#include "LA.h"

#include <cmath>


using namespace std;

/* Create Graphite lattice according to
 * http://arxiv.org/pdf/1001.2681.pdf
 */
int main(int argc, char *argv[]) {
    double a[3] = {0.3816, 0.0, 0.0};
    double b[3] = {0.3816*0.5, 0.3816*(sqrt(3)/2.0), 0.0};
    //double c[3] = {0.0, 0.0, 0.8996};
    double c[3] = {0.0, 0.0, 2*0.8996};
    Lattice lat;
    lat.init(hexagonal, primitive, a, b, c);

    Basis base;

    molecule_t molecule;

    molecule.cid = 0;
    molecule.r[0] = 0.0;
    molecule.r[1] = 0.0;
    molecule.r[2] = 0.0;
    base.addMolecule(molecule);

    molecule.r[0] = a[0];
    molecule.r[1] = 0.0;
    molecule.r[2] = c[2]/2;
    base.addMolecule(molecule);

    Generator gen;
    double bBoxMin[3] = {0, 0, 0};
    //double bBoxMax[3] = {3, 3, c[2]};
    double bBoxMax[3] = {20, 20, 4*0.8996};
    double origin[3] = {0.0, 0.0, 0.0};
    gen.init(lat, base, origin, bBoxMin, bBoxMax);

    while(gen.getMolecule(&molecule)) {
        cout << molecule.cid << "\t" << molecule.r[0] << ", " << molecule.r[1] << ", " << molecule.r[2] << endl;
    }
    return 0;
}

