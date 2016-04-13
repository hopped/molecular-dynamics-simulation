#include <iostream>

#include "Lattice.h"


using namespace std;

int main(int argc, char *argv[]) {
    double a[3] = {1.0, 0.0, 0.0};
    double b[3] = {0.0, 1.0, 0.0};
    double c[3] = {0.0, 0.0, 1.0};
    long startdims[3] = {0,0,0};
    long dims[3] = {2, 2, 2};
    Lattice lat;
    lat.initDims(cubic, face, a, b, c, startdims, dims);
    double r[3];
    while(lat.getPoint(r)) {
        cout << r[0] << ", " << r[1] << ", " << r[2] << endl;
    }
    return 0;
}

