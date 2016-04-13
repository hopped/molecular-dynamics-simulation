/*
 * Copyright (c) 2013-2014 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "Generator.h"

#include "LA.h"
#include "LesSolver.h"
#include "Objects.h"

#include <iostream>
#include <cmath>


void Generator::init(Lattice& lattice, Basis& basis, double origin[3], Object *object) {
    _lattice = lattice;
    _basis = basis;
    _object = object;

    object->getBboxMin(_bBoxMin);
    object->getBboxMax(_bBoxMax);
    for(int d = 0; d < 3; d++) {
        _origin[d] = origin[d];
    }
    _baseCount = 0;

    /* transpose (a,b,c) */
    double *A[3];
    for(int i = 0; i < 3; i++) {
        A[i] = new double[3];
        A[i][0] = _lattice.a()[i];
        A[i][1] = _lattice.b()[i];
        A[i][2] = _lattice.c()[i];
    }

    /* As we use grid coordinates the origin has to be subtracted.
     * This displacement will be added later on again in getMolecule. */
    double boxMin[3], boxMax[3];
    for(int d = 0; d < 3; d++) {
        boxMin[d] = _bBoxMin[d] - _origin[d];
        boxMax[d] = _bBoxMax[d] - _origin[d];
    }
    double x0[3];
    double x1[3];
    LesSolve(A, (double*) boxMin, 3, x0);
    LesSolve(A, (double*) boxMax, 3, x1);

    long startdims[3];
    long enddims[3];
    for(int d = 0; d < 3; d++) {
        startdims[d] = floor(x0[d]);
        enddims[d] = ceil(x1[d]);
    }

    double rtmp[3];
    double x[3];
    for(int j = 0; j < 3; j++) {
        for(int d = 0; d < 3; d++) {
            rtmp[d] = boxMin[d];
        }
        rtmp[j] = boxMax[j];
        LesSolve(A, rtmp, 3, x);
        for(int d = 0; d < 3; d++) {
            long coord = floor(x[d]);
            if(coord < startdims[d]) { startdims[d] = coord; }
            coord = ceil(x[d]);
            if(coord > enddims[d]) { enddims[d] = coord; }
        }
    }

    _lattice.setDimsMin(startdims);
    _lattice.setDimsMax(enddims);
}

int Generator::getMolecule(generator_molecule_t* molecule) {
    for(;;) {
        if(_baseCount == 0) {
            if(_lattice.getPoint(_lattice_point) == 0) {
                return 0;
            }
        }

        generator_molecule_t molecule_base;
        molecule_base = _basis.getMolecule(_baseCount);

        for(int d = 0; d < 3; d++) {
            molecule->r[d] = _origin[d] + _lattice_point[d] + molecule_base.r[d];
        }
        molecule->cid = molecule_base.cid;

        _baseCount = (_baseCount + 1) % _basis.numMolecules();

        /* _bBoxMin[0] != _bBoxMax[0] means "no bounding box mode */
        if(_bBoxMin[0] != _bBoxMax[0] && _object->isInside(molecule->r)) {
            break;
        }
    }

	return 1;
}

generator_t* generator_create() {
    return new generator_t();
}

void generator_destroy(generator_t* generator) {
    delete generator;
}

void generator_init_box(generator_t* generator, lattice_t* lattice, basis_t* basis, double origin[3], double bBoxMin[3], double bBoxMax[3]) {
    Cuboid *box = new Cuboid(bBoxMin, bBoxMax);
    generator->init(*lattice, *basis, origin, box);
}

void generator_init_sphere(generator_t* generator, lattice_t* lattice, basis_t* basis, double origin[3], double center[3], double radius) {
    Sphere *sphere = new Sphere(center, radius);
    generator->init(*lattice, *basis, origin, sphere);
}

int generator_getMolecule(generator_t* generator, generator_molecule_t* molecule) {
    return generator->getMolecule(molecule);
}
