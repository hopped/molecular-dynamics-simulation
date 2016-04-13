/*
 * Copyright (c) 2012-2014 Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "Lattice.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
using namespace std;

/** List of the names of the 7 Bravais lattices */
static const char* LatticeSystemNames[] = {
	"triclinic",
	"monoclinic",
	"orthorombic",
	"tetragonal",
	"rhomboedral",
	"hexagonal",
	"cubic"
};

/** List with the number of cetneres in the different lattice centering types */
static int LatticeCenteringNums[6] = { 1, 2, 4, 2, 2, 2 };

/** List of the names of the centerings */
static const char* LatticeCenteringNames[] = {
    "primitive",
    "body",
    "face",
    "base A",
    "base B",
    "base C"
};

/** Array holding the relative coordinates of the lattice centers in the a,b,c system */
static const double LatticeCenteringCoords[6][4][3] = {
    { /* primitive */
        {0.0, 0.0, 0.0}
    },
    { /* body */
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    },
    { /* face */
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}
    },
    { /* base A */
        {0.0, 0.0, 0.0},
        {0.0, 0.5, 0.5}
    },
    { /* base B */
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.5}
    },
    { /* base C */
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0}
    }
};

void Lattice::init(LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3]) {
	_system    = system;
	_centering = centering;
	for(int d = 0; d < 3; d++) {
		_a[d] = a[d];
		_b[d] = b[d];
		_c[d] = c[d];
	}
	_centeringCounter = 0;
}

void Lattice::initDims(LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3], long dimsMin[3], long dimsMax[3]) {
    init(system, centering, a, b, c);
    setDimsMin(dimsMin);
    setDimsMax(dimsMax);
}

int Lattice::getPoint(double* r) {

    /* after iterating over all centering positions for current cell, move to next cell */
    if(_centeringCounter >= LatticeCenteringNums[_centering]) {
        _centeringCounter = 0;
        _pos[0]++;
        /* for hexagonal lattic we have to skip mid points of hexagons */
        if( (_system == hexagonal) && ( (_pos[0] - _pos[1] + 2) % 3 == 0) ) {
            _pos[0]++;
        }
        if(_pos[0] >= _dimsMax[0]) {
            _pos[0] = _dimsMin[0];
            _pos[1]++;
            /* for hexagonal lattic we have to skip mid points of hexagons */
            if( (_system == hexagonal) && ( (_pos[0] - _pos[1] + 2) % 3 == 0) ) {
                _pos[0]++;
            }
            if(_pos[1] >= _dimsMax[1]) {
                _pos[1] = _dimsMin[1];
                _pos[2]++;
                if(_pos[2] >= _dimsMax[2]) {
                    return 0;
                }
            }
        }
    }

	double ia = _pos[0] + LatticeCenteringCoords[_centering][_centeringCounter][0];
	double ib = _pos[1] + LatticeCenteringCoords[_centering][_centeringCounter][1];
	double ic = _pos[2] + LatticeCenteringCoords[_centering][_centeringCounter][2];
	for(int d = 0; d < 3; d++) {
		r[d] = ia * _a[d] + ib * _b[d] + ic * _c[d];
	}
	_centeringCounter++;
	return 1;
}

int Lattice::checkValidity(){
	/* Checks for validity of system centering combination */
	switch(_system) {
		case triclinic:
			switch(_centering) {
				case primitive:
					break;
				default:
					return 0;
			}
			break;
		case monoclinic:
			switch(_centering) {
				case primitive:	case base_A: case base_B: case base_C:
					break;
				default:
					return 0;
			}
			break;
		case orthorombic:
			switch(_centering) {
				case primitive: case base_A: case base_B: case base_C: case body: case face:
					break;
				default:
					return 0;
			}
			break;
		case tetragonal:
			switch(_centering) {
				case primitive: case body:
					break;
				default:
					return 0;
			}
			break;
		case rhomboedral:
			switch(_centering) {
				case primitive:
					break;
				default:
					return 0;
			}
			break;
		case hexagonal:
			switch(_centering) {
				case primitive:
					break;
				default:
					return 0;
			}
			break;
		case cubic:
			switch(_centering) {
				case primitive: case body: case face:
					break;
				default:
					return 0;
			}
			break;
		default:
			return 0;
	}
	/* TODO: Check validity of lattice vectors for the given system. */
	return 1;
}

void Lattice::setDimsMin(long dimsMin[3]) {
    for(int d = 0; d < 3; d++) {
        _dimsMin[d] = dimsMin[d];
        _pos[d] = dimsMin[d];
    }
}

void Lattice::setDimsMax(long dimsMax[3]) {
    for(int d = 0; d < 3; d++) {
        _dimsMax[d] = dimsMax[d];
    }
}

const char* Lattice::systemName() {
    return LatticeSystemNames[_system];
}
const char* Lattice::centeringName() {
    return LatticeCenteringNames[_centering];
}

int Lattice::numCenters(LatticeCentering centering) {
    return LatticeCenteringNums[centering];
}
LatticeSystem Lattice::system(char* str) {
    if( strcmp( "triclinic", str ) == 0 ) { return triclinic; }
    else if( strcmp( "monoclinic", str ) == 0 ) { return monoclinic; }
    else if( strcmp( "orthorombic", str ) == 0 ) { return orthorombic; }
    else if( strcmp( "tetragonal", str ) == 0 ) { return tetragonal; }
    else if( strcmp( "rhomboedral", str ) == 0 ) { return rhomboedral; }
    else if( strcmp( "hexagonal", str ) == 0 ) { return hexagonal; }
    else if( strcmp( "cubic", str ) == 0 ) { return cubic; }
    else { return unknownSystem; }
}
LatticeCentering Lattice::centering(char* str) {
    if( strcmp( "primitive", str ) == 0 ) { return primitive; }
    else if( strcmp( "body", str ) == 0 ) { return body; }
    else if( strcmp( "face", str ) == 0 ) { return face; }
    else if( strcmp( "base_A", str ) == 0 ) { return base_A; }
    else if( strcmp( "base_B", str ) == 0 ) { return base_B; }
    else if( strcmp( "base_C", str ) == 0 ) { return base_C; }
    else { return unknownCentering; }
}

const double* Lattice::a() {
    return _a;
}

const double* Lattice::b() {
    return _b;
}

const double* Lattice::c() {
    return _c;
}

lattice_t* lattice_create() {
    return new lattice_t();
}

void lattice_destroy(lattice_t* lattice) {
    delete lattice;
}

void lattice_init(lattice_t* lattice, LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3]) {
    lattice->init(system, centering, a, b, c);
}

void lattice_initDims(lattice_t* lattice, LatticeSystem system, LatticeCentering centering, double a[3], double b[3], double c[3], long dimsMin[3], long dimsMax[3]) {
    lattice->initDims(system, centering, a, b, c, dimsMin, dimsMax);
}

int lattice_getPoint(lattice_t* lattice, double r[3]) {
    return lattice->getPoint(r);
}

int lattice_checkValidits(lattice_t *lattice) {
    return lattice->checkValidity();
}

const char* lattice_systemName(lattice_t* lattice) {
    return lattice->systemName();
}

const char* lattice_centeringName(lattice_t* lattice) {
    return lattice->centeringName();
}

int lattice_numCenters(LatticeCentering centering) {
    return Lattice::numCenters(centering);
}

LatticeSystem lattice_system(char* str) {
    return Lattice::system(str);
}

LatticeCentering lattice_centering(char* str) {
    return Lattice::centering(str);
}

const double* lattice_a(lattice_t* lattice) {
    return lattice->a();
}

const double* lattice_b(lattice_t* lattice) {
    return lattice->b();
}

const double* lattice_c(lattice_t* lattice) {
    return lattice->c();
}
