/*
 * Copyright (c) 2012-2014 Christoph Niethammer <christoph.niethammer@googlemail.com>
 * Copyright (c) 2012-2014 Mohamad Amer Wafai <amerwafai@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef CONFIG_H
#define CONFIG_H

#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_STRING "0.1"

#define MOLECULE_BLOCK_CAPACITY 1024  /**< allocation unit for container */
#define CELL_CAPACITY 8

/* Parameters for the velocity scaling thermostat */
#define VELOCITY_SCALING_FACTOR_MAX 1.05 /**< maximal allowed global velocity scaling factor */
#define VELOCITY_SCALING_FACTOR_MIN 0.95 /**< minimal allowed global velocity scaling factor */
#define VELOCITY_SCALING_E_KIN_LIMIT_FACTOR 10.0 /**< In case of an out of range scaling factor limit for kinetic energy per molecule as factor times target temperature. */
/* #define VELOCITY_SCALING_E_ROT_LIMIT_FACTOR   may be necessary in future with RDF for molecules */

#include "typedefs.h"
#include "utils/generator/Lattice.h"

typedef struct {
    int povray_output;
    int ascii_output;
    real cutoff_radius;
    real cutoff_radius_sq;
    size_t count;
    real m; /**< mass of molecules */
    unsigned int molecule_container_type;
    LatticeSystem lattice_system;
    LatticeCentering lattice_centering;
    int gridgenerator_random_seed;
    int gridgenerator_dist_boltzmann;
    double gridgenerator_random_shift;
    double externalForceField[3];
    double generator_dropRadius;
    double density;
} config_t;

#ifdef MPI
enum {
    X=0, Y, Z
};
enum {
    p=0, n
};
#endif

#define DISABLED 0
#define ENABLED  1

#define USE_STREAMS

#ifndef CONFIG_SRC
extern config_t config;
#endif

#endif
