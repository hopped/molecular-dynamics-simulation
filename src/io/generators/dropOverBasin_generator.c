/*
 * Copyright (c) 2014      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "io/generators/dropOverBasin_generator.h"

#include "constants.h"
#include "utils/gaussrand.h"
#include "utils/generator/Lattice.h"
#include "utils/generator/Basis.h"
#include "utils/generator/Generator.h"
#include "utils/generator/Molecule.h"


#include "typedefs.h"
#include "molecule.h"
#include "basis.h"
#include "ensemble.h"
#include "domain.h"
#include "phasespace.h"
#include "unit_cell.h"
#include "utils/logger.h"

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>

void dropOverBasin_generator(ensemble_t* ensemble, double hBasin, double dropCenter[3], double dropRadius, domain_t* domain, phasespace_t* psp) {
    real E_avg = 0.;
    if (ensemble->type == NVE) {
        E_avg = ensemble->E / ensemble->N;
    } else
    if (ensemble->type == NVT) {
        /** @todo TODO: Take care about actual dof. */
        E_avg = 3./2. * kB * ensemble->T;
    } else {
        LOG(Error, "Can only handle NVE or NVT ensembles.");
        exit(1);
    }

    if(hBasin > dropCenter[2] - dropRadius) {
        LOG(Error, "Drop inside basin\n");
        exit(1);
    }
    real v_avg  = sqrt(2*E_avg / config.m);
    real v_avg1 = sqrt(1/3.)*v_avg; /* 1/3rd of energy in single direction*/


    /* Determine the number and diameter of unit cells. */
    basis_t *basis;
    basis = basis_create();
    molecule_t molecule;
    generator_molecule_t gmolecule;
    gmolecule.cid = 0;
    gmolecule.r[0] = 0.0;
    gmolecule.r[1] = 0.0;
    gmolecule.r[2] = 0.0;
    basis_addMolecule(basis, gmolecule);

    LatticeCentering centering = face;
    LatticeSystem system = cubic;
    LOG(Info, "Drop over basin generator\n");
    LOG(Info, "Box lengths: %"PRIreal" x %"PRIreal" x %"PRIreal"\n", domain->L[0], domain->L[1], domain->L[2]);

    /* create molecules in the basin */
    double rho = config.density;
    LOG(Info, "rho: %"PRIreal"\n", rho);
    long num_molecules_per_crystal_cell = basis_numMolecules(basis) * lattice_numCenters(centering);
    LOG(Debug, "Number of molecules per crystall cell: %ld\n", num_molecules_per_crystal_cell);
    double crystal_cell_volume = num_molecules_per_crystal_cell / rho;
    double l = pow(crystal_cell_volume, 1./3.);
    double a[3], b[3], c[3];
    a[0] = l; a[1] = 0; a[2] = 0;
    b[0] = 0; b[1] = l; b[2] = 0;
    c[0] = 0; c[1] = 0; c[2] = l;

    double VBasin = domain->L[0]*domain->L[0]*hBasin;
    double VSphere = 4./3. * M_PI * pow(dropRadius, 3.);
    /* approximate number of molecules, add 100 extra places to prevent segfaults. */
    ensemble->N = ceil((VBasin + VSphere) * rho) + 100;
    psp_init(psp, ensemble, domain); // Allocate memory for the molecules

    long count = 0;
    long countDrop = 0;
    real global_momentum[3] = {0., 0., 0.};
    real global_mass = 0.;

    double origin[3] = {l, l, l};
    double boxMax[3];
    boxMax[0] = domain->L[0] - l;
    boxMax[1] = domain->L[1] - l;
    boxMax[2] = hBasin;
    lattice_t *lattice;
    lattice = lattice_create();
    lattice_init(lattice, system, centering, a, b, c);

    generator_t *basinGenerator;
    basinGenerator = generator_create();
    generator_init_box(basinGenerator, lattice, basis, origin, origin, boxMax);

    generator_t *dropGenerator;
    dropGenerator = generator_create();
    generator_init_sphere(dropGenerator, lattice, basis, dropCenter, dropCenter, dropRadius);

    srand(config.gridgenerator_random_seed);
    /* Insert molecules to the domain */
    int rc = generator_getMolecule(basinGenerator, &gmolecule);
    /* If all molecules in the basin are created generate molecules for the drop. */
    if (rc == 0) {
        rc = generator_getMolecule(dropGenerator, &gmolecule);
        countDrop++;
    }
    while(rc != 0) {
        int d;
        if( config.gridgenerator_random_shift >= 0 ) {
            real rnd_a, rnd_b, rnd_c;
            real k = config.gridgenerator_random_shift;
            rnd_a = k * random() / RAND_MAX;
            rnd_b = k * random() / RAND_MAX;
            rnd_c = k * random() / RAND_MAX;
            for(d = 0; d < 3; d++) {
                gmolecule.r[d] += rnd_a*a[d] + rnd_b*b[d] + rnd_c*c[d];
            }
        }
        for(d = 0; d < 3; d++) {
            molecule.r[d] = gmolecule.r[d];
        }

        if( config.gridgenerator_dist_boltzmann == 1 ) {
            /* maxwell boltzmann velocity distribution
             */
            molecule.v[0] = v_avg1*gaussrand();
            molecule.v[1] = v_avg1*gaussrand();
            molecule.v[2] = v_avg1*gaussrand();
        }
        else {
            /* uniform distribution of velocity orientations
             * see e.g. http://mathworld.wolfram.com/SpherePointPicking.html
             */
            double omega = 2*M_PI * (double)random() / RAND_MAX;
            double phi = acos( (2.0 * (double)random() / RAND_MAX) - 1);
            molecule.v[0] = v_avg*sin(phi)*cos(omega);
            molecule.v[1] = v_avg*sin(phi)*sin(omega);
            molecule.v[2] = v_avg*cos(phi);
        }
        molecule.id = count;
        global_momentum[0] += config.m * molecule.v[0];
        global_momentum[1] += config.m * molecule.v[1];
        global_momentum[2] += config.m * molecule.v[2];
        global_mass += config.m;

        psp_add_molecule(psp, &molecule);
        count++;
        rc = generator_getMolecule(basinGenerator, &gmolecule);
        /* If all molecules in the basin are created generate molecules for the drop. */
        if (rc == 0) {
            rc = generator_getMolecule(dropGenerator, &gmolecule);
            countDrop += rc;
        }
    }
    LOG(Info, "Created molecules total %ld, drop: %ld\n", count, countDrop);
    if( config.gridgenerator_random_shift > 0 ) {
        psp_update(psp, domain);
    }

    generator_destroy(basinGenerator);
    basis_destroy(basis);
    lattice_destroy(lattice);

    real center_of_mass_velocity[3] = {0., 0., 0.};
    center_of_mass_velocity[0] = global_momentum[0] / global_mass;
    center_of_mass_velocity[1] = global_momentum[1] / global_mass;
    center_of_mass_velocity[2] = global_momentum[2] / global_mass;

    LOG(Info, "center of mass velocity: V = %"PRIreal", %"PRIreal", %"PRIreal"\n", center_of_mass_velocity[0], center_of_mass_velocity[1], center_of_mass_velocity[2]);
    mc_shift_velocity(psp->mc, center_of_mass_velocity);

    mc_reset_forces_and_momenta(psp->mc);
}
