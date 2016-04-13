/*
 * Copyright (c) 2012-2014 Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "io/generators/grid_generator.h"

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
#include "utils/logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>



void grid_generator(ensemble_t *ensemble, LatticeSystem system, LatticeCentering centering, domain_t *domain, phasespace_t *psp) {

    real E_avg = 0.;
    if (ensemble->type == NVE) {
        E_avg = ensemble->E / ensemble->N;
    } else
    if (ensemble->type == NVT) {
        // TODO: Take care about actual dof.
        E_avg = 3./2. * kB * ensemble->T;
    } else {
        LOG(Error, "Can only handle NVE or NVT ensembles.");
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

    LOG(Info, "Grid generator\n");
    LOG(Debug, "Grid system: %d\n", system);
    LOG(Debug, "Grid centering: %d\n", centering);

    long num_molecules_per_crystal_cell = basis_numMolecules(basis) * lattice_numCenters(centering);
    LOG(Debug, "Number of molecules per crystall cell: %ld\n", num_molecules_per_crystal_cell);
    size_t num_crystal_cells = ensemble->N / num_molecules_per_crystal_cell; /* number of crystal cells */
    long num_crystal_cells_per_dim[3];
    num_crystal_cells_per_dim[0] = round(pow(num_crystal_cells, 1./3.));
    num_crystal_cells_per_dim[1] = num_crystal_cells_per_dim[0];
    num_crystal_cells_per_dim[2] = num_crystal_cells_per_dim[0];
    num_crystal_cells = num_crystal_cells_per_dim[0] * num_crystal_cells_per_dim[1] * num_crystal_cells_per_dim[2];
    LOG(Debug, "Number of crystall cells: %ld (%ld x %ld x %ld)\n", num_crystal_cells, num_crystal_cells_per_dim[0], num_crystal_cells_per_dim[1], num_crystal_cells_per_dim[2]);
    real rho = ensemble->N / domain_get_volume(domain);
    LOG(Info, "rho: %"PRIreal"\n", rho);

    size_t num_molecules = num_molecules_per_crystal_cell * num_crystal_cells_per_dim[0] * num_crystal_cells_per_dim[1] * num_crystal_cells_per_dim[2];
    if( num_molecules != ensemble->N ) {
        LOG(Warning, "Changing number of molecules from N=%lu to N=%lu\n", ensemble->N, num_molecules);
        ensemble->N = num_molecules;
        real old_V = domain_get_volume(domain);
        real new_V = num_molecules / rho;
        real box_length[3];
        box_length[0] = pow(new_V, 1./3.);
        box_length[1] = box_length[0];
        box_length[2] = box_length[0];
        domain_set_size(domain, box_length);

        LOG(Warning, "Changing volume from V=%"PRIreal" to V=%"PRIreal"\n", old_V, domain_get_volume(domain));
    }
    psp_init(psp, ensemble, domain); // Allocate memory for the molecules

    long count = 0;
    real global_momentum[3] = {0., 0., 0.};
    real global_mass = 0.;

    double a[3], b[3], c[3];
    LOG(Info, "Box lengths: %"PRIreal" x %"PRIreal" x %"PRIreal"\n", domain->L[0], domain->L[1], domain->L[2]);
    a[0] = domain->L[0] / num_crystal_cells_per_dim[0];
    a[1] = 0;
    a[2] = 0;
    b[0] = 0;
    b[1] = domain->L[1] / num_crystal_cells_per_dim[1];
    b[2] = 0;
    c[0] = 0;
    c[1] = 0;
    c[2] = domain->L[2] / num_crystal_cells_per_dim[2];
    lattice_t *lattice;
    lattice = lattice_create();
    lattice_init(lattice, system, centering, a, b, c);

    double origin[3] = {0.0, 0.0, 0.0};
    double boxMax[3];
    /* Calculate the max corner to circumvent problems due to numerical errors in lattice vector computation instead of using domain->L. */
    boxMax[0] = a[0] * num_crystal_cells_per_dim[0] + b[0] * num_crystal_cells_per_dim[1] + c[0] * num_crystal_cells_per_dim[2] - 0.1;
    boxMax[1] = a[1] * num_crystal_cells_per_dim[0] + b[1] * num_crystal_cells_per_dim[1] + c[1] * num_crystal_cells_per_dim[2] - 0.1;
    boxMax[2] = a[2] * num_crystal_cells_per_dim[0] + b[2] * num_crystal_cells_per_dim[1] + c[2] * num_crystal_cells_per_dim[2] - 0.1;

    generator_t *generator;
    generator = generator_create();
    generator_init_box(generator, lattice, basis, origin, origin, boxMax);

    srand(config.gridgenerator_random_seed);
    /* Insert molecules to the domain */
    while(generator_getMolecule(generator, &gmolecule)) {
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
    }

    if( config.gridgenerator_random_shift > 0 ) {
        psp_update(psp, domain);
    }

    generator_destroy(generator);
    basis_destroy(basis);
    lattice_destroy(lattice);

    assert(num_molecules == count);

    real center_of_mass_velocity[3] = {0., 0., 0.};
    center_of_mass_velocity[0] = global_momentum[0] / global_mass;
    center_of_mass_velocity[1] = global_momentum[1] / global_mass;
    center_of_mass_velocity[2] = global_momentum[2] / global_mass;

    LOG(Info, "center of mass velocity: V = %"PRIreal", %"PRIreal", %"PRIreal"\n", center_of_mass_velocity[0], center_of_mass_velocity[1], center_of_mass_velocity[2]);
    mc_shift_velocity(psp->mc, center_of_mass_velocity);
}

