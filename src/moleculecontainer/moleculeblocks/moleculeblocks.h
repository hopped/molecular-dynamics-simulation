/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef MOLECULE_CONTAINER_H
#define MOLECULE_CONTAINER_H

#include "config.h"
#include "typedefs.h"

#include "../moleculecontainer.h"
#include "molecule.h"
#include "molecule_cell.h"
#include "ensemble.h"
#include "domain.h"

#include <stdio.h>

// TODO: Molecules should be stored by component later on to ease vecroization.

/** molecule container 
 *
 * consists out of cells with diameter >= cutoff radius 
 * used in the link cell agorithm.
 *
 * cell coordinates start with (0,0,0) in the halo and with (1,1,1) in the
 * real domain.
 */
typedef struct mb_t {
    real cutoff_radius;
    size_t num_cells;
    long num_molecules; /**< Number of molecules stored in this container (including halos, ...). */
    long num_cells_per_dim[3];
    real cell_size[3];
    molecule_cell_t *cells;
} mb_t;

void* mb_alloc(); /**< allocate a new molecule container */
void mb_free(void *container); /**< free molecule container */
void mb_decompose(void *container, domain_t *domain);
void mb_destroy(void * container); /**< destroy and deallocate help memory */
void mb_estimate_halos( void * container, domain_t * domain);
void mb_init(void *container, ensemble_t *ensemble, domain_t * domain); 
void mb_update(void *container, domain_t *domain);
void mb_print_ascii(void *container, char *filename);
void mb_print_pov(void *container, char *filename);

void mb_save_to_file(char *filename, void *container);
void mb_read_from_file(char *filename, void *container);

void mb_add_molecule(void *container, molecule_t *m);

long mb_get_num_molecules(void *container);

void mb_calc_forces(void *container, real *U_pot);

void mb_print_stats(void *container, ensemble_t *ensemble);
void mb_clear_halo(void *container);
void mb_reset_forces_and_momenta(void *container);
void mb_shift_velocity(void *container, real v[3]);

void mb_update_halo(void *container, domain_t* domain);
void mb_integrate_pref(void *container, simulation_t *simulation);
void mb_integrate_postf(void *container, simulation_t *simulation);
void mb_copyBack(void *container);
real mb_calc_E_kin(void *container);
void mb_apply_thermostat(void *container, real T_cur, real T_target);

void mb_calc_center_of_mass_velocity(void *container, real v[3]);

#endif /* MOLECULE_CONTAINER_H */
