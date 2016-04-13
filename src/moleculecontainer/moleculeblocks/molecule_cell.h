/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef MOLECULE_CELL_H
#define MOLECULE_CELL_H

#include "config.h"
#include "typedefs.h"

#include "molecule_block.h"
#include "molecule.h"
#include "ensemble.h"
#include "simulation.h"

#include <stdio.h>

/** molecule cell 
 *
 * A molecule cell corresponds to a cell used in the link-cell algorithm. [1]
 *
 * [1] 
 */
typedef struct molecule_cell_t molecule_cell_t;
struct molecule_cell_t {
    molecule_cell_t *neighbours[26];
    molecule_block_t *data;
};

void mcell_alloc(molecule_cell_t *cell);
void mcell_free(molecule_cell_t *cell);

void mcell_add_molecule(molecule_t *m, molecule_cell_t *cell);
static inline
size_t mcell_get_num_molecules(molecule_cell_t *cell) {
    size_t count = 0;
    molecule_block_t *mblock = cell->data;
    while(mblock != NULL) {
        count += mblock->size;
        mblock = mblock->next;
    }
    return count;
}
void mcell_print_ascii(molecule_cell_t *cell, FILE *fh);
void mcell_print_pov(molecule_cell_t *cell, FILE *fh);

void mcell_calc_forces(molecule_cell_t *cell,  real *U_pot);

void mcell_calc_intra_forces(molecule_cell_t *cell,  real *U_pot);
void mcell_calc_inter_forces(molecule_cell_t *cell1, molecule_cell_t *cell2, real *U_pot);

void mcell_clear(molecule_cell_t *cell);
void mcell_copy(molecule_cell_t *source, molecule_cell_t *dest);
void mcell_move_particles(real dr[3], molecule_cell_t *cell);
void mcell_reset_forces_and_momenta(molecule_cell_t *cell);
void mcell_shift_velocity(molecule_cell_t *cell, real v[3]);

void mcell_integrate_pref(simulation_t *simulation, molecule_cell_t *cell);
void mcell_integrate_postf(simulation_t *simulation, molecule_cell_t *cell);
real mcell_calc_E_kin(molecule_cell_t *cell);

void mcell_scale_velocity(molecule_cell_t *cell, real f, real T_target);
#endif /* MOLECULE_CELL_H */
