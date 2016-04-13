/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef MOLECULE_BLOCK_H
#define MOLECULE_BLOCK_H

#include "config.h"
#include "typedefs.h"
#include "ensemble.h"
#include "simulation.h"

#include <stdio.h>

typedef struct molecule_block_t molecule_block_t;
struct molecule_block_t {
    size_t size;
    molecule_block_t *next;
    unsigned char flags[CELL_CAPACITY];  /**< flag (empty/deleted = 0) */
    
    /*** should be aligned here ***/

    real   r[3*CELL_CAPACITY];  /**< positions */
    real   F[3*CELL_CAPACITY];  /**< forces */
    real   v[3*CELL_CAPACITY];  /**< velocities */
    size_t id[CELL_CAPACITY];   /**< molecule ids */
};

void molecule_block_alloc(molecule_block_t **block);
void molecule_block_free(molecule_block_t *block);

void molecule_block_clear(molecule_block_t *block);
void molecule_block_move(real dr[3], molecule_block_t *block);
void molecule_block_copy(molecule_block_t * source, molecule_block_t *dest);

void molecule_block_reset_forces_and_momenta(molecule_block_t *block);
void molecule_block_shift_velocity(molecule_block_t *block, real v[3]);

/** first step of the verlet integration scheme. */
void molecule_block_integrate_pref(simulation_t *simulation, molecule_block_t *block);
/** second step of the verlet integration scheme. */
void molecule_block_integrate_postf(simulation_t *simulation, molecule_block_t *block);
real molecule_block_calc_E_kin(molecule_block_t *block);

void molecule_block_print_ascii(molecule_block_t *block, FILE *fh);
void molecule_block_print_povray(molecule_block_t *block, FILE *fh);

void molecule_block_scale_velocity(molecule_block_t *block, real f, real T_target);
#endif /* MOLECULE_BLOCK_H */
