/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "molecule_block.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "utils/logger.h"


void molecule_block_alloc(molecule_block_t **block) {
    (*block) = (molecule_block_t *) malloc(sizeof(molecule_block_t));
    (*block)->size = 0;
    (*block)->next = NULL;
    memset((*block)->flags, 0, CELL_CAPACITY * sizeof(unsigned char));
}

void molecule_block_free(molecule_block_t *block) {
	assert(block->next == NULL);
    free(block);
}

void molecule_block_clear(molecule_block_t *block) {
    memset(block->flags, 0, CELL_CAPACITY * sizeof(unsigned char));
    block->size = 0;
    block->next = NULL;
}

void molecule_block_copy(molecule_block_t * source, molecule_block_t *dest) {
	assert(dest->size == 0);
    molecule_block_t *nextBak = dest->next;
    memcpy(dest, source, sizeof(molecule_block_t));
    dest->next = nextBak;;
}

void molecule_block_move(real dr[3], molecule_block_t *block) {
    int d, i;
    for(d = 0; d < 3; d++) {
        for(i = 0; i < CELL_CAPACITY; i++) {
	    real flag = block->flags[i];
            block->r[d*CELL_CAPACITY + i] = flag * (block->r[d*CELL_CAPACITY + i] + dr[d]);
        }
    }
}

void molecule_block_reset_forces_and_momenta(molecule_block_t *block) {
    memset(block->F, 0, 3*CELL_CAPACITY * sizeof(real));
}

void molecule_block_integrate_pref(simulation_t *simulation, molecule_block_t *block) {
	int d, i;
	real a = 0.5 * simulation->dt / config.m;
	for(d = 0; d < 3; d++) {
		for(i = 0; i < CELL_CAPACITY; i++) {
			real flag = block->flags[i];
			block->v[d*CELL_CAPACITY + i] += flag * a * block->F[d*CELL_CAPACITY + i];
			block->r[d*CELL_CAPACITY + i] += flag * simulation->dt * block->v[d*CELL_CAPACITY + i];
		}
	}
}

void molecule_block_integrate_postf(simulation_t *simulation, molecule_block_t *block) {
	int d, i;
	real a = 0.5 * simulation->dt / config.m;
	for(d = 0; d < 3; d++) {
		for(i = 0; i < CELL_CAPACITY; i++) {
            real flag = block->flags[i];
			block->v[d*CELL_CAPACITY + i] += flag * a * block->F[d*CELL_CAPACITY + i];
		}
	}
}

real molecule_block_calc_E_kin(molecule_block_t *block) {
    int d, i;
    real e_kin = 0.0;
	real e_tmp[CELL_CAPACITY];
	memset(e_tmp, 0, CELL_CAPACITY * sizeof(real));
    for(d = 0; d < 3; d++) {
        for(i = 0; i < CELL_CAPACITY; i++) {
            real flag = block->flags[i];
            e_tmp[i] += flag * 0.5 * config.m * block->v[d*CELL_CAPACITY + i]*block->v[d*CELL_CAPACITY + i];
        }
    }
    for(i = 0; i < CELL_CAPACITY; i++) {
		e_kin += e_tmp[i];
	}
    return e_kin;
}

void molecule_block_print_ascii(molecule_block_t *block, FILE *fh) {
    int i;
    fprintf(fh, "# size: %lu\n", block->size);
    fprintf(fh, "# ID flag r, v, F\n");
    for(i = 0; i < CELL_CAPACITY; i++) {
        fprintf(fh, "%4.1d %4.0"PRIreal
                " [%8.4"PRIreal", %8.4"PRIreal", %8.4"PRIreal"]"
                " [%8.4"PRIreal", %8.4"PRIreal", %8.4"PRIreal"]"
                " [%8.4"PRIreal", %8.4"PRIreal", %8.4"PRIreal"]\n", 
                i, (double) block->flags[i],
                block->r[0*CELL_CAPACITY + i], block->r[1*CELL_CAPACITY + i], block->r[2*CELL_CAPACITY + i],
                block->v[0*CELL_CAPACITY + i], block->v[1*CELL_CAPACITY + i], block->v[2*CELL_CAPACITY + i],
                block->F[0*CELL_CAPACITY + i], block->F[1*CELL_CAPACITY + i], block->F[2*CELL_CAPACITY + i]
               );
    }
}


void molecule_block_print_povray(molecule_block_t *block, FILE *fh) {
	int i;
	for( i = 0; i < CELL_CAPACITY; i++) {
		if(block->flags[i] > 0) {
			fprintf(fh, "sphere { <%"PRIreal", %"PRIreal", %"PRIreal">, 0.15   pigment { color rgb <1, 0, 0> } }\n",
					block->r[0*CELL_CAPACITY + i], block->r[1*CELL_CAPACITY + i], block->r[2*CELL_CAPACITY + i]);
		}
	}
}

void molecule_block_shift_velocity(molecule_block_t *block, real v[3]) {
    int i;
    for(i = 0; i < CELL_CAPACITY; i++) {
        block->v[0*CELL_CAPACITY + i] -= v[0];
        block->v[1*CELL_CAPACITY + i] -= v[1];
        block->v[2*CELL_CAPACITY + i] -= v[2];
    }
}

void molecule_block_scale_velocity(molecule_block_t *block, real f, real T_target) {
    int i;
    if(f > VELOCITY_SCALING_FACTOR_MAX || f < VELOCITY_SCALING_FACTOR_MIN) {
        LOG(Warning, "Thermostat scaling factor outside of range [%"PRIreal",%"PRIreal"]: %"PRIreal"\n",
            VELOCITY_SCALING_FACTOR_MIN, VELOCITY_SCALING_FACTOR_MAX, f);
        real E_kin_limit = VELOCITY_SCALING_E_KIN_LIMIT_FACTOR * T_target;
        for(i = 0; i < CELL_CAPACITY; i++) {
            real E_kin = 0.5*config.m* (block->v[0*CELL_CAPACITY + i]*block->v[0*CELL_CAPACITY + i]+block->v[1*CELL_CAPACITY + i]*block->v[1*CELL_CAPACITY + i] + block->v[2*CELL_CAPACITY + i]*block->v[2*CELL_CAPACITY + i]);
            if(E_kin > E_kin_limit) {
                LOG(Warning, "E_kin: %"PRIreal"/%"PRIreal"\n", E_kin, E_kin_limit);
                real ff = sqrt(E_kin_limit/E_kin);
                block->v[0*CELL_CAPACITY + i] *= ff;
                block->v[1*CELL_CAPACITY + i] *= ff;
                block->v[2*CELL_CAPACITY + i] *= ff;
            }
        }
    }
    else {
        for(i = 0; i < CELL_CAPACITY; i++) {
            block->v[0*CELL_CAPACITY + i]*= f;
            block->v[1*CELL_CAPACITY + i]*= f;
            block->v[2*CELL_CAPACITY + i]*= f;
        }
    }
}
