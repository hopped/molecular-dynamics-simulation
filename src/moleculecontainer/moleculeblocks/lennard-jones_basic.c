/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "config.h"
#include "typedefs.h"

#include "lj_parameter.h"
#include "lennard-jones_basic.h"
#include "ensemble.h"
#include "moleculeblocks.h"

#include <math.h>
#include <assert.h>

void lj_basic_block_calc(molecule_block_t *block1, molecule_block_t *block2, real *U_pot) {
	size_t i;
    real U_pot_tmp = 0.;

	for(i = 0; i < CELL_CAPACITY; i++) {
        if( block1->flags[i] == 0) 
            continue;
        size_t j;
		for(j = 0; j < CELL_CAPACITY; j++) {
            if( block2->flags[j] == 0) 
                continue;
			
			real dr[3];
			real dr2 = 0.;
            int d;
			for( d = 0; d < 3; d++ ) {
				dr[d] = block2->r[d*CELL_CAPACITY + j] - block1->r[d*CELL_CAPACITY + i];
				dr2 += dr[d] * dr[d];
			}
			if( dr2 > config.cutoff_radius_sq )
                continue;
			real invdr2 = 1. / dr2;
            config.count++;


			/* Lennard Jones interaction forces */
            // TODO: Check epsilon24 factor usage here
            // TODO: Improvements by to remove the - as lj12 is the square of lj6 ... ?
			real lj6 = sigma2 * invdr2;
			lj6 = lj6 * lj6 * lj6;
			real lj12 = lj6 * lj6;
			real lj12m6 = lj12 - lj6;
			real u_pot  = epsilon24 * lj12m6;
			real factor = epsilon24 * (lj12 + lj12m6) * invdr2;
#ifndef NDEBUG
            if(isnan(factor)) {
                fprintf(stderr, "Blocks\n");
                fprintf(stderr, "i: %lu, j: %lu, dr2: %"PRIreal", invdr2: %"PRIreal"\n", i, j, dr2, invdr2);
                molecule_block_print_ascii(block1, stderr);
                molecule_block_print_ascii(block2, stderr);
            }
#endif

			for( d = 0; d < 3; d++ ) {
				block1->F[d*CELL_CAPACITY + i] -= factor * dr[d];
			}
			U_pot_tmp += u_pot;
		}
	} 
    U_pot_tmp /= 6.0*2.0; /* we use epsilon24 instead of epsilon4 in potential calculation and do not use Newton 3 */
    *U_pot += U_pot_tmp;
}

void lj_basic_diagblock_calc(molecule_block_t *block, real *U_pot) {
	size_t i;
    real U_pot_tmp = 0.;

	for(i = 0; i < CELL_CAPACITY; i++) {
        if( block->flags[i] == 0) 
            continue;
        size_t j;
		for(j = i + 1; j < CELL_CAPACITY; j++) {
            if( block->flags[j] == 0) 
                continue;
			real dr[3];
			real dr2 = 0.;
            int d;
			for( d = 0; d < 3; d++ ) {
				dr[d] = block->r[d*CELL_CAPACITY + j] - block->r[d*CELL_CAPACITY + i];
				dr2 += dr[d] * dr[d];
			}
			if( dr2 > config.cutoff_radius_sq )
                continue;
			real invdr2 = 1. / dr2;
            config.count++;

			/* Lennard Jones interaction forces */
			real lj6 = sigma2 * invdr2;
			lj6 = lj6 * lj6 * lj6;
			real lj12 = lj6 * lj6;
			real lj12m6 = lj12 - lj6;
			real u_pot = epsilon24 * lj12m6;
			real factor = epsilon24 * (lj12 + lj12m6) * invdr2;
#ifndef NDEBUG
            if(isnan(factor)) {
                fprintf(stderr, "Blocks\n");
                fprintf(stderr, "i: %lu, j: %lu, dr2: %"PRIreal", invdr2: %"PRIreal"\n", i, j, dr2, invdr2);
                molecule_block_print_ascii(block, stderr);
            }
#endif

			for( d = 0; d < 3; d++ ) {
				block->F[d*CELL_CAPACITY + i] -= factor * dr[d];
				block->F[d*CELL_CAPACITY + j] += factor * dr[d];
				//block->F[d*CELL_CAPACITY + i] -= block->F[d*CELL_CAPACITY + i];
			}
			U_pot_tmp += u_pot;
		}
	}
    U_pot_tmp /= 6.0; /* we use epsilon24 instead of epsilon4 in potential calculation */
    *U_pot += U_pot_tmp;
}


