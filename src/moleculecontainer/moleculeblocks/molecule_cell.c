/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#include "molecule_cell.h"
#include "molecule_block.h"
#include "lennard-jones_basic.h"

#include <stdlib.h>

#include <assert.h>

/** @todo init */
void mcell_alloc(molecule_cell_t *cell) {
    molecule_block_alloc(&cell->data);
    assert(cell->data != NULL);
}
void mcell_free(molecule_cell_t *cell) {
    mcell_clear(cell);
    molecule_block_free(cell->data);
}

void mcell_clear(molecule_cell_t *cell) {
    
    molecule_block_t *mblock = cell->data->next;
    molecule_block_t *next_mblock;

    /* delete additional blocks */
    while(mblock != NULL) {
        next_mblock = mblock->next;
        mblock->next=NULL;
        molecule_block_free(mblock);
        mblock = next_mblock;
    }
    /* clear original block */
    molecule_block_clear(cell->data);
}

void mcell_add_molecule(molecule_t *m, molecule_cell_t *cell) {
    molecule_block_t *mblock;
    mblock = cell->data;
    /* TODO: Lock cell. (Locking on insert position may cause cache line conflicts.) */
    /* TODO: We might remember the cell if there is room to insert another molecule */
    while (mblock->size == CELL_CAPACITY) {
        if (mblock->next == NULL)
            molecule_block_alloc(&mblock->next);
        mblock = mblock->next;
    } 
    int i;
    for(i = 0; i < CELL_CAPACITY; i++) {
		/** @todo may go into molecule_block */
        if(mblock->flags[i] == 0) {
            mblock->flags[i] = 1;
            mblock->size++;
	    mblock->id[i] = m->id;
	    mblock->r[0*CELL_CAPACITY + i] = m->r[0];
	    mblock->r[1*CELL_CAPACITY + i] = m->r[1];
	    mblock->r[2*CELL_CAPACITY + i] = m->r[2];
	    mblock->v[0*CELL_CAPACITY + i] = m->v[0];
	    mblock->v[1*CELL_CAPACITY + i] = m->v[1];
	    mblock->v[2*CELL_CAPACITY + i] = m->v[2];
	    /** @todo Do we need to copy forces here?
	    mblock->F[0*CELL_CAPACITY + i] = m->F[0];
	    mblock->F[1*CELL_CAPACITY + i] = m->F[1];
	    mblock->F[2*CELL_CAPACITY + i] = m->F[2];
	    */
            return;
        }
    }
}

void mcell_print_ascii(molecule_cell_t *cell, FILE *fh) {
    molecule_block_t *mblock = cell->data;
    while (mblock != NULL) {
		molecule_block_print_ascii(mblock, fh);
        mblock = mblock->next;
    }
}

void mcell_print_pov(molecule_cell_t *cell, FILE *fh) {
    molecule_block_t *mblock = cell->data;
    while (mblock != NULL) {
		molecule_block_print_povray(mblock, fh);
        mblock = mblock->next;
    }
}

void mcell_calc_forces(molecule_cell_t *cell, real *U_pot) {
	mcell_calc_intra_forces(cell, U_pot);
	molecule_cell_t *neighbour_cell;
	int i;
	for (i = 0; i < 26; i++) {
		neighbour_cell = cell->neighbours[i];
        assert(neighbour_cell != NULL);
		mcell_calc_inter_forces(cell, neighbour_cell, U_pot);
	}
}

void mcell_calc_intra_forces(molecule_cell_t *cell, real *U_pot) {
	molecule_block_t *block1;
        molecule_block_t *start;
	block1 = cell->data;
	start = cell->data;
        while(block1 != NULL) {
           lj_basic_diagblock_calc(block1, U_pot);
           molecule_block_t *block2;
           block2 = start;
  	   while(block2 != NULL) {
              if(block1 != block2)
                lj_basic_block_calc(block1, block2, U_pot);
   	      block2 = block2->next;
	   }
	   block1 = block1->next;
	}
}

void mcell_calc_inter_forces(molecule_cell_t *cell1, molecule_cell_t *cell2, real *U_pot) {
	molecule_block_t *block1;
	block1 = cell1->data;
        while(block1 != NULL) {
           molecule_block_t *block2;
           block2 = cell2->data;
   	   while(block2 != NULL) {
             lj_basic_block_calc(block1, block2, U_pot);
	     block2 = block2->next;
	   }   
	   block1 = block1->next;
        }
}

void mcell_copy(molecule_cell_t *source, molecule_cell_t *dest) {
    molecule_block_t *block_source = source->data;
    molecule_block_t *block_dest;
	if(dest->data == NULL) {
		mcell_alloc(dest);
	}
    assert(dest != NULL);
    assert(dest->data != NULL);
    mcell_clear(dest);
    assert(mcell_get_num_molecules(dest) == 0);
    block_dest = dest->data;

    while(block_source != NULL) {
        molecule_block_copy(block_source, block_dest);
		if(block_source->next != NULL) {
			molecule_block_alloc(&block_dest->next);
		}
        block_source = block_source->next;
        block_dest = block_dest->next;
    }
}

void mcell_move_particles(real dr[3], molecule_cell_t *cell) {
    molecule_block_t *block;
    block = cell->data;
    while(block != NULL) {
        molecule_block_move(dr, block);
        block = block->next;
    }
}

void mcell_reset_forces_and_momenta(molecule_cell_t *cell) {
    molecule_block_t *block;
    block = cell->data;
    while(block != NULL) {
        molecule_block_reset_forces_and_momenta(block);
        block = block->next;
    }
}

void mcell_integrate_pref(simulation_t *simulation, molecule_cell_t *cell) {
    molecule_block_t *block;
    block = cell->data;
    while(block != NULL) {
        molecule_block_integrate_pref(simulation, block);
        block = block->next;
    }
}

void mcell_integrate_postf(simulation_t *simulation, molecule_cell_t *cell) {
    molecule_block_t *block;
    block = cell->data;
    while(block != NULL) {
        molecule_block_integrate_postf(simulation, block);
        block = block->next;
    }
}
real mcell_calc_E_kin(molecule_cell_t *cell) {
    molecule_block_t *block;
    real E_kin = 0.0;
    block = cell->data;
    while(block != NULL) {
        E_kin += molecule_block_calc_E_kin(block);
        block = block->next;
    }
    return E_kin;
}

void mcell_shift_velocity(molecule_cell_t *cell, real v[3]) {
    molecule_block_t *block;
    block = cell->data;
    while(block != NULL) {
        molecule_block_shift_velocity(block, v);
        block = block->next;
    }
}

void mcell_scale_velocity(molecule_cell_t *cell, real f, real T_target) {
    molecule_block_t *block;
    block = cell->data;
    while(block != NULL) {
        molecule_block_scale_velocity(block, f, T_target);
        block = block->next;
    }
}
