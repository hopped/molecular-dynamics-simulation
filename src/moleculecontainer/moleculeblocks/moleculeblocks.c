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
#include "constants.h"

#include "moleculeblocks.h"
#include "molecule_cell.h"

#include "ensemble.h"
#include "domain.h"

#include "utils/logger.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

static size_t get_cell_id(long ix, long iy, long iz, mb_t *mc) {
    size_t index;
    if (ix < 0)
        ix += mc->num_cells_per_dim[0];
    else if (ix >= mc->num_cells_per_dim[0])
        ix -= mc->num_cells_per_dim[0];
    if (iy < 0)
        iy += mc->num_cells_per_dim[1];
    else if (iy >= mc->num_cells_per_dim[1])
        iy -= mc->num_cells_per_dim[1];
    if (iz < 0)
        iz += mc->num_cells_per_dim[2];
    else if (iz >= mc->num_cells_per_dim[2])
        iz -= mc->num_cells_per_dim[2];

    index = ( mc->num_cells_per_dim[1] * iz + iy) * mc->num_cells_per_dim[0] + ix;
    return index;
}

static void get_cell_coords_from_coordinate(real x, real y, real z, mb_t *mc, long coords[3]) {
    coords[0] = x / mc->cell_size[0] + 1;
    coords[1] = y / mc->cell_size[1] + 1;
    coords[2] = z / mc->cell_size[2] + 1;
}

static void mb_get_cell_coords(mb_t *mc, size_t cell_id, long coords[3]) {
    coords[0] = cell_id % mc->num_cells_per_dim[0];
    cell_id /= mc->num_cells_per_dim[0];
    coords[1] = cell_id % mc->num_cells_per_dim[1];
    cell_id /= mc->num_cells_per_dim[1];
    coords[2] = cell_id % mc->num_cells_per_dim[2];
    assert(cell_id < mc->num_cells_per_dim[2]);
}

static size_t get_cell_index_from_coordinate(real x, real y, real z, mb_t *mc) {
    size_t index;
    long coords[3];
    get_cell_coords_from_coordinate(x, y, z, mc, coords);
    index = get_cell_id(coords[0], coords[1], coords[2], mc);
    return index;
}

void mb_estimate_halos( void * container, domain_t * domain){}

void* mb_alloc() {
    mb_t *mc;
    mc = (mb_t *) malloc(sizeof(mb_t));
    return (mc_t*) mc;
}

void mb_init(void *container, ensemble_t *ensemble, domain_t * domain) {
    mb_t *mc = (mb_t *) container;
    mc->num_molecules = 0;
    mc->num_cells = 1;
    mc->cutoff_radius = config.cutoff_radius;
    assert(mc->cutoff_radius > 0);
    int d;
    for(d = 0; d < 3; d++) {
        mc->num_cells_per_dim[d] = floor(domain->L[d] / mc->cutoff_radius) + 2; /* 2 more cells in each direction to handle boundary conditinos and halo */
        assert(mc->num_cells_per_dim[d] > 2);
        mc->num_cells *= mc->num_cells_per_dim[d];
        mc->cell_size[d] = domain->L[d] / (mc->num_cells_per_dim[d] - 2);
        assert(mc->cell_size[d] <= domain->L[d]);
    }
    assert(mc->num_cells>0);
    mc->cells = (molecule_cell_t *) malloc(mc->num_cells * sizeof(molecule_cell_t));
    long i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                mcell_alloc(cell);
                if (i != 0 && j != 0 && k != 0 &&
                    i != mc->num_cells_per_dim[0] - 1 && j != mc->num_cells_per_dim[1] - 1 && k != mc->num_cells_per_dim[2] - 1){
		    int nId = 0;
		    int di, dj, dk;
		    for(dk = -1; dk < 2; dk++) {
		    	for(dj = -1; dj < 2; dj++) {
		    		for(di = -1 ; di < 2; di++) {
		    			size_t neighbour_cell_id = get_cell_id( i + di , j + dj , k + dk, mc);
		    			if(neighbour_cell_id == cell_id)
		    				continue;
		    			cell->neighbours[nId++] = &mc->cells[neighbour_cell_id];
// 		    			printf("cell id %lu, nId %d, neighbour id%lu\n", cell_id, nId, neighbour_cell_id);
		    		}
		    	}
		    }
                }
            }
        }
    }
}

void mb_update(void *container, domain_t *domain) {
    mb_t *mc = (mb_t *) container;
    long i, j, k;
#pragma omp parallel for default(shared) private(i,j,k) collapse(3)
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;

                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                molecule_block_t *mblock;

                mblock = cell->data;
                while (mblock != NULL) {
                    int w;
                    for( w = 0; w < CELL_CAPACITY; w++) {
                        if(mblock->flags[w] != 0) {
                            /* apply periodic boundary condtions */
                            int d;
                            for(d = 0; d < 3; d ++) {
                               if(mblock->r[d*CELL_CAPACITY + w] < 0) mblock->r[d*CELL_CAPACITY + w] += domain->L[d];
                               else if (mblock->r[d*CELL_CAPACITY + w] >= domain->L[d]) mblock->r[d*CELL_CAPACITY + w] -= domain->L[d];
                            }
                            size_t dest_cell_id = get_cell_index_from_coordinate(mblock->r[0*CELL_CAPACITY + w], mblock->r[1*CELL_CAPACITY + w], mblock->r[2*CELL_CAPACITY + w], mc);
                            assert(dest_cell_id < mc->num_cells);
                            if( cell_id != dest_cell_id ) {
                                /* safe molecule data to struct */
                                molecule_t m;
                                m.id   = mblock->id[w];
                                m.r[0] = mblock->r[0*CELL_CAPACITY + w];
                                m.r[1] = mblock->r[1*CELL_CAPACITY + w];
                                m.r[2] = mblock->r[2*CELL_CAPACITY + w];
                                m.v[0] = mblock->v[0*CELL_CAPACITY + w];
                                m.v[1] = mblock->v[1*CELL_CAPACITY + w];
                                m.v[2] = mblock->v[2*CELL_CAPACITY + w];
                                m.F[0] = mblock->F[0*CELL_CAPACITY + w];
                                m.F[1] = mblock->F[1*CELL_CAPACITY + w];
                                m.F[2] = mblock->F[2*CELL_CAPACITY + w];
                                /* delete molecule */
                                mblock->flags[w] = 0;
                                mblock->size--;
                                /* calculate target cell and add molecule */    
                                //DEBUG("Molecule coords [dest]: %"PRIreal", %"PRIreal", %"PRIreal"\n", m.r[0], m.r[1], m.r[2]);
                                //DEBUG("Domain size: %"PRIreal", %"PRIreal", %"PRIreal"\n", mc->cell_size[0], mc->cell_size[1], mc->cell_size[2]);
                                molecule_cell_t *dest_cell = &mc->cells[dest_cell_id];
                                #pragma omp critical 
                                mcell_add_molecule(&m, dest_cell);
                                /* debug output */
                                long coords_src[3], coords_dest[3];
                                mb_get_cell_coords(mc, cell_id, coords_src);
                                mb_get_cell_coords(mc, dest_cell_id, coords_dest);
                                DEBUG("Moving Molecule %lu from cell %lu (%ld,%ld,%ld) to cell %lu (%ld,%ld,%ld)\n", m.id, cell_id, coords_src[0], coords_src[1], coords_src[2], dest_cell_id, coords_dest[0], coords_dest[1], coords_dest[2]);
                            }
                        }
                    }
                    mblock = mblock->next;
                }

            }
        }
    }
}

void mb_free(void *container){
    mb_t *mc = (mb_t *) container;
    size_t i;
    for(i = 0; i < mc->num_cells; i++) {
        mcell_free(&mc->cells[i]);
    }
    free(mc->cells);
    free(mc);
}

void mb_clear_halo(void *container) {
    mb_t *mc = (mb_t *) container;
    long i, j, k;
    for(k = 0; k < mc->num_cells_per_dim[2]; k++) {
        for(j = 0; j < mc->num_cells_per_dim[1]; j++) {
            for(i = 0; i < mc->num_cells_per_dim[0]; i++) {
                if( i == 0 || i == mc->num_cells_per_dim[0] - 1 || 
                        j == 0 || j == mc->num_cells_per_dim[1] - 1 || 
                        k == 0 || k == mc->num_cells_per_dim[2] - 1 ) 
                {
                    size_t cell_id = get_cell_id(i, j, k, mc);
                    mcell_clear(&mc->cells[cell_id]);
                }
            }
        }
    }
}


void mb_print_ascii(void *container, char *filename){
    int print_halo=1;
    mb_t *mc = (mb_t *) container;
    FILE *fh;
    fh = fopen(filename, "w+");

    long i, j, k;
    if( print_halo == 0 ) {
        for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
            for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
                for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
                    size_t cell_id = get_cell_id( i, j, k ,mc);
                    molecule_cell_t *cell;
                    cell = &mc->cells[cell_id];
                        printf( "Cell (%ld,%ld,%ld)\n", i, j, k);
                    mcell_print_ascii(cell, fh);
                }
            }
        }
    }
    else {
        for(k = 0; k < mc->num_cells_per_dim[2]; k++) {
            for(j = 0; j < mc->num_cells_per_dim[1]; j++) {
                for(i = 0; i < mc->num_cells_per_dim[0]; i++) {
                    size_t cell_id = get_cell_id( i, j, k ,mc);
                    molecule_cell_t *cell;
                    cell = &mc->cells[cell_id];
                        printf( "Cell (%ld,%ld,%ld)\n", i, j, k);
                    mcell_print_ascii(cell, fh);
                }
            }
        }
    }
    fclose(fh);
}

void mb_print_pov(void *container, char *filename){
    mb_t *mc = (mb_t *) container;
    FILE *fh;
    fh = fopen(filename, "w+");

    size_t i;
    fprintf(fh, "camera {\n location <-2, -2, -2>\n look_at <0, 0, 0>\n}\n");
    fprintf(fh, "light_source {\n <2,-5, 0>,  color rgb <1, 1, 1>\n shadowless\n}\n");
    for( i = 0; i < mc->num_cells; i++) {
        molecule_cell_t *cell;
        cell = &mc->cells[i];
        mcell_print_pov(cell, fh);
    }
    fclose(fh);
}

void mb_save_to_file(char *filename, void *container){
}
void mb_read_from_file(char *filename, void *container){
}

void mb_add_molecule(void *container, molecule_t *m){
    mb_t *mc = (mb_t *) container;
    size_t cell_id = get_cell_index_from_coordinate(m->r[0], m->r[1], m->r[2], mc);
    //long coords[3];
    //get_cell_coords_from_coordinate(m->r[0], m->r[1], m->r[2], mc, coords);
    //DEBUG("ADDING MOLECULE TO CELL %ld, %ld, %ld\n", coords[0], coords[1], coords[2]);
    mcell_add_molecule(m, &mc->cells[cell_id]);
    mc->num_molecules++;
}

long mb_get_num_molecules(void *container) {
    mb_t *mc = (mb_t *) container;
    size_t count = 0;
    size_t i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                size_t n = mcell_get_num_molecules(cell);
                count += n;
            }
        }
    }
    return count;
}

static void mb_update_corners(mb_t *mc, domain_t *domain) {
    size_t corner_coords[8][3] = {
        { 0, 0, 0 },
        { 0, 0, mc->num_cells_per_dim[2] - 1},
        { 0, mc->num_cells_per_dim[1] - 1, 0 },
        { 0, mc->num_cells_per_dim[1] - 1, mc->num_cells_per_dim[2] - 1 },
        { mc->num_cells_per_dim[0] - 1, 0, 0 },
        { mc->num_cells_per_dim[0] - 1, 0, mc->num_cells_per_dim[2] - 1 },
        { mc->num_cells_per_dim[0] - 1, mc->num_cells_per_dim[1] - 1, 0 },
        { mc->num_cells_per_dim[0] - 1, mc->num_cells_per_dim[1] - 1, mc->num_cells_per_dim[2] - 1 }
    };
    int offsets[8][3] = {
        { -1, -1, -1 },
        { -1, -1,  1 },
        { -1,  1, -1 },
        { -1,  1,  1 },
        {  1, -1, -1 },
        {  1, -1,  1 },
        {  1,  1, -1 },
        {  1,  1,  1 }
    };
    /* l: lower, h: higher */
    int i;
    for( i = 0; i < 8; i++ ) {
        /* 2*offset because we have a halo on the oposite side, too. */
        size_t source_id =  get_cell_id( corner_coords[i][0] + 2*offsets[i][0], corner_coords[i][1] + 2*offsets[i][1], corner_coords[i][2] + 2*offsets[i][2], mc);
        size_t dest_id   =  get_cell_id( corner_coords[i][0], corner_coords[i][1], corner_coords[i][2], mc);
        mcell_copy(&mc->cells[source_id], &mc->cells[dest_id]);
        int d;
        real dr[3];
        for( d = 0; d < 3; d++)
            dr[d] = offsets[i][d] * domain->L[d];
        mcell_move_particles(dr, &mc->cells[dest_id]);
        //DEBUG("CORNER %d\n", i);
        //DEBUG("dest: [%ld, %ld, %ld]\n", corner_coords[i][0], corner_coords[i][1], corner_coords[i][2]); 
        //DEBUG("src:  [%ld, %ld, %ld]\n", corner_coords[i][0] + 2*offsets[i][0], corner_coords[i][1] + 2*offsets[i][1], corner_coords[i][2] + 2*offsets[i][2]);
    }
}

static void mb_update_edges(mb_t *mc, domain_t *domain) {
    long offsets[3][4][3] = {
        {
            {0, -1, -1},
            {0, -1,  1},
            {0,  1, -1},
            {0,  1,  1}
        },
        {
            {-1, 0, -1},
            {-1, 0,  1},
            { 1, 0, -1},
            { 1, 0,  1}
        },
        {
            {-1, -1, 0},
            {-1,  1, 0},
            { 1, -1, 0},
            { 1,  1, 0}
        }
    };
    long edges[3][4][3] = {
        { /* x direction */
            {0, 0, 0},
            {0, 0, mc->num_cells_per_dim[2] - 1},
            {0, mc->num_cells_per_dim[1] - 1, 0},
            {0, mc->num_cells_per_dim[1] - 1, mc->num_cells_per_dim[2] - 1}
        },
        { /* y direction */
            {0, 0, 0},
            {0, 0, mc->num_cells_per_dim[2] - 1},
            {mc->num_cells_per_dim[0] - 1, 0, 0},
            {mc->num_cells_per_dim[0] - 1, 0, mc->num_cells_per_dim[2] - 1}
        },
        { /* z direction */
            {0, 0, 0},
            {0, mc->num_cells_per_dim[1] - 1, 0},
            {mc->num_cells_per_dim[0] - 1, 0, 0},
            {mc->num_cells_per_dim[0] - 1, mc->num_cells_per_dim[1] - 1, 0}
        }
    }; /* fixed coordinates for the halo edges in this direction */    

    int k;
    int direction;
    for(direction = 0; direction < 3; direction++) {
        for(k = 0; k < 4; k++) {
            long i;
            long coords[3];
            coords[0] = edges[direction][k][0];
            coords[1] = edges[direction][k][1];
            coords[2] = edges[direction][k][2];
            real dr[3] = {offsets[direction][k][0] * domain->L[0], offsets[direction][k][1] * domain->L[1], offsets[direction][k][2] * domain->L[2]};
            dr[direction] = 0;

            for(i = 1; i < mc->num_cells_per_dim[direction] - 1; i++) {
                coords[direction] = i;
                size_t dest_id   =  get_cell_id(coords[0], coords[1], coords[2], mc);
                size_t source_id =  get_cell_id(coords[0] + 2*offsets[direction][k][0], coords[1] + 2*offsets[direction][k][1], coords[2] + 2*offsets[direction][k][2], mc);
                //DEBUG("EDGE direction=%d, i=%d\n", direction, i);
                //DEBUG("dest: [%ld, %ld, %ld]\n", coords[0], coords[1], coords[2]); 
                //DEBUG("src:  [%ld, %ld, %ld]\n", coords[0] + 2*offsets[direction][k][0], coords[1] + 2*offsets[direction][k][1], coords[2] + 2*offsets[direction][k][2]);
                mcell_copy(&mc->cells[source_id], &mc->cells[dest_id]);
                mcell_move_particles(dr, &mc->cells[dest_id]);
            }
        }
    }
}

static void mb_update_surfaces(mb_t *mc, domain_t * domain) {
    int direction;
    int offsets[2] = { -1, 1 };

    /* xy surfaces */
    size_t z_planes[2] = {0, mc->num_cells_per_dim[2] - 1}; /* halo plane to fill */
    for(direction = 0; direction < 2; direction++) {
        size_t iz = z_planes[direction];
        size_t ix;
        real dr[3] = {0.0, 0.0, 0.0};
        dr[2] = offsets[direction] * domain->L[2];
        for(ix = 1; ix < mc->num_cells_per_dim[0] - 1; ix++) {
            size_t iy;
            for(iy = 1; iy < mc->num_cells_per_dim[1] - 1; iy++) {
                size_t dest_direction   =  get_cell_id(ix, iy, iz, mc);
                size_t source_direction =  get_cell_id(ix, iy, iz + 2*offsets[direction], mc);
                //DEBUG("SURFACE direction=%d, ix=%d, iy=%d\n", direction, ix, iy);
                //DEBUG("dest: [%ld, %ld, %ld]\n", ix, iy, iz); 
                //DEBUG("src:  [%ld, %ld, %ld]\n", ix, iy, iz + 2*offsets[direction]);
                mcell_copy(&mc->cells[source_direction], &mc->cells[dest_direction]);
                mcell_move_particles(dr, &mc->cells[dest_direction]);
            }
        }
    }
    /* xz surfaces */
    size_t y_planes[2] = {0, mc->num_cells_per_dim[1] - 1}; /* halo plane to fill */
    for(direction = 0; direction < 2; direction++) {
        size_t iy = y_planes[direction];
        size_t ix;
        real dr[3] = {0.0, 0.0, 0.0};
        dr[1] = offsets[direction] * domain->L[1];
        for(ix = 1; ix < mc->num_cells_per_dim[0] - 1; ix++) {
            size_t iz;
            for(iz = 1; iz < mc->num_cells_per_dim[2] - 1; iz++) {
                size_t dest_direction   =  get_cell_id(ix, iy, iz, mc);
                size_t source_direction =  get_cell_id(ix, iy + 2*offsets[direction], iz, mc);
                mcell_copy(&mc->cells[source_direction], &mc->cells[dest_direction]);
                mcell_move_particles(dr, &mc->cells[dest_direction]);
            }
        }
    }
    /* yz surfaces */
    size_t x_planes[2] = {0, mc->num_cells_per_dim[0] - 1}; /* halo plane to fill */
    for(direction = 0; direction < 2; direction++) {
        size_t ix = x_planes[direction];
        size_t iy;
        real dr[3] = {0.0, 0.0, 0.0};
        dr[0] = offsets[direction] * domain->L[0];
        for(iy = 1; iy < mc->num_cells_per_dim[1] - 1; iy++) {
            size_t iz;
            for(iz = 1; iz < mc->num_cells_per_dim[2] - 1; iz++) {
                size_t dest_direction   =  get_cell_id(ix, iy, iz, mc);
                size_t source_direction =  get_cell_id(ix + 2*offsets[direction], iy, iz, mc);
                //DEBUG("SURFACE direction=%d, ix=%d, iy=%d\n", direction, ix, iy);
                //DEBUG("dest: [%ld, %ld, %ld]\n", ix, iy, iz); 
                //DEBUG("src:  [%ld, %ld, %ld]\n", ix + 2*offsets[direction], iy, iz);
                mcell_copy(&mc->cells[source_direction], &mc->cells[dest_direction]);
                mcell_move_particles(dr, &mc->cells[dest_direction]);
            }
        }
    }

}

void mb_update_halo(void *container, domain_t* domain) {
    mb_t *mc = (mb_t *) container;
    /* copy cells to boundary to fullfill the periodic boundary condition. */
    /* we have 26 directions */

    /* 8 corners */ 
    mb_update_corners(mc, domain);
    /* 12 edges */
    mb_update_edges(mc, domain);
    /* 6 lateral surfaces */
    mb_update_surfaces(mc, domain);
}


void mb_calc_forces(void *container, real *U_pot) {
    mb_t *mc = (mb_t *) container;
    *U_pot = 0.0;
    real U_pot_local=0.0;
    size_t i, j, k;
#pragma omp parallel for default(shared) private(i,j,k) collapse(3) reduction(+:U_pot_local)
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                size_t n = mcell_get_num_molecules(cell);
                if(n > 0) {
                    mcell_calc_forces(cell, &U_pot_local);
                }
            }
        }
    }
    *U_pot=U_pot_local;
    DEBUG("Number of computed interactions: %lu\n", config.count);
}

void mb_print_stats(void *container, ensemble_t *ensemble) {
    mb_t *mc = (mb_t *) container;
    size_t min_molecules_in_cell = ensemble->N;
    size_t max_molecules_in_cell = 0;

    size_t i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                size_t n = mcell_get_num_molecules(cell);
                if(n > max_molecules_in_cell) 
                    max_molecules_in_cell = n;
                else if (n < min_molecules_in_cell) 
                    min_molecules_in_cell = n;
                //if(n == 0)
                //    DEBUG("Warning: Cell (%lu, %lu, %lu) has 0 molecules!\n", i, j ,k);
            }
        }
    }
    long num_real_cells = 1;
    int d;
    for(d = 0; d < 3; d++) {
        num_real_cells *= (mc->num_cells_per_dim[d] - 2);
    }
    LOG(Info, "Domain size: %"PRIreal", %"PRIreal", %"PRIreal"\n", domain.L[0], domain.L[1], domain.L[2]);
    LOG(Info, "Cells (# cells, #cells per dim(x,y,z), cell size(x,y,z)): %lu, (%lu, %lu, %lu), (%"PRIreal", %"PRIreal", %"PRIreal")\n", 
            mc->num_cells, mc->num_cells_per_dim[0], mc->num_cells_per_dim[1], mc->num_cells_per_dim[2],
            mc->cell_size[0], mc->cell_size[1], mc->cell_size[2]);
    LOG(Info, "Molecules per cell (min, max, avg): %lu, %lu, %lf\n", 
            min_molecules_in_cell, max_molecules_in_cell, (double) ensemble->N / (double) num_real_cells);
}

void mb_reset_forces_and_momenta(void *container) {
    mb_t *mc = (mb_t *) container;
    size_t i, j, k;
    for(i = 0; i < mc->num_cells_per_dim[0] - 0; i++) {
        for(j = 0; j < mc->num_cells_per_dim[1] - 0; j++) {
            for(k = 0; k < mc->num_cells_per_dim[2] - 0; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                mcell_reset_forces_and_momenta(cell);
            }
        }
    }
}

void mb_integrate_pref(void *container, simulation_t *simulation) {
    mb_t *mc = (mb_t *) container;
    size_t i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                mcell_integrate_pref(simulation, cell);
            }
        }
    }
}

void mb_integrate_postf(void *container, simulation_t *simulation) {
    mb_t *mc = (mb_t *) container;
    size_t i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                mcell_integrate_postf(simulation, cell);
            }
        }
    }
}

real mb_calc_E_kin(void *container) {
    real E_kin = 0.0;
    mb_t *mc = (mb_t *) container;
    size_t i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                E_kin += mcell_calc_E_kin(cell);
            }
        }
    }
    return E_kin;
}

void mb_shift_velocity(void *container, real v[3]) {
    mb_t *mc = (mb_t *) container;
    size_t i, j, k;
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                mcell_shift_velocity(cell, v);
            }
        }
    }
}

void mb_apply_thermostat(void *container, real T_cur, real T_target) {
    assert(container != NULL);
    assert(T_cur > 0.0);
    assert(T_target > 0.0);
    
    mb_t *mc = (mb_t *) container;
    
    size_t i, j, k;
    real f = sqrt(T_target / T_cur);
    LOG(Warning, "Thermorstat scaling factor: f = %"PRIreal"\n", f);
    
    for(i = 1; i < mc->num_cells_per_dim[0] - 1; i++) {
        for(j = 1; j < mc->num_cells_per_dim[1] - 1; j++) {
            for(k = 1; k < mc->num_cells_per_dim[2] - 1; k++) {
                size_t cell_id;
                cell_id = get_cell_id( i, j, k ,mc);
                molecule_cell_t * cell = &mc->cells[cell_id];
                mcell_scale_velocity(cell, f, T_target);
            }
        }
    }
}
void mb_decompose(void *container, domain_t *domain){
}

void mb_destroy(void * container){
}
void mb_copyBack(void *container) { 
}


