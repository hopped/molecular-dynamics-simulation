/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

/**
 * @file
 *
 * Molecule container basicN2.
 */

#ifndef BASICN2_H
#define BASICN2_H

#include "config.h"
#include "typedefs.h"

#include "ensemble.h"
#include "../moleculecontainer.h"
#include "molecule.h"
#include "domain.h"

/**
 * Basic molecule container iterating over all $N^2$ molecule pairs.
 * Apply Lennard Jones potential only to molecules pairs with distance
 * less or equal to the cutoff radius.
 */

#ifdef BASICN2_C
int Sx, Sy, Sz, NSD;
#else
extern int Sx, Sy, Sz, NSD;
#endif

#ifdef MPI
typedef struct mpi_buffer_t{
  long id;
  real r[3];
  real v[3];
  int  t[4];
}mpi_buffer_t;
#endif

typedef struct basicN2 {
#ifdef MPI
    int *t; // travel flag in X,Y,Z direction (+1 or -1) ==> from one psp to next one
    mpi_buffer_t * send_buffer;
    mpi_buffer_t * recv_buffer;
    long Max_transfer_sz; // Maximum transfer size
    long Max_receive_sz; // Maximum receive size
    long TXs[2]; // Count of send molecules in  X_Positive and negatives direction
    long TXr[2]; // Count of reveive molecules in  X_Positive and negatives direction
    long TYs[2]; // Count of send molecules in Y_Positive and negatives direction
    long TYr[2]; // Count of reveive molecules in  Y_Positive and negatives direction
    long TZs[2]; // Count of send molecules in  Z_Positive and negatives direction
    long TZr[2]; // Count of reveive molecules in Z_Positive and negatives direction
    int Nx[2];   // Neighbours in X direction
    int Ny[2];   // Neighbours in Y direction
    int Nz[2];   // Neighbours in Z direction
    int myrank;
#endif
    long capacity; /**< capacity of container */
    long N; /**< number of molecules in container */
    real *r; /**< positions */
    real *F; /**< forces */
    real *v; /**< velocities */
    long *id; /**< molecule ids */
    real cutoff_radius; /**< Cutoff radius for the Lennard Jones potential. */
    struct basicN2 *halo; /**< Molecule container used to store halo moleucles. */
} basicN2;

#define X(mc, member, i) (mc->member[i])
#define Y(mc, member, i) (mc->member[i+mc->capacity])
#define Z(mc, member, i) (mc->member[i+2*mc->capacity])
#define ID(mc, i) (mc->id[i])

void* basicN2_alloc();
void basicN2_init(void *mc, ensemble_t *ensemble, domain_t * domain);
void basicN2_free(void *mc);
void basicN2_destroy(void * container);
void basicN2_decompose(void * container, domain_t *domain);
void basicN2_estimate_halos( void * container, domain_t * domain);

void basicN2_add_molecule(void *mc, molecule_t *m);
void basicN2_del_molecule_by_id(void *mc, long id);

void basicN2_reset_forces_and_momenta(void *mc);
void basicN2_calc_forces(void *mc, real *U_pot);
void basicN2_integrate_pref(void *mc, simulation_t *simulation);
void basicN2_integrate_postf(void *mc, simulation_t *simulation);
void basicN2_copyBack(void *mc);
void basicN2_update(void *mc, domain_t *domain);
void basicN2_update_halo(void *mc, domain_t* domain);

void basicN2_print_stats(void *mc, ensemble_t *ensemble);
void basicN2_print_ascii(void *mc, char *filename);
void basicN2_print_ascii_xyz(void *mc, char *filename);
void basicN2_print_pov(void *mc, char *filename);

real basicN2_calc_E_kin(void *mc);
void basicN2_apply_thermostat(void *mc, real T_cur, real T_target);

long basicN2_get_num_molecules(void *mc);
void basicN2_clear_halo(void *mc);

void basicN2_shift_velocity(void *container, real v[3]);


#endif /* BASICN2_H */
