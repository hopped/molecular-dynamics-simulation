/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 * Copyright (c) 2012-2014 Mohamad Amer Wafai <amerwafai@gmail.com>
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
 * Molecul container interface definition.
 */
#ifndef MOLECULECONTAINER_H
#define MOLECULECONTAINER_H

#include "config.h"
#include "typedefs.h"

#include "simulation.h"
#include "ensemble.h"
#include "domain.h"
#include "molecule.h"

typedef enum {
    MOLECULEBLOCKS = 0,
    BASICN2,
    mc_type_MAX
} mc_type_t;



typedef struct mc_t {
    mc_type_t type;
    void *container;
} mc_t;

typedef struct mc_methods_list_t {
    void* (* alloc)();
    void (* destroy )(void *mc);
    void (* estimate_halos)( void *mc, domain_t * domain);
    void (* free)(void *mc);
    void (* init)(void *mc, ensemble_t *ensemble, domain_t * domain); 
    void (* decompose )(void *mc, domain_t *domain);
    void (* update)(void *mc, domain_t *domain);
    void (* print_ascii)(void *mc, char *filename);
    void (* print_pov)(void *mc, char *filename);
    void (* add_molecule)(void *mc, molecule_t *m);
    long (* get_num_molecules)(void *mc);
    void (* calc_forces)(void *mc, real *U_pot);
    void (* print_stats)(void *mc, ensemble_t *ensemble);
    void (* clear_halo)(void *mc);
    void (* reset_forces_and_momenta)(void *mc);
    void (* shift_velocity)(void *mc, real v[3]);
    void (* update_halo)(void *mc, domain_t *domain);
    void (* integrate_pref)(void *mc, simulation_t *simulation);
    void (* integrate_postf)(void *mc, simulation_t *simulation);
    void (* copyBack)(void *mc);
    real (* calc_E_kin)(void *mc);
    void (* apply_thermostat)(void *mc, real T_cur, real T_target);
    void (* calc_center_of_mass_velocity)(void *mc, real v[3]);
} mc_methods_list_t;


mc_t *mc_alloc(mc_type_t type); /**< allocate a new molecule container */
void mc_free(mc_t *mc); /**< free molecule container */
void mc_destroy(mc_t *mc); /**< destroy molecule container */

void mc_init(mc_t *mc, ensemble_t *ensemble, domain_t * domain); 
void mc_update(mc_t *mc, domain_t *domain);
void mc_print_ascii(mc_t *mc, char *filename);
void mc_print_pov(mc_t *mc, char *filename);
void mc_decompose(mc_t * mc, domain_t * domain);
void mc_estimate_halos( mc_t * mc, domain_t * domain);

void mc_add_molecule(mc_t *mc, molecule_t *m);
long mc_get_num_molecules(mc_t *mc);

void mc_calc_forces(mc_t *mc, real *U_pot);

void mc_print_stats(mc_t *mc, ensemble_t *ensemble);
void mc_clear_halo(mc_t *mc);
void mc_reset_forces_and_momenta(mc_t *mc);
void mc_shift_velocity(mc_t *mc, real v[3]);

void mc_update_halo(mc_t *mc, domain_t *domain);
void mc_integrate_pref(mc_t *mc, simulation_t *simulation);
void mc_integrate_postf(mc_t *mc, simulation_t *simulation);
void mc_copyBack(mc_t *mc);
real mc_calc_E_kin(mc_t *mc);
void mc_apply_thermostat(mc_t *mc, real T_cur, real T_target);

void mc_calc_center_of_mass_velocity(mc_t *mc, real v[3]);

mc_type_t get_mc_type(char *typename);
#endif /* MOLECULECONTAINER_H */
