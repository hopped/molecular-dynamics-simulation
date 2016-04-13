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

#include "moleculecontainer.h"

#include "moleculecontainer/moleculeblocks/moleculeblocks.h"
#include "moleculecontainer/basicN2/basicN2.h"
#include "utils/logger.h"

#include <stdlib.h>
#include <string.h>

static mc_methods_list_t mc_method_list[] = { 
[MOLECULEBLOCKS] = {
    .alloc = mb_alloc,
    .free = mb_free,
    .destroy = mb_destroy,
    .estimate_halos = mb_estimate_halos,
    .decompose = mb_decompose,
    .init = mb_init,
    .copyBack = mb_copyBack,
    .update = mb_update,
    .print_ascii = mb_print_ascii,
    .print_pov = mb_print_pov,
    .add_molecule = mb_add_molecule,
    .get_num_molecules = mb_get_num_molecules,
    .calc_forces = mb_calc_forces,
    .print_stats = mb_print_stats,
    .clear_halo = mb_clear_halo,
    .reset_forces_and_momenta = mb_reset_forces_and_momenta,
    .shift_velocity = mb_shift_velocity,
    .update_halo = mb_update_halo,
    .integrate_pref = mb_integrate_pref,
    .integrate_postf = mb_integrate_postf,
    .calc_E_kin = mb_calc_E_kin,
    .apply_thermostat = mb_apply_thermostat,
},
[BASICN2] = {
    .alloc = basicN2_alloc,
    .destroy = basicN2_destroy,
    .estimate_halos = basicN2_estimate_halos,
    .decompose = basicN2_decompose,
    .free = basicN2_free,
    .init = basicN2_init,
    .copyBack = basicN2_copyBack,
    .update = basicN2_update,
    .print_ascii = basicN2_print_ascii,
    .print_pov = basicN2_print_pov,
    .add_molecule = basicN2_add_molecule,
    .get_num_molecules = basicN2_get_num_molecules,
    .calc_forces = basicN2_calc_forces,
    .print_stats = basicN2_print_stats,
    .clear_halo = basicN2_clear_halo,
    .reset_forces_and_momenta = basicN2_reset_forces_and_momenta,
    .shift_velocity = basicN2_shift_velocity,
    .update_halo = basicN2_update_halo,
    .integrate_pref = basicN2_integrate_pref,
    .integrate_postf = basicN2_integrate_postf,
    .calc_E_kin = basicN2_calc_E_kin,
    .apply_thermostat = basicN2_apply_thermostat,
}
};


mc_type_t get_mc_type(char *typename) {
    mc_type_t type;
    if( strcmp("MOLECULEBLOCKS", typename) == 0 ) {
        type = MOLECULEBLOCKS;
    } 
    else if( strcmp("BASICN2", typename) == 0 ) {
        type = BASICN2;
    }
    else {
        LOG(Error, "Unknown molecule container type '%s'\n", typename);
        type = mc_type_MAX;
    }
    return type;
}

mc_t *mc_alloc(mc_type_t type) { 
    mc_t * mc = (mc_t *) malloc (sizeof(mc_t));
    mc->type = type;
    mc->container = mc_method_list[type].alloc(); 
    return mc;
}

void mc_free(mc_t *mc) { 
    mc_method_list[mc->type].free(mc->container);
    free(mc);
}
void mc_estimate_halos( mc_t *mc, domain_t * domain){mc_method_list[mc->type].estimate_halos(mc->container, domain);}
void mc_destroy(mc_t *mc) {mc_method_list[mc->type].destroy(mc->container);}
void mc_decompose(mc_t *mc, domain_t * domain) {mc_method_list[mc->type].decompose(mc->container, domain);}
void mc_init(mc_t *mc, ensemble_t *ensemble, domain_t * domain) { mc_method_list[mc->type].init(mc->container, ensemble, domain); }
void mc_update(mc_t *mc, domain_t *domain) { mc_method_list[mc->type].update(mc->container, domain); }
void mc_print_ascii(mc_t *mc , char *filename) { mc_method_list[mc->type].print_ascii(mc->container, filename); }
void mc_print_pov(mc_t *mc , char *filename) { mc_method_list[mc->type].print_pov(mc->container, filename); }

void mc_add_molecule(mc_t *mc , molecule_t *m) { mc_method_list[mc->type].add_molecule(mc->container, m); }
long mc_get_num_molecules(mc_t *mc) { return mc_method_list[mc->type].get_num_molecules(mc->container); }

void mc_calc_forces(mc_t *mc , real *U_pot) { mc_method_list[mc->type].calc_forces(mc->container, U_pot); }
real mc_calc_E_kin(mc_t *mc) { return mc_method_list[mc->type].calc_E_kin(mc->container); }

void mc_print_stats(mc_t *mc , ensemble_t *ensemble) { mc_method_list[mc->type].print_stats(mc->container, ensemble); }
void mc_clear_halo(mc_t *mc) { mc_method_list[mc->type].clear_halo(mc->container); }
void mc_reset_forces_and_momenta(mc_t *mc) { mc_method_list[mc->type].reset_forces_and_momenta(mc->container); }
void mc_shift_velocity(mc_t *mc , real v[3]) { mc_method_list[mc->type].shift_velocity(mc->container, v); }

void mc_update_halo(mc_t *mc, domain_t *domain) { mc_method_list[mc->type].update_halo(mc->container,domain); }
void mc_integrate_pref(mc_t *mc , simulation_t *simulation) { mc_method_list[mc->type].integrate_pref(mc->container, simulation); }
void mc_integrate_postf(mc_t *mc , simulation_t *simulation) { mc_method_list[mc->type].integrate_postf(mc->container, simulation); }
void mc_copyBack(mc_t *mc) { mc_method_list[mc->type].copyBack(mc->container); }
void mc_apply_thermostat(mc_t *mc , real T_cur, real T_target) { mc_method_list[mc->type].apply_thermostat(mc->container, T_cur, T_target); }
