/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 * Copyright (c) 2012-2014      Mohamad Amer Wafai <amerwafai@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#define PHASESPACE_SRC

#include "config.h"
#include "typedefs.h"
#include "constants.h"
#include "lj_parameter.h"

#include "phasespace.h"
#include "moleculecontainer/moleculecontainer.h"
#include "ensemble.h"
#include "utils/logger.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

phasespace_t* psp_alloc() {
    phasespace_t* psp = (phasespace_t *) malloc(sizeof(phasespace_t));
    assert(psp != NULL);
    psp->mc = mc_alloc(config.molecule_container_type);
    return psp;
}

void psp_free(phasespace_t *psp) {
    mc_free(psp->mc);
    free(psp);
}

void psp_destroy(phasespace_t * psp){
    mc_destroy(psp->mc);
}

void psp_decompose(phasespace_t *psp, domain_t * domain){
    mc_decompose(psp->mc, domain);
}

void psp_init(phasespace_t *psp, ensemble_t *ensemble, domain_t * domain) {
    mc_init(psp->mc, ensemble, domain);
}

void psp_estimate_halos( phasespace_t * psp, domain_t * domain){
     mc_estimate_halos(psp->mc, domain);
}

void psp_generate(phasespace_t *psp){
}

void psp_update(phasespace_t *psp, domain_t* domain){
    /* Now update the particle container as molecules moved to their final position. */
    mc_update(psp->mc, domain);
}

void psp_update_halo(phasespace_t *psp, domain_t* domain){
    /* Now update the halos. This may include exchanging molecules with neigbours. */
    mc_update_halo(psp->mc, domain);
}

void psp_integrate_pref(phasespace_t *psp, simulation_t *simulation, domain_t *domain){
    mc_integrate_pref(psp->mc, simulation);
}

void psp_integrate_postf(phasespace_t *psp, simulation_t *simulation, domain_t *domain){
    /* Complete the integration step. */
    mc_integrate_postf(psp->mc, simulation);
}
void psp_copyBack(phasespace_t *psp){
    mc_copyBack(psp->mc);
}
void psp_calc_forces(phasespace_t *psp, real *U_pot){
    mc_calc_forces(psp->mc, U_pot);
}

void psp_add_molecule(phasespace_t *psp, molecule_t *m){
    mc_add_molecule(psp->mc, m);
}

void psp_calc_ensemble_values(phasespace_t *psp, domain_t *domain, ensemble_t *ensemble) {
    ensemble->N = mc_get_num_molecules(psp->mc);
    ensemble->V = domain_get_volume(domain);
    ensemble->E_kin = mc_calc_E_kin(psp->mc);
    ensemble->T = ensemble->E_kin / (3./2. * kB * ensemble->N) ;
    ensemble->E = ensemble->E_kin + ensemble->U_pot;
    /* Long range corrections */
    /* standard LJ long range correction: */
    real rc = config.cutoff_radius;
    if( rc <= 0 ) {
        LOG(Warning, "Cutoff radius <= 0 (rc = %"PRIreal"), cannot determin long range cutoff correction.\n", rc);
    }
    else {
        real rho = ensemble->N / ensemble->V;
        real epsilon = epsilon24 / 24.;
        real sigma = sqrt(sigma2);
        long N = ensemble->N;
        real UpotCorrLJ = N * 8./3. * M_PI * rho * epsilon * pow(sigma,3.0) * (pow(sigma/rc,9)/3.0 - pow(sigma/rc,3));
        LOG(Info, "U_pot LJ long range correction: %"PRIreal"\n", UpotCorrLJ);
        ensemble->U_pot += UpotCorrLJ;
    }
}

void psp_reset_forces_and_momenta(phasespace_t *psp) {
    mc_reset_forces_and_momenta(psp->mc);
}

void psp_apply_thermostat(phasespace_t *psp, real T_cur, real T_target) {
    mc_apply_thermostat(psp->mc, T_cur, T_target);
}
