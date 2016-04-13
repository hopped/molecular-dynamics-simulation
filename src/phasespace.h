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

#ifndef PHASESPACE_H
#define PHASESPACE_H

#include "typedefs.h"
#include "molecule.h"
#include "moleculecontainer/moleculecontainer.h"
#include "ensemble.h"
#include "simulation.h"
#include "domain.h"

/**
 * Phase spcae object
 */
typedef struct {
    mc_t *mc; /**< Molecule container holding the moleucles of the phase spce. */
} phasespace_t;

phasespace_t* psp_alloc(); /**< Allocate storage for phase space object. */
void psp_free(phasespace_t *psp); /**< Free storage of a phase space object. */
void psp_init(phasespace_t *psp, ensemble_t *ensemble, domain_t * domain); /**< Initialize phase space to fit ensemble */
void psp_update_halo(phasespace_t *psp, domain_t* domain); /**< Update halos of phase space*/
void psp_update(phasespace_t *psp, domain_t* domain); /**< Update of phase space*/
void psp_decompose(phasespace_t *psp, domain_t * domain);
void psp_destroy(phasespace_t * psp);
void psp_estimate_halos( phasespace_t * psp, domain_t * domain); 

void psp_generate(phasespace_t *psp);
void psp_integrate_pref(phasespace_t *psp, simulation_t *simulation, domain_t *domain); /**< Execute stuff before force calculation. */
void psp_integrate_postf(phasespace_t *psp, simulation_t *simulation, domain_t *domain); /**< Execute stuff after force calculation. */
void psp_copyBack(phasespace_t *psp);
void psp_calc_forces(phasespace_t *psp, real *U_pot); /**< Calculate foreces. */

/** Add moleucle to phase space 
 * @param[in,out] psp        phase spcae container to which molecule shall be added
 * @param[in]     m          molecule to be added
 */
void psp_add_molecule(phasespace_t *psp, molecule_t *m); /**< Add a molecule to the phase space. */

void psp_calc_ensemble_values(phasespace_t *psp, domain_t *domain, ensemble_t *ensemble); /**< Calculate ensemble values for the given phase space. */

void psp_reset_forces_and_momenta(phasespace_t *psp); /**< Reset forces and momenta of molecules in the phase space to zero. */

/** Apply thermostat
 * @param[in,out] psp        phase spcae container which shall be thermalized
 * @param[in]     T_cur      current temperature
 * @param[in]     T_target   target temperature
 */
void psp_apply_thermostat(phasespace_t *psp, real T_cur, real T_target);
#endif /* PHASESPACE_H */
