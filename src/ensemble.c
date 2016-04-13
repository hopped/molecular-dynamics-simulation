/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#define ENSEMBLE_SRC
#include "ensemble.h"

#include "typedefs.h"
#include "utils/logger.h"

#include <stdio.h>


//ensemble_t ensemble;

/* global energies */
real U_pot = 0.; /**< global potential energy */
real U_kin = 0.; /**< global kinetic energy */
real U = 0.; /**< global energy, sum of potential and kinetic energy */

static const char * const ensemble_names[_ENSEMBLE_TYPE_MAX] = {
    [NVE] = "microcanonical ensemble",
    [NVT] = "canonical ensemble", 
    [NPE] = "", 
    [NPT] = "isothermal-isobaric ensemble", 
    [MVE] = "", 
    [MVT] = "grand canonical ensemble", 
    [MPE] = "", 
    [MPT] = ""
};
	
static const char * const ensemble_acronyms[_ENSEMBLE_TYPE_MAX] = {
    [NVE] = "NVE", 
    [NVT] = "NVT", 
    [NPE] = "NPE", 
    [NPT] = "NPT",
    [MVE] = "MVE",
    [MVT] = "MVT",
    [MPE] = "MPE",
    [MPT] = "MPT"
};

const char * ensemble_get_name(ensemble_type_t ensemble) {
    return ensemble_names[ensemble];
}

const char * ensemble_get_acronym(ensemble_type_t ensemble) {
    return ensemble_acronyms[ensemble];
}

void ensemble_init(ensemble_t *ensemble) {
    ensemble->N  = 0;
    ensemble->NPS  = 0;
    ensemble->mu = 0.;
    ensemble->V = 0.;
    ensemble->p = 0.;
    ensemble->E = 0.;
    ensemble->T = 0.;
    ensemble->U_pot = 0.;
    ensemble->E_kin = 0.;
    ensemble->virial = 0.;
    ensemble->CPU_PART = 0;
    ensemble->type = NVE;
}

void ensemble_print_info(ensemble_t *ensemble) {
    LOG(Info, "N=%lu, mu=%"PRIreal", V=%"PRIreal", p=%"PRIreal", E=%"PRIreal", T=%"PRIreal", U_pot=%"PRIreal", E_kin=%"PRIreal", NPS=%lu, CPU_PART=%lu\n",
            ensemble->N, ensemble->mu, ensemble->V, ensemble->p, ensemble->E, ensemble->T, ensemble->U_pot/ensemble->N, ensemble->E_kin, ensemble->NPS, ensemble->CPU_PART);
    printf("N=%lu, mu=%"PRIreal", V=%"PRIreal", p=%"PRIreal", E=%"PRIreal", T=%"PRIreal", U_pot=%"PRIreal", E_kin=%"PRIreal", NPS=%lu, CPU_PART=%lu\n",
            ensemble->N, ensemble->mu, ensemble->V, ensemble->p, ensemble->E, ensemble->T, ensemble->U_pot/ensemble->N, ensemble->E_kin, ensemble->NPS, ensemble->CPU_PART);
}

void ensemble_print_config(ensemble_t *ensemble) {
    LOG(Info, "Ensemlble type: %s (%s)\n", ensemble_get_acronym(ensemble->type), ensemble_get_name(ensemble->type));
    LOG(Info, "N=%lu, mu=%"PRIreal", V=%"PRIreal", p=%"PRIreal", E=%"PRIreal", T=%"PRIreal"\n",
            ensemble->N, ensemble->mu, ensemble->V, ensemble->p, ensemble->E, ensemble->T);
}








