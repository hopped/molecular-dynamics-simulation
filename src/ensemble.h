/*
 * Copyright (c) 2012      Christoph Niethammer <christoph.niethammer@googlemail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include "config.h"
#include "typedefs.h"


/** ensemble type
 *
 * The ensemble type is encoded as a bit field where an extensive variable is
 * encoded by an unset bit and an intensive variables by a set bit.
 */
typedef enum ensemble_type_t { 
    NVE=0x0, 
    NVT=0x1, 
    NPE=0x2, 
    NPT=0x3, 
    MVE=0x4, 
    MVT=0x5, 
    MPE=0x6, 
    MPT=0x7,
    _ENSEMBLE_TYPE_MAX
} ensemble_type_t;

const char * ensemble_get_name(ensemble_type_t ensemble); /**< Get the common physical name of ensemble. */
const char * ensemble_get_acronym(ensemble_type_t ensemble); /** Get the 3 letter acronym of the ensemble. */


typedef struct ensemble_t {
    size_t N; /**< number of molecules */
    size_t NPS; /**< number of molecules per subdomain */
    size_t CPU_PART; /**< number of subdomain running on CPU */
    real mu;  /**< chemical potential*/
    real V;   /**< volume */
    real p;   /**< pressure */
    real E;   /**< energy */
    real T;   /**< temperature */
    real U_pot;  /**< potential energy of the interactions (cache) */
    real E_kin;  /**< kinetic energy of molecules (cache) */
	/* 64 byte cache line boundary if type of real is double */
    real virial; /**< virial */
    enum ensemble_type_t  type;
} ensemble_t;

void ensemble_init(ensemble_t *ensemble);
void ensemble_print_info(ensemble_t *ensemble); /**< Print ensemble info to stdout and global logfile. */

void ensemble_print_config(ensemble_t *ensemble); /**< Print ensemble config parameters to logfile */
#endif /* ENSEMBLE_H */
